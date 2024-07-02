import logging
import numpy
logger = logging.getLogger(__name__)
from metadataviewer.dao.model import IDAO
from metadataviewer.model import Table, Column, StrRenderer, IntRenderer, FloatRenderer,Page


class NumpyDao(IDAO):
    """ DAO for the metadata viewer to read numpy files """
    def __init__(self, filename: str):
        self.file = filename
        self._data = None

    def getData(self)-> numpy.ndarray:
        if self._data is None:
            self._data = numpy.load(self.file)

        return self._data

    def fillPage(self, page: Page, actualColumn: str, orderAsc: bool) -> None:

        # moving to the first row of the page
        pageNumber = page.getPageNumber()
        pageSize = page.getPageSize()
        firstRow = pageNumber * pageSize - pageSize
        data = self.getData()
        limit = min(firstRow+pageSize, data.size)
        table = page.getTable()

        # If sorting changed
        if table.hasSortingChanged():
            # do the sorting
            sortingColIdx = table.getColumnIndexFromLabel(table.getSortingColumn())
            sortingCol = data.dtype.descr[sortingColIdx][0]
            data.sort(order=sortingCol)

            if table.isSortingAsc():
                data = numpy.flipud(data)

        for i in range(firstRow, limit):
            row = data[i]
            values = [value for value in row]
            id = i
            page.addRow((int(id), values))

    def fillTable(self, table: Table, objectManager) -> None:

        for descr in self.getData().dtype.descr:
            field, fType = descr[0], descr[1]
            if len(descr) == 3:
                fType = "S"
            newCol = Column(name=field, renderer=self.getRenderer(fType))
            table.addColumn(newCol)

    def getRenderer(self, fType):
        if fType.startswith("<u"):
            return IntRenderer()
        elif fType.startswith("<f"):
            return FloatRenderer()
        else:
            return StrRenderer()

    def getTables(self):
        return {"data": Table("data")}

    @staticmethod
    def getCompatibleFileTypes() -> list:
        return ["cs"]

    def getTableRowCount(self, tableName: str) -> int:
        return self.getData().size

    def getTableWithAdditionalInfo(self):
        pass

    def getSelectedRangeRowsIds(self, tableName, startRow, numberOfRows, column, reverse=True) -> list:
        pass
