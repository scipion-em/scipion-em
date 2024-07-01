import logging
from functools import lru_cache

import numpy

from pwem.emlib.image.image_handler import ImageReadersRegistry

logger = logging.getLogger(__name__)

import metadataviewer
from metadataviewer.dao.model import IDAO
from metadataviewer.model import Table, Column, StrRenderer, IntRenderer, FloatRenderer,Page
from metadataviewer.model.renderers import ImageReader
from pwem.viewers.mdviewer.star_dao import StarFile
from pwem.viewers.mdviewer.sqlite_dao import ScipionSetsDAO
class ScipionImageReader(ImageReader):
    @classmethod
    def getCompatibleFileTypes(cls) -> list:
        return list(ImageReadersRegistry._readers.keys())

    @classmethod
    @lru_cache
    def open(cls, path):
        ext = path.split('.')[-1]
        imageReader = ImageReadersRegistry._readers[ext]
        return imageReader.open(path)






class NumpyDao(IDAO):
    """ DAO for the metadata viewer to read numpy files """
    def __init__(self, filename: str):
        self.file = filename
        self._data = None

    def getData(self)-> numpy.ndarray:
        if self._data is None:
            self._data = numpy.load(self.file)

        return self._data

    def fillPage(self, page: Page, actualColumn: int, orderAsc: bool) -> None:

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
            sortingColIdx= table.getSortingColumnIndex()
            sortingCol = data.dtype.descr[sortingColIdx][0]
            sortingAsc = table.isSortingAsc()
            data.sort(order=sortingCol)

            if table.isSortingAsc():
                data=numpy.flipud(data)


        for i in range(firstRow,limit):
            row = data[i]
            values = [value for value in row]
            id = i

            page.addRow((int(id), values))

    def fillTable(self, table: Table, objectManager) -> None:

        for descr in self.getData().dtype.descr:
            field, ftype = descr[0], descr[1]

            if len(descr)==3:
                ftype="S"
            newCol = Column(name=field, renderer=self.getRenderer(ftype))
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



def extendMDViewer(om: metadataviewer.model.ObjectManager):
    """ Function to extend the object manager with DAOs and readers"""
    om.registerDAO(ScipionSetsDAO)
    om.registerReader(ScipionImageReader)
    om.registerDAO(NumpyDao)
    om.registerDAO(StarFile)