# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *          Pablo Conesa Mingo         (pconesa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import logging

from metadataviewer.model import Table

logger = logging.getLogger()

from metadataviewer.dao.model import IDAO


class StarFile(IDAO):
    """
    Class to handle STAR or XMD files.
    """
    def __init__(self, inputFile):
        self._file = self.__loadFile(inputFile)
        self._tableCount = {}
        self._tableData = {}
        self._labels = {}
        self._tables = {}
        self._labelsTypes = {}

    def __loadFile(self, inputFile):
        try:
            return open(inputFile, 'r')
        except Exception as e:
            logger.error("The file could not be opened. Make sure the path is "
                         "correct: \n %s" % e)
            return None

    def fillTable(self, table, objectManager):
        """Create the table structure"""
        self._loadStarFileInfo(table)

    def fillPage(self, page, actualColumn, orderAsc):
        """
        Fill a page taking into account the page number and size.
        Variables actualColumn and orderAsc control the table order
        """
        table = page.getTable()
        tableName = table.getName()
        if table.hasSortingChanged():
            columnIndex = table.getColumnIndexFromLabel(table.getSortingColumn())
            self.sort(tableName, columnIndex, orderAsc)

        # moving to the first row of the page
        pageNumber = page.getPageNumber()
        pageSize = page.getPageSize()
        firstRow = pageNumber * pageSize - pageSize
        endRow = pageNumber * pageSize  # getting more rows

        logger.debug("Creating the page with rows: %d -> %d " % (firstRow, endRow))
        for row in self._iterRowLines(tableName, firstRow, endRow):
            if row:
                page.addRow(row)

    def getTableRowCount(self, tableName):
        """Get the number of rows of a given table"""
        return self._tableCount[tableName]

    def getTables(self):
        """ Return all the names of the data_ blocks found in the file and
            fill all labels and data of every table name """
        if not self._tables:
            f = self._file
            line = f.readline()
            logger.debug("Reading the star file completely and filling "
                         "the model().")
            while line:
                # While searching for a data line, we will store the offsets
                # for any data_ line that we find
                if line.startswith('data_'):
                    table = line.strip()
                    alias = table.replace('data_', '')
                    tbl = Table(table)
                    tbl.setAlias(alias)
                    self._tables[table] = tbl
                    data = self._getData()
                    self._tableCount[table] = data[0]
                    self._labels[table] = ['id'] + data[1]
                    self._labelsTypes[table] = data[2]
                    self._tableData[table] = data[3]
                line = f.readline()

        return self._tables

    def _getData(self):
        """ Method to get all information of the table (labels, data,...)"""
        self._findLabelLine()
        data = []
        line, labels = self._getLabels()
        count = 0
        f = self._file
        firstRow = line.split()
        labelsTypes = [int]
        for i in range(len(firstRow)):
            labelsTypes.append(_guessType(firstRow[i]))

        while line and not line.startswith('\n'):
            line = str(count+1) + ' ' + line
            data.append(line.split())
            count += 1
            line = f.readline().strip()
        return count, labels, labelsTypes, data

    def _getLabels(self):
        """Get the table labels"""
        logger.debug("Getting the table labels...")
        line = self._line
        labels = []
        while line.startswith('\n'):
            line = self._file.readline()
        while line.startswith('_'):
            parts = line.split()
            labels.append(parts[0][1:])
            line = self._file.readline()
        while line.startswith('\n'):
            line = self._file.readline()
        return line, labels

    def _loadStarFileInfo(self, table):
        """Create the table structure"""
        logger.debug("Creating the table columns...")
        colNames = self._labels[table.getName()]
        values = self._tableData[table.getName()][0]
        table.createColumns(colNames, values)
        table.setAlias(table.getAlias())

    def _findLabelLine(self):
        """Find the first labels line in the star file"""
        line = ''
        foundLoop = False

        rawLine = self._file.readline()
        while rawLine:
            if rawLine.startswith('_'):
                line = rawLine
                break
            elif rawLine.startswith('loop_'):
                foundLoop = True
            rawLine = self._file.readline()

        self._line = line.strip()
        self._foundLoop = foundLoop

    def _iterRowLines(self, tableName, firstRow, endRow):
        """Iter over the table in a range of rows """
        if self._tableCount[tableName] == 1:
            yield 1, self._tableData[tableName][0]
            return
        if firstRow + endRow > self._tableCount[tableName]:
            endRow = self._tableCount[tableName]
        for i in range(firstRow, endRow):
            values = self._tableData[tableName][i]
            yield int(values[0]), values

    def close(self):
        if getattr(self, '_file', None):
            self._file.close()
            self._file = None

    @classmethod
    def getCompatibleFileTypes(cls):
        """Return a list of compatible extension of files"""
        logger.debug("Selected StarFile DAO")
        return ['star', 'xmd']

    def sort(self, tableName, column, sortAsc=True):
        """ Sort the table in place using the provided column.
            :param column is a number, it is the index of one column. """
        _columType = self._labelsTypes[tableName][column]
        orderList = sorted(self._tableData[tableName],
                           key=lambda x: _columType(x[column]),
                           reverse=not sortAsc)
        self._tableData[tableName] = orderList

    def getSelectedRangeRowsIds(self, tableName, startRow, numberOfRows, column, reverse=True):
        """Return a range of rows starting at 'startRow' an amount of
           'numberOfRows' """
        logger.debug("Reading the table %s and selected a range of rows %d - %d" % (tableName, startRow+1, numberOfRows + 1))
        col = 0
        for i in range(len(self._labels[tableName])):
            if self._labels[tableName][i] == 'id':
                col = i
                break
        table = self._tableData[tableName]
        rowsIds = [int(table[row][col]) for row in range(startRow-1, startRow+numberOfRows)]
        return rowsIds

    def getTableWithAdditionalInfo(self):
        """Return a tuple with the table that need to show additional info and
        the column that we need to show"""
        return None

# -----------------------------------Utils methods -----------------------


def _guessType(strValue):
    try:
        int(strValue)
        return int
    except ValueError:
        try:
            float(strValue)
            return float
        except ValueError:
            return str





