# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
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
from datetime import datetime

logger = logging.getLogger()

import sqlite3
from metadataviewer.dao.model import IDAO, Table

EXTENDED_COLUMN_NAME = '_representative'
ALLOWED_COLUMNS_TYPES = ['String', 'Float', 'Integer', 'Boolean', 'Matrix']
ADITIONAL_INFO_DISPLAY_COLUMN_LIST = ['_size', 'id']
EXCLUDED_COLUMNS = ['label', 'comment', 'creation', '_streamState']
CLASS_OBJECT = 1
REPRESENTATIVE_OBJECT = 2
CLASS_ELEMENTS = 3


class SqliteFile(IDAO):
    """  Class to manipulate Scipion Sqlite files. """
    userSubsetCreationCallback = None

    def __init__(self, sqliteFile):
        self._names = []
        self._file = sqliteFile
        self._con = self.__loadDB(sqliteFile)
        self._con.row_factory = self._dictFactory
        self._tableCount = {}
        self._tables = {}
        self._labels = {}
        self._labelsTypes = {}
        self._aliases = {}
        self._columnsMap = {}
        self._extendedColumn = None
        self._tableWithAdditionalInfo = None
        self._objectsType = {}

    def __loadDB(self, sqliteFile):
        """Load a sqlite file"""
        try:
            return sqlite3.connect(f"file:{sqliteFile}?mode=ro", uri=True)
        except Exception as e:
            logger.error("The file could not be opened. Make sure the path is "
                         "correct: \n %s" % e)
            return None

    def hasExtendedColumn(self):
        """Return if the table need to extend a column. That column is used to
        renderer an image that is composed by other two columns"""
        return self._extendedColumn is not None

    def composeDataTables(self, tablesNames):
        """This method is used to generate a dictionary with the principal
           tables mapping the dependencies with other tables"""
        tablesNames = sorted(tablesNames)
        for tableName in tablesNames:
            divTable = tableName.split('_')
            if len(divTable) > 1:
                if divTable[0].startswith('Class') and divTable[1].startswith('Class') and tableName not in self._tables:
                    objectTable = divTable[0] + '_Objects'
                    self._tables[objectTable] = tableName
                    self._names.append(objectTable)

    def composeObjectType(self):
        """Define the different objects types"""
        # General type defined into Properties table
        firstRow = self.getTableRow('Properties', 0)
        objectType = firstRow['value']
        self._objectsType[self._aliases['objects']] = objectType

        for alias in self._aliases.values():
            objectTypeAux = alias.split('_')
            if len(objectTypeAux) == 2:
                sufix = 's' if objectTypeAux[0][-1] in "aeiouAEIOU" or objectTypeAux[0][-1] != 's' else ''
                objectType = 'SetOf%s%s' % (objectTypeAux[1], sufix)

                if objectTypeAux[1] not in self._objectsType:
                    self._objectsType[objectTypeAux[1]] = objectType

    def composeTableAlias(self, tableName):
        """Create an alias for the given table"""
        firstRow = self.getTableRow(tableName, 0)
        className = firstRow['class_name']
        if tableName.__contains__('_'):
            alias = tableName.split('_')[0] + '_' + className
        else:
            alias = className
        return alias

    def getTableNames(self):
        """ Return all the table names found in the database. """
        if not self._names:
            self._tables = {'objects': 'classes'}
            self._names = ['objects']

            res = self._con.execute("SELECT name FROM sqlite_master WHERE type='table'")
            tablesNames = [row['name'] for row in res.fetchall()]
            self.composeDataTables(tablesNames)

            for tableName in self._names:
                # Getting the first row to identify the labels and theirs type
                firstRow = self.getTableRow(tableName, 0, classes=self._tables[tableName])
                self._labels[tableName] = [key for key in firstRow.keys() if key not in EXCLUDED_COLUMNS]
                alias = self.composeTableAlias(self._tables[tableName])
                self._aliases[tableName] = alias

                labelsTypes = []
                self._tableCount[tableName] = self.getRowsCount(tableName)
                for key, value in firstRow.items():
                    if key not in EXCLUDED_COLUMNS:
                        labelsTypes.append(_guessType(value))
                self._labelsTypes[tableName] = labelsTypes

            self.composeObjectType()

        if len(self._tables) > 1:
            self._tableWithAdditionalInfo = 'objects'
        return self._names

    def findColbyName(self, colNames, colName):
        """Return a column index given a column name"""
        for i, col in enumerate(colNames):
            if colName == col:
                return i
        return None

    def updateExtendColumn(self, table):
        """Find the columns that need to extend and keep the indexes"""
        tableName = table.getName()
        colNames = self._labels[tableName]
        indexCol = self.findColbyName(colNames, '_index')
        representativeCol = self.findColbyName(colNames, '_filename')

        if indexCol and representativeCol:
            logger.debug("The columns _index and _filename have been found. "
                         "We will proceed to create a new column with the "
                         "values of these columns.")
            self._extendedColumn = indexCol, representativeCol
        else:
            indexCol = self.findColbyName(colNames, '_representative._index')
            representativeCol = self.findColbyName(colNames,
                                                   '_representative._filename')
            if indexCol and representativeCol:
                logger.debug("The columns _representative._index and "
                             "_representative._filename have been found. "
                             "We will proceed to create a new column with the "
                             "values of these columns.")
                self._extendedColumn = indexCol, representativeCol

    def generateTableActions(self, table, objectManager):
        """Generate actions for a given table in order to create subsets"""
        alias = table.getAlias()
        labels = list(self._objectsType.keys())
        objectTypes = list(self._objectsType.values())
        aliasSplit = alias.split('_')
        if alias.startswith('Class') and len(aliasSplit) == 1:

            table.addAction(labels[0], lambda: self.createSubsetCallback(table, objectTypes[0], objectManager))
            table.addAction(labels[1], lambda: self.createSubsetCallback(table, objectTypes[1], objectManager))
            if alias == 'Class2D':
                table.addAction('Averages', lambda: self.createSubsetCallback(table, 'SetOfAverages', objectManager))
            else:
                table.addAction('Volumes', lambda: self.createSubsetCallback(table, 'SetOfVolumes',  objectManager))

        else:
            if len(aliasSplit) == 2:
                table.addAction(aliasSplit[1], lambda: self.createSubsetCallback(table, self._objectsType[aliasSplit[1]], objectManager))
            else:
                table.addAction(aliasSplit[0], lambda: self.createSubsetCallback(table, self._objectsType[aliasSplit[0]], objectManager))

    def fillTable(self, table, objectManager):
        """Create the table structure (columns) and set the table alias"""
        tableName = table.getName()
        colNames = self._labels[tableName]
        self.updateExtendColumn(table)

        row = self.getTableRow(tableName, 0, classes=self._tables[tableName])
        values = [value for key, value in row.items() if key not in EXCLUDED_COLUMNS]
        if self._extendedColumn:
            logger.debug("Creating an extended column: %s" % EXTENDED_COLUMN_NAME)
            colNames.insert(self._extendedColumn[1] + 1, EXTENDED_COLUMN_NAME)
            values.insert(self._extendedColumn[1] + 1, str(values[self._extendedColumn[0]]) + '@' + str(values[self._extendedColumn[1]]))
        table.createColumns(colNames, values)
        table.setAlias(self._aliases[tableName])
        self.generateTableActions(table, objectManager)

    def fillPage(self, page, actualColumn=0, orderAsc=True):
        """
        Read the given table from the sqlite and fill the page(add rows)
        """
        tableName = page.getTable().getName()
        # moving to the first row of the page
        pageNumber = page.getPageNumber()
        pageSize = page.getPageSize()
        firstRow = pageNumber * pageSize - pageSize
        limit = pageSize

        column = self._labels[tableName][actualColumn]
        mode = 'ASC' if orderAsc else 'DESC'
        self.updateExtendColumn(page.getTable())

        for row in self.iterTable(tableName, start=firstRow, limit=limit,
                                  classes=self._tables[tableName],
                                  orderBy=column, mode=mode):
            if row:
                if 'id' in row.keys():
                    id = row['id']
                    values = [value for key, value in row.items() if key not in EXCLUDED_COLUMNS]
                    # Checking if exists an extended column
                    if self.hasExtendedColumn():
                        values.insert(self._extendedColumn[1] + 1,
                                      str(values[self._extendedColumn[0]]) + '@' + values[self._extendedColumn[1]])
                    page.addRow((int(id), values))

    def getRowsCount(self, tableName):
        """ Return the number of elements in the given table. """
        logger.debug("Reading the table %s" %tableName)
        return self._con.execute(f"SELECT COUNT(*) FROM {tableName}").fetchone()['COUNT(*)']

    def getTableRowCount(self, tableName):
        return self._tableCount[tableName]

    def getSelectedRangeRowsIds(self, tableName, startRow, numberOfRows, column, reverse=True):
        """Return a range of rows starting at 'startRow' an amount
           of 'numberOfRows' """

        logger.debug("Reading the table %s and selected a range of rows %d - %d" % (tableName,startRow, numberOfRows + 1))
        mode = 'ASC' if reverse else 'DESC'
        col = self._getColumnMap(tableName, column)
        if col == None:
            col = column
        query = "SELECT id FROM %s ORDER BY %s %s LIMIT %d , %d" % (tableName, col, mode, startRow - 1, numberOfRows + 1)
        rowsList = self._con.execute(query).fetchall()
        rowsIds = [row['id'] for row in rowsList]
        return rowsIds

    def iterTable(self, tableName, **kwargs):
        """
        Method to iterate over the table's rows
        :param tableName: the name of the table
        :param kwargs:
                limit: integer value to limit the number of elements
                start: start from a given element
                classes: read column names from a 'classes' table
                orderBy: clause to sort given a column name
                mode: sort direction ASC or DESC
        """
        query = f"SELECT * FROM {tableName}"

        if 'mode' in kwargs:
            if 'orderBy' in kwargs:
                if kwargs['orderBy']:
                    column = self._getColumnMap(tableName, kwargs['orderBy'])
                    if not column:
                        column = kwargs['orderBy']

                    query += f" ORDER BY {column}"

            if kwargs['mode']:
                query += f" {kwargs['mode']}"

        if 'start' in kwargs and 'limit' not in kwargs:
            kwargs['limit'] = -1

        if 'limit' in kwargs:
            query += f" LIMIT {kwargs['limit']}"

        if 'start' in kwargs:
            query += f" OFFSET {kwargs['start']}"

        if 'classes' not in kwargs:
            res = self._con.execute(query)
            while row := res.fetchone():
                yield row
        else:  # Mapping the column names and  including only the allowed columns
            self._columnsMap[tableName] = {row['column_name']: row['label_property']
                          for row in self.iterTable(kwargs['classes']) if row['class_name'] in ALLOWED_COLUMNS_TYPES}
            self._excludedColumns = {row['column_name']: row['label_property']
                          for row in self.iterTable(kwargs['classes']) if row['class_name'] not in ALLOWED_COLUMNS_TYPES}

            def _row_factory(cursor, row):
                fields = [column[0] for column in cursor.description]
                rowFact = {self._columnsMap[tableName].get(k, k): v for k, v in zip(fields, row) if k not in self._excludedColumns}
                return rowFact

            # Modify row factory to modify column names
            self._con.row_factory = _row_factory
            res = self._con.execute(query)
            while row := res.fetchone():
                yield row
            # Restore row factory
            self._con.row_factory = self._dictFactory

    def getTableAliases(self):
        """Return the tables aliases"""
        return self._aliases

    def _getColumnMap(self, tableName, column):
        """Get the column name that has been mapped"""
        for key, value in self._columnsMap[tableName].items():
            if value == column:
                return key
        return None

    def getTableRow(self, tableName, rowIndex, **kwargs):
        """ Get a given row by index. Extra args are passed to iterTable. """
        kwargs['start'] = rowIndex
        kwargs['limit'] = 1
        for row in self.iterTable(tableName, **kwargs):
            return row

    def getTableWithAdditionalInfo(self):
        """Return a tuple with the table that need to show additional info and
        the column that we need to show"""
        return self._tableWithAdditionalInfo, ADITIONAL_INFO_DISPLAY_COLUMN_LIST

    def createSubsetCallback(self, table: Table, objectType: str, objectManager):
        """Create a subset"""
        selection = table.getSelection().getSelection()
        tableName = table.getName()
        elementsCount = len(selection)
        if not elementsCount:
            elementsCount = self._tableCount[tableName]
        subsetName = objectManager.getGui().getSubsetName(objectType, elementsCount)
        if subsetName:
            format = '%Y%m%d%H%M%S'
            now = datetime.now()
            timestamp = now.strftime(format)
            path = 'Logs/selection_%s.txt' % timestamp
            self.writeSelection(table, path)
            if tableName != 'objects':
                path += ',%s' % tableName.split('Objects')[0]
            self.userSubsetCreationCallback(subsetName, path, objectType)

    def writeSelection(self, table: Table, path):
        """ Create a file with the selected rows ids"""
        tableName = table.getName()
        rowsIds = table.getSelection().getSelection().keys()
        if not rowsIds:
            rowsIds = [i+1 for i in range(self._tableCount[tableName])]
        try:
            with open(path, 'w') as file:
                for rowId in rowsIds:
                    file.write(str(rowId) + ' ')
                file.close()
            logger.debug(f"The file: {path} was created correctly.")
        except Exception as e:
            logger.error(f"Error creating the file: {e}")

    def close(self):
        if getattr(self, '_con', None):
            self._con.close()
            self._con = None

    def _dictFactory(self, cursor, row):
        fields = [column[0] for column in cursor.description]
        return {key: value for key, value in zip(fields, row)}

    def getCompatibleFileTypes(self):
        """Return a list of compatible extension of files"""
        logger.debug("Selected SqliteFile DAO")
        return ['sqlite']

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()


# --------- Helper functions  ------------------------
def _guessType(strValue):
    if strValue is None:
        return str('None')
    try:
        int(strValue)
        return int
    except ValueError:
        try:
            float(strValue)
            return float
        except ValueError:
            return str

