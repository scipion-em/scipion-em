import logging
import os
import sys
from functools import lru_cache
from subprocess import Popen

from pwem.emlib.image.image_handler import ImageReadersRegistry
from pwem.viewers.mdviewer.volumeViewer import launchViewer

logger = logging.getLogger(__name__)
from datetime import datetime

import sqlite3
import tempfile
import time

import pyworkflow as pw
import metadataviewer
from metadataviewer.dao.model import IDAO
from metadataviewer.model import Table, Column, BoolRenderer, ImageRenderer, StrRenderer
from metadataviewer.model.renderers import ImageReader, Action


SCIPION_OBJECT_ID = "SCIPION_OBJECT_ID"
SCIPION_PORT = "SCIPION_PORT"


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


ALLOWED_COLUMNS_TYPES = ['String', 'Float', 'Integer', 'Boolean', 'Matrix',
                         'CsvList']
ADITIONAL_INFO_DISPLAY_COLUMN_LIST = ['_size', 'id']
EXCLUDED_COLUMNS = ['label', 'comment', 'creation', '_streamState']
PERMANENT_COLUMNS = ['id', 'enabled']
CLASS_OBJECT = 1
REPRESENTATIVE_OBJECT = 2
CLASS_ELEMENTS = 3
EXTENDED_COLUMN_NAME = 'stack'
ENABLED_COLUMN = 'enabled'
PROPERTIES_TABLE = 'Properties'
OBJECT_TABLE = 'objects'


class SqliteFile(IDAO):
    """  Class to manipulate Scipion Sqlite files. """

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
                if divTable[-1].startswith('Class') and tableName not in self._tables:
                    objectTable = tableName.replace(divTable[-1], '') + 'Objects'
                    self._tables[objectTable] = tableName
                    self._names.append(objectTable)

        self._tables[PROPERTIES_TABLE] = PROPERTIES_TABLE
        self._names.append(PROPERTIES_TABLE)

    def composeObjectType(self):
        """Define the different objects types"""
        # General type defined into Properties table
        firstRow = self.getTableRow(PROPERTIES_TABLE, 0)
        objectType = firstRow['value']
        self._objectsType[self._aliases[OBJECT_TABLE]] = objectType

        for alias in self._aliases.values():
            objectTypeAux = alias.split('_')
            if len(objectTypeAux) == 2:
                sufix = 's' if objectTypeAux[0][-1] in "aeiouAEIOU" or objectTypeAux[0][-1] != 's' else ''
                objectType = 'SetOf%s%s' % (objectTypeAux[1], sufix)

                if objectTypeAux[1] not in self._objectsType:
                    self._objectsType[objectTypeAux[1]] = objectType

    def composeTableAlias(self, tableName):
        """Create an alias for the given table"""
        if tableName != PROPERTIES_TABLE:
            firstRow = self.getTableRow(tableName, 0)
            className = firstRow['class_name']
            if tableName.__contains__('_'):
                tableSplit = tableName.split('_')
                lenTableSplit = len(tableSplit)
                if lenTableSplit > 1:
                    alias = tableName.replace(tableSplit[-1],'') + className
            else:
                alias = className
        else:
            alias = PROPERTIES_TABLE
        return alias

    def getTableNames(self):
        """ Return all the table names found in the database. """
        initTime = time.time()
        if not self._names:
            self._tables = {OBJECT_TABLE: 'classes'}
            self._names = [OBJECT_TABLE]

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
                for key, value in firstRow.items():
                    if key not in EXCLUDED_COLUMNS:
                        labelsTypes.append(_guessType(value))
                self._labelsTypes[tableName] = labelsTypes

            self.composeObjectType()

        if len(self._tables) > 2: # Assuming that there are more tables than just Object and Properties
            self._tableWithAdditionalInfo = OBJECT_TABLE
        endTime = time.time()
        logging.info("GetTableNames: %f" % (endTime - initTime))
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
        fileNameCol = self.findColbyName(colNames, '_filename')

        if indexCol and fileNameCol:
            logger.debug("The columns _index and _filename have been found. "
                         "We will proceed to create a new column with the "
                         "values of these columns.")
            self._extendedColumn = indexCol, fileNameCol
        else:
            indexCol = self.findColbyName(colNames, '_representative._index')
            fileNameCol = self.findColbyName(colNames,
                                                   '_representative._filename')
            if indexCol and fileNameCol:
                logger.debug("The columns _representative._index and "
                             "_representative._filename have been found. "
                             "We will proceed to create a new column with the "
                             "values of these columns.")
                self._extendedColumn = indexCol, fileNameCol

    def generateTableActions(self, table, objectManager):
        """Generate actions for a given table in order to create subsets"""
        if self.getScipionPort() and table.getName() != PROPERTIES_TABLE:
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
            elif alias.startswith('Class'):
                    table.addAction(aliasSplit[1], lambda: self.createSubsetCallback(table, self._objectsType[aliasSplit[1]], objectManager))
            elif alias in self._objectsType and self._objectsType[alias].startswith('SetOf'):
                table.addAction(alias, lambda: self.createSubsetCallback(table, self._objectsType[alias], objectManager))

    def fillTable(self, table, objectManager):
        """Create the table structure (columns) and set the table alias"""
        initTime = time.time()
        tableName = table.getName()
        colNames = self._labels[tableName]
        self.updateExtendColumn(table)

        row = self.getTableRow(tableName, 0, classes=self._tables[tableName])
        if 'id' not in row:
            table.setHasColumnId(False)
        values = [value for key, value in row.items() if key not in EXCLUDED_COLUMNS]

        for index, colName in enumerate(colNames):

            isFileNameCol = self.hasExtendedColumn() and index == self._extendedColumn[1]

            if colName == ENABLED_COLUMN:
                renderer = BoolRenderer()
            elif isFileNameCol:
                renderer = StrRenderer()
            else:
                renderer = table.guessRenderer(str(values[index]))
                if isinstance(renderer, ImageRenderer):
                    imageExt = str(values[index]).split('.')[-1]
                    self.addExternalProgram(renderer, imageExt)

            newCol = Column(colName, renderer)
            newCol.setIsSorteable(True)
            newCol.setIsVisible(objectManager.isLabelVisible(colName))
            table.addColumn(newCol)

            if isFileNameCol:
                logger.debug("Creating an extended column: %s" % EXTENDED_COLUMN_NAME)
                imageExt = str(values[index]).split('.')[-1]
                if values[index] is not None:
                    renderer = ImageRenderer()
                self.addExternalProgram(renderer, imageExt)
                extraCol = Column(colName, renderer)
                extraCol.setIsVisible(newCol.isVisible())
                extraCol.setIsSorteable(False)
                table.addColumn(extraCol)
                newCol.setIsVisible(False)
                newCol.setName(EXTENDED_COLUMN_NAME)

        table.setAlias(self._aliases[tableName])
        self.generateTableActions(table, objectManager)
        endTime = time.time()
        logging.info("Create the table structure: %f" % (endTime - initTime))

    def addExternalProgram(self, renderer: ImageRenderer, imageExt: str):
        self.addChimera(renderer, imageExt)
        self.addImageJ(renderer)
        self.addImageViewer(renderer, imageExt)

    def addChimera(self, renderer: ImageRenderer, imageExt: str):
        chimeraPath = os.environ.get('CHIMERA_HOME', None)
        if chimeraPath is not None:
            if imageExt not in ['st', 'stk', 'mrcs']:
                icon = pw.findResource('chimera.png')

                def openChimeraCallback(imagePath):
                    path = imagePath.split('@')[-1]
                    program = os.path.join(chimeraPath, 'bin', 'ChimeraX')
                    cmd = program + ' "%s"' % path
                    Popen(cmd, shell=True, cwd=os.getcwd())

                renderer.addAction(Action('ChX', icon, 'Open with ChimeraX',
                                          openChimeraCallback))

    def addImageJ(self, renderer: ImageRenderer):
        imageJPath = os.environ.get('IMAGEJ_BINARY_PATH', None)
        if imageJPath is not None:
            icon = pw.findResource('Imagej.png')

            def openImageJCallback(imagePath):
                imageSplit = imagePath.split('@')
                selectedSlice = int(imageSplit[0])
                path = imageSplit[-1]
                program = os.path.join(imageJPath)
                macro = r"""
path = "%s";
open(path);
slice = %d; 
Stack.setSlice(slice);
                """ % (os.path.abspath(path), selectedSlice)

                with tempfile.NamedTemporaryFile(delete=False) as tempFile:
                    tempFile.write(macro.encode('utf-8'))
                    macroPath = tempFile.name

                cmd = program + ' -macro %s' % macroPath
                Popen(cmd, shell=True, cwd=os.getcwd())

            renderer.addAction(Action('IJ', icon, 'Open with ImageJ',
                                      openImageJCallback))

    def addImageViewer(self, renderer: ImageRenderer, imageExt: str):
        icon = pw.findResource('file_vol.png')

        def openImageViewerCallback(imagePath):
            path = imagePath.split('@')[-1]

            pythonPath = sys.executable
            program = '%s -m pwem.viewers.mdviewer.volumeViewer' % pythonPath
            cmd = program + ' "%s"' % path
            Popen(cmd, shell=True, cwd=os.getcwd())

        if imageExt in ['mrc', 'stk']:
            renderer.addAction(Action('IV', icon, 'Open with ImageViewer',
                                      openImageViewerCallback))

    def fillPage(self, page, actualColumn=0, orderAsc=True):
        """
        Read the given table from the sqlite and fill the page(add rows)
        """
        initTime = time.time()
        tableName = page.getTable().getName()
        # moving to the first row of the page
        pageNumber = page.getPageNumber()
        pageSize = page.getPageSize()
        firstRow = pageNumber * pageSize - pageSize
        limit = pageSize

        column = self._labels[tableName][actualColumn]
        mode = 'ASC' if orderAsc else 'DESC'
        self.updateExtendColumn(page.getTable())

        for rowcount, row in enumerate(self.iterTable(tableName, start=firstRow, limit=limit,
                                       classes=self._tables[tableName],
                                       orderBy=column, mode=mode)):
            if row:
                values = [value for key, value in row.items() if key not in EXCLUDED_COLUMNS]
                if 'id' in row.keys():
                    id = row['id']
                else:
                    id = rowcount

                # Checking if exists an extended column
                if self.hasExtendedColumn() and tableName != PROPERTIES_TABLE:
                    if values[self._extendedColumn[1]] is not None:
                        values.insert(self._extendedColumn[1] + 1,
                                      str(values[self._extendedColumn[0]]) + '@' + values[self._extendedColumn[1]])
                    else:
                        values.insert(self._extendedColumn[1] + 1, None)

                page.addRow((int(id), values))
        endTime = time.time()
        logging.info("Fill the page: %f" % (endTime - initTime))

    def getRowsCount(self, tableName):
        """ Return the number of elements in the given table. """
        logger.debug("Reading the table %s" %tableName)
        if tableName == OBJECT_TABLE:
            size = int(self._con.execute(f"SELECT * FROM {PROPERTIES_TABLE} WHERE key='_size'").fetchall()[0]['value'])
            return size
        return self._con.execute(f"SELECT COUNT(*) FROM {tableName}").fetchone()['COUNT(*)']

    def getTableRowCount(self, tableName):
        if tableName not in self._tableCount:
            initTime = time.time()
            self._tableCount[tableName] = self.getRowsCount(tableName)
            endTime = time.time()
            logging.info("Table count: %f" % (endTime - initTime))
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

    def getColumnsValues(self, tableName, columns, xAxis, selection, limit,
                         useSelection, reverse=True):
        """Get the values of the selected columns in order to plot them"""

        logger.debug("Reading the table %s and selected some columns values...")
        cols = columns
        if xAxis and xAxis not in cols:
            cols.append(xAxis)
        columnNames = []
        for column in cols:
            col = self._getColumnMap(tableName, column) or column
            columnNames.append(col)

        columnNames = ", ".join(columnNames)

        col = self._getColumnMap(tableName, xAxis)
        if col is not None:
            xAxis = col

        mode = 'ASC' if reverse else 'DESC'
        orderBy = ' ORDER BY %s %s' % (xAxis, mode) if xAxis else ''
        limit = ' LIMIT %d' % limit if limit is not None else ''
        where = f" WHERE id in ({', '.join(map(str, selection.getSelection().keys()))})" if selection.getCount() > 1 and useSelection else ''

        query = "SELECT %s FROM %s %s %s %s" % (columnNames, tableName, where,
                                                 orderBy, limit)
        selectedColumns = self._con.execute(query).fetchall()

        columnsValues = {}

        firstValue = selectedColumns[0]
        for colName in columns:
            col = self._getColumnMap(tableName, colName) or colName
            columnsValues[colName] = [firstValue[col]]

        for pos, value in enumerate(selectedColumns):
            if pos > 0:
                for colName in columns:
                    col = self._getColumnMap(tableName, colName) or colName
                    columnsValues[colName].append((value[col]))

        return columnsValues

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

        if 'classes' not in kwargs or kwargs['classes'] == PROPERTIES_TABLE:
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
        if tableName in self._columnsMap:
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
            path +="," # Always add a comma, it is expected by the user subset protocol
            if tableName != OBJECT_TABLE:
                path += tableName.split(OBJECT_TABLE)[0]

        self.sendSubsetCreationMessage(path, objectType, subsetName)

    def sendSubsetCreationMessage(self, selectionFile, outputClassName, label):

        import socket

        # Create a client socket
        clientSocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM);

        # Connect to the local server where Scipion GUI is listening
        clientSocket.connect(("127.0.0.1", int(self.getScipionPort())));

        # We should create this message:
        # run protocol ProtUserSubSet inputObject=380 sqliteFile='...','' outputClassName=SetOfTiltSeries other='' label='create subset'
        data = f"run protocol ProtUserSubSet inputObject={self.getScipionObjectId()} " \
               f"sqliteFile='{selectionFile}' outputClassName='{outputClassName}' other='' label='{label}'"

        clientSocket.send(data.encode());


    def getScipionPort(self):
        """ Returns Scipion port or None if not in the environment"""
        return os.getenv(SCIPION_PORT, '1300')

    def getScipionObjectId(self):
        """ Returns Scipion object id"""
        return os.getenv(SCIPION_OBJECT_ID)

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

    @classmethod
    def getCompatibleFileTypes(cls):
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


def extendMDViewer(om:metadataviewer.model.ObjectManager):
    """ Function to extend the object manager with DAOs and readers"""
    om.registerDAO(SqliteFile)
    om.registerReader(ScipionImageReader)