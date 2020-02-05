# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from collections import OrderedDict

from pyworkflow.object import ObjectWrap

from ..lib import *


LABEL_TYPES = {
    LABEL_SIZET: float,
    LABEL_DOUBLE: float,
    LABEL_INT: int,
    LABEL_BOOL: bool
}


class Row:
    """ Support Xmipp class to store label and value pairs
    corresponding to a Metadata row.
    """

    def __init__(self):
        self._labelDict = OrderedDict()  # Dictionary containing labels and values
        self._objId = None  # Set this id when reading from a metadata

        # Add set and get method as alias to setValue and getValue
        # in this way the Row will behave more like a dict
        self.set = self.setValue
        self.get = self.getValue

    def getObjId(self):
        return self._objId

    def hasLabel(self, label):
        return self.containsLabel(label)

    def containsLabel(self, label):
        # Allow getValue using the label string
        if isinstance(label, str):
            label = str2Label(label)
        return label in self._labelDict

    def removeLabel(self, label):
        if self.hasLabel(label):
            del self._labelDict[label]

    def setValue(self, label, value):
        """args: this list should contains tuples with
        MetaData Label and the desired value"""
        # Allow setValue using the label string
        if isinstance(label, str):
            label = str2Label(label)
        self._labelDict[label] = value

    def getValue(self, label, default=None):
        """ Return the value of the row for a given label. """
        # Allow getValue using the label string
        if isinstance(label, str):
            label = str2Label(label)
        return self._labelDict.get(label, default)

    def getValueAsObject(self, label, default=None):
        """ Same as getValue, but making an Object wrapping. """
        return ObjectWrap(self.getValue(label, default))

    def readFromMd(self, md, objId):
        """ Get all row values from a given id of a metadata. """
        self._labelDict.clear()
        self._objId = objId

        for label in md.getActiveLabels():
            self._labelDict[label] = md.getValue(label, objId)

    def addToMd(self, md):
        self.writeToMd(md, md.addObject())

    def writeToMd(self, md, objId):
        """ Set back row values to a metadata row. """
        for label, value in self._labelDict.items():
            # TODO: Check how to handle correctly unicode type
            # in Xmipp and Scipion
            if type(value) is str:
                value = str(value)
            try:
                md.setValue(label, value, objId)
            except Exception as ex:
                import sys
                print("XmippMdRow.writeToMd: Error writing value to metadata.", file=sys.stderr)
                print("                     label: %s, value: %s, type(value): %s" % (
                    label2Str(label), value, type(value)), file=sys.stderr)
                raise ex

    def readFromFile(self, fn):
        md = MetaData(fn)
        self.readFromMd(md, md.firstObject())

    def writeToFile(self, fn):
        md = MetaData()
        self.writeToMd(md, md.addObject())
        md.write(fn)

    def copyFromRow(self, other):
        for label, value in other._labelDict.items():
            self.setValue(label, value)

    def clone(self):
        """ Return another Row that have exactly the same
        values as self.
        """
        row = Row()
        row.copyFromRow(self)
        row._objId = self._objId

        return row

    def __str__(self):
        s = '{'
        for k, v in self._labelDict.items():
            s += '  %s = %s\n' % (label2Str(k), v)
        return s + '}'

    def __iter__(self):
        return self._labelDict.items()

    def __getitem__(self, item):
        return self.getValue(item)

    def __setitem__(self, key, value):
        return self.setValue(key, value)

    def containsAll(self, labels):
        """ Check if all labels are present in the row.
        Params:
            row: the Row object.
            labels: either a dict or list object containing the labels
                (in the case of dicts, label are the dict.values())
        """
        values = labels.values() if isinstance(labels, dict) else labels
        return all(self.containsLabel(l) for l in values)

    def containsAny(self, labels):
        """ Check if at least one of labels is present in the row.
        Params:
            row: the Row object.
            labels: either a dict or list object containing the labels
                (in the case of dicts, label are the dict.values())
        """
        values = labels.values() if isinstance(labels, dict) else labels
        return any(self.containsLabel(l) for l in values)

    def printDict(self):
        """ Fancy printing of the row, mainly for debugging. """
        print(str(self))


class RowMetaData:
    """ This class is a wrapper for MetaData in row mode.
    Where only one object is used.
    """

    def __init__(self, filename=None):
        self._md = MetaData()
        self._md.setColumnFormat(False)
        self._id = self._md.addObject()

        if filename:
            self.read(filename)

    def setValue(self, label, value):
        self._md.setValue(label, value, self._id)

    def getValue(self, label):
        return self._md.getValue(label, self._id)

    def write(self, filename, mode=MD_APPEND):
        self._md.write(filename, mode)

    def read(self, filename):
        self._md.read(filename)
        self._md.setColumnFormat(False)
        self._id = self._md.firstObject()

    def containsLabel(self, label):
        return self._md.containsLabel(label)

    def __str__(self):
        return str(self._md)


def label2Python(label):
    """ Return the Python type (int, float, bool) for a given 
    metadata label (LABEL_INT, LABEL_DOUBLE..etc)
    """
    return LABEL_TYPES.get(labelType(label), str)


def getLabel(value):
    """ Return the label value either from an int value or an string. """
    if isinstance(value, int):
        return value
    elif isinstance(value, str):
        return str2Label(value)
    else:
        raise Exception("Invalid value type (%s) for label. " % type(value))


def getFirstRow(mdOrFn):
    """ Return the first object of a metadata.
    Params:
        mdOrFn: you can pass a metadata or a filename as argument.
    """

    if isinstance(mdOrFn, str):
        md = MetaData()
        md.read(mdOrFn, 1)
    else:  # mdOrFn is MetaData
        md = mdOrFn
        
    if md.getParsedLines():
        row = Row()
        row.readFromMd(md, md.firstObject())
    else:
        row = None
    
    return row


def getSize(filename):
    """ Return the metadata size without parsing entirely. """
    md = MetaData()
    md.read(filename, 1)
    return md.getParsedLines()


def isEmpty(filename):
    """ Use getMdSize to check if metadata is empty. """
    return getSize(filename) == 0


def iterRows(md, sortByLabel=None):
    """ Iterate over the rows of the given metadata.
    Params:
        md: a MetaData object or a filename (MetaData will be read)
        sortByLabel: a label to sort the metadata before iterate.
    """
    # If md is string, take as filename and create the metadata

    if isinstance(md, str):
        md = MetaData(md)

    if sortByLabel is not None:
        md.sort(sortByLabel)

    row = Row()
    
    for objId in md:
        row.readFromMd(md, objId)
        yield row


def dropColumns(mdObj, *labels):
    """ Drop all columns from a given metadata.
    Labels can be either string or int.
    """
    for l in labels:
        mdObj.removeLabel(getLabel(l))


def keepColumns(mdObj, *labels):
    """ Drop all columns from mdObj that are not in labels.
    Labels can be either string or int.
    """
    # Handle string or int labels input
    keepLabels = {getLabel(l) for l in labels}

    for l in mdObj.getActiveLabels():
        if l not in keepLabels:
            mdObj.removeLabel(l)


def joinBlocks(inputMd, blockPrefix=None):
    mdImages = MetaData()
    mdAll = MetaData()
    mdBlocks = getBlocksInMetaDataFile(inputMd)
    
    for mdBlock in mdBlocks:
        if blockPrefix is not None:
            if mdBlock.startswith(blockPrefix):
                mdImages.read(mdBlock + "@" + inputMd)
                mdAll.unionAll(mdImages)
        else:
            mdImages.read(mdBlock + "@" + inputMd)
            mdAll.unionAll(mdImages)
    return mdAll


class SetMdIterator:
    """ Class to iterate over an input set and skip
    elements not present in metadata.
    This class can be used in copyItems when the number
    of elements in the set is higher that in metadata.
    """
    def __init__(self, md, sortByLabel=None, 
                 keyLabel=MDL_ITEM_ID,
                 updateItemCallback=None,
                 skipDisabled=False):
        
        if updateItemCallback is None:
            raise Exception('Set an updateItemCallback')
        
        self.iterMd = iterRows(md, sortByLabel) 
        self.keyLabel = keyLabel
        self.updateItemCallback = updateItemCallback
        self.skipDisabled = skipDisabled
        self.__nextRow()
        
    def __nextRow(self):
        try:
            self.lastRow = next(self.iterMd)
        except StopIteration:
            self.lastRow = None
            
    def updateItem(self, item, row):
        """ This function should be passed to copyItems
        as callback and it will filter the items
        not present in the metadata.
        """
        row = self.lastRow
        
        if row is not None:
            if row.hasLabel(MDL_ENABLED):
                enabled = row.getValue(MDL_ENABLED)
            else:
                enabled = 1

        if (row is None or
                item.getObjId() != row.getValue(self.keyLabel)):
            item._appendItem = False
            
        elif enabled == -1 and self.skipDisabled:
            item._appendItem = False
            self.__nextRow()
            
        else:
            item._appendItem = True
            self.updateItemCallback(item, row)
            self.__nextRow()
