# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from datetime import datetime

import emtable

import pyworkflow.tests as pwtests

import pwem.emlib.metadata as md


class TestMetaData(pwtests.unittest.TestCase):
    
    _labels = [pwtests.WEEKLY]

    def _newMd(self, n=5):
        md0 = md.MetaData()
        xcoor = range(n)
        ycoor = [x*x for x in xcoor]
        for i in range(n):
            self._addRow(md0, '%02d@proj.stk' % i, xcoor[i], ycoor[i])
        return md0

    def _addRow(self, md0, imageFn, xcoor, ycoor):
        objId = md0.addObject()
        md0.setValue(md.MDL_IMAGE, imageFn, objId)
        md0.setValue(md.MDL_XCOOR, xcoor, objId)
        md0.setValue(md.MDL_YCOOR, ycoor, objId)

    def test_removeDuplicates(self):
        md0 = self._newMd()
        md1 = self._newMd()

        # If removing without labels, this metadata should remain the same
        md1.removeDuplicates()
        self.assertEqual(md0, md1)

        # We can use labels for removeDuplicates
        self._addRow(md1, '00@proj.stk', 0, 0)
        md1.removeDuplicates()
        self.assertEqual(md0, md1)

        self._addRow(md1, '06@proj.stk', 0, 0)
        self.assertNotEqual(md0, md1)
        md1.removeDuplicates()
        self.assertNotEqual(md0, md1)
        md1.removeDuplicates(md.MDL_XCOOR)
        self.assertEqual(md0, md1)

        md1.clear()
        self._addRow(md1, '00@proj.stk', 0, 1)
        self._addRow(md1, '00@proj.stk', 0, 2)
        self._addRow(md1, '00@proj.stk', 0, 3)
        md1.removeDuplicates(md.MDL_IMAGE)
        self.assertEqual(md1.size(), 1)
        objId = md1.firstObject()
        self.assertEqual(md1.getValue(md.MDL_YCOOR, objId), 1)

    def test_dropKeepColumns(self):
        md0 = self._newMd()
        md1 = self._newMd()

        self.assertEqual(md0, md1)

        md.dropColumns(md0, md.MDL_XCOOR, md.MDL_YCOOR)
        self.assertEqual(md0.getActiveLabels(), [md.MDL_IMAGE])

        md.keepColumns(md1, "image")
        self.assertEqual(md1.getActiveLabels(), [md.MDL_IMAGE])

        self.assertEqual(md0, md1)

    def iterRowBinding(self, n):
        md0 = self._newMd(n)
        count = 0
        _dt = datetime.now()
        for row in enumerate(md.iterRows(md0)):
            self.assertIsNotNone(row)
            rowValues = row[1]
            image = rowValues.getValue(label='image')
            xcorr = rowValues.getValue(label='xcoor')
            ycorr = rowValues.getValue(label='ycoor')
            self.assertEqual(xcorr*xcorr, ycorr)
            count += 1
        elapsedTime = datetime.now() - _dt
        self.assertTrue(elapsedTime.seconds < 8,
                        msg="Metadata iteration is too slow")
        self.assertEqual(count, n)
        return elapsedTime

    def iterRowsEMTable(self, n):
        cols = ['image', 'xcoor', 'ycoor']
        table = emtable.Table(columns=cols)
        xcoor = range(n)
        ycoor = [x * x for x in xcoor]
        _dict = dict
        for i in range(n):
            _dict = {'image': '%02d@proj.stk' % i,
                     'xcoor': xcoor[i],
                     'ycoor': ycoor[i]}
            table.addRow(*_dict.values())

        starFile = '/tmp/test.star'
        with open(starFile, 'w') as f:
            f.write("# Star file generated with Scipion\n")
            f.write("# version 30001\n")
            table.writeStar(f, tableName='test')

        count = 0
        _dt = datetime.now()

        for row in enumerate(table.iterRows(fileName=starFile, tableName='test')):
            self.assertIsNotNone(row)
            rowValues = row[1]
            image = rowValues.get('image')
            xcorr = rowValues.get('xcoor')
            ycorr = rowValues.get('ycoor')
            self.assertEqual(xcorr * xcorr, ycorr)
            count += 1

        elapsedTime = datetime.now() - _dt
        self.assertTrue(elapsedTime.seconds < 1,
                        msg="Metadata iteration is too slow")
        self.assertEqual(count, n)
        return elapsedTime

    def test_iterRowMetadata(self):
        """
        Check if iterating over N rows using the binding is faster than emtable
        """
        n = 10000  # number of .star file rows
        iterRowBindingElapsedTime = self.iterRowBinding(n)
        iterRowsEMTableElapsedTime = self.iterRowsEMTable(n)

        print("Binding iteration : %s" % iterRowBindingElapsedTime)
        print("Emtable iteration : %s" % iterRowsEMTableElapsedTime)

        self.assertTrue(iterRowsEMTableElapsedTime < iterRowBindingElapsedTime,
                        msg="Iterating over %d rows using the binding is "
                            "faster than using emtable" % n)

