# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software: you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation, either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program.  If not, see <https://www.gnu.org/licenses/>.
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

import pyworkflow.viewer as pwviewer
import pyworkflow.utils as pwutils
import pyworkflow.object as pwobj

# FIXME: rename data by objects
import pwem.objects as emobj
import pwem.emlib.metadata as md
from pwem.viewers import ObjectView


class EmProtocolViewer(pwviewer.ProtocolViewer):
    """ Subclass of ProtocolViewer that have specific functions
    for EM data visualization.
    """
    def objectView(self, filenameOrObject, **kwargs):
        """ This is a wrapper around the ObjectView constructor, just to
        avoid passing the project and protocol, since both are know
        here in the ProtocolViewer.
        Params:
            filenameOrObject: This parameter can be either a filename or an
                object that has 'getFileName' method.
            **kwargs: Can receive extra keyword-arguments that will be passed
                to the ObjectView constructor
        """
        fn = None

        if isinstance(filenameOrObject, str):
            # If the input is a string filename, we should take the object id
            # from the protocol. This assumes that self.protocol have been
            # previously set
            fn = filenameOrObject
            strId = self.getProtocolId()

        elif isinstance(filenameOrObject, pwobj.Object):

            if filenameOrObject.getObjId() is None:
                strId = self.getProtocolId()
            else:
                strId = filenameOrObject.strId()

            if hasattr(filenameOrObject, 'getLocation'):
                # In this case fn will be a location tuple that will be
                # correctly handled by the showj DataView
                fn = filenameOrObject.getLocation()
            elif hasattr(filenameOrObject, 'getFileName'):
                # If the input is an object, we can take the id from it
                fn = filenameOrObject.getFileName()

        if fn is None:
            raise Exception("Incorrect input object, it should be 'string' or "
                            "'Object' (with 'getLocation' or 'getFileName' "
                            "methods).")

        return ObjectView(self._project, strId, fn, **kwargs)

    def createVolumesSqlite(self, files, path, samplingRate,
                            updateItemCallback=None):
        pwutils.cleanPath(path)
        volSet = emobj.SetOfVolumes(filename=path)
        volSet.setSamplingRate(samplingRate)

        for volFn in files:
            vol = emobj.Volume()
            vol.setFileName(volFn)
            if updateItemCallback:
                updateItemCallback(vol)
            volSet.append(vol)
        volSet.write()
        volSet.close()

        return volSet

    def createAngDistributionSqlite(self, sqliteFn, numberOfParticles,
                                    itemDataIterator):
        if not os.path.exists(sqliteFn):
            # List of list of 3 elements containing angleTilt, anglePsi, weight
            projectionList = []

            def getCloseProjection(angleRot, angleTilt):
                """ Get an existing projection close to angleRot, angleTilt.
                Return None if not found close enough.
                """
                for projection in projectionList:
                    if (abs(projection[0] - angleRot) <= 0.5 and
                            abs(projection[1] - angleTilt) <= 0.5):
                        return projection
                return None

            weight = 1. / numberOfParticles

            for angleRot, angleTilt in itemDataIterator:
                projection = getCloseProjection(angleRot, angleTilt)
                if projection is None:
                    projectionList.append([angleRot, angleTilt, weight])
                else:
                    projection[2] = projection[2] + weight

            mdProj = md.MetaData()

            for projection in projectionList:
                mdRow = md.Row()
                mdRow.setValue(md.MDL_ANGLE_ROT, projection[0])
                mdRow.setValue(md.MDL_ANGLE_TILT, projection[1])
                mdRow.setValue(md.MDL_WEIGHT, projection[2])
                mdRow.writeToMd(mdProj, mdProj.addObject())
            mdProj.write(sqliteFn)
