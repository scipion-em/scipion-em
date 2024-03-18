# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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


import pyworkflow as pw
from pyworkflow.mapper.sqlite_db import SqliteDb
from pyworkflow.utils.properties import Message

import pwem.objects.data_flexhub as pwobj


# FIXME: Do not use this methods and remove in the future
class ProtFlexBase:
    def _createSet(self, SetClass, template, suffix, **kwargs):
        """ Create a set and set the filename using the suffix.
        If the file exists, it will be deleted. """
        setFn = self._getPath(template % suffix)
        # Close the connection to the database if
        # it is open before deleting the file
        pw.utils.cleanPath(setFn)

        SqliteDb.closeConnection(setFn)
        setObj = SetClass(filename=setFn, **kwargs)
        return setObj

    def _createSetOfParticlesFlex(self, suffix='', **kwargs):
        return self._createSet(pwobj.SetOfParticlesFlex,
                               'flexparticles%s.sqlite', suffix, **kwargs)

    def _createSetOfClassesFlex(self, flexParticles, suffix='', **kwargs):
        classes = self._createSet(pwobj.SetOfClassesFlex,
                                  'flexClasses%s.sqlite', suffix, **kwargs)
        classes.setImages(flexParticles)
        return classes

    def _createSetOfClassesStructFlex(self, flexParticles, suffix='', **kwargs):
        classes = self._createSet(pwobj.SetOfClassesStructFlex,
                                  'flexClassesStruct%s.sqlite', suffix, **kwargs)
        classes.setImages(flexParticles)
        return classes

    def _createSetOfVolumesFlex(self, suffix='', **kwargs):
        return self._createSet(pwobj.SetOfVolumesFlex,
                               'flexvolumes%s.sqlite', suffix, **kwargs)

    def _createSetOfAtomStructFlex(self, suffix='', **kwargs):
        return self._createSet(pwobj.SetOfAtomStructFlex,
                               'flexstructs%s.sqlite', suffix, **kwargs)

    def _getOutputSuffix(self, cls):
        """ Get the name to be used for a new output.
        For example: output3DCoordinates7.
        It should take into account previous outputs
        and number with a higher value.
        """
        maxCounter = -1
        for attrName, _ in self.iterOutputAttributes(cls):
            suffix = attrName.replace(self.OUTPUT_PREFIX, '')
            try:
                counter = int(suffix)
            except:
                counter = 1  # when there is not number assume 1
            maxCounter = max(counter, maxCounter)

        return str(maxCounter+1) if maxCounter > 0 else ''  # empty if not output
