# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
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

import os

import pyworkflow.viewer as pwviewer
import pyworkflow.utils as pwutils

from pwem import Domain
from pwem import emlib
import pwem.objects as emobj
import pwem.protocols as emprot

from .views import (ObjectView, MicrographsView, CoordinatesObjectView,
                    ClassesView, Classes3DView, CtfView, DataView)
from .showj import (RENDER, SAMPLINGRATE, ORDER, VISIBLE, MODE, MODE_MD,
                    SORT_BY, getJvmMaxMemory, launchTiltPairPickerGUI, LABELS)
from ..convert.headers import Ccp4Header


class RegistryViewerConfig:
    config = {}

    @classmethod
    def getConfig(cls, type):
        return cls.config.get(type, None)

    @classmethod
    def registerConfig(cls, type, config=None):
        """
        Allow registration of other objects
        :param type: Class to register --> SetOfXXX
        :param config: Optional, dictionary with the fields configuration, otherwise all fields will be shown.

        :return: Nothing

        """
        if config is None:
            config = {MODE: MODE_MD}

        cls.config[type] = config

# Registering viewer config
RegistryViewerConfig.registerConfig(emobj.SetOfPDBs, {ORDER: 'id _filename ',
                                   VISIBLE: 'id _filename ',
                                   MODE: MODE_MD,
                                   RENDER: "no"})


#SetOfParticles
labels = ('id enabled _index _filename _xmipp_zScore _xmipp_cumulativeSSNR '
                      '_sampling _xmipp_scoreByVariance _xmipp_scoreEmptiness '
                      '_ctfModel._defocusU _ctfModel._defocusV _ctfModel._defocusAngle '
                      '_transform._matrix')
RegistryViewerConfig.registerConfig(emobj.SetOfParticles,
                                   {ORDER: labels,
                                   VISIBLE: labels,
                                    SORT_BY: '_xmipp_zScore asc',
                                    RENDER: '_filename'})


# SetOfVolumes
labels = 'id enabled comment _filename '
RegistryViewerConfig.registerConfig(emobj.SetOfVolumes,
                                   {ORDER: labels,
                                   VISIBLE: labels,
                                    SORT_BY: '_xmipp_zScore asc',
                                    RENDER: '_filename'})

# SetOfMovies
labels = 'id _filename _samplingRate _acquisition._dosePerFrame _acquisition._doseInitial '
RegistryViewerConfig.registerConfig(emobj.SetOfMovies,
                                   {ORDER: labels,
                                   VISIBLE: labels,
                                    MODE: MODE_MD,
                                    RENDER: 'no'})

#MicrographsTiltPair
labels = 'id enabled _untilted._filename _tilted._filename'
renderLabels = '_untilted._filename _tilted._filename'

RegistryViewerConfig.registerConfig(emobj.MicrographsTiltPair,
                                   {ORDER: labels,
                                   VISIBLE: labels,
                                    MODE: MODE_MD,
                                    RENDER: renderLabels})


#ParticlesTiltPair
labels = 'id enabled _untilted._filename _tilted._filename'
renderLabels = '_untilted._filename _tilted._filename'
RegistryViewerConfig.registerConfig(emobj.ParticlesTiltPair,
                                   {ORDER: labels,
                                   VISIBLE: labels,
                                    MODE: MODE_MD,
                                    RENDER: renderLabels})

# SetOfClasses2D
labels = 'enabled id _size _representative._filename _filename'
RegistryViewerConfig.registerConfig(emobj.SetOfClasses2D,
                                   {ORDER: labels,
                                   VISIBLE: labels,
                                    RENDER: '_representative._filename',
                                    SORT_BY: '_size desc',
                                    LABELS: 'id _size'})



class DataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _name = "Xmipp"
    _environments = [pwviewer.DESKTOP_TKINTER, pwviewer.WEB_DJANGO]
    _targets = [
        emobj.Image,
        emobj.SetOfClasses2D,
        emobj.SetOfClasses3D,
        emobj.SetOfCoordinates,
        emobj.SetOfCTF,
        emobj.SetOfParticles,
        emobj.SetOfMovies,
        emobj.SetOfNormalModes,
        emobj.SetOfPrincipalComponents,
        emobj.SetOfPDBs,
        emobj.SetOfAtomStructs,
        emprot.ProtParticlePicking,
        emprot.ProtImportMovies,
        # TiltPairs related data
        emobj.CoordinatesTiltPair,
        emobj.MicrographsTiltPair,
        emobj.ParticlesTiltPair,
        emobj.SetOfImages,
    ]

    @classmethod
    def registerConfig(cls, type, config=None):
        """
        Allow registration of other objects
        :param type: Class to register --> SetOfXXX
        :param config: Optional, dictionary with the fields configuration, otherwise all fields will be shown.

        :return: Nothing

        """
        cls._targets.append(type)
        RegistryViewerConfig.registerConfig(type, config)

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _addObjView(self, obj, fn, viewParams={}):
        objView = ObjectView(self._project, obj.strId(), fn,
                             viewParams=viewParams)
        self._views.append(objView)
        return objView

    def _visualize(self, obj, **kwargs):
        cls = type(obj)

        if issubclass(cls, emobj.SetOfClasses2D):
            return [ClassesView(self._project, obj.strId(), obj.getFileName(), **kwargs)]


        # Try registry first
        config = RegistryViewerConfig.getConfig(cls)

        if config is not None:
            fn = obj.getFileName()
            objView = self._addObjView(obj, fn, config)
            if issubclass(cls, emobj.SetOfMovies):
                # For movies increase the JVM memory by 1 GB, just in case
                objView.setMemory(getJvmMaxMemory() + 1)
            return self._views

        if issubclass(cls, emobj.Volume):
            fn = emlib.image.ImageHandler.locationToXmipp(obj)

            if fn.endswith(':mrc'):
                # fn may came in the form of 000001@Runs/..../vol.mrc". We need to pass an actual file to Ccp4Header
                actualFile = fn.split("@")[-1]
                ccp4header = Ccp4Header(actualFile, readHeader=True)
                if ccp4header.getISPG() == 1:
                    fn = fn.replace(':mrc', '')
                else:
                    print('Developers warning: %s headers do not match the standard for volumes. '
                          'You can use fixVolume(pwem.convert.headers.py) '
                          'method to fix the headers' % fn)
            self._addObjView(obj, fn,
                             {RENDER: 'image',
                              SAMPLINGRATE: obj.getSamplingRate()})

        elif issubclass(cls, emobj.Image):
            fn = emlib.image.ImageHandler.locationToXmipp(obj)
            self._addObjView(obj, fn)

        # elif issubclass(cls, emobj.SetOfPDBs):
        #     fn = obj.getFileName()
        #     labels = 'id _filename '
        #     self._addObjView(obj, fn, {ORDER: labels,
        #                                VISIBLE: labels,
        #                                MODE: MODE_MD,
        #                                RENDER: "no"})

        # elif issubclass(cls, emobj.SetOfMovies):
        #     fn = obj.getFileName()
        #     # Enabled for the future has to be available
        #     labels = ('id _filename _samplingRate _acquisition._dosePerFrame '
        #               '_acquisition._doseInitial ')
        #     moviesView = self._addObjView(obj, fn, {ORDER: labels,
        #                                             VISIBLE: labels,
        #                                             MODE: MODE_MD,
        #                                             RENDER: "no"})
        #     # For movies increase the JVM memory by 1 GB, just in case
        #     moviesView.setMemory(getJvmMaxMemory() + 1)

        elif issubclass(cls, emobj.SetOfMicrographs):
            self._views.append(MicrographsView(self._project, obj, **kwargs))

        # elif issubclass(cls, emobj.MicrographsTiltPair):
        #     labels = 'id enabled _untilted._filename _tilted._filename'
        #     renderLabels = '_untilted._filename _tilted._filename'
        #     self._addObjView(obj, obj.getFileName(), {ORDER: labels,
        #                                               VISIBLE: labels,
        #                                               MODE: MODE_MD,
        #                                               RENDER: renderLabels})
        #
        # elif issubclass(cls, emobj.ParticlesTiltPair):
        #     labels = 'id enabled _untilted._filename _tilted._filename'
        #     renderLabels = '_untilted._filename _tilted._filename'
        #     self._addObjView(obj, obj.getFileName(), {ORDER: labels,
        #                                               VISIBLE: labels,
        #                                               RENDER: renderLabels,
        #                                               MODE: MODE_MD})

        elif issubclass(cls, emobj.SetOfCoordinates):
            # FIXME: Remove dependency on xmipp3 plugin to visualize coordinates
            xmipp3 = Domain.importFromPlugin('xmipp3',
                                             errorMsg="xmipp3 plugin is required "
                                                      "now to visualize coordinates.")
            micSet = obj.getMicrographs()  # accessing mics to provide metadata file
            if micSet is None:
                raise Exception('visualize: SetOfCoordinates has no micrographs set.')

            mdFn = getattr(micSet, '_xmippMd', None)
            if mdFn:
                fn = mdFn.get()
            else:  # happens if protocol is not an xmipp one
                fn = self._getTmpPath(micSet.getName() + '_micrographs.xmd')
                xmipp3.convert.writeSetOfMicrographs(micSet, fn)
            tmpDir = self._getTmpPath(obj.getName())
            pwutils.cleanPath(tmpDir)
            pwutils.makePath(tmpDir)
            # FIXME: (JMRT) We are always writing the SetOfCoordinates and removing
            # the tmpDir, we need to take into account if the user have pick
            # some particles in the tmpDir and have not save them, that now
            # will loose all picked particles.
            # A possible solution could be to alert that changes have not been
            # written during modification of tmpDir or create a new Xmipp picking
            # protocol to continue picking later without loosing the coordinates.
            xmipp3.convert.writeSetOfCoordinates(tmpDir, obj)
            self._views.append(CoordinatesObjectView(self._project, fn,
                                                     tmpDir, self.protocol,
                                                     inTmpFolder=True))

        # elif issubclass(cls, emobj.SetOfParticles):
        #     fn = obj.getFileName()
        #     labels = ('id enabled _index _filename _xmipp_zScore _xmipp_cumulativeSSNR '
        #               '_sampling _xmipp_scoreByVariance _xmipp_scoreEmptiness '
        #               '_ctfModel._defocusU _ctfModel._defocusV _ctfModel._defocusAngle '
        #               '_transform._matrix')
        #     self._addObjView(obj, fn, {ORDER: labels,
        #                                VISIBLE: labels,
        #                                SORT_BY: '_xmipp_zScore asc',
        #                                RENDER: '_filename'})

        # elif issubclass(cls, emobj.SetOfVolumes):
        #     fn = obj.getFileName()
        #     labels = 'id enabled comment _filename '
        #     self._addObjView(obj, fn, {MODE: MODE_MD,
        #                                ORDER: labels,
        #                                VISIBLE: labels,
        #                                RENDER: '_filename'})

        # elif issubclass(cls, emobj.SetOfClasses2D):
        #     self._views.append(ClassesView(self._project, obj.strId(),
        #                                    obj.getFileName(), **kwargs))
        #
        elif issubclass(cls, emobj.SetOfClasses3D):
            self._views.append(Classes3DView(self._project, obj.strId(),
                                             obj.getFileName()))

        elif issubclass(cls, emobj.SetOfImages):
            self._views.append(ObjectView(self._project, obj.strId(),
                                          obj.getFileName(), **kwargs))

        elif issubclass(cls, emobj.SetOfCTF):
            self._views.append(CtfView(self._project, obj))

        elif issubclass(cls, emobj.CoordinatesTiltPair):
            # FIXME: Remove dependency on xmipp3 plugin to visualize coordinates
            xmipp3 = Domain.importFromPlugin('xmipp3',
                                             errorMsg="xmipp3 plugin is required "
                                                      "now to visualize coordinates.")
            tmpDir = self._getTmpPath(obj.getName())
            pwutils.makePath(tmpDir)

            mdFn = os.path.join(tmpDir, 'input_micrographs.xmd')
            xmipp3.convert.writeSetOfMicrographsPairs(
                obj.getUntilted().getMicrographs(),
                obj.getTilted().getMicrographs(), mdFn)

            # TODO: Review this if ever a non Xmipp CoordinatesTiltPair is available
            xmipp3.convert.writeSetOfCoordinates(tmpDir, obj.getUntilted())
            xmipp3.convert.writeSetOfCoordinates(tmpDir, obj.getTilted())
            launchTiltPairPickerGUI(mdFn, tmpDir, self.protocol)

        elif issubclass(cls, emprot.ProtParticlePicking):
            if obj.getOutputsSize() >= 1:
                self._visualize(obj.getCoords())

        elif issubclass(cls, emprot.ProtImportMovies):
            movs = obj.outputMovies
            self._visualize(movs)
            gainFn = movs.getGain()
            if gainFn is not None:
                if os.path.exists(gainFn):
                    self._views.append(DataView(gainFn))

        elif issubclass(cls, emobj.SetOfNormalModes):
            self._views.append(DataView(obj.getFileName()))

        elif issubclass(cls, emobj.EMSet):
            self._views.append(DataView(obj.getFileName()))

        return self._views


class BasicMDViewer(DataViewer):
    """ Viewer to show sets in the metadataviewer to allow manual subsets when
    DataViewer shows it in a different way: e.g: Coordinates2D."""
    _name = "Basic MD viewer"
    _targets = [emobj.SetOfCoordinates, emobj.SetOfFSCs]

    def _visualize(self, obj, **kwargs):
        fn = obj.getFileName()
        self._addObjView(obj, fn, {MODE: MODE_MD})
        return self._views
