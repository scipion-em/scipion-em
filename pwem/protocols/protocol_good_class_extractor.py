# **************************************************************************
# *
# * Authors:     Daniel Marchan [1]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import os
import time
import sys
import matplotlib.pyplot as plt
from pyworkflow.utils import prettyTime
import pyworkflow.protocol.params as params
from pyworkflow.object import Set
from pyworkflow.protocol import ProtStreamingBase, STEPS_PARALLEL
from pwem.protocols import EMProtocol
from pwem.objects import SetOfParticles
#import logging
#logger = logging.getLogger(__file__)


OUTPUT_PARTICLES = "outputParticles"
OUTPUT_DISCARDED_PARTICLES = "outputParticlesDiscarded"
LAST_DONE_FILE = "last_done.txt"


class ProtGoodClassesExtractor(EMProtocol, ProtStreamingBase):
    """ Extracts items from a SetOfClasses based on a list of IDs or a set of given good averages/classes
    """

    _label = "good classes extractor"
    outputsToDefine = {}
    _possibleOutputs = {OUTPUT_PARTICLES: SetOfParticles,
                        OUTPUT_DISCARDED_PARTICLES: SetOfParticles}
    # Mode
    LIST_CLASSES = 0
    LIST_IDS = 1

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputClasses', params.PointerParam,
                      pointerClass='SetOfClasses',
                      label='Input classes',
                      help='Set of classes to extract items from.')
        form.addParam('mode', params.EnumParam, choices=['list_classes', 'list_ids'],
                      label="Select the source from which to extract the good references", default=self.LIST_CLASSES,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='This option allows for either get the good classes from a set of classes '
                           'or from a list of ids.')
        form.addParam('inputGoodClasses', params.PointerParam,
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      label='Good references',
                      condition="mode==%d" % self.LIST_CLASSES,
                      help='Set of good reference to extract particles from the inputClasses.')
        form.addParam('inputGoodListIds', params.StringParam,
                      label='Good references',
                      condition="mode==%d" % self.LIST_IDS,
                      help='List of good reference IDs, separated by commas, '
                           'to extract particles from the inputClasses.')

        form.addParallelSection(threads=3, mpi=1)

    # -------------------------- INSERT steps functions ---------------------------
    def stepsGeneratorStep(self) -> None:
        """
        This step should be implemented by any streaming protocol.
        It should check its input and when ready conditions are met
        call the self._insertFunctionStep method.
        """
        self.newDeps = []
        self.initialStep()

        while not self.finish:
            if not self._newParticlesToProcess():
                 self.info('No new particles')
            else:
                with self._lock:
                    classSet = self._loadInputClassesSet()

                self.isStreamClosed = classSet.getStreamState()

                if self.selectGood:  # Only happens once
                    selectStep = self._insertFunctionStep(self.selectGoodClasses,
                                                          prerequisites=[])

                extractStep = self._insertFunctionStep(self.extractElements, classSet,
                                                       prerequisites=selectStep)

                self.newDeps.append(extractStep)

            if self.isStreamClosed == Set.STREAM_CLOSED:
                self.info('Stream closed')
                # Finish everything and close output sets
                self._insertFunctionStep(self.closeOutputStep,
                                         prerequisites=self.newDeps)
                self.finish = True

            sys.stdout.flush()
            time.sleep(60)

    # --------------------------- STEPS functions -------------------------------
    def initialStep(self):
        self.finish = False
        self.selectGood = True
        self.isStreamClosed = False
        self.goodParticles = []
        self.badParticles = []
        self.particlesDistribution = {'good': [], 'bad': []}
        self.goodClassesIDs = []
        self.dictsTimes = {}

    def extractElements(self, inputClasses):
        """
        Method to extract the particles from the selected classes, this method generates two output sets:
            - accepted particles
            - discarded particles
        """
        output = self._loadOutputSet(OUTPUT_PARTICLES, "")
        outputDiscarded = self._loadOutputSet(OUTPUT_DISCARDED_PARTICLES, "discarded")

        with self._lock:
            # For each class (order by number of items)
            for clazz in inputClasses.iterItems(orderBy="_size", direction="DESC"):
                # Make the query to load only the new particles
                where = None
                if str(clazz.getObjId()) in self.dictsTimes:
                    lastTime = str(self.dictsTimes[str(clazz.getObjId())])
                    where = 'creation>"' + lastTime + '"'
                    self.debug('Last creation time in class %d: %s'
                               % (clazz.getObjId(), lastTime))

                # Two sets of particles:
                if clazz.getObjId() in self.goodClassesIDs:  # Accepted particles
                    for image in clazz.iterItems(orderBy='creation', direction='ASC', where=where):
                        tmp_accepted = image.getObjCreation()
                        newImage = image.clone()
                        output.append(newImage)
                        self.goodParticles.append(image.getObjId())
                    self.dictsTimes[str(clazz.getObjId())] = tmp_accepted  # Store the latest time
                else:  # Discarded particles
                    for image in clazz.iterItems(orderBy='creation', direction='ASC', where=where):
                        tmp_discarded = image.getObjCreation()
                        newImageDiscarded = image.clone()
                        outputDiscarded.append(newImageDiscarded)
                        self.badParticles.append(image.getObjId())
                    self.dictsTimes[str(clazz.getObjId())] = tmp_discarded  # Store the latest time

        self.info('Size output %d and size discarded output %d' % (len(output), len(outputDiscarded)))
        self.debug(str(self.dictsTimes))

        if len(output) > 0:
            self._updateOutputSet(OUTPUT_PARTICLES, output, self.isStreamClosed)
        if len(outputDiscarded) > 0:
            self._updateOutputSet(OUTPUT_DISCARDED_PARTICLES, outputDiscarded, self.isStreamClosed)

        self._writeLastDone(self.dictsTimes)
        self._createPlots()

    def selectGoodClasses(self):
        """
        Select only the good Classes from the Averages or from the IDs list
        """
        if self.mode == self.LIST_CLASSES:
            self.goodClassesIDs = self.inputGoodClasses.get().getIdSet()
        else:
            self.goodClassesIDs = self._getGoodIds()

        self.info('Good classes IDs:')
        self.info(self.goodClassesIDs)
        self.selectGood = False
        # Change to list of filenames? For the moment it matched so not needed
        # inputGoodClasses.get().getUniqueValues('filename'), clazz.getRepresentative()._filename, clazz.getObjName()

    def closeOutputStep(self):
        self.info("Size of good particles output: %d" % len(self.goodParticles))
        self.info("Size of bad particles rejected: %d" % len(self.badParticles))
        self._closeOutputSet()

# --------------------------- UTILS functions -----------------------------
    def _loadOutputSet(self, outputName, suffix):
        """
        Load the output set if it exists or create a new one.
        """
        outputSet = getattr(self, outputName, None)
        if outputSet is None:
            outputSet = self._createSetOfParticles(suffix)
            images = self.inputClasses.get().getImages()
            outputSet.copyInfo(images)
            outputSet.setStreamState(Set.STREAM_OPEN)
        else:
            outputSet.enableAppend()

        return outputSet

    def _newParticlesToProcess(self):
        classesFile = self.inputClasses.get().getFileName()
        now = datetime.now()
        self.lastCheck = getattr(self, 'lastCheck', now)
        mTime = datetime.fromtimestamp(os.path.getmtime(classesFile))
        self.debug('Last check: %s, modification: %s'
                   % (prettyTime(self.lastCheck),
                      prettyTime(mTime)))
        # If the input have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and self.dictsTimes:
            newParticlesBool = False
        else:
            newParticlesBool = True

        self.lastCheck = now
        return newParticlesBool

    def _loadInputClassesSet(self):
        """ Returns te input set of particles"""
        classSet = self.inputClasses.get()
        classSet.loadAllProperties()

        return classSet

    def _getGoodIds(self):
        ids = self.inputGoodListIds.get().split(',')
        listIDs = [int(id) for id in ids]
        return listIDs

    def _writeLastDone(self, creationTimeDict):
        """ Write to a text file the last item creation time done. """
        dictStr = str(creationTimeDict)

        with open(self._getExtraPath(LAST_DONE_FILE), 'w') as f:
            f.write(dictStr)

    def _getLastDone(self):
        # Open the file in read mode and read the number
        with open(self._getExtraPath(LAST_DONE_FILE), "r") as file:
            content = file.read()
        dictTimes = eval(content)

        return dictTimes

    def _createPlots(self):
        balancePlot(len(self.goodParticles), len(self.badParticles),
                    self._getExtraPath('particle_distribution.png'))
        self.particlesDistribution['good'].append(len(self.goodParticles))
        self.particlesDistribution['bad'].append(len(self.badParticles))
        if len(self.particlesDistribution['good']) >= 3:
            balanceOverTimePlot(self.particlesDistribution['good'], self.particlesDistribution['bad'],
                                self._getExtraPath('cumulative_distribution.png'))

    def getDistributionPlot(self):
        return self._getExtraPath('particle_distribution.png')

    def getDistributionTimePlot(self):
        return self._getExtraPath('cumulative_distribution.png')

# --------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        errors = []
        return errors


# --------------------------- EXTRA functions --------------------------------
def balancePlot(good_particles, bad_particles, fileName):
    # Labels for the classes
    classes = ['Good', 'Bad']
    # Values for the classes
    particle_counts = [good_particles, bad_particles]
    # Define colors for the bars
    colors = ['#007ACC', '#FF585D']
    # Create a bar plot with custom colors and formal style
    fig, ax = plt.subplots(figsize=(8, 6))  # Adjust the figure size
    ax.bar(classes, particle_counts, color=colors, edgecolor='black', linewidth=1.2)
    # Customize axis labels and title
    ax.set_xlabel('Classes', fontsize=14)
    ax.set_ylabel('Number of Particles', fontsize=14)
    ax.set_title('Particle Distribution in Classes', fontsize=16)
    # Add grid lines
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    # Customize font size for tick labels
    ax.tick_params(axis='both', which='major', labelsize=12)
    # Save the figure as an image (e.g., PNG)
    plt.savefig(fileName, dpi=300, bbox_inches='tight')

def balanceOverTimePlot(cumulative_good, cumulative_bad, fileName):
    time_points = range(1, len(cumulative_good)+1)
    # Create a plot for the cumulative distributions over time
    plt.figure(figsize=(10, 6))  # Adjust the figure size
    plt.plot(time_points, cumulative_good, marker='o', label='Good Classes', color='blue')
    plt.plot(time_points, cumulative_bad, marker='o', label='Bad Classes', color='red')
    plt.xlabel('Time (Updates in time)', fontsize=14)
    plt.ylabel('Cumulative Distribution (number of particles)', fontsize=14)
    plt.title('Cumulative Distribution Over Time', fontsize=16)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    # Save the figure as an image (e.g., PNG)
    plt.savefig(fileName, dpi=300, bbox_inches='tight')