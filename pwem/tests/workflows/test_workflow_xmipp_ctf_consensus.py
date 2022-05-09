# ***************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
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
# ***************************************************************************/

import pyworkflow.tests as pwtests

from pwem import Domain
import pwem.protocols as emprot


class TestCtfConsensus(pwtests.BaseTest):
    """ Check if the Xmipp-CTFconsensus rejects CTFs (and the corresponding mics)
        when two CTF estimations give different results,
        and accept when the two estimations give similar results.
    """

    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.dataset = pwtests.DataSet.getDataSet('xmipp_tutorial')
        cls.micsFn = cls.dataset.getFile('allMics')

    def checkCTFs(self, protConsensus, refCTFs, refMics, label='',
                  avgCTF=None, MDmerging=False):
        outputCTF = getattr(protConsensus, "outputCTF" + label, None)
        outputMicrographs = getattr(protConsensus, "outputMicrographs" + label, None)

        self.assertIsNotNone(outputCTF,
                             "There was a problem with the CTF-Consensus. "
                             "No outputCTF%s is created." % label)
        self.assertIsNotNone(outputMicrographs,
                             "There was a problem with the CTF-Consensus. "
                             "No outputMicrographs%s is created." % label)

        self.assertEqual(outputCTF.getSize(), refCTFs.getSize(),
                         "The outputCTF%s size is wrong." % label)
        self.assertEqual(outputMicrographs.getSize(), refCTFs.getSize(),
                         "The outputMicrographs%s size is wrong." % label)
        self.assertTupleEqual(outputMicrographs.getDim(), refMics.getDim(),
                              "The outputMicrographs%s dimension is wrong." % label)

        # we will check that the first CTF in the set make sense.
        firstCTF = outputCTF.getFirstItem()

        if not avgCTF:
            refCTF = refCTFs.getFirstItem()
            self.assertTrue(firstCTF.equalAttributes(refCTF),
                            "The outputCTF%s has different attributes "
                            "than the input." % label)
        else:
            self.assertAlmostEqual(avgCTF['defocusU'], firstCTF.getDefocusU(), delta=200,
                                   msg="DefocusU doesn't match when defocus averaging.")
            self.assertAlmostEqual(avgCTF['defocusV'], firstCTF.getDefocusV(), delta=100,
                                   msg="DefocusV doesn't match when defocus averaging.")
            self.assertAlmostEqual(avgCTF['defocusAngle'], firstCTF.getDefocusAngle(), delta=1,
                                   msg="DefocusAngle doesn't match when defocus averaging.")

        if MDmerging:
            for label in ['CritCorr13', 'CritIceness', 'CritCtfMargin',
                          'CritPsdCorr90', 'Q0', 'CritNonAstigmaticValidty',
                          'CritfirstZeroRatio', 'CritFirstZero']:
                MDlabel = '_xmipp_ctf' + label
                self.assertTrue(MDlabel in firstCTF.getObjDict().keys(),
                                "'%s' metadata not found in the result. "
                                "Bad merging" % MDlabel)

    def test1(self):
        # Import a set of micrographs
        protImport = self.newProtocol(emprot.ProtImportMicrographs,
                                      filesPath=self.micsFn,
                                      samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs,
                             "There was a problem with the import")

        # Create pure noise micrographs (to force a Discarded consensus set)
        protStream = self.newProtocol(emprot.ProtCreateStreamData,
                                      xDim=9216,
                                      yDim=9441,
                                      nDim=3,
                                      samplingRate=1.237,
                                      setof=2,  # 2 -> SetOfRandomMicrographs
                                      creationInterval=1)
        self.proj.launchProtocol(protStream, wait=False)

        # Computes the CTF with Xmipp
        XmippProtCTFMicrographs = Domain.importFromPlugin('xmipp3.protocols',
                                                          'XmippProtCTFMicrographs',
                                                          doRaise=True)
        protCTF1 = self.newProtocol(XmippProtCTFMicrographs)
        protCTF1.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protCTF1, wait=False)

        # Computes the CTF with CTFFind4
        ProtCTFFind = Domain.importFromPlugin('cistem.protocols',
                                              'CistemProtCTFFind', doRaise=True)
        protCTF2 = self.newProtocol(ProtCTFFind)
        protCTF2.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protCTF2, wait=False)

        self._waitOutput(protStream, "outputMicrographs")
        # Computes the CTF with CTFFind4 for the noise mics
        protCTF3 = self.newProtocol(ProtCTFFind)
        protCTF3.inputMicrographs.set(protStream.outputMicrographs)
        self.proj.launchProtocol(protCTF3, wait=False)

        # Computes the Consensus of GOOD CTFs
        self._waitOutput(protCTF1, "outputCTF")
        self._waitOutput(protCTF2, "outputCTF")
        XmippProtCTFConsensus = Domain.importFromPlugin('xmipp3.protocols',
                                                        'XmippProtCTFConsensus',
                                                        doRaise=True)
        protCTFcons = self.newProtocol(XmippProtCTFConsensus,
                                       objLabel='default (pass)',
                                       useDefocus=False,
                                       useAstigmatism=False,
                                       useResolution=False,
                                       calculateConsensus=True,
                                       averageDefocus=False,
                                       includeSecondary=False)
        protCTFcons.inputCTF.set(protCTF1.outputCTF)
        protCTFcons.inputCTF2.set(protCTF2.outputCTF)
        self.launchProtocol(protCTFcons)

        protCTF1.outputCTF.load()  # Needed to update the setSize
        self.checkCTFs(protCTFcons,
                       refMics=protImport.outputMicrographs,
                       refCTFs=protCTF1.outputCTF)

        # Computes the Consensus comparing a good CTF to a RANDOM one
        self._waitOutput(protCTF3, "outputCTF")
        protCTFcons2 = self.newProtocol(XmippProtCTFConsensus,
                                        objLabel='default (block)',
                                        useDefocus=False,
                                        useAstigmatism=False,
                                        useResolution=False,
                                        calculateConsensus=True,
                                        averageDefocus=False,
                                        includeSecondary=False)

        protCTFcons2.inputCTF.set(protCTF1.outputCTF)
        protCTFcons2.inputCTF2.set(protCTF3.outputCTF)
        self.launchProtocol(protCTFcons2)
        self.checkCTFs(protCTFcons2,
                       refMics=protImport.outputMicrographs,
                       refCTFs=protCTF1.outputCTF,
                       label="Discarded")

        # Averaging CTF parameters
        protCTFcons3 = self.newProtocol(XmippProtCTFConsensus,
                                        objLabel='defocus average',
                                        useDefocus=False,
                                        useAstigmatism=False,
                                        useResolution=False,
                                        calculateConsensus=True,
                                        averageDefocus=True)
        protCTFcons3.inputCTF.set(protCTF1.outputCTF)
        protCTFcons3.inputCTF2.set(protCTF2.outputCTF)
        self.launchProtocol(protCTFcons3)

        protCTF1.outputCTF.load()  # Needed to update the set
        protCTF2.outputCTF.load()  # Needed to update the set
        ctfAveraged = {'defocusU': 24025.6729,
                       'defocusV': 23610.2071,
                       'defocusAngle': 57.1943}
        self.checkCTFs(protCTFcons3,
                       refMics=protImport.outputMicrographs,
                       refCTFs=protCTF1.outputCTF,
                       avgCTF=ctfAveraged)

        # merging Metadata columns
        protCTFcons4 = self.newProtocol(XmippProtCTFConsensus,
                                        objLabel='metadata merge',
                                        useDefocus=False,
                                        useAstigmatism=False,
                                        useResolution=False,
                                        calculateConsensus=True,
                                        includeSecondary=True)
        protCTFcons4.inputCTF.set(protCTF2.outputCTF)
        protCTFcons4.inputCTF2.set(protCTF1.outputCTF)
        self.launchProtocol(protCTFcons4)

        protCTF1.outputCTF.load()  # Needed to update the set
        protCTF2.outputCTF.load()  # Needed to update the set
        self.checkCTFs(protCTFcons4,
                       refMics=protImport.outputMicrographs,
                       refCTFs=protCTF2.outputCTF,
                       MDmerging=True)
