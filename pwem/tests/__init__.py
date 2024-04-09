#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
#                Laura del Cano         (ldelcano@cnb.csic.es)
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
from pyworkflow.tests import DataSet

def defineDatasets():
    """ Defines datasets for testing."""
    DataSet(name='igbmc_gempicker', folder='igbmc_gempicker',
            files={
                'allMics': 'micrographs/*.mrc',
                'templates': 'templates/*.mrc'})

    DataSet(name='xmipp_tutorial', folder='xmipp_tutorial',
            files={
                'aligned_particles': 'gold/aligned_particles.sqlite',
                'allMics': 'micrographs/*.mrc',
                'boxingDir': 'pickingEman/info',
                'ctfGold': 'gold/xmipp_ctf.ctfparam',
                'ctf': 'ctfs/ctfs.sqlite',
                'images10': 'gold/images10.xmd',
                'mic1': 'micrographs/BPV_1386.mrc',
                'mic2': 'micrographs/BPV_1387.mrc',
                'mic3': 'micrographs/BPV_1388.mrc',
                'micsGoldSqlite': 'gold/micrographs_gold.sqlite',
                'micsGoldSqlite2': 'gold/micrographs2_gold.sqlite',
                'micsSqlite': 'micrographs/micrographs.sqlite',
                'particles': 'particles/*.hdf',
                'particles1': 'particles/BPV_1386_ptcls.hdf',
                'particles2': 'particles/BPV_particles.sqlite',
                'particles3': 'particles/BPV_particles_aligned.sqlite',
                'posAllDir': 'pickingXmipp/pickedAll',
                'posSupervisedDir': 'pickingXmipp/pickedSupervised',
                'vol1': 'volumes/BPV_scale_filtered_windowed_64.vol',
                'vol2': 'volumes/volume_1_iter_002.mrc',
                'vol3': 'volumes/volume_2_iter_002.mrc',
                'vol4': 'volumes/helix_59_4__6_7.vol',
                'vol110': 'volumes/BPV_scale_filtered_windowed_110.vol',
                'volumes': 'volumes/*.mrc',
                'mask2d': 'masks/mask2d.xmp',
                'mask3d': 'masks/mask.vol',
                'ctfsDir': 'ctfs'})

    DataSet(name='mda', folder='hemoglobin_mda',
            files={
                'particles': 'particles/*.spi',
                'volumes': 'volumes/*.spi',
                'averages': 'averages/averages.stk'})

    DataSet(name='nma', folder='nma',
            files={'pdb': 'pdb/AK.pdb'})

    DataSet(name='relion_tutorial', folder='relion_tutorial',
            files={'allMics': 'micrographs/*.mrc',
                   'boxingDir': 'pickingEman',
                   'posAllDir': 'pickingXmipp',
                   'volume': 'volumes/reference.mrc'})

    DataSet(name='relion30_tutorial', folder='relion30_tutorial', files={})

    DataSet(name='ribo_movies', folder='ribo_movies',
            files={
                'movies': 'movies/1??_*.mrcs',
                'posAllDir': 'pickingXmipp',
                'volume': 'volumes/reference.mrc'})

    DataSet(name='model', folder='model',
            files={'classesSelection': 'gold/classes_selection.sqlite'})

    DataSet(name='rct', folder='rct',
            files={
                'classes': 'classes/classes2D_stable_core.sqlite',
                'positions': 'positions',
                'tilted': 'micrographs/F_rct_t*.tif',
                'untilted': 'micrographs/F_rct_u*.tif'})

    DataSet(name='groel', folder='groel',
            files={
                'averages': 'classes/level_classes.stk',
                'mic1': 'micrographs/groel-stain.mrc',
                'pickingDogpicker': 'pickingDogpicker'})

    # TODO: Files to be deleted after checking no plugin or old scipion
    #  installations (1.2, ...) are using them:
    #  reconstRotandShiftFlip_Gold_output.vol (moved to reconstruction/gold)
    #  reconstRotandShift_Gold_output.vol (moved to reconstruction/gold)
    #  reconstRotOnly_Gold_output.vol (moved to reconstruction/gold)

    DataSet(name='emx', folder='emx',
            files={
                'alignFlip': 'alignment/alignFlip.mrcs',
                'alignFlip2': 'alignment/alignFlip2.mrcs',
                'alignShiftRotExp': 'alignment/Test1/images.mrc',
                'alignShiftRotExpGold': 'alignment/gold/align_shift_rot_exp_gold.mrcs',

                'alignReconst': 'alignment/recons.mrcs',
                'alignReconstShift': 'alignment/recons_sh.mrcs',
                'alignRotOnly': 'alignment/align_rot_only.mrcs',
                'alignRotOnlyGold': 'alignment/gold/align_rot_only_gold.mrcs',
                # TODO: transform geometry not working in stacks of volumes
                'alignRotOnly3D': 'alignment/align_rot_only_3d.mrcs',
                # 'alignRotShift': 'alignment/align_rot_shift.mrcs',
                'alignShiftOnly': 'alignment/align_shift_only.mrcs',
                'alignShiftOnlyGold':
                    'alignment/gold/align_shift_only_gold.mrcs',
                # TODO: transform geometry not working in stacks of volumes
                'alignShiftOnly3D': 'alignment/align_shift_only_3d.mrcs',
                'alignShiftRot': 'alignment/align_shift_rot.mrcs',
                'alignShiftRotGold': 'alignment/gold/align_shift_rot_gold.mrcs',
                # TODO: transform geometry not working in stacks of volumes
                'alignShiftRot3D': 'alignment/align_shift_rot_3d.mrcs',
                'coordinatesGoldT1':
                    'coordinates/Test1/coordinates_gold.sqlite',
                'coordinatesT1': 'coordinates/Test1/coordinates.emx',
                'defocusParticleT2': 'defocusParticle/particles.emx',
                'emxMicrographCtf1': 'MicrographsCTF/ctfindCTFEstimation.emx',
                'emxMicrographCtf1Gold':
                    'MicrographsCTF/ctfindCTFEstimation.sqlite',
                'micrographsGoldT2': 'defocusParticle/micrographs_gold.sqlite',
                'particlesGoldT2': 'defocusParticle/particles_gold.sqlite',
                'reconstRotOnly': 'reconstruction/reconstRotOnly.mrcs',
                'reconstRotOnlyGold':
                    'reconstruction/gold/reconstRotOnly_output.vol',
                'reconstRotandShift': 'reconstruction/reconstRotandShift.mrcs',
                'reconstRotandShiftGold':
                    'reconstruction/gold/reconstRotandShift_output.vol',
                'reconstRotandShiftFlip':
                    'reconstruction/reconstRotandShiftFlip.mrcs',
                'reconstRotandShiftFlipGold':
                    'reconstruction/gold/reconstRotandShiftFlip_output.vol',
                'reconstRotOnlyGoldRln':
                    'reconstruction/gold/reconstRotOnly_rln_output.vol',
                'reconstRotandShiftGoldRln':
                    'reconstruction/gold/reconstRotandShift_rln_output.vol',
                'reconstRotandShiftFlipGoldRln':
                    'reconstruction/gold/reconstRotandShiftFlip_rln_output.vol'
            })

    DataSet(name='movies', folder='movies',
            files={'eer': 'FoilHole_29355520_Data_29330528_29330530_20200329_234551_Fractions.mrc.eer'})

    DataSet(name='jmbFalconMovies', folder='jmbFalconMovies',
            files={'movie1': 'Falcon_2014_07_16-16_32_43_0.mrcs'})

    DataSet(name='grigorieff', folder='grigorieff', files={})

    DataSet(name='initial_volume', folder='initial_volume',
            files={'groel': 'groel_64.stk',
                   'ribosome': 'ribosomeFrank_64_preprocessed.stk',
                   'bpv': 'bpv_100.stk'})

    DataSet(name='general', folder='general',
            files={'vol_helix': 'volumes/helix_59_4__6_7.vol'})

    DataSet(name='SymVirus', folder='SymVirus',
            files={'whole_vol_half': 'whole_vol_half.vol',
                   'large_vol_half_th': 'large_vol_half_th.vol',
                   'small_vol_half_th': 'small_vol_half_th.vol'})

    DataSet(name='xmipp_programs', folder='xmipp_programs', files={})

    DataSet(name='resmap', folder='resmap',
            files={'betagal': 'betaGal.mrc',
                   'betagal_half1': 'betaGal_half01.mrc',
                   'betagal_half2': 'betaGal_half02.mrc',
                   'betagal_mask': 'betaGal_mask.mrc'})

    DataSet(name='10010', folder='10010',
            files={'initialVolume': 'initialVolume.vol',
                   'particles': 'particles.sqlite'})

    DataSet(name='model_building_tutorial', folder='model_building_tutorial',
            files={'vol1': 'volumes/1ake_4-5A.mrc',
                   'pdb1': 'PDBx_mmCIF/1ake_start.pdb'})

    DataSet(name='eman', folder='eman',
            files={
                'micU': 'mics/ip3r10252011-0005_0-2.hdf',
                'micT': 'mics/ip3r10252011-0005_10.hdf',
                'volume': '3DmapIP3R1_clip_376apix.mrc',
                'coords': 'coords/*.json',
                'ctfs': 'ctfs/'
            })

    DataSet(name='deepConsensusPicking', folder='deepConsensusPicking',
            files={
                'allMics': 'micrographs/*.mrc',
                'autoPickingCoordinates': 'autoPickingCoordinates'
            })
