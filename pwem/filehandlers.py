# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
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
"""
This modules contains file handlers to be registered in the scipion file browser
"""
import os
from os.path import dirname

from pyworkflow.gui.browser import FileHandler, isStandardImage
from pyworkflow import gui
import pyworkflow.utils as pwutils

class ImageFileHandler(FileHandler):
    import xmippLib
    _image = xmippLib.Image()
    _index = ''

    def _getImageString(self, filename):
        if isStandardImage(filename):
            return "Image file."
        x, y, z, n = xmippLib.getImageSize(filename)
        objType = 'Image'
        dimMsg = "*%(objType)s file*\n  dimensions: %(x)d x %(y)d"
        expMsg = "Columns x Rows "
        if z > 1:
            dimMsg += " x %(z)d"
            expMsg += " x Slices"
            objType = 'Volume'
        if n > 1:
            dimMsg += " x %(n)d"
            expMsg += " x Objects"
            objType = 'Stack'
        return (dimMsg + "\n" + expMsg) % locals()

    def _getImagePreview(self, filename):
        dim = 128

        if isStandardImage(filename):
            self.tkImg = gui.getImage(os.path.abspath(filename),
                                      tkImage=True, maxheight=dim)
        else:
            fn = self._index + filename
            self.tkImg = gui.getTkImage(self._image, fn, dim)

        return self.tkImg

    def getFilePreview(self, objFile):
        fn = objFile.getPath()
        return self._getImagePreview(fn), self._getImageString(fn)

    def getFileActions(self, objFile):
        from viewers import DataView
        fn = objFile.getPath()
        return [('Open with Xmipp viewer', lambda: DataView(fn).show(),
                 pwutils.Icon.ACTION_VISUALIZE)]


class ParticleFileHandler(ImageFileHandler):
    def getFileIcon(self, objFile):
        return 'file_image.gif'


class VolFileHandler(ImageFileHandler):
    def getFileIcon(self, objFile):
        return 'file_vol.gif'


class StackHandler(ImageFileHandler):
    _index = '1@'

    def getFileIcon(self, objFile):
        return 'file_stack.gif'


class ChimeraHandler(FileHandler):

    def getFileActions(self, objFile):
        from viewers import ChimeraView
        fn = objFile.getPath()
        return [('Open with Chimera', lambda: ChimeraView(fn).show(),
                 pwutils.Icon.ACTION_VISUALIZE)]

    def getFileIcon(self, objFile):
        return 'file_text.gif'


class MdFileHandler(ImageFileHandler):
    def getFileIcon(self, objFile):
        return 'file_md.gif'

    def _getImgPath(self, mdFn, imgFn):
        """ Get ups and ups until finding the relative location to images. """
        path = dirname(mdFn)
        import xmippLib
        index, fn = xmippLib.FileName(imgFn).decompose()

        while path and path != '/':
            newFn = os.path.join(path, fn)
            if os.path.exists(newFn):
                if index:
                    newFn = '%d@%s' % (index, newFn)
                return newFn
            path = dirname(path)

        return None

    def _getMdString(self, filename, block=None):
        md = xmippLib.MetaData()
        if block:
            md.read(block + '@' + filename)
        else:
            md.read(filename, 1)
        labels = md.getActiveLabels()
        msg = "Metadata items: *%d*\n" % md.getParsedLines()
        msg += "Metadata labels: " + ''.join(
            ["\n   - %s" % xmippLib.label2Str(l)
             for l in labels])

        imgPath = None
        for label in labels:
            if xmippLib.labelIsImage(label):
                imgPath = self._getImgPath(filename,
                                           md.getValue(label, md.firstObject()))
                break
        if imgPath:
            self._imgPreview = self._getImagePreview(imgPath)
            self._imgInfo = self._getImageString(imgPath)
        return msg

    def getFilePreview(self, objFile):
        self._imgPreview = None
        self._imgInfo = None
        filename = objFile.getPath()
        ext = pwutils.getExt(filename)

        if ext == '.xmd' or ext == '.ctfparam' or ext == '.pos' or ext == '.doc':
            msg = "*Metadata File* "
            blocks = xmippLib.getBlocksInMetaDataFile(filename)
            nblocks = len(blocks)
            if nblocks <= 1:
                mdStr = self._getMdString(filename)
                msg += "  (single block)\n"
                if self._imgInfo:
                    msg += "\nFirst item: \n" + self._imgInfo
                msg += '\n' + mdStr
            else:
                mdStr = self._getMdString(filename, blocks[0])
                msg += "  (%d blocks) " % nblocks
                if self._imgInfo:
                    msg += "\nFirst item: \n" + self._imgInfo
                msg += "\nFirst block: \n" + mdStr
                msg += "\nAll blocks:" + ''.join(
                    ["\n  - %s" % b for b in blocks])
        elif ext == '.star':
            msg = "*Relion STAR file* \n"
            msg += self._getMdString(filename)

        return self._imgPreview, msg