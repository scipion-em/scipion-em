# **************************************************************************
# *
# * Authors: David Herreros Calero    (dherreros@cnb.csic.es)
# *
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


import numpy as np
import re

from pwem.convert.transformations import euler_matrix


class DraggablePoint:
    '''Class to create a dragabble point in Matplotlib using a scatter plot (in 3D)'''
    lock = None  # Only one can be animated at a time (probably useful in the future)
    def __init__(self, point, figure, axes, plot, M):
        self.point = point
        self.figure = figure
        self.axes = axes
        self.plot = plot
        self.M = M
        self.press = None
        self.background = None

    def getxyz(self, event):
        s = self.axes.format_coord(event.xdata, event.ydata)
        out = [float(re.sub(r'[^\x00-\x7F]+', '-', x.split('=')[1].strip())) for x in s.split(',')]
        return out

    def remove_projection_direction(self, coords, prev_point):
        prev_point = np.copy(prev_point)
        rot_matrix = euler_matrix(self.M[0], self.M[1], self.M[2], 'szyx')
        prev_point_rot = np.dot(rot_matrix, np.hstack([prev_point, 1]))
        new_point = np.dot(rot_matrix, np.hstack([coords, 1]))
        new_point[0] = prev_point_rot[0]
        coords = np.dot(np.linalg.inv(rot_matrix), new_point)
        return coords

    def connect(self):
        'connect to all the events we need'
        self.axes.disable_mouse_rotation()
        self.cidpress = self.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = self.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        if event.inaxes != self.axes:
            return
        if DraggablePoint.lock is not None:
            return
        coords = self.getxyz(event)
        coords = self.remove_projection_direction(coords, self.point)
        self.point = np.array((coords[0], coords[1], coords[2]))
        self.press = self.point, coords[0], coords[1], coords[2]
        DraggablePoint.lock = self

        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.figure.canvas
        axes = self.axes
        self.plot.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.axes.bbox)

        # now redraw just the rectangle
        self.plot.set_offsets([coords[0], coords[1]])
        self.plot.set_3d_properties(coords[2], 'z')

        # and blit just the redrawn area
        self.plot.do_3d_projection(renderer=self.figure._cachedRenderer)
        canvas.restore_region(self.background)
        self.axes.draw_artist(self.plot)
        canvas.flush_events()
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        if DraggablePoint.lock is not self:
            return
        if event.inaxes != self.axes:
            return
        self.point, xpress, ypress, zpress = self.press
        coords = self.getxyz(event)
        coords = self.remove_projection_direction(coords, self.point)
        dx = coords[0] - xpress
        dy = coords[1] - ypress
        dz = coords[2] - zpress
        self.point = np.array((self.point[0]+dx, self.point[1]+dy, self.point[2]+dz))

        canvas = self.figure.canvas
        axes = self.axes
        # restore the background region
        canvas.restore_region(self.background)

        # now redraw just the rectangle
        self.plot.set_offsets([self.point[0], self.point[1]])
        self.plot.set_3d_properties(self.point[2], 'z')

        # and blit just the redrawn area
        self.plot.do_3d_projection(renderer=self.figure._cachedRenderer)
        canvas.restore_region(self.background)
        self.axes.draw_artist(self.plot)
        canvas.flush_events()
        canvas.blit(axes.bbox)

    def on_release(self, event):
        'on release we reset the press data'
        if DraggablePoint.lock is not self:
            return

        self.press = None
        DraggablePoint.lock = None

        # turn off the rect animation property and reset the background
        self.plot.set_animated(False)
        self.background = None

        # redraw the full figure
        self.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.axes.mouse_init()
        self.figure.canvas.mpl_disconnect(self.cidpress)
        self.figure.canvas.mpl_disconnect(self.cidrelease)
        self.figure.canvas.mpl_disconnect(self.cidmotion)
