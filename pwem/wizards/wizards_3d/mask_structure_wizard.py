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
import tkinter as tk

from matplotlib.widgets import RadioButtons, Slider
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from pyworkflow.gui.plotter import plt

from ...convert import AtomicStructHandler
from .callbacks import DraggablePoint


class MaskStructureWizard(object):
    '''Create a mask for an atomic structure interactively. Masks currently implemented:
            - Spherical mask
    '''
    def __init__(self, filename):
        ah = AtomicStructHandler()
        ah.read(filename)
        atomIterator = ah.getStructure().get_atoms()
        coords = np.asarray([np.append(atom.get_coord(), 1) for atom in atomIterator])
        self.shift = np.mean(coords, axis=0)
        self.coords = np.copy(coords) - self.shift
        self.coordsDownsampled = np.copy(coords) - self.shift
        self.origin = np.array([0, 0, 0])
        self.pressed = False
        self.radio = None
        self.cb = None
        self.sphere_artist = None
        self.radius = 0
        self.running = True
        self.root = tk.Tk()
        self.root.resizable(False, False)  # For matplotlib <= 3.3.x
        self.fig = plt.Figure(figsize=plt.figaspect(1)*1.5)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.ax_3d = self.fig.add_subplot(projection='3d')
        # plt.style.use('dark_background')

    def get_sphere_params(self):
        origin_shifted = self.origin + self.shift
        return np.hstack([origin_shifted, self.radius])

    def is_window_closed(self):
        self.running = False

    def set_axes_equal(self, ax: plt.Axes):
        """Set 3D plot axes to equal scale.

        Make axes of 3D plot have equal scale so that spheres appear as
        spheres and cubes as cubes.  Required since `ax.axis('equal')`
        and `ax.set_aspect('equal')` don't work on 3D.
        """
        limits = np.array([
            ax.get_xlim3d(),
            ax.get_ylim3d(),
            ax.get_zlim3d(),
        ])
        origin = np.mean(limits, axis=1)
        radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
        self._set_axes_radius(ax, origin, radius)

    def _set_axes_radius(self, ax, origin, radius):
        x, y, z = origin
        ax.set_xlim3d([x - radius, x + radius])
        ax.set_ylim3d([y - radius, y + radius])
        ax.set_zlim3d([z - radius, z + radius])

    def plotScatter(self):
        self.ax_3d.clear()
        xi = self.coordsDownsampled[:, 0]
        yi = self.coordsDownsampled[:, 1]
        zi = self.coordsDownsampled[:, 2]
        ori_x, ori_y, ori_z = self.origin[0], self.origin[1], self.origin[2]
        plt.ion()
        self.ax_3d.scatter(xi, yi, zi, s=12, c='purple', edgecolors='k', alpha=0.3)
        self.M = [-self.ax_3d.azim * np.pi / 180, self.ax_3d.elev * np.pi / 180, 0]
        scatter_origin = self.ax_3d.scatter(ori_x, ori_y, ori_z, s=100, c='cyan', edgecolors='k')
        self.plot_sphere(self.radius)
        self.dr = DraggablePoint(self.origin, self.fig, self.ax_3d, scatter_origin, self.M)
        self.ax_3d.set_axis_off()
        # self.ax_3d.set_box_aspect([1, 1, 1]) # For matplotlib => 3.3.x
        self.set_axes_equal(self.ax_3d)

    def plot_sphere(self, radius):
        self.radius = radius
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        x = radius * np.outer(np.cos(u), np.sin(v)) + self.origin[0]
        y = radius * np.outer(np.sin(u), np.sin(v)) + self.origin[1]
        z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + self.origin[2]
        if self.sphere_artist == None:
            self.sphere_artist = self.ax_3d.plot_surface(x, y, z, rstride=4, cstride=4, color='r',
                                                         linewidth=0, alpha=0.1)
        else:
            self.sphere_artist.remove()
            self.sphere_artist = self.ax_3d.plot_surface(x, y, z, rstride=4, cstride=4, color='r',
                                                         linewidth=0, alpha=0.1)

    def downsamplingPC(self, voxel_size):
        mode = 'barycenter'
        non_empty_voxel_keys, inverse, nb_pts_per_voxel = np.unique(((self.coords - np.min(self.coords, axis=0))
                                                                     // voxel_size).astype(int), axis=0,
                                                                    return_inverse=True,
                                                                    return_counts=True)
        idx_pts_vox_sorted = np.argsort(inverse)
        voxel_grid = {}
        grid_barycenter, grid_candidate_center = [], []
        last_seen = 0
        for idx, vox in enumerate(non_empty_voxel_keys):
            voxel_grid[tuple(vox)] = self.coords[idx_pts_vox_sorted[last_seen:last_seen + nb_pts_per_voxel[idx]]]
            grid_barycenter.append(np.mean(voxel_grid[tuple(vox)], axis=0))
            last_seen += nb_pts_per_voxel[idx]
        if mode == 'barycenter':
            self.coordsDownsampled = np.asarray(grid_barycenter)
        self.plotScatter()

    def press_shift(self, event):
        if event.key == 'shift':
            self.pressed = True
            self.dr.connect()

    def release_shift(self, event):
        if self.pressed and event.key == 'shift':
            self.pressed = False
            self.dr.disconnect()
            self.fig.canvas.mpl_connect('button_release_event', self.on_release)
            self.origin = self.dr.point
            self.plot_sphere(self.radius)
            self.fig.canvas.draw()

    def on_release(self, event):
        # self.M = event.inaxes.get_proj()
        self.M = [-self.ax_3d.azim * np.pi / 180, self.ax_3d.elev * np.pi / 180, 0]
        self.dr.M = self.M

    def change_view(self, event):
        if event == "X":
            self.ax_3d.view_init(elev=0., azim=0.)
            self.M = [-0., 0., 0]
            self.dr.M = self.M
        elif event == "Y":
            self.ax_3d.view_init(elev=0., azim=90.)
            self.M = [-np.pi / 2., 0., 0]
            self.dr.M = self.M
        elif event == "Z":
            self.ax_3d.view_init(elev=90., azim=0.)
            self.M = [0, np.pi / 2., 0]
            self.dr.M = self.M
        self.fig.canvas.draw()

    def initializePlot(self):
        self.downsamplingPC(5.01)

        # Buttons
        axcolor = 'grey'
        # rax = self.fig.add_axes([0.1, 0.4, 0.12, 0.25], facecolor=axcolor)  # For matplotlib => 3.3.x
        rax = self.fig.add_axes([0.05, 0.5, 0.12, 0.15], facecolor=axcolor)  # For matplotlib <= 3.3.x
        self.radio = RadioButtons(rax, ('X', 'Y', 'Z'), activecolor='navy')
        self.radio.on_clicked(self.change_view)
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)

        # Slider
        sax = self.fig.add_axes([0.2, 0.02, 0.65, 0.03], facecolor=axcolor)
        self.slider = Slider(sax, 'Downsampling', 0.01, 10, valinit=5.01, valstep=0.2, color='navy')
        self.slider.on_changed(self.downsamplingPC)

        srax = self.fig.add_axes([0.2, 0.06, 0.65, 0.03], facecolor=axcolor)
        max_radius = np.round(1.5 * np.amax(np.linalg.norm(self.coords, axis=1)))
        self.slider_radius = Slider(srax, 'Radius', 0, max_radius, valinit=0, valstep=1, color='navy')
        self.slider_radius.on_changed(self.plot_sphere)

        # Toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, self.root)
        toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)
        self.root.protocol("WM_DELETE_WINDOW", self.is_window_closed)

        self.fig.canvas.mpl_connect('key_press_event', self.press_shift)
        self.canvas.mpl_connect('key_release_event', self.release_shift)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)

        # GUI Running Loop
        while self.running:
            self.root.update_idletasks()
            self.root.update()
        self.root.destroy()
