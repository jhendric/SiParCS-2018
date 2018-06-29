'''#!/Users/wbd1/anaconda3/bin/python3'''

from tkinter import *
from tkinter import ttk

import mkl
import xarray as xa

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.collections import LineCollection, PolyCollection
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy.crs as ccrs
import cartopy.feature
from cartopy.mpl.patch import geos_to_path

import numpy as np
import pandas as pd
import math
from plot_2D_obs_initial import plot_2D_obs
np.set_printoptions(threshold=np.nan) #without this setting, self.levels will be incomplete
import time
import itertools



'''

Incorporates plot_2D_obs_initial.py into a GUI with slider bars
and drop down menus.

'''



class GUI_3D_obs_initial:

    
    
    def __init__(self, window, grid_col, grid_row):
        
        self.plotter = plot_2D_obs('../obs_series/obs_epoch_001.nc')
        self.original_data = self.plotter.data
        #print('initial test: ', np.unique(self.plotter.data.where(
        #    self.plotter.data.obs_types == 5, drop = True).z.values))
        #print('initial test: ')
        #x = self.plotter.data.where(self.plotter.data.obs_types == 5, drop = True)
        #print(x.where(abs(x.z - 460.0) < 1e-10, drop = True))
        print(np.unique(self.plotter.obs_types.values).size)
        self.window = window
        self.window.grid_columnconfigure(0, weight = 1)
        self.window.grid_rowconfigure(0, weight = 1)

        #a mainframe
        self.main_frame = ttk.Frame(self.window, padding = "8")
        self.main_frame.grid(column = grid_col, row = grid_row, sticky = "N, S, E, W") 
        self.main_frame.grid_columnconfigure(0, weight = 1) #weights for whole grid
        self.main_frame.grid_rowconfigure(0, weight = 1) #weights for whole grid

        self.style = ttk.Style()

        #obs parameter variables

        #get counts for each obs_type
        unique, counts = np.unique(self.plotter.data.obs_types.values, return_counts = True)
        count_dict = dict(zip(unique, counts))
        
        #strip useless characters from string interpretation of obs types, and add counts
        obs_type_dict_sparse = [x.replace('[', '').replace(']', '').
                                replace(',', '').replace('\'', '').
                                replace('dict_keys(', '')
                                for x in self.plotter.obs_type_dict.keys()]
        
        self.obs_type_names = StringVar(value = [str(count_dict[self.plotter.obs_type_dict[x]]) +
                                                 " : " + x for x in obs_type_dict_sparse])
                                                 

        #GUI config

        #observation selection
        self.obs_frame = ttk.Frame(self.main_frame, padding = "2")
        self.obs_frame.grid(column = 2, row = 1, sticky = "N, S, E, W")
        ttk.Label(self.obs_frame, text = "Observation Type Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.obs_menu = Listbox(self.obs_frame, listvariable = self.obs_type_names,
                                height = 18, width = 40, exportselection = False)
        self.obs_menu.grid(column = 1, row = 2, rowspan = 1, sticky = "N, S, E, W")
        
        self.obs_menu.bind('<Return>', lambda event : self.populate('qc', self.qc_menu, event))
        
        self.obs_menu.selection_set(0)
        self.obs_menu.event_generate('<<ListboxSelect>>')
        #print('size of obs box: ', self.obs_menu.size())

        #obs scrollbar
        self.obs_bar = ttk.Scrollbar(self.obs_frame, orient = VERTICAL, command = self.obs_menu.yview)
        self.obs_menu.configure(yscrollcommand = self.obs_bar.set)
        self.obs_bar.grid(column = 2, row = 2, rowspan = 2,  sticky = "N, S, E")
         
        self.qc = StringVar()
        
        #qc selection
        
        self.qc_frame = ttk.Frame(self.main_frame, padding = "2")
        self.qc_frame.grid(column=2, row = 2, sticky = "N, S, E, W")
        ttk.Label(self.qc_frame, text = "DART QC Value Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.qc_menu = Listbox(self.qc_frame, listvariable = self.qc,
                               height = 8, width = 40, selectmode = "extended", exportselection = False)
        self.qc_menu.grid(column = 1, row = 2, sticky ="N, S, E, W")

        #for use in populating and clearing menus (in populate function)
        self.data_obs_types = 1
        self.data_qc = 2
        
        self.data_request_dict = {
            'data_qc' : 'obs_types'
        }
        
        self.menu_hierarchy = [self.obs_menu, self.qc_menu]
        self.data_hierarchy = ['original_data', 'data_qc']
        
        #populate levels

        #populate qc
        self.qc_key = {0 : '0 - Assimilated O.K.',
                       1 : '1 - Evaulated O.K., not assimilated because namelist specified evaluate only',
                       2 : '2 - Assimilated, but posterior forward operator failed',
                       3 : '3 - Evaluated, but posterior forward operator failed',
                       4 : '4 - Prior forward operator failed',
                       5 : '5 - Not used because of namelist control',
                       6 : '6 - Rejected because incoming data QC higher than namelist control',
                       7 : '7 - Rejected because of outlier threshold test',
                       8 : '8 - Failed vertical conversion'}
        
        self.populate('qc', self.qc_menu)

        #current plotting occurs only with press of enter from qc menu
        self.qc_menu.bind('<Return>', self.plot_3D)
        for i in range(len(self.qc.get())):
            self.qc_menu.selection_set(i)
        self.qc_menu.event_generate('<<ListboxSelect>>')

        
        #qc scrollbar
        self.qc_bar = ttk.Scrollbar(self.qc_frame, orient = HORIZONTAL, command = self.qc_menu.xview)
        self.qc_menu.configure(xscrollcommand = self.qc_bar.set)
        self.qc_bar.grid(column = 1, row = 3, rowspan = 1, sticky ="N, S, E, W")

        #selection of what value type to plot
        self.val_type = StringVar()
        self.val_type.set('QC')
        self.fill_frame = ttk.Frame(self.main_frame, padding = "2")
        self.fill_frame.grid(column = 2, row = 3, sticky = "N, S, E, W")
        ttk.Label(self.fill_frame, text = "Please select type of value to plot").grid(column = 1, row = 1, sticky= "E, W")
        self.qc_button = ttk.Radiobutton(self.fill_frame, text = 'QC', variable = self.val_type, value = 'QC')
        self.qc_button.grid(column = 1, row = 2, sticky = "N, S, E, W")
        self.val_button = ttk.Radiobutton(self.fill_frame, text = 'Observation value', variable = self.val_type,
                                          value = "Observation value")
        self.val_button.grid(column = 1, row = 3, sticky = "N, S, E, W")

        #plot button
        self.plot_button = ttk.Button(self.main_frame, text = "Plot", command = self.plot_3D, padding = "2")
        self.plot_button.grid(column = 2, row = 4, sticky = "N, S, E, W")

        
        #for plotting later
        self.markers = ['o', 'v', 'H', 'D', '^', '<', '8',
                        's', 'p', '>',  '*', 'h', 'd']

        #these markers do not seem to have border color capabilities
        #'x', '_', '|'

        s = ttk.Style()
        s.theme_use('clam')

    def populate(self, variable_name, menu, event = None):
        #event arg is passed by menu events, variable is the data variable to be manipulated, menu is the menu to change
        print('populating ' + variable_name)
        #clear lower level menus
        for i in range(self.menu_hierarchy.index(menu), len(self.menu_hierarchy)):
            self.menu_hierarchy[i].delete('0', 'end')

        #get currently selected values
        
        indices = None

        #used to dynamically access object variables
        var = 'data_' + variable_name

        #TODO: this differs from the 2D obs and should be more modularized at some point to match hierarchy idea
        if var == 'data_qc':
            indices = [self.plotter.obs_type_dict[self.obs_menu.get(val).split(" : ", 1)[1]]
                       for val in self.obs_menu.curselection()]

        #retrieve relevant data for this level of the hierarchy
        setattr(self, var,
                self.plotter.filter_test(getattr(self, self.data_hierarchy[self.data_hierarchy[1:].index(var)]),
                                         (self.data_request_dict[var], indices)))

        #set corresponding menu variables
        if var == 'data_levels':
            self.levels.set(value = np.unique(getattr(self, var).z.values))

        elif var == 'data_qc':
            unique, counts = np.unique(getattr(self, var).qc_DART.values, return_counts = True)
            count_dict = dict(zip(unique, counts))
            print(count_dict)
            self.qc.set(value = [str(count_dict[val]) + " : " + str(self.qc_key[val]) for val in unique])
            print(np.unique(getattr(self, var).qc_DART.values).size)
            
        #should work in class scope since menu is a self variable
        if (menu.get(0) == '['):
            menu.delete('0')
            
        if (menu.get(0)[0] == '['):
            first = menu.get(0)[1:]
            menu.delete('0')
            menu.insert(0, first)
            
        if (menu.get('end')[-1] == ']'):
            last = menu.get('end')[:-1]
            menu.delete('end')
            menu.insert(END, last)    
        
    def plot_3D(self, event = None):
        a = time.time()
        #event arg is passed by menu events
        print('plotting')
        #print(self.plotter.obs_type_dict.values())
        print('currently selected ob types: ', self.obs_menu.curselection())

        qc = [np.int64(self.qc_menu.get(val).split(": ", 1)[1][0]) for val in self.qc_menu.curselection()]

        print(qc)
        
        #make figure and canvas to draw on
        fig = Figure(figsize = (12,8))
        #ax = fig.add_axes([0.01, 0.01, 0.98, 0.98], projection = ccrs.PlateCarree())
        canvas = FigureCanvasTkAgg(fig, master = self.main_frame)
        canvas.get_tk_widget().grid(column = 1, row = 1, rowspan = 3, sticky = "N, S, E, W")
        self.main_frame.grid_columnconfigure(1, weight = 1)
        self.main_frame.grid_rowconfigure(1, weight = 1)

        
        #have to set up a separate toolbar frame because toolbar doesn't like gridding with others
        self.toolbar_frame = ttk.Frame(self.main_frame)
        self.toolbar = NavigationToolbar2TkAgg(canvas, self.toolbar_frame)
        self.toolbar_frame.grid(column = 1, row = 4, sticky = "N, S, E, W")
        #disable part of the coordinate display functionality, else everything flickers (may need for smaller window)
        #ax.format_coord = lambda x, y: ''
        
        data = self.plotter.filter_test(self.data_qc, ('qc_DART', qc))
        
        target_projection = ccrs.PlateCarree()

        #weird redundant line that works for transformation
        transform = target_projection.transform_points(target_projection, data.lons.values, data.lats.values)

        lons, lats = transform[:, 0], transform[:, 1]

        concat = lambda iterable: list(itertools.chain.from_iterable(iterable))

        feature = cartopy.feature.NaturalEarthFeature('Physical', 'land', '110m')
        
        geoms = feature.geometries()
        
        geoms = [target_projection.project_geometry(geom, feature.crs) for geom in geoms]
        
        paths = concat(geos_to_path(geom) for geom in geoms)
        
        polys = concat(path.to_polygons() for path in paths)

        ax = Axes3D(fig, xlim = [-180, 180], ylim = [-90, 90])
        z_max = max(data.z.values)

        #deal with orientation of graph for different vert types
        if data.vert_types.values[0] < 0 or data.vert_types.values[0] == 3:
            #surface, heights, or unknown units where up is positive direction
            lc = Poly3DCollection(polys, edgecolor = 'black', facecolor = 'green', closed = False)
            ax.set_zlim(bottom = 0, top = z_max)
            ax.add_collection3d(lc)
        else:
            #level or pressure units where down is positive direction
            polys = [[(point[0], point[1], z_max) for point in shape] for shape in polys]
            lc = Poly3DCollection(polys, edgecolor = 'black', facecolor = 'green', closed = False)
            ax.set_zlim(bottom = z_max, top = 0)
            ax.add_collection3d(lc)
        
        #print(data.obs_types.values)
        print(data.obs_types.values.size)
        print(np.unique(data.qc_DART.values))
        print(np.unique(data.obs_types.values))

        '''
        #get indices where obs_types change (array is sorted in filter_disjoint)
        indices = np.where(data.obs_types.values[:-1] != data.obs_types.values[1:])[0]
        indices[0:indices.size] += 1
        indices = np.insert(indices, 0, 0)
        indices = np.append(indices, data.obs_types.values.size)
        print(indices)'''
        
        plot_values = None
        cmap = None
        max_value = None
        print('num times: ', np.unique(data.times.values).size)
        #colormap for QC values
        if self.val_type.get() == 'QC':
            cmap = plt.get_cmap('gist_ncar', 9)
            plot_values = data.qc_DART
            max_value = 8
        elif self.val_type.get() == 'Observation value':
            cmap = plt.get_cmap('jet', np.unique(data.values).size)
            plot_values = data.values.ravel()
            max_value = max(plot_values)
        #cmap = plt.get_cmap('gist_ncar', 9)
        #ecmap = plt.get_cmap('jet', 90)

        print(plot_values.shape)
        #plot each observation type separately to get differing edge colors and markers
        
        #ax.scatter(lons, lats, data.z, c = data.qc_DART,
        #           cmap = cmap, vmin = 0, vmax = 9, s = 35, alpha = 0.5)
        ax.scatter(lons, lats, data.z, c = plot_values, cmap = cmap, s = 35, alpha = 0.5)
        #label = self.plotter.obs_type_inverse.get(data.obs_types.values[start]))
        #marker = self.markers[i % len(self.markers)], transform = ccrs.PlateCarree())
        
        '''
        #legend positioning
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend( bbox_to_anchor = (1, 1),
                  fontsize = 7, framealpha = 0.25)'''

        #make color bar
        sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(0,max_value))
        sm._A = []
        cbar = plt.colorbar(sm, ax=ax, orientation = 'horizontal', pad = 0.05)
        #cbar.ax.get_xaxis().labelpad = 15

        #set up color bar title according to plot type
        if self.val_type.get() == 'QC':
            cbar.ax.set_xlabel('DART QC Value')
        elif self.val_type.get() == 'Observation value':
            #get observation name
            cbar.ax.set_xlabel(self.obs_menu.get(self.obs_menu.curselection()).split(" : ", 1)[1])

        
        #TODO: make fill colors in legend transparent to avoid confusion
        #leg = ax.get_legend()
        #for handle in leg.legendHandles:
        #    handle.set_fill_color('None')
        ax.set_aspect('auto')
        #fig.tight_layout()
        print(time.time()-a)

        s= ttk.Style()
        s.theme_use('clam')
        
root = Tk()       
widg = GUI_3D_obs_initial(root, 0, 0)
widg.plot_3D()
print('on to mainloop')
root.style = ttk.Style()
print(root.style.theme_use())
root.style.theme_use('clam')
print(root.style.theme_use())
root.mainloop()
