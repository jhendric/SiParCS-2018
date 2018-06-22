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
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import math
from plot_2D_obs_initial import plot_2D_obs
np.set_printoptions(threshold=np.nan) #without this setting, self.levels will be incomplete
'''

Incorporates plot_2D_obs_initial.py into a GUI with slider bars
and drop down menus.

'''



class GUI_2D_obs_initial:

    
    
    def __init__(self, window, grid_col, grid_row):
        
        self.plotter = plot_2D_obs('../obs_series/obs_epoch_001.nc')
        #print('initial test: ', np.unique(self.plotter.data.where(
        #    self.plotter.data.obs_types == 5, drop = True).z.values))
        #print('initial test: ')
        #x = self.plotter.data.where(self.plotter.data.obs_types == 5, drop = True)
        #print(x.where(abs(x.z - 460.0) < 1e-10, drop = True))
        print(np.unique(self.plotter.data.obs_types.values))
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

        #need: dict of obs_type names to obs_type values, list of obs_type names


        #strip useless characters from string interpretation of obs types
        self.obs_type_names = StringVar(value = [x.replace('[', '').replace(']', '').
                                                 replace(',', '').replace('\'', '').
                                                 replace('dict_keys(', '')
                                                 for x in self.plotter.obs_type_dict.keys()])
                                                 
        
        #self.obs_type_names = StringVar(value = self.plotter.obs_type_dict.keys())
        
        #print('obs_type_names: ', self.obs_type_names.get())

        #GUI config

        #observation selection
        self.obs_frame = ttk.Frame(self.main_frame, padding = "2")
        self.obs_frame.grid(column = 2, row = 1, sticky = "N, S, E, W")
        ttk.Label(self.obs_frame, text = "Observation Type Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.obs_menu = Listbox(self.obs_frame, listvariable = self.obs_type_names,
                                height = 10, selectmode = "extended", exportselection = False)
        self.obs_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")
        self.obs_menu.bind('<Return>', self.populate_levels)
        for i in range(len(self.obs_type_names.get())):
            self.obs_menu.selection_set(i)
        self.obs_menu.event_generate('<<ListboxSelect>>')
        #print('size of obs box: ', self.obs_menu.size())

        #current obs types
        obs_indices = [val+1 for val in self.obs_menu.curselection()]

        #print('indices of ob types: ', obs_indices)
        #print(type(obs_indices))

        
        #self.data = self.plotter.filter_disjoint(self.plotter.data, ('obs_types', obs_indices))
        #self.levels = StringVar(value = np.unique(self.data.z.values))

        self.levels = StringVar()
        
        #level selection
        self.level_frame = ttk.Frame(self.main_frame, padding = "2")
        self.level_frame.grid(column = 2, row = 2, sticky = "N, S, E, W")
        ttk.Label(self.level_frame, text = "Observation Level Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.level_menu = Listbox(self.level_frame, listvariable = self.levels,
                                  height = 10, selectmode = "extended", exportselection = False)
        self.level_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")
        
        self.populate_levels()
        #current plotting occurs only with press of enter from level menu
        self.level_menu.bind('<Return>', self.plot_2D)
        self.level_menu.selection_set(1)
        self.level_menu.event_generate('<<ListboxSelect>>')

        #for plotting later
        self.markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h',
                           'H', 'x', 'D', 'd', '|', '_']

        
    def populate_levels(self, event = None):
        #event arg is passed by menu events
        print('populating levels')
        self.level_menu.delete('0', 'end')

        #get indices of currently selected observation types
        obs_indices = [self.plotter.obs_type_dict[self.obs_menu.get(val)]
                       for val in self.obs_menu.curselection()]

        print(obs_indices)
        #print('indices of ob types: ', obs_indices)
        #print(type(obs_indices))

        self.data = self.plotter.filter_test(self.plotter.data, ('obs_types', obs_indices))
        
        self.levels.set(value = np.unique(self.data.z.values))
        print(np.unique(self.data.z.values).size)
        
        if (self.level_menu.get(0) == '['):
            self.level_menu.delete('0')
            
        if (self.level_menu.get(0)[0] == '['):
            first = self.level_menu.get(0)[1:]
            self.level_menu.delete('0')
            self.level_menu.insert(0, first)
            
        if (self.level_menu.get('end')[-1] == ']'):
            last = self.level_menu.get('end')[:-1]
            self.level_menu.delete('end')
            self.level_menu.insert(END, last)
        
    def plot_2D(self, event = None):
        #event arg is passed by menu events
        print('plotting')
        #print(self.plotter.obs_type_dict.values())
        print('currently selected ob types: ', self.obs_menu.curselection())


        #current levels
        #level_indices = [val + 1 for val in self.level_menu.curselection()]

        levels = [np.float64(self.level_menu.get(val)) for val in self.level_menu.curselection()]
        
        #print('indices of levels :', level_indices)
        #print(type(level_indices))
        
        #print('levels: ', levels)

        #print('test line: ', self.data.where(self.data.z == levels[0], drop = True))

        #make figure and canvas to draw on
        fig = Figure(figsize = (9,6))
        ax = fig.add_axes([0.01, 0.01, 0.98, 0.98], projection = ccrs.PlateCarree())
        canvas = FigureCanvasTkAgg(fig, master = self.main_frame)
        canvas.get_tk_widget().grid(column = 1, row = 1, rowspan = 2, sticky = "N, S, E, W")
        self.main_frame.grid_columnconfigure(1, weight = 1)
        self.main_frame.grid_rowconfigure(1, weight = 1)

        #have to set up a separate toolbar frame because toolbar doesn't like gridding with others
        self.toolbar_frame = ttk.Frame(self.main_frame)
        self.toolbar = NavigationToolbar2TkAgg(canvas, self.toolbar_frame)
        self.toolbar_frame.grid(column = 1, row = 3, sticky = "N, S, E, W")
        #disable part of the coordinate display functionality, else everything flickers
        #ax.format_coord = lambda x, y: ''
        
        data = self.plotter.filter_test(self.data, ('z', levels))

        #print(data.obs_types.values)
        print(data.obs_types.values.size)
        print(np.unique(data.qc_DART.values))
        print(np.unique(data.obs_types.values))

        #get indices where obs_types change (array is sorted in filter_disjoint)
        indices = np.where(data.obs_types.values[:-1] != data.obs_types.values[1:])[0]
        indices[0:indices.size] += 1
        indices = np.insert(indices, 0, 0)
        indices = np.append(indices, data.obs_types.values.size)
        print(indices)
        
        ax.stock_img()
        ax.gridlines()
        ax.coastlines()
        
        #colormap for QC values
        cmap = plt.get_cmap('gist_ncar', 9)
        ecmap = plt.get_cmap('jet', 90)
        #image = ax.scatter(data.lons, data.lats, c = data.qc_DART.values, cmap = cmap, vmin = 0, vmax = 9,
        #           s = 100, marker = "+", transform = ccrs.PlateCarree())

        '''brute force solution test

        for obs_type in np.unique(data.obs_types.values):
            d = data.where(data.obs_types == obs_type)
            ax.scatter(d.lons, d.lats, c = d.qc_DART.values, cmap = cmap, vmin = 0, vmax = 9, s = 100,
                       marker = "+", label = obs_type, transform = ccrs.PlateCarree())'''

        #plot each observation type separately to get differing edge colors and markers
        '''
        for i in range(indices.size - 1):
            start = indices[i]
            end = indices[i+1]
            ax.scatter(data.lons[start:end], data.lats[start:end], c = data.qc_DART.values[start:end],
                       cmap = cmap, vmin = 0, vmax = 9, s = 50, edgecolors = ecmap(1-float(i/indices.size)),
                       label = self.plotter.obs_type_inverse.get(data.obs_types.values[start]),
                       marker = self.markers[i % len(self.markers)], transform = ccrs.PlateCarree())'''
        for i in range(indices.size - 1):
            start = indices[i]
            end = indices[i+1]
            ax.scatter(data.lons[start:end], data.lats[start:end], c = ecmap(1-float(i/indices.size)),
                       cmap = ecmap, vmin = 0, vmax = 9, s = 50, edgecolors = cmap(data.qc_DART.values),
                       label = self.plotter.obs_type_inverse.get(data.obs_types.values[start]),
                       marker = self.markers[i % len(self.markers)], transform = ccrs.PlateCarree())

        #legend positioning
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc = 'top left', bbox_to_anchor = (1, 1),
                  fontsize = 6, framealpha = 0.25)

        #make color bar
        sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(0,8))
        sm._A = []
        cbar = plt.colorbar(sm, ax=ax, orientation = 'horizontal', pad = 0.05)
        #cbar.ax.get_xaxis().labelpad = 15
        cbar.ax.set_xlabel('DART QC Value')

        
        #TODO: make fill colors in legend transparent to avoid confusion
        #leg = ax.get_legend()
        #for handle in leg.legendHandles:
        #    handle.set_fill_color('None')
        ax.set_aspect('auto')
        fig.tight_layout()


root = Tk()
style = ttk.Style()
#style.theme_use('classic')
widg = GUI_2D_obs_initial(root, 0, 0)
widg.plot_2D()
print('on to mainloop')
root.mainloop()
