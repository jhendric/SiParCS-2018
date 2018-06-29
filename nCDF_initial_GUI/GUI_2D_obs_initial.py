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
import time
'''

Incorporates plot_2D_obs_initial.py into a GUI with slider bars
and drop down menus.

'''



class GUI_2D_obs_initial:

    
    
    def __init__(self, window, grid_col, grid_row):
        
        self.plotter = plot_2D_obs('../obs_series/obs_epoch_001.nc')
        self.original_data = self.plotter.data
        #print('initial test: ', np.unique(self.plotter.data.where(
        #    self.plotter.data.obs_types == 5, drop = True).z.values))
        #print('initial test: ')
        #x = self.plotter.data.where(self.plotter.data.obs_types == 5, drop = True)
        #print(x.where(abs(x.z - 460.0) < 1e-10, drop = True))
        #print(np.unique(self.plotter.obs_types.values).size)
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
                                                 
        
        #self.obs_type_names = StringVar(value = self.plotter.obs_type_dict.keys())
        
        #print('obs_type_names: ', self.obs_type_names.get())

        #GUI config

        #observation selection
        self.obs_frame = ttk.Frame(self.main_frame, padding = "2")
        self.obs_frame.grid(column = 2, row = 1, sticky = "N, S, E, W")
        ttk.Label(self.obs_frame, text = "Observation Type Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.obs_menu = Listbox(self.obs_frame, listvariable = self.obs_type_names,
                                height = 18, width = 40, selectmode = "extended", exportselection = False)
        self.obs_menu.grid(column = 1, row = 2, rowspan = 2, sticky = "N, S, E, W")
        
        #self.obs_menu.bind('<Return>', self.populate_levels)
        self.obs_menu.bind('<Return>', lambda event : self.populate('levels', self.level_menu, event))
        
        for i in range(len(self.obs_type_names.get())):
            self.obs_menu.selection_set(i)
        self.obs_menu.event_generate('<<ListboxSelect>>')
        #print('size of obs box: ', self.obs_menu.size())

        #current obs types
        #obs_indices = [val+1 for val in self.obs_menu.curselection()]

        #obs scrollbar
        self.obs_bar = ttk.Scrollbar(self.obs_frame, orient = VERTICAL, command = self.obs_menu.yview)
        self.obs_menu.configure(yscrollcommand = self.obs_bar.set)
        self.obs_bar.grid(column = 2, row = 2, rowspan = 2,  sticky = "N, S, E")
        
        #print('indices of ob types: ', obs_indices)
        #print(type(obs_indices))

        
        #self.data = self.plotter.filter_disjoint(self.plotter.data, ('obs_types', obs_indices))
        #self.levels = StringVar(value = np.unique(self.data.z.values))

        self.levels = StringVar()
        
        #level selection
        
        self.level_frame = ttk.Frame(self.main_frame, padding = "2")
        self.level_frame.grid(column = 2, row = 2, sticky = "N, S, E, W")
        ttk.Label(self.level_frame, text = "Observation Level Selection").grid(column = 1,
                                                                               row = 1, sticky = "E, W")
        self.level_menu = Listbox(self.level_frame, listvariable = self.levels,
                                  height = 18, width = 40, selectmode = "extended", exportselection = False)
        self.level_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")

        #self.level_menu.bind('<Return>', self.populate_qc)
        self.level_menu.bind('<Return>', lambda event : self.populate('qc', self.qc_menu, event))
        
        #level scrollbar
        self.level_bar = ttk.Scrollbar(self.level_frame, orient = VERTICAL, command = self.level_menu.yview)
        self.level_menu.configure(yscrollcommand = self.level_bar.set)
        self.level_bar.grid(column = 2, row = 2, rowspan = 2, sticky = "N, S, E")
        
        self.qc = StringVar()
        
        #qc selection
        
        self.qc_frame = ttk.Frame(self.main_frame, padding = "2")
        self.qc_frame.grid(column=2, row = 3, sticky = "N, S, E, W")
        ttk.Label(self.qc_frame, text = "DART QC Value Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.qc_menu = Listbox(self.qc_frame, listvariable = self.qc,
                               height = 8, width = 40, selectmode = "extended", exportselection = False)
        self.qc_menu.grid(column = 1, row = 2, sticky ="N, S, E, W")

        #for use in populating and clearing menus (in populate function)
        self.data_obs_types = 1
        self.data_levels = 2
        self.data_qc = 3
        '''self.data_dict = {
            'obs_types' : self.data_obs_types,
            'levels' : self.data_levels,
            'qc' : self.data_qc
        }'''
        self.data_request_dict = {
            'data_levels' : 'obs_types',
            'data_qc' : 'z'
        }
        self.menu_hierarchy = [self.obs_menu, self.level_menu, self.qc_menu]
        #self.data_hierarchy = [self.plotter.data, self.data_obs_types, self.data_levels, self.data_qc]
        self.data_hierarchy = ['original_data', 'data_levels', 'data_qc']
        
        #populate levels
        
        self.populate('levels', self.level_menu)
        self.level_menu.selection_set(1)
        self.level_menu.event_generate('<<ListboxSelect>>')

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
        self.qc_menu.bind('<Return>', self.plot_2D)
        for i in range(len(self.qc.get())):
            self.qc_menu.selection_set(i)
        self.qc_menu.event_generate('<<ListboxSelect>>')

        
        #qc scrollbar
        self.qc_bar = ttk.Scrollbar(self.qc_frame, orient = HORIZONTAL, command = self.qc_menu.xview)
        self.qc_menu.configure(xscrollcommand = self.qc_bar.set)
        self.qc_bar.grid(column = 1, row = 3, rowspan = 1, sticky ="N, S, E, W")
        
        
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
        
        if var == 'data_levels':
            indices = [self.plotter.obs_type_dict[self.obs_menu.get(val).split(" : ", 1)[1]]
                       for val in self.obs_menu.curselection()]

        elif var == 'data_qc':

            indices = [np.float64(self.level_menu.get(val)) for val in self.level_menu.curselection()]

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
        
    def plot_2D(self, event = None):
        a = time.time()
        #event arg is passed by menu events
        print('plotting')
        #print(self.plotter.obs_type_dict.values())
        print('currently selected ob types: ', self.obs_menu.curselection())


        #current levels
        #level_indices = [val + 1 for val in self.level_menu.curselection()]

        #levels = [np.float64(self.level_menu.get(val)) for val in self.level_menu.curselection()]
        #print([self.qc_menu.get(val).split(":", 1) for val in self.qc_menu.curselection()])
        qc = [np.int64(self.qc_menu.get(val).split(": ", 1)[1][0]) for val in self.qc_menu.curselection()]

        print(qc)
        #print('indices of levels :', level_indices)
        #print(type(level_indices))
        
        #print('levels: ', levels)

        #print('test line: ', self.data.where(self.data.z == levels[0], drop = True))

        #make figure and canvas to draw on
        fig = Figure(figsize = (12,8))
        ax = fig.add_axes([0.01, 0.01, 0.98, 0.98], projection = ccrs.PlateCarree())
        canvas = FigureCanvasTkAgg(fig, master = self.main_frame)
        canvas.get_tk_widget().grid(column = 1, row = 1, rowspan = 3, sticky = "N, S, E, W")
        self.main_frame.grid_columnconfigure(1, weight = 1)
        self.main_frame.grid_rowconfigure(1, weight = 1)

        #have to set up a separate toolbar frame because toolbar doesn't like gridding with others
        self.toolbar_frame = ttk.Frame(self.main_frame)
        self.toolbar = NavigationToolbar2TkAgg(canvas, self.toolbar_frame)
        self.toolbar_frame.grid(column = 1, row = 4, sticky = "N, S, E, W")
        #disable part of the coordinate display functionality, else everything flickers
        #ax.format_coord = lambda x, y: ''
        
        data = self.plotter.filter_test(self.data_qc, ('qc_DART', qc))
        
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
        ax.legend( bbox_to_anchor = (1, 1),
                  fontsize = 7, framealpha = 0.25)

        #make color bar
        sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(0,9))
        sm._A = []
        cbar = plt.colorbar(sm, ax=ax, orientation = 'horizontal', pad = 0.05)
        cbar.ax.set_xlabel('DART QC Value')

        #center colorbar ticks and labels
        labels = np.arange(0, 9, 1)
        loc = labels + 0.5
        cbar.set_ticks(loc)
        cbar.set_ticklabels(labels)
        
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
widg = GUI_2D_obs_initial(root, 0, 0)
widg.plot_2D()
print('on to mainloop')
root.style = ttk.Style()
print(root.style.theme_use())
root.style.theme_use('clam')
print(root.style.theme_use())
root.mainloop()
