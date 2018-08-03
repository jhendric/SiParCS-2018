'''#!/Users/wbd1/anaconda3/bin/python3'''

from tkinter import *
from tkinter import ttk

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
from plot_2D_obs import plot_2D_obs
np.set_printoptions(threshold=np.nan) #without this setting, self.levels will be incomplete
import itertools


class GUI3DObs:

    '''

    Incorporates plot_2D_obs_initial.py into a GUI for plotting observation values
    and QC values in 3D.

    '''
    
    def __init__(self, window, grid_col, grid_row, obs_sequence):

        '''Initialize GUI for plotting observation values and QC values in 3D

        Keyword arguments:
        window -- the root window holding all GUI elements
        grid_col -- the column in the root window that will contain the main tkinter frame
        grid_row -- the row in the root window that will contain the main tkinter frame
        obs_sequence -- path to a DART obs sequence file

        '''
        
        self.plotter = plot_2D_obs(obs_sequence)
        self.original_data = self.plotter.data

        self.window = window
        self.window.grid_columnconfigure(0, weight = 1)
        self.window.grid_rowconfigure(0, weight = 1)

        #a mainframe
        self.main_frame = ttk.Frame(self.window, padding = "8")
        self.main_frame.grid(column = grid_col, row = grid_row, sticky = "N, S, E, W")

        #resizing
        self.main_frame.grid_columnconfigure(0, weight = 1) #weights for whole grid
        self.main_frame.grid_rowconfigure(0, weight = 1) #weights for whole grid
        self.main_frame.grid_columnconfigure(1, weight = 20)
        self.main_frame.grid_columnconfigure(2, weight = 1)
        self.main_frame.grid_rowconfigure(1, weight = 1) #weights for whole grid
        self.main_frame.grid_rowconfigure(2, weight = 1)
        self.main_frame.grid_rowconfigure(3, weight = 1)
        self.main_frame.grid_rowconfigure(4, weight = 1)
        self.main_frame.grid_rowconfigure(5, weight = 1)
        
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
        self.obs_menu = Listbox(self.obs_frame, listvariable = self.obs_type_names, #height = 18, width = 40,
                                exportselection = False)
        self.obs_menu.grid(column = 1, row = 2, rowspan = 1, sticky = "N, S, E, W")
        
        self.obs_menu.bind('<Return>', lambda event : self.populate('times', self.times_menu, event))
        
        self.obs_menu.selection_set(0)
        self.obs_menu.event_generate('<<ListboxSelect>>')

        #obs scrollbar
        self.obs_bar = ttk.Scrollbar(self.obs_frame, orient = VERTICAL, command = self.obs_menu.yview)
        self.obs_menu.configure(yscrollcommand = self.obs_bar.set)
        self.obs_bar.grid(column = 2, row = 2, rowspan = 2,  sticky = "N, S, W")

        #resizing
        self.obs_frame.grid_columnconfigure(1, weight = 1)
        self.obs_frame.grid_columnconfigure(2, weight = 1)
        self.obs_frame.grid_rowconfigure(1, weight = 1)
        self.obs_frame.grid_rowconfigure(2, weight = 1)

        #time selection
        self.times = StringVar()
        
        self.times_frame = ttk.Frame(self.main_frame, padding = "2")
        self.times_frame.grid(column = 2, row = 2, sticky = "N, S, E, W")
        ttk.Label(self.times_frame, text = "Observation Times Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.times_menu = Listbox(self.times_frame, listvariable = self.times, #height = 18, width = 40,
                                  selectmode = "extended", exportselection = False)
        self.times_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")
        self.times_menu.bind('<Return>', lambda event : self.populate('qc', self.qc_menu, event))

        self.times_bar = ttk.Scrollbar(self.times_frame, orient = VERTICAL, command = self.times_menu.yview)
        self.times_menu.configure(yscrollcommand = self.times_bar.set)
        self.times_bar.grid(column = 2, row = 2, rowspan = 2, sticky = "N, S, W")

        #resizing
        self.times_frame.grid_columnconfigure(1, weight = 1)
        self.times_frame.grid_columnconfigure(2, weight = 1)
        self.times_frame.grid_rowconfigure(1, weight = 1)
        self.times_frame.grid_rowconfigure(2, weight = 1)
        
        #qc selection

        self.qc = StringVar()
        
        self.qc_frame = ttk.Frame(self.main_frame, padding = "2")
        self.qc_frame.grid(column=2, row = 3, sticky = "N, S, E, W")
        ttk.Label(self.qc_frame, text = "DART QC Value Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.qc_menu = Listbox(self.qc_frame, listvariable = self.qc, #height = 8, width = 40,
                               selectmode = "extended", exportselection = False)
        self.qc_menu.grid(column = 1, row = 2, sticky ="N, S, E, W")
                
        #qc scrollbar
        self.qc_bar = ttk.Scrollbar(self.qc_frame, orient = HORIZONTAL, command = self.qc_menu.xview)
        self.qc_menu.configure(xscrollcommand = self.qc_bar.set)
        self.qc_bar.grid(column = 1, row = 3, rowspan = 1, sticky ="N, S, E, W")

        #resizing
        self.qc_frame.grid_rowconfigure(1, weight = 1)
        self.qc_frame.grid_rowconfigure(2, weight = 1)
        self.qc_frame.grid_rowconfigure(3, weight = 1)
        self.qc_frame.grid_columnconfigure(1, weight = 1)
        self.qc_frame.grid_columnconfigure(2, weight = 1)
        
        #for use in populating and clearing menus (in populate function)
        self.data_obs_types = 1
        self.data_times = 2
        self.data_qc = 3
        
        self.data_request_dict = {
            'data_times' : 'obs_types',
            'data_qc' : 'times'
        }
        
        self.menu_hierarchy = [self.obs_menu, self.times_menu, self.qc_menu]
        self.data_hierarchy = ['original_data', 'data_times', 'data_qc']

        #populate times
        self.populate('times', self.times_menu)

        #plot all times by default
        for i in range(len(self.times.get())):
            self.times_menu.selection_set(i)
        self.times_menu.event_generate('<<ListboxSelect>>')
        
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

        #current plotting occurs with press of enter from qc menu or click of plot button
        self.qc_menu.bind('<Return>', self.plot_3D)
        for i in range(len(self.qc.get())):
            self.qc_menu.selection_set(i)
        self.qc_menu.event_generate('<<ListboxSelect>>')

        
        #selection of what value type to plot
        self.val_type = StringVar()
        self.val_type.set('QC')
        self.fill_frame = ttk.Frame(self.main_frame, padding = "2")
        self.fill_frame.grid(column = 2, row = 4, sticky = "N, S, E, W")
        ttk.Label(self.fill_frame, text = "Please select type of value to plot").grid(column = 1, row = 1, sticky= "E, W")
        self.qc_button = ttk.Radiobutton(self.fill_frame, text = 'QC', variable = self.val_type, value = 'QC')
        self.qc_button.grid(column = 1, row = 2, sticky = "N, S, E, W")
        self.val_button = ttk.Radiobutton(self.fill_frame, text = 'Observation value', variable = self.val_type,
                                          value = "Observation value")
        self.val_button.grid(column = 1, row = 3, sticky = "N, S, E, W")

        #resizing
        self.fill_frame.rowconfigure(1, weight = 1)
        self.fill_frame.rowconfigure(2, weight = 1)
        self.fill_frame.rowconfigure(3, weight = 1)
        self.fill_frame.columnconfigure(1, weight = 1)
        
        #plot button
        self.plot_button = ttk.Button(self.main_frame, text = "Plot", command = self.plot_3D, padding = "2")
        self.plot_button.grid(column = 2, row = 5, sticky = "N, S, E, W")

        
        
        #for plotting later
        self.markers = ['o', 'v', 'H', 'D', '^', '<', '8',
                        's', 'p', '>',  '*', 'h', 'd']

        #these markers do not seem to have border color capabilities
        #'x', '_', '|'

        s = ttk.Style()
        s.theme_use('clam')

    def populate(self, variable_name, menu, event = None):
        
        '''Populate levels, time, and QC menus based on which selections in a menu have
        been modified.
        
        Keyword arguments:
        variable_name -- data variable to be populated
        menu -- corresponding menu to change
        event -- argument passed automatically by any tkinter menu event
        
        '''
        
        #clear lower level menus
        for i in range(self.menu_hierarchy.index(menu), len(self.menu_hierarchy)):
            self.menu_hierarchy[i].delete('0', 'end')

        #get currently selected values
        
        indices = None

        #used to dynamically access object variables
        var = 'data_' + variable_name

        #TODO: this differs from the 2D obs and should be more modularized at some point to match hierarchy idea

        if var == 'data_times':
            indices = [self.plotter.obs_type_dict[self.obs_menu.get(val).split(" : ", 1)[1]]
                       for val in self.obs_menu.curselection()]
        
        elif var == 'data_qc':
            indices = [np.datetime64(self.times_menu.get(val).replace('\'', '')) for val in self.times_menu.curselection()]
            
        #retrieve relevant data for this level of the hierarchy
        setattr(self, var,
                self.plotter.filter(getattr(self, self.data_hierarchy[self.data_hierarchy[1:].index(var)]),
                                         (self.data_request_dict[var], indices)))

        #set corresponding menu variables
        if var == 'data_times':
            self.times.set(value = np.unique(getattr(self, var).times.values))
            
        elif var == 'data_qc':
            unique, counts = np.unique(getattr(self, var).qc_DART.values, return_counts = True)
            count_dict = dict(zip(unique, counts))
            self.qc.set(value = [str(count_dict[val]) + " : " + str(self.qc_key[val]) for val in unique])
            
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

        '''Plot observation QC values on a global 2D map

        Keyword arguments:
        event -- an argument passed by any tkiner menu event. Has no influence on output but
        tkinter requires it to be passed to any method called by a menu event

        '''

        qc = [np.int64(self.qc_menu.get(val).split(": ", 1)[1][0]) for val in self.qc_menu.curselection()]
        
        #make figure and canvas to draw on
        fig = Figure(figsize = (12, 8))
        #ax = fig.add_axes([0.01, 0.01, 0.98, 0.98], projection = ccrs.PlateCarree())
        canvas = FigureCanvasTkAgg(fig, master = self.main_frame)
        canvas.get_tk_widget().grid(column = 1, row = 1, rowspan = 4, sticky = "N, S, E, W")
        
        #have to set up a separate toolbar frame because toolbar doesn't like gridding with others
        self.toolbar_frame = ttk.Frame(self.main_frame)
        self.toolbar = NavigationToolbar2TkAgg(canvas, self.toolbar_frame)
        self.toolbar_frame.grid(column = 1, row = 5, sticky = "N, S, E, W")
        self.toolbar.grid(column = 1, row = 1, sticky = "N, S, E, W")
        #resizing
        self.toolbar_frame.grid_columnconfigure(1, weight = 1)
        self.toolbar_frame.grid_rowconfigure(1, weight = 1)
        
        #disable part of the coordinate display functionality, else everything flickers (may need for smaller window)
        #ax.format_coord = lambda x, y: ''
        
        data = self.plotter.filter(self.data_qc, ('qc_DART', qc))
        
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

        ax = Axes3D(fig, [0.02, 0.02, 0.97, 0.97],  xlim = [-180, 180], ylim = [-90, 90])
        
        z_max = max(data.z.values)

        #deal with orientation of graph for different vert types
        if data.vert_types.values[0] < 0 or data.vert_types.values[0] == 3:
            #surface, heights, or unknown units where up is positive direction
            #don't need to do 3D polygons since it will plot 2D polygons at z = 0 by default
            lc = PolyCollection(polys, edgecolor = 'black', facecolor = 'green', closed = False)
            ax.set_zlim(bottom = 0, top = z_max)
            ax.add_collection3d(lc)
        else:
            #level or pressure units where down is positive direction
            polys = [[(point[0], point[1], z_max) for point in shape] for shape in polys]
            lc = Poly3DCollection(polys, edgecolor = 'black', facecolor = 'green', closed = False)
            ax.set_zlim(bottom = z_max, top = 0)
            ax.add_collection3d(lc)
        
        plot_values = None
        cmap = None
        max_value = None
        min_value = None
        
        #colormap for QC values
        if self.val_type.get() == 'QC':
            cmap = plt.get_cmap('gist_ncar', 9)
            plot_values = data.qc_DART
            max_value = 8
            min_value = 0
        elif self.val_type.get() == 'Observation value':
            cmap = plt.get_cmap('jet', np.unique(data.values).size)
            plot_values = data.values.ravel()
            max_value = max(plot_values)
            min_value = min(plot_values)

        #plot each observation type separately to get differing edge colors and markers
        
        ax.scatter(lons, lats, data.z, c = plot_values, cmap = cmap, vmin = min_value, vmax = max_value, s = 35, alpha = 0.5)
        

        #make color bar
        sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(0,max_value+1))
        sm._A = []
        cbar = plt.colorbar(sm, ax=ax, orientation = 'horizontal', pad = 0.04)

        #set up color bar title according to plot type
        if self.val_type.get() == 'QC':
            
            cbar.ax.set_xlabel('DART QC Value')
            #also center colorbar ticks and labels
            labels = np.arange(0, max_value + 1, 1)
            loc = labels + 0.5
            cbar.set_ticks(loc)
            cbar.set_ticklabels(labels)
            
        elif self.val_type.get() == 'Observation value':
            
            #get observation name
            cbar.ax.set_xlabel(self.obs_menu.get(self.obs_menu.curselection()).split(" : ", 1)[1])
            
        ax.set_aspect('auto')
        #ax.margins(0.5, tight = False)
        s= ttk.Style()
        s.theme_use('clam')


def main(obs_sequence):

    '''create a tkinter GUI for plotting observation QC values in 2D

    Keyword arguments:
    obs_sequence -- path to a DART obs sequence file

    '''
    
    root = Tk()
    root.title("3D Observation Plotter")
    widg = GUI3DObs(root, 0, 0, obs_sequence)
    widg.plot_3D()
    root.style = ttk.Style()
    root.style.theme_use('clam')
    root.mainloop()
        
if __name__ == '__main__':
    #only cmd line argument is obs sequence file name
    main(sys.argv[1])

    
