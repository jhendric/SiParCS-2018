'''#!/Users/wbd1/anaconda3/bin/python3'''

from tkinter import *
from tkinter import ttk

import xarray as xa

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy.crs as ccrs

import numpy as np

from plot_2D_obs import plot_2D_obs
np.set_printoptions(threshold=np.nan) #without this setting, self.levels will be incomplete


class GUI2DObs:
    '''
    
    Incorporates plot_2D_obs_initial.py into a GUI for plotting observation QC values in 2D.
    
    '''    
    
    def __init__(self, window, grid_col, grid_row, obs_sequence):

        '''Initialize GUI for plotting observation QC values in 2D
        
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
        self.main_frame.grid_columnconfigure(1, weight = 20) 
        self.main_frame.grid_columnconfigure(2, weight = 1)
        self.main_frame.grid_rowconfigure(0, weight = 1)
        self.main_frame.grid_rowconfigure(1, weight = 1) 
        self.main_frame.grid_rowconfigure(2, weight = 1)
        self.main_frame.grid_rowconfigure(3, weight = 1)
        self.main_frame.grid_rowconfigure(4, weight = 1)
        
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
                                selectmode = "extended", exportselection = False)
        self.obs_menu.grid(column = 1, row = 2, rowspan = 2, sticky = "N, S, E, W")

        self.obs_menu.bind('<Return>', lambda event : self.populate('levels', self.level_menu, event))
        
        for i in range(len(self.obs_type_names.get())):
            self.obs_menu.selection_set(i)
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

        
        self.levels = StringVar()
        
        #level selection
        
        self.level_frame = ttk.Frame(self.main_frame, padding = "2")
        self.level_frame.grid(column = 2, row = 2, sticky = "N, S, E, W")
        ttk.Label(self.level_frame, text = "Observation Level Selection").grid(column = 1,
                                                                               row = 1, sticky = "E, W")
        self.level_menu = Listbox(self.level_frame, listvariable = self.levels, #height = 18, width = 40,
                                  selectmode = "extended", exportselection = False)
        self.level_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")

        self.level_menu.bind('<Return>', lambda event : self.populate('qc', self.qc_menu, event))
        
        #level scrollbar
        self.level_bar = ttk.Scrollbar(self.level_frame, orient = VERTICAL, command = self.level_menu.yview)
        self.level_menu.configure(yscrollcommand = self.level_bar.set)
        self.level_bar.grid(column = 2, row = 2, rowspan = 2, sticky = "N, S, W")

        #resizing
        self.level_frame.grid_rowconfigure(1, weight = 1)
        self.level_frame.grid_rowconfigure(2, weight = 1)
        self.level_frame.grid_columnconfigure(1, weight = 1)
        self.level_frame.grid_columnconfigure(2, weight = 1)
        
        self.qc = StringVar()
        
        #qc selection
        
        self.qc_frame = ttk.Frame(self.main_frame, padding = "2")
        self.qc_frame.grid(column=2, row = 3, sticky = "N, S, E, W")
        ttk.Label(self.qc_frame, text = "DART QC Value Selection").grid (column = 1, row = 1, sticky = "E, W")
        self.qc_menu = Listbox(self.qc_frame, listvariable = self.qc, #height = 8, width = 40,
                               selectmode = "extended", exportselection = False)
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
        self.qc_bar.grid(column = 1, row = 3, rowspan = 1, sticky = "N, S, E, W")

        #resizing
        self.qc_frame.grid_rowconfigure(1, weight = 1)
        self.qc_frame.grid_rowconfigure(2, weight = 1)
        self.qc_frame.grid_columnconfigure(1, weight = 1)
        self.qc_frame.grid_columnconfigure(2, weight = 1)
        
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
        
        if var == 'data_levels':
            indices = [self.plotter.obs_type_dict[self.obs_menu.get(val).split(" : ", 1)[1]]
                       for val in self.obs_menu.curselection()]

        elif var == 'data_qc':

            indices = [np.float64(self.level_menu.get(val)) for val in self.level_menu.curselection()]

        #retrieve relevant data for this level of the hierarchy
        setattr(self, var,
                self.plotter.filter(getattr(self, self.data_hierarchy[self.data_hierarchy[1:].index(var)]),
                                         (self.data_request_dict[var], indices)))

        #set corresponding menu variables
        if var == 'data_levels':
            self.levels.set(value = np.unique(getattr(self, var).z.values))

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
        
    def plot_2D(self, event = None):

        '''Plot observation QC values on a global 2D map

        Keyword arguments:
        event -- an argument passed by any tkiner menu event. Has no influence on output but
        tkinter requires it to be passed to any method called by a menu event

        '''
        
        qc = [np.int64(self.qc_menu.get(val).split(": ", 1)[1][0]) for val in self.qc_menu.curselection()]

        #make figure and canvas to draw on
        fig = Figure(figsize = (12,8))
        ax = fig.add_axes([0.01, 0.01, 0.98, 0.98], projection = ccrs.PlateCarree())
        canvas = FigureCanvasTkAgg(fig, master = self.main_frame)
        canvas.get_tk_widget().grid(column = 1, row = 1, rowspan = 2, sticky = "N, S, E, W")

        #have to set up a separate toolbar frame because toolbar doesn't like gridding with others
        self.toolbar_frame = ttk.Frame(self.main_frame)
        self.toolbar = NavigationToolbar2TkAgg(canvas, self.toolbar_frame)
        self.toolbar.grid(column = 1, row = 1, sticky = "S, E, W")
        self.toolbar_frame.grid(column = 1, row = 3, sticky = "S, E, W")

        #resizing
        self.toolbar_frame.grid_columnconfigure(1, weight = 1)
        self.toolbar_frame.grid_rowconfigure(1, weight = 1)
        
        #disable part of the coordinate display functionality, else everything flickers
        #ax.format_coord = lambda x, y: ''
        
        data = self.plotter.filter(self.data_qc, ('qc_DART', qc))
        

        #get indices where obs_types change (array is sorted in filter_disjoint)
        indices = np.where(data.obs_types.values[:-1] != data.obs_types.values[1:])[0]
        indices[0:indices.size] += 1
        indices = np.insert(indices, 0, 0)
        indices = np.append(indices, data.obs_types.values.size)
        
        ax.stock_img()
        ax.gridlines()
        ax.coastlines()
        
        #colormap for QC values
        cmap = plt.get_cmap('gist_ncar', 9)
        ecmap = plt.get_cmap('jet', 90)

        
        #plot each observation type separately to get differing edge colors and markers
        
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

        ax.set_aspect('auto')

        s= ttk.Style()
        s.theme_use('clam')


def main(obs_sequence):

    '''create a tkinter GUI for plotting observation QC values in 2D

    Keyword arguments:
    obs_sequence -- path to a DART obs sequence file

    '''
    
    root = Tk()
    root.title("2D Observation Plotter")
    widg = GUI2DObs(root, 0, 0, obs_sequence)
    #widg.plot_2D()
    root.style = ttk.Style()
    root.style.theme_use('clam')
    root.mainloop()
        
if __name__ == '__main__':
    #only cmd line argument is obs sequence file name
    main(sys.argv[1])


