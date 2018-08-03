from tkinter import *
from tkinter import ttk

import xarray as xa

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy.crs as ccrs

import numpy as np
from read_innov import GenerateInnov
np.set_printoptions(threshold = np.nan) #without this setting, self.levels will be incomplete


class GUIInnov:
    '''
    
    Generates innovation plots in a GUI using read_innov.py
    
    '''
    
    def __init__(self, window, grid_col, grid_row, initial, final):

        '''Initialize GUI with available observation types and levels
        
        Keyword arguments:
        window -- the root window holding all GUI elements
        grid_col -- the column in the root window that will contain the main tkinter frame
        grid_row -- the row in the root window that will contain the main tkinter frame
        initial -- path of pre-DART model restart file
        final -- path of post-DART model restart file
        
        '''
        
        self.gen = GenerateInnov(initial, final)

        self.initial = self.gen.initial
        self.final = self.gen.final

        #used to hold DataArray of a specific variable type
        self.initial_var = None
        self.final_var = None

        #used to hold innovation of a specific variable type at a specific level
        self.innov = None

        #No convention for z coordinate among different files. Whatever is found is stored here
        self.level_type = None

        #resizing
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

        data_types = [name for name in self.final.data_vars]

        data_types_sparse = [name.replace('[', '').replace(']', '').
                            replace(',', '').replace('\'', '')
                             for name in data_types]

        self.data_type_names = StringVar(value = data_types_sparse)

        #GUI config

        #data type selection

        self.data_frame = ttk.Frame(self.main_frame, padding = "2")
        self.data_frame.grid(column = 2, row = 1, sticky = "N, S, E, W")
        ttk.Label(self.data_frame, text = "State Variable Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.data_menu = Listbox(self.data_frame, listvariable = self.data_type_names,
                                 #height = 18, width = 40,
                                 exportselection = False)
        self.data_menu.grid(column = 1, row = 2, rowspan = 1, sticky = "N, S, E, W")
        
        self.data_menu.bind('<Return>', lambda event : self.populate('levels', self.level_menu, event))
        
        self.data_menu.selection_set(0)
        self.data_menu.event_generate('<<ListboxSelect>>')

        #data type scrollbar
        self.data_bar = ttk.Scrollbar(self.data_frame, orient = VERTICAL, command = self.data_menu.yview)
        self.data_menu.configure(yscrollcommand = self.data_bar.set)
        self.data_bar.grid(column = 2, row = 2,  sticky = "N, S, W")

        #resizing
        self.data_frame.grid_columnconfigure(1, weight = 1)
        self.data_frame.grid_columnconfigure(2, weight = 1)
        self.data_frame.grid_rowconfigure(1, weight = 1)
        self.data_frame.grid_rowconfigure(2, weight = 1)
        

        #level selection

        self.levels = StringVar()

        self.level_frame = ttk.Frame(self.main_frame, padding = "2")
        self.level_frame.grid(column = 2, row = 2, sticky = "N, S, E, W")
        ttk.Label(self.level_frame, text = "Data Level Selection").grid(column = 1,
                                                                               row = 1, sticky = "E, W")
        self.level_menu = Listbox(self.level_frame, listvariable = self.levels,
                                  #height = 18, width = 40,
                                  exportselection = False)
        self.level_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")
        self.populate('levels', self.level_menu)
        self.level_menu.selection_set(0)
        self.level_menu.bind('<Return>', self.plot)

        #resizing
        self.level_frame.grid_columnconfigure(1, weight = 1)
        self.level_frame.grid_columnconfigure(2, weight = 1)
        self.level_frame.grid_rowconfigure(1, weight = 1)
        self.level_frame.grid_rowconfigure(2, weight = 1)
        
        #level scrollbar
        self.level_bar = ttk.Scrollbar(self.level_frame, orient = VERTICAL, command = self.level_menu.yview)
        self.level_menu.configure(yscrollcommand = self.level_bar.set)
        self.level_bar.grid(column = 2, row = 2, rowspan = 2, sticky = "N, S, W")

        #plot button
        self.plot_button = ttk.Button(self.main_frame, text = "Plot", command = self.plot, padding = "2")
        self.plot_button.grid(column = 2, row = 3, sticky = "N, S, E, W")
        
        
    def populate(self, variable_name, menu, event = None):
        #event arg is passed by menu events, variable is the data variable to be manipulated, menu is the menu to change

        '''Populate levels menu based on observation type selected
        
        Keyword arguments:
        variable_name -- data variable to be populated (levels in this GUI)
        menu -- corresponding menu to change
        event -- argument passed automatically by any tkinter menu event
        
        '''
        
        #clear lower level menus

        self.level_menu.delete('0', 'end')
            
        #get currently selected values
        
        indices = None

        #used to dynamically access object variables
        
        var = 'data_' + variable_name

        if var == 'data_levels':
            
            self.initial_var = self.initial[self.data_menu.get(self.data_menu.curselection())]
            self.final_var = self.final[self.data_menu.get(self.data_menu.curselection())]
            
            #find what coordinates levels are in 

            self.level_type = None

            for dim in self.final_var.dims:
                #get dimension type
                if dim.lower() in ('plevel', 'hlevel', 'surface', 'level', 'lev',
                                     'height', 'elevation', 'elev', 'z', 'h', 'depth'):
                    self.level_type = dim
                    break
            if self.level_type != None:    
                self.levels.set(value = np.unique(self.final_var[self.level_type].values))
            else:
                self.levels.set('N/A')
        
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

    def plot(self, event = None):

        '''Plot difference in model state data for a selected observation type and level
        
        Keyword arguments:
        event -- an argument passed by any tkinter menu event. Has no influence on output but
        tkinter requires it to be passed to any method called by a menu event.

        '''
        
        fig = Figure(figsize = (12, 8))
        ax = fig.add_axes([0.03, 0.03, 0.95, 0.95], projection = ccrs.PlateCarree())
        #cax = fig.add_axes([0.955, 0.03, 0.99, 0.99])
        
        canvas = FigureCanvasTkAgg(fig, master = self.main_frame)
        canvas.get_tk_widget().grid(column = 1, row = 1, rowspan = 2, sticky = "N, S, E, W")

        #have to set up a separate toolbar frame because toolbar doesn't like gridding with others
        self.toolbar_frame = ttk.Frame(self.main_frame)
        self.toolbar = NavigationToolbar2TkAgg(canvas, self.toolbar_frame)
        self.toolbar_frame.grid(column = 1, row = 3, sticky = "N, S, E, W")
        self.toolbar.grid(column = 1, row = 1, sticky = "N, S ,E , W")

        #resizing
        self.toolbar_frame.grid_columnconfigure(1, weight = 1)
        self.toolbar_frame.grid_rowconfigure(1, weight = 1)
        
        
        var_name = self.data_menu.get(self.data_menu.curselection())
        if self.level_menu.get(self.level_menu.curselection()).split(':')[0] == 'N/A':
            #don't try and index by level if no levels actually exist for that variable type
            self.final_var = self.final[var_name]
            self.initial_var = self.initial[var_name]
            self.innov = self.final_var - self.initial_var
            self.innov.attrs = self.final_var.attrs

            if 'long_name' in self.innov.attrs:
                self.innov.attrs['long_name'] = 'Difference in ' + self.innov.attrs['long_name']

            #divider = make_axes_locatable(ax)
            #cax = divider.append_axes("right", size = "5%", pad = 0.05)
            test = self.innov.plot(ax = ax, transform = ccrs.PlateCarree(), cmap = 'bwr_r',
                                   cbar_kwargs = {'orientation' : 'vertical',
                                                  'fraction' : 0.1})
            
            ax.set_title('Innovation for ' + var_name + '\n \n')
        else:
            #get current level
            level = np.float64(float(self.level_menu.get(self.level_menu.curselection()).split(':')[0]))
            #narrow down data to correct data variable, correct level
            self.innov = self.gen.diff(var_name, level, self.level_type)
            self.innov.plot(ax = ax, transform = ccrs.PlateCarree(), cmap = 'bwr_r')
            ax.set_title('Innovation for ' + var_name +  ' @ ' + str(level) + '\n \n')

        ax.margins(0.5)    
        ax.gridlines(draw_labels = True)
        ax.coastlines()

def main(initial, final):

    '''create a tkinter GUI for plotting innovation data

    Keyword arguments:
    initial -- path of pre-DART model restart file
    final -- path of post-DART model restart file
    
    '''
    
    root = Tk()
    root.title("Innovation Plotter")
    widg = GUIInnov(root, 0, 0, initial, final)
    root.mainloop()
        
if __name__ == '__main__':
    #cmd line arguments are initial conditions file name and final conditions file name
    main(sys.argv[1], sys.argv[2])


                                   
        
