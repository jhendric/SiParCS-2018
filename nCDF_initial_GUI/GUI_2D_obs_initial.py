'''#!/Users/wbd1/anaconda3/bin/python3'''

from tkinter import *
from tkinter import ttk
import mkl
import xarray as xa
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import math
from plot_2D_obs_initial import plot_2D_obs

'''

Incorporates plot_2D_obs_initial.py into a GUI with slider bars
and drop down menus.

'''

class GUI_2D_obs_initial:

    def plot_2D(self):

        print(self.obs_menu.curselection())
        obs_indices = [self.plotter.obs_type_dict.get(val) for val in self.obs_menu.curselection()]

        #may be the wrong datatype

        print(obs_indices)
        print(type(obs_indices))

        level_indices = self.level_menu.curselection()

        print(level_indices)
        print(type(level_indices))
        
        
        fig = Figure(figsize = (6,6))
        #this ax line may not function correctly
        ax = plt.axes(projection = ccrs.PlateCarree())
        canvas = FigureCanvasTkAgg(fig, master = self.main_frame)
        canvas.get_tk_widget().grid(column = 1, row = 1, rowspan = 4, sticky = "N, S, E, W")
        self.main_frame.grid_columnconfigure(1, weight = 1)
        self.main_frame.grid_rowconfigure(1, weight = 1)

        data = self.plotter.filter_disjoint(('obs_types', obs_indices), ('z', level_indices))

        ax.stock_img()
        ax.gridlines()
        ax.coastlines()
        ax.scatter(data.lons, data.lats, c = data.qc_DART, s = 100, marker = "+", transform = ccrs.PlateCarree())
        ax.legend(loc = 'upper left')
        fig.tight_layout()
    
    def __init__(self, window, grid_col, grid_row):
        self.plotter = plot_2D_obs('../obs_series/obs_epoch_001.nc')
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

        print(list(self.plotter.obs_type_dict.keys()))
    
        self.obs_type_names = StringVar(value = list(self.plotter.obs_type_dict.keys()))

        print(self.obs_type_names)

        self.levels = StringVar(value = self.plotter.z.values)

        #GUI config

        #observation selection
        self.obs_frame = ttk.Frame(self.main_frame, padding = "2")
        self.obs_frame.grid(column = 2, row = 1, sticky = "N, S, E, W")
        ttk.Label(self.obs_frame, text = "Observation Type Selection").grid(column = 1, row = 1, sticky = "N, S, E, W")
        self.obs_menu = Listbox(self.obs_frame, listvariable = self.obs_type_names,
                                height = 10, selectmode = "extended")
        self.obs_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")
        self.obs_menu.bind('<Return>', self.plot_2D)
        self.obs_menu.selection_set(4)

        #level selection
        self.level_frame = ttk.Frame(self.main_frame, padding = "2")
        self.level_frame.grid(column = 3, row = 1, sticky = "N, S, E, W")
        ttk.Label(self.level_frame, text = "Observation Level Selection").grid(column = 1, row = 1, sticky = "N, S, E, W")
        self.level_menu = Listbox(self.level_frame, listvariable = self.levels,
                                  height = 10, selectmode = "extended")
        self.level_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")
        self.level_menu.bind('<Return>', self.plot_2D)
        self.level_menu.selection_set(1)        

        self.plot_2D


root = Tk()
style = ttk.Style()
#style.theme_use('classic')
widg = GUI_2D_obs_initial(root, 0, 0)
widg.plot_2D()
root.mainloop()
