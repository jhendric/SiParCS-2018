from tkinter import *
from tkinter import ttk

import xarray as xa

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import cartopy.crs as ccrs

import numpy as np
import pandas as pd
import math
from read_innov import GenerateInnov
np.set_printoptions(threshold = np.nan) #without this setting, self.levels will be incomplete
import time

'''

Generates innovation plots in a GUI using read_innov.py

'''

class GUIInnov:

    def __init__(self, window, grid_col, grid_row, initial, final):

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
        
        self.window = window
        self.window.grid_columnconfigure(0, weight = 1)
        self.window.grid_rowconfigure(0, weight = 1)

        #a mainframe
        self.main_frame = ttk.Frame(self.window, padding = "8")
        self.main_frame.grid(column = grid_col, row = grid_row, sticky = "N, S, E, W") 
        self.main_frame.grid_columnconfigure(0, weight = 1) #weights for whole grid
        self.main_frame.grid_rowconfigure(0, weight = 1) #weights for whole grid

        data_types = [name for name in self.final.data_vars]

        data_types_sparse = [name.replace('[', '').replace(']', '').
                            replace(',', '').replace('\'', '')
                             for name in data_types]

        self.data_type_names = StringVar(value = data_types_sparse)

        #GUI config

        #data type selection

        self.data_frame = ttk.Frame(self.main_frame, padding = "2")
        self.data_frame.grid(column = 2, row = 1, sticky = "N, S, E, W")
        ttk.Label(self.data_frame, text = "Observation Type Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.data_menu = Listbox(self.data_frame, listvariable = self.data_type_names,
                                height = 18, width = 40, exportselection = False)
        self.data_menu.grid(column = 1, row = 2, rowspan = 1, sticky = "N, S, E, W")
        
        self.data_menu.bind('<Return>', lambda event : self.populate('levels', self.level_menu, event))
        
        self.data_menu.selection_set(0)
        self.data_menu.event_generate('<<ListboxSelect>>')
        #print('size of data box: ', self.data_menu.size())

        #data type scrollbar
        self.data_bar = ttk.Scrollbar(self.data_frame, orient = VERTICAL, command = self.data_menu.yview)
        self.data_menu.configure(yscrollcommand = self.data_bar.set)
        self.data_bar.grid(column = 2, row = 2, rowspan = 2,  sticky = "N, S, E")

        #level selection

        self.levels = StringVar()

        self.level_frame = ttk.Frame(self.main_frame, padding = "2")
        self.level_frame.grid(column = 2, row = 2, sticky = "N, S, E, W")
        ttk.Label(self.level_frame, text = "Data Level Selection").grid(column = 1,
                                                                               row = 1, sticky = "E, W")
        self.level_menu = Listbox(self.level_frame, listvariable = self.levels,
                                  height = 18, width = 40, exportselection = False)
        self.level_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")
        self.populate('levels', self.level_menu)
        self.level_menu.selection_set(0)
        self.level_menu.bind('<Return>', self.plot)
        
        #level scrollbar
        self.level_bar = ttk.Scrollbar(self.level_frame, orient = VERTICAL, command = self.level_menu.yview)
        self.level_menu.configure(yscrollcommand = self.level_bar.set)
        self.level_bar.grid(column = 2, row = 2, rowspan = 2, sticky = "N, S, E")

        #plot button
        self.plot_button = ttk.Button(self.main_frame, text = "Plot", command = self.plot, padding = "2")
        self.plot_button.grid(column = 2, row = 5, sticky = "N, S, E, W")
        
        
    def populate(self, variable_name, menu, event = None):
        #event arg is passed by menu events, variable is the data variable to be manipulated, menu is the menu to change
        
        print('populating' + variable_name)
        
        #clear lower level menus

        self.level_menu.delete('0', 'end')
            
        #get currently selected values
        
        indices = None

        #used to dynamically access object variables
        
        var = 'data_' + variable_name

        #TODO: this differs from the 3D obs and should be more modularized at some point to match hierarchy idea

        if var == 'data_levels':
            
            self.initial_var = self.initial[self.data_menu.get(self.data_menu.curselection())]
            self.final_var = self.final[self.data_menu.get(self.data_menu.curselection())]
            print(self.final_var)                  
            #find what coordinates levels are in 

            self.level_type = None

            for dim in self.final_var.dims:
                #get dimension type
                print(dim)
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
            
        print(self.levels.get())

    def plot(self, event = None):
        #a = time.time()
        #event arg is passed by menu events
        print('plotting')
        print('currently selected data type: ', self.data_menu.get(self.data_menu.curselection()))
        print('currently selected level: ', self.level_menu.get(self.level_menu.curselection()))
        
        
        
        fig = Figure(figsize = (12, 8))
        ax = fig.add_axes([0.01, 0.01, 0.98, 0.98], projection = ccrs.PlateCarree())
        
        canvas = FigureCanvasTkAgg(fig, master = self.main_frame)
        canvas.get_tk_widget().grid(column = 1, row = 1, rowspan = 3, sticky = "N, S, E, W")
        self.main_frame.grid_columnconfigure(1, weight = 1)
        self.main_frame.grid_rowconfigure(1, weight = 1)

        #have to set up a separate toolbar frame because toolbar doesn't like gridding with others
        self.toolbar_frame = ttk.Frame(self.main_frame)
        self.toolbar = NavigationToolbar2TkAgg(canvas, self.toolbar_frame)
        self.toolbar_frame.grid(column = 1, row = 4, sticky = "N, S, E, W")
        
        var_name = self.data_menu.get(self.data_menu.curselection())
        if self.level_menu.get(self.level_menu.curselection()).split(':')[0] == 'N/A':
            #don't try and index by level if no levels actually exist for that variable type
            self.final_var = self.final[var_name]
            self.initial_var = self.initial[var_name]
            self.innov = self.final_var - self.initial_var
            self.innov.attrs = self.final_var.attrs
            if 'long_name' in self.innov.attrs:
                self.innov.attrs['long_name'] = 'Difference in ' + self.innov.attrs['long_name']
            print(self.innov)
            print(np.mean(self.innov.values))
            self.innov.plot(ax = ax, transform = ccrs.PlateCarree(), cmap = 'bwr_r')
            #ax.set_ytitle('Difference between final and initial ' + var_name)
            ax.set_title('Innovation for ' + var_name)
        else:
            #get current level
            level = np.float64(float(self.level_menu.get(self.level_menu.curselection()).split(':')[0]))
            print('level: ', level)
            #narrow down data to correct data variable, correct level
            self.innov = self.gen.diff(var_name, level, self.level_type)
            print(self.innov)
            print(np.mean(self.innov.values))
            self.innov.plot(ax = ax, transform = ccrs.PlateCarree(), cmap = 'bwr_r')            
            ax.set_title('Innovation for ' + var_name +  ' @ ' + str(level))
            #ax.set_ytitle('Difference between final and initial ' + var_name + ' @ ' + str(level))
        plt.title = 'test'
        
        ax.gridlines(draw_labels = True)
        ax.coastlines()

        #ax.set_title('test')



root = Tk()
widg = GUIInnov(root, 0, 0, 'cam_input.0001.nc', 'cam_output.0001.nc')
#widg.plot()
print('on to mainloop')
root.mainloop()
                                   
        
