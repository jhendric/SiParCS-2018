from tkinter import *
from tkinter import ttk

import xarray as xa

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import numpy as np
import pandas as pd
import math
from read_diag import ReadDiag
np.set_printoptions(threshold = np.nan) #without this setting, self.levels will be incomplete
import time
import itertools
from datetime import datetime

'''

Uses read_diag.py for plotting obs diag info in a GUI

'''


class GUIVertProf:

    def __init__(self, window, grid_col, grid_row):

        self.reader = ReadDiag('obs_diag_output.nc')
        self.original_data = self.reader.full_data

        self.window = window
        self.window.grid_columnconfigure(0, weight = 1)
        self.window.grid_rowconfigure(0, weight = 1)

        #a mainframe
        self.main_frame = ttk.Frame(self.window, padding = "8")
        self.main_frame.grid(column = grid_col, row = grid_row, sticky = "N, S, E, W") 
        self.main_frame.grid_columnconfigure(0, weight = 1) #weights for whole grid
        self.main_frame.grid_rowconfigure(0, weight = 1) #weights for whole grid

        #get observation type names
        #'comment' marks the beginning of the list of observation type names
        start_index = list(self.original_data.attrs).index('comment') + 1

        obs_types = [name for name in list(self.original_data.attrs)[start_index:]]

        #strip potential useless characters from string interpretation of obs types

        obs_types_sparse = [name.replace('[', '').replace(']', '').
                            replace(',', '').replace('\'', '')
                            for name in obs_types]

        self.obs_type_names = StringVar(value = obs_types_sparse)

        regions = [self.reader.bytes_to_string(name) for name in self.original_data.region_names.values]

        regions_sparse = [name.replace('[', '').replace(']', '').
                            replace(',', '').replace('\'', '')
                            for name in regions]

        self.region_names = StringVar(value = regions_sparse)
        
        self.forecast = None
        self.analysis = None
        
        '''
        #dictionary of obs types to their indices (shouldn't need this since should just need the names)
        self.obs_type_dict = {obs_type : index for obs_type, index in
                              zip(obs_types_sparse, list(self.original_data.attrs.values())[start_index])}'''
        

        #GUI config

        #observation selection

        self.obs_frame = ttk.Frame(self.main_frame, padding = "2")
        self.obs_frame.grid(column = 2, row = 1, sticky = "N, S, E, W")
        ttk.Label(self.obs_frame, text = "Observation Type Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.obs_menu = Listbox(self.obs_frame, listvariable = self.obs_type_names,
                                height = 18, width = 40, exportselection = False)
        self.obs_menu.grid(column = 1, row = 2, rowspan = 1, sticky = "N, S, E, W")
        
        self.obs_menu.bind('<Return>', self.plot)
        
        self.obs_menu.selection_set(0)
        self.obs_menu.event_generate('<<ListboxSelect>>')
        #print('size of obs box: ', self.obs_menu.size())

        #obs scrollbar
        self.obs_bar = ttk.Scrollbar(self.obs_frame, orient = VERTICAL, command = self.obs_menu.yview)
        self.obs_menu.configure(yscrollcommand = self.obs_bar.set)
        self.obs_bar.grid(column = 2, row = 2, rowspan = 2,  sticky = "N, S, E")

        #region selection

        self.region_frame = ttk.Frame(self.main_frame, padding = "2")
        self.region_frame.grid(column = 2, row = 2, sticky = "N, S, E, W")
        ttk.Label(self.region_frame, text = "Observation Region Selection").grid(column = 1,
                                                                               row = 1, sticky = "E, W")
        self.region_menu = Listbox(self.region_frame, listvariable = self.region_names,
                                  height = 18, width = 40, exportselection = False)
        self.region_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")
        self.region_menu.selection_set(0)
        self.region_menu.bind('<Return>', self.plot)
        
        #region scrollbar
        self.region_bar = ttk.Scrollbar(self.region_frame, orient = VERTICAL, command = self.region_menu.yview)
        self.region_menu.configure(yscrollcommand = self.region_bar.set)
        self.region_bar.grid(column = 2, row = 2, rowspan = 2, sticky = "N, S, E")
        
        #plot button
        self.plot_button = ttk.Button(self.main_frame, text = "Plot", command = self.plot, padding = "2")
        self.plot_button.grid(column = 2, row = 5, sticky = "N, S, E, W")
        
    def plot(self, event = None):
        #a = time.time()
        #event arg is passed by menu events
        print('plotting')
        print('currently selected ob type: ', self.obs_menu.curselection())
        obs_type = self.obs_menu.get(self.obs_menu.curselection())
        print('obs type: ', obs_type)

        forecast = self.reader.get_variable(self.obs_menu.get(self.obs_menu.curselection()),
                                            'vertical profile', 'forecast', self.original_data)
        analysis = self.reader.get_variable(self.obs_menu.get(self.obs_menu.curselection()),
                                            'vertical profile', 'analysis', self.original_data)

        #narrow to one region
        forecast_region = self.reader.filter_single(forecast, ('region', self.region_menu.curselection()))
        analysis_region = self.reader.filter_single(analysis, ('region', self.region_menu.curselection()))
        
        print('forecast: ', forecast_region)
        print('analysis: ', analysis_region)

        possible_obs = self.reader.filter_single(forecast_region, ('copy', 'Nposs'))
        used_obs = self.reader.filter_single(forecast_region, ('copy', 'Nused'))

        #only need rmse from forecast and analysis

        forecast_region = self.reader.filter_single(forecast_region, ('copy', 'rmse'))
        analysis_region = self.reader.filter_single(analysis_region, ('copy', 'rmse'))
        print('forecast filtered to rmse: ', forecast_region)
        print('analysis filtered to rmse: ', analysis_region)

        #get level type
        level_type = None
        for coord in forecast_region.coords:
            if coord != 'copy' and coord != 'region':
                level_type = coord

        forecast_levels = forecast_region[level_type]
        analysis_levels = analysis_region[level_type]
        
        fig, ax = plt.subplots(1, 1, figsize = (8, 8))
        canvas = FigureCanvasTkAgg(fig, master = self.main_frame)
        canvas.get_tk_widget().grid(column = 1, row = 1, rowspan = 3, sticky = "N, S, E, W")
        self.main_frame.grid_columnconfigure(1, weight = 1)
        self.main_frame.grid_rowconfigure(1, weight = 1)

        #have to set up a separate toolbar frame because toolbar doesn't like gridding with others
        self.toolbar_frame = ttk.Frame(self.main_frame)
        self.toolbar = NavigationToolbar2TkAgg(canvas, self.toolbar_frame)
        self.toolbar_frame.grid(column = 1, row = 4, sticky = "N, S, E, W")
        
        
        #print('forecast region: ', forecast_region)
        #print('analysis region: ', analysis_region)
        #get rid of nan values by getting masks of only valid values, then indexing into them during plotting
        print('forecast region with nans: ', forecast_region.values)
        forecast_mask = np.array(list(filter(lambda v: v == v, forecast_region.values.flatten())))
        forecast_mask = ~np.isnan(forecast_region.values)
        #print(forecast_no_nans.flatten())
        forecast_no_nans = forecast_region.values[forecast_mask]
        #print(forecast_mask.size, forecast_region.time.values.size)
        forecast_levels_no_nans = forecast_levels.values[forecast_mask.flatten()]
        analysis_mask = ~np.isnan(analysis_region.values)
        analysis_no_nans = analysis_region.values[analysis_mask]
        analysis_levels_no_nans = analysis_levels.values[analysis_mask.flatten()]
        #end = forecast_no_nans.size
        #forecast_levels_no_nans = forecast_region.time.values[:end]
        #analysis_no_nans = np.array(list(filter(lambda v: v == v, analysis_region.values)))
        #analysis_levels_no_nans = analysis_region.time.values[:end]
        print(forecast_no_nans.size)
        print('forecast cleaned of nans: ', forecast_no_nans)
        print('analysis cleaned of nans: ', analysis_no_nans)
        #do not plot regions with no values
        if forecast_no_nans.size == 0:
            ax.text(0.5, 0.5, 'No valid rmse data in this region')
            ax.set_title(self.region_menu.get(self.region_menu.curselection()) + '     ' +
                 'forecast: mean = ' + str(np.nanmean(forecast_region.values.flatten())) + '     ' +
                 'analysis: mean = ' + str(np.nanmean(analysis_region.values.flatten())))
            return

        #plot both scatter and line for forecast and analysis to achieve connected appearance
        ax.scatter(y = forecast_levels_no_nans, x = forecast_no_nans,
                   edgecolors = 'black', marker = 'x', s = 15)
        ax.plot(forecast_no_nans.flatten(), forecast_levels_no_nans,
                'kx-', label = 'forecast')
        ax.scatter(y = analysis_levels_no_nans,
                   x = analysis_no_nans,
                   edgecolors = 'red', marker = 'o', s = 15, facecolors = 'none')
        ax.plot(analysis_no_nans.flatten(), analysis_levels_no_nans,
                 'ro-', mfc = 'none', label = 'analysis')

        #set min and max y limits to be wider than actual limits for a nicer plot
        #pad y axis by 10% on both sides
        pad_y_axis = .10 * (max(forecast_levels.values) - min(forecast_levels.values))
        ax.set_ylim((forecast_levels_no_nans[0] - pad_y_axis),
                    (forecast_levels_no_nans[-1] + pad_y_axis))
        #ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%t'))
        #ax.xaxis.set_minor_formatter(mdates.DateFormatter('%m/%d/%t'))

        print(ax.get_yticks())
        #pad x axis by 20% at right
        x_max = max(np.nanmax(forecast_region.values.flatten()), np.nanmax(analysis_region.values.flatten()))
        x_max = x_max + .20 * x_max
        ax.set_xlim(0, x_max)

        #add horizontal and vertical lines
        for i in range(1, int(x_max)):
            ax.axvline(x = i, ls = ':')

        for lev in ax.get_yticks():
            ax.axhline(y = lev, ls = ':')


        #subplot title
        #print('forecast: ', forecast_region.values.flatten())
        #print('forecast mean: ', np.nanmean(forecast_region.values.flatten()))
        #print('analysis: ', analysis_region.values.flatten())
        #print('analysis mean: ', np.nanmean(analysis_region.values.flatten()))            
        '''ax.set_title(self.region_menu.get(self.region_menu.curselection()) + '     ' +
                     'forecast: mean = ' + str(round(np.nanmean(forecast_region.values.flatten()), 5)) + '     ' +
                     'analysis: mean = ' + str(round(np.nanmean(analysis_region.values.flatten()), 5)))'''


        ax.set_xlabel(self.region_menu.get(self.region_menu.curselection()) + '\n' + 'rmse')
        ax.legend(loc = 'upper left', framealpha = 0.25)


        #need to basically plot two plots on top of each other to get 2 y scales
        ax_twin = ax.twiny()
        ax_twin.scatter(y = possible_obs[level_type].values, x = possible_obs.values,
                        color = 'blue', marker = 'o', s = 15, facecolors = 'none')
        ax_twin.scatter(y = used_obs[level_type].values, x = used_obs.values,
                        color = 'blue', marker = 'x', s = 15)
        '''
        print('possible obs: ', possible_obs_region.values.flatten())
        print('used obs: ', used_obs_region.values.flatten())
        print('difference in obs: ', possible_obs_region.values.flatten()-used_obs_region.values.flatten())'''

        x_max = max(max(possible_obs.values.flatten()), max(used_obs.values.flatten()))
        x_max = x_max + .20 * x_max

        ax_twin.set_xlim(0, x_max)
        ax_twin.set_xlabel('# of obs: o = poss, x = used', color = 'blue')
        #ax.fmt_xdata = mdates.DateFormatter('%m/%d')
        #ax.autofmt_xdate()
        '''for tick in ax.get_xticklabels():
            tick.set_rotation(45)
        ax.tick_params(labelsize = 8)'''
            
        #need to add plot title and spacing adjustments
        
        fig.suptitle(str(obs_type) + '\n' + self.region_menu.get(self.region_menu.curselection()) + '     ' +
                     'forecast: mean = ' + str(round(np.nanmean(forecast_region.values.flatten()), 5)) + '     ' +
                     'analysis: mean = ' + str(round(np.nanmean(analysis_region.values.flatten()), 5)))
        #fig.subplots_adjust(hspace = 0.8)
        
root = Tk()
widg = GUIVertProf(root, 0, 0)
#widg.plot()
print('on to mainloop')
root.mainloop()
        
