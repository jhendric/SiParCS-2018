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


class GUIObsDiagInitial:

    def __init__(self, window, grid_col, grid_row, diag):

        self.reader = ReadDiag(diag)
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

        self.forecast = None
        self.analysis = None
        self.level_type = None
        
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
        
        self.obs_menu.bind('<Return>', lambda event : self.populate('levels', self.level_menu, event))
        
        self.obs_menu.selection_set(0)
        self.obs_menu.event_generate('<<ListboxSelect>>')
        #print('size of obs box: ', self.obs_menu.size())

        #obs scrollbar
        self.obs_bar = ttk.Scrollbar(self.obs_frame, orient = VERTICAL, command = self.obs_menu.yview)
        self.obs_menu.configure(yscrollcommand = self.obs_bar.set)
        self.obs_bar.grid(column = 2, row = 2, rowspan = 2,  sticky = "N, S, E")

        #level selection

        self.levels = StringVar()

        self.level_frame = ttk.Frame(self.main_frame, padding = "2")
        self.level_frame.grid(column = 2, row = 2, sticky = "N, S, E, W")
        ttk.Label(self.level_frame, text = "Observation Level Selection").grid(column = 1,
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
        
        #for i in range(self.menu_hierarchy.index(menu), len(self.menu_hierarchy)):
        #    self.menu_hierarchy[i].delete('0', 'end')

        self.level_menu.delete('0', 'end')
            
        #get currently selected values
        
        indices = None

        #used to dynamically access object variables
        
        var = 'data_' + variable_name

        #TODO: this differs from the 3D obs and should be more modularized at some point to match hierarchy idea

        if var == 'data_levels':
            
            self.forecast = self.reader.get_variable(self.obs_menu.get(self.obs_menu.curselection()),
                                            'time series', 'forecast', self.original_data)
            self.analysis = self.reader.get_variable(self.obs_menu.get(self.obs_menu.curselection()),
                                                     'time series', 'analysis', self.original_data)
            #find whether level types are plevel, hlevel, or surface

            self.level_type = None

            for coord in self.forecast.coords:
                #get coordinate type
                print(coord)
                if coord.lower() in ('plevel', 'hlevel', 'surface'):
                    self.level_type = coord
                    break

            #surface may not show up in coords because it has no dimensions, so this is a workaround
            if self.level_type == None:
                self.level_type = 'surface'
            #self.levels.set(value = np.unique(self.forecast[self.level_type].values))
            for level in np.unique(self.forecast[self.level_type].values):
                #add warning for levels that are filled entirely with nan's
                forecast = self.reader.filter_single(self.forecast, (self.level_type, level))
                nanmean = np.nanmean(self.reader.filter_single(forecast, ('copy', 'rmse')))
                #a nan value will never equal any value, so if the mean outputs nan this will be false
                print(level, nanmean, nanmean == nanmean)
                if nanmean != nanmean:
                    menu.insert(END, str(level) + ': No data found for this level, will not plot')
                else:
                    menu.insert(END, level)
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
        print('currently selected ob type: ', self.obs_menu.curselection())
        print(self.level_menu.get(self.level_menu.curselection()))
        level = np.float64(float(self.level_menu.get(self.level_menu.curselection()).split(':')[0]))
        obs_type = self.obs_menu.get(self.obs_menu.curselection())
        print('level: ', level)
        print('obs type: ', obs_type)
        print('level type: ', self.level_type)
        
        print('forecast before level filtering: ', self.forecast)
        print('analysis before level filtering: ', self.analysis)
        #print('forecast test: ', self.forecast.where(self.forecast.plevel == 687.5))
        forecast = self.reader.filter_single(self.forecast, (self.level_type, level))
        analysis = self.reader.filter_single(self.analysis, (self.level_type, level))
        
        print('forecast: ', forecast)
        print('analysis: ', analysis)

        possible_obs = self.reader.filter_single(forecast, ('copy', 'Nposs'))
        used_obs = self.reader.filter_single(forecast, ('copy', 'Nused'))

        #only need rmse from forecast and analysis

        forecast = self.reader.filter_single(forecast, ('copy', 'rmse'))
        analysis = self.reader.filter_single(analysis, ('copy', 'rmse'))
        print('forecast filtered to rmse: ', forecast)
        print('analysis filtered to rmse: ', analysis)

        
        #need to change this to tolerate different numbers of subplots (regions)
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize = (8, 8))
        canvas = FigureCanvasTkAgg(fig, master = self.main_frame)
        canvas.get_tk_widget().grid(column = 1, row = 1, rowspan = 3, sticky = "N, S, E, W")
        self.main_frame.grid_columnconfigure(1, weight = 1)
        self.main_frame.grid_rowconfigure(1, weight = 1)

        #have to set up a separate toolbar frame because toolbar doesn't like gridding with others
        self.toolbar_frame = ttk.Frame(self.main_frame)
        self.toolbar = NavigationToolbar2TkAgg(canvas, self.toolbar_frame)
        self.toolbar_frame.grid(column = 1, row = 4, sticky = "N, S, E, W")

        print('region 1 mean: ', np.nanmean(forecast.where(forecast.region == 1, drop = True)))
        print('region 2 mean: ', np.nanmean(forecast.where(forecast.region == 2, drop = True)))
        print('region 3 mean: ', np.nanmean(forecast.where(forecast.region == 3, drop = True)))
        
        for index, ax in enumerate((ax1, ax2, ax3)):

            forecast_region = forecast.where(forecast.region == index + 1, drop = True)
            analysis_region = analysis.where(analysis.region == index + 1, drop = True)
            possible_obs_region = possible_obs.where(possible_obs.region == index + 1, drop = True)
            used_obs_region = used_obs.where(used_obs.region == index + 1, drop = True)
            
            #print('forecast region: ', forecast_region)
            #print('analysis region: ', analysis_region)
            #get rid of nan values by getting masks of only valid values, then indexing into them during plotting
            print('forecast region with nans: ', forecast_region.values)
            forecast_mask = np.array(list(filter(lambda v: v == v, forecast_region.values)))
            forecast_mask = ~np.isnan(forecast_region.values)
            #print(forecast_no_nans.flatten())
            forecast_no_nans = forecast_region.values[forecast_mask]
            print(forecast_mask.size, forecast_region.time.values.size)
            forecast_times_no_nans = forecast_region.time.values[forecast_mask.flatten()]
            analysis_mask = ~np.isnan(analysis_region.values)
            analysis_no_nans = analysis_region.values[analysis_mask]
            analysis_times_no_nans = analysis_region.time.values[analysis_mask.flatten()]
            #end = forecast_no_nans.size
            #forecast_times_no_nans = forecast_region.time.values[:end]
            #analysis_no_nans = np.array(list(filter(lambda v: v == v, analysis_region.values)))
            #analysis_times_no_nans = analysis_region.time.values[:end]
            print(forecast_no_nans.size)
            print('forecast cleaned of nans: ', forecast_no_nans)
            print('analysis cleaned of nans: ', analysis_no_nans)
            #do not plot regions with no values
            if forecast_no_nans.size == 0:
                ax.text(0.5, 0.5, 'No valid rmse data in this region')
                ax.set_title(str(self.reader.bytes_to_string(self.original_data['region_names'].values[index])) + '     ' +
                     'forecast: mean = ' + str(np.nanmean(forecast_region.values.flatten())) + '     ' +
                     'analysis: mean = ' + str(np.nanmean(analysis_region.values.flatten())))
                continue
            
            #plot both scatter and line for forecast and analysis to achieve connected appearance
            ax.scatter(x = forecast_times_no_nans, y = forecast_no_nans,
                       edgecolors = 'black', marker = 'x', s = 15)
            ax.plot(forecast_times_no_nans,
                    forecast_no_nans.flatten(), 'kx-', label = 'forecast')
            ax.scatter(x = analysis_times_no_nans,
                       y = analysis_no_nans,
                       edgecolors = 'red', marker = 'o', s = 15, facecolors = 'none')
            ax.plot(analysis_times_no_nans,
                    analysis_no_nans.flatten(), 'ro-', mfc = 'none', label = 'analysis')

            #set min and max x limits to be wider than actual limits for a nicer plot
            #pad x axis by 10% on both sides
            pad_x_axis = .10 * (max(forecast_region.time.values) - min(forecast_region.time.values))
            ax.set_xlim(((forecast_times_no_nans[0] - pad_x_axis).astype("M8[ms]")),
                        (forecast_times_no_nans[-1] + pad_x_axis).astype("M8[ms]"))
            #ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%t'))
            #ax.xaxis.set_minor_formatter(mdates.DateFormatter('%m/%d/%t'))
            
            print(ax.get_xticks())
            #pad y axis by 20% at top
            y_max = max(np.nanmax(forecast_region.values.flatten()), np.nanmax(analysis_region.values.flatten()))
            y_max = y_max + .20 * y_max
            ax.set_ylim(0, y_max)

            #add horizontal and vertical lines
            for i in range(1, int(y_max)):
                ax.axhline(y = i, ls = ':')

            for time in ax.get_xticks():
                ax.axvline(x = time, ls = ':')

            
            #subplot title
            #print('forecast: ', forecast_region.values.flatten())
            #print('forecast mean: ', np.nanmean(forecast_region.values.flatten()))
            #print('analysis: ', analysis_region.values.flatten())
            #print('analysis mean: ', np.nanmean(analysis_region.values.flatten()))            
            ax.set_title(str(self.reader.bytes_to_string(self.original_data['region_names'].values[index])) + '     ' +
                     'forecast: mean = ' + str(round(np.nanmean(forecast_region.values.flatten()), 5)) + '     ' +
                         'analysis: mean = ' + str(round(np.nanmean(analysis_region.values.flatten()), 5)))
            
            
            ax.set_ylabel(str(self.reader.bytes_to_string(self.original_data['region_names'].values[index])) + '\n' + 'rmse')
            ax.legend(loc = 'upper left', framealpha = 0.25)
            
            
            #need to basically plot two plots on top of each other to get 2 y scales
            ax_twin = ax.twinx()
            ax_twin.scatter(x = possible_obs_region.time.values, y = possible_obs_region.values,
                            color = 'blue', marker = 'o', s = 15, facecolors = 'none')
            ax_twin.scatter(x = used_obs_region.time.values, y = used_obs_region.values,
                            color = 'blue', marker = 'x', s = 15)
            '''
            print('possible obs: ', possible_obs_region.values.flatten())
            print('used obs: ', used_obs_region.values.flatten())
            print('difference in obs: ', possible_obs_region.values.flatten()-used_obs_region.values.flatten())'''

            y_max = max(max(possible_obs_region.values), max(used_obs_region.values))
            y_max = y_max + .20 * y_max
            
            ax_twin.set_ylim(0, y_max)
            ax_twin.set_ylabel('# of obs: o = poss, x = used', color = 'blue')
            #ax.fmt_xdata = mdates.DateFormatter('%m/%d')
            #ax.autofmt_xdate()
            for tick in ax.get_xticklabels():
                tick.set_rotation(45)
            ax.tick_params(labelsize = 8)
            
        #need to add plot title and spacing adjustments
        fig.suptitle(str(obs_type) + ' @ ' + str(level))
        fig.subplots_adjust(hspace = 0.8)


def main(diag):
    
    root = Tk()
    widg = GUIObsDiagInitial(root, 0, 0, diag)
    widg.plot()
    print('on to mainloop')
    root.mainloop()

if __name__ == '__main__':
    #only cmd line argument is obs diag file name
    main(sys.argv[1])
        
