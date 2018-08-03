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
from read_diag import ReadDiag
np.set_printoptions(threshold = np.nan) #without this setting, self.levels will be incomplete


class GUIObsDiag:

    '''
    
    Uses read_diag.py for plotting evolution of RMSE in a GUI
    
    '''
    
    def __init__(self, window, grid_col, grid_row, diag):

        '''Initialize GUI for plotting evolution of RMSE
        
        Keyword arguments:
        window -- the root window holding all GUI elements
        grid_col -- the column in the root window that will contain the main tkinter frame
        grid_row -- the row in the root window that will contain the main tkinter frame
        diag -- an obs_diag output netCDF file
        
        '''
        
        self.reader = ReadDiag(diag)
        self.original_data = self.reader.full_data

        self.window = window
        
        #resizing
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
        self.main_frame.grid_rowconfigure(1, weight = 1)
        self.main_frame.grid_rowconfigure(2, weight = 1)
        self.main_frame.grid_rowconfigure(3, weight = 1)
        self.main_frame.grid_rowconfigure(4, weight = 1)
        self.main_frame.grid_rowconfigure(5, weight = 1)
        
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
        self.level_type = None

        #GUI config

        #observation selection

        self.obs_frame = ttk.Frame(self.main_frame, padding = "2")
        self.obs_frame.grid(column = 2, row = 1, sticky = "N, S, E, W")
        ttk.Label(self.obs_frame, text = "Observation Type Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.obs_menu = Listbox(self.obs_frame, listvariable = self.obs_type_names,
                                #width = 40,
                                exportselection = False)
        self.obs_menu.grid(column = 1, row = 2, rowspan = 1, sticky = "N, S, E, W")
        
        self.obs_menu.bind('<Return>', lambda event : self.populate('levels', self.level_menu, event))
        
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
        
        #level selection

        self.levels = StringVar()

        self.level_frame = ttk.Frame(self.main_frame, padding = "2")
        self.level_frame.grid(column = 2, row = 2, sticky = "N, S, E, W")
        ttk.Label(self.level_frame, text = "Observation Level Selection").grid(column = 1,
                                                                               row = 1, sticky = "E, W")
        self.level_menu = Listbox(self.level_frame, listvariable = self.levels, #width = 40,
                                  exportselection = False)
        self.level_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")
        self.populate('levels', self.level_menu)
        self.level_menu.selection_set(0)
        self.level_menu.bind('<Return>', self.plot)
        
        #level scrollbar
        self.level_bar = ttk.Scrollbar(self.level_frame, orient = VERTICAL, command = self.level_menu.yview)
        self.level_menu.configure(yscrollcommand = self.level_bar.set)
        self.level_bar.grid(column = 2, row = 2, rowspan = 2, sticky = "N, S, W")

        #resizing
        self.level_frame.grid_columnconfigure(1, weight = 1)
        self.level_frame.grid_columnconfigure(2, weight = 1)
        self.level_frame.grid_rowconfigure(1, weight = 1)
        self.level_frame.grid_rowconfigure(2, weight = 1)
        
        #region selection

        self.region_frame = ttk.Frame(self.main_frame, padding = "2")
        self.region_frame.grid(column = 2, row = 3, sticky = "N, S, E, W")
        ttk.Label(self.region_frame, text = "Observation Region Selection").grid(column = 1,
                                                                               row = 1, sticky = "E, W")
        self.region_menu = Listbox(self.region_frame, listvariable = self.region_names, #width = 40,
                                   exportselection = False)
        self.region_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")
        self.region_menu.selection_set(0)
        self.region_menu.bind('<Return>', self.plot)
        
        #region scrollbar
        self.region_bar = ttk.Scrollbar(self.region_frame, orient = VERTICAL, command = self.region_menu.yview)
        self.region_menu.configure(yscrollcommand = self.region_bar.set)
        self.region_bar.grid(column = 2, row = 2, rowspan = 2, sticky = "N, S, E, W")

        #resizing
        self.region_frame.grid_rowconfigure(1, weight = 1)
        self.region_frame.grid_rowconfigure(2, weight = 1)
        self.region_frame.grid_columnconfigure(1, weight = 1)
        self.region_frame.grid_columnconfigure(2, weight = 1)
        
        
        #plot button
        self.plot_button = ttk.Button(self.main_frame, text = "Plot", command = self.plot, padding = "2")
        self.plot_button.grid(column = 2, row = 4, sticky = "N, S, E, W")
        
    def populate(self, variable_name, menu, event = None):

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
            
            self.forecast = self.reader.get_variable(self.obs_menu.get(self.obs_menu.curselection()),
                                            'time series', 'forecast', self.original_data)
            self.analysis = self.reader.get_variable(self.obs_menu.get(self.obs_menu.curselection()),
                                                     'time series', 'analysis', self.original_data)
            #find whether level types are plevel, hlevel, or surface

            self.level_type = None

            for coord in self.forecast.coords:
                #get coordinate type
                if coord.lower() in ('plevel', 'hlevel', 'surface'):
                    self.level_type = coord
                    break

            #surface may not show up in coords because it has no dimensions, so this is a workaround
            if self.level_type == None:
                self.level_type = 'surface'
                
            for level in np.unique(self.forecast[self.level_type].values):
                #add warning for levels that are filled entirely with nan's
                forecast = self.reader.filter_single(self.forecast, (self.level_type, level))
                nanmean = np.nanmean(self.reader.filter_single(forecast, ('copy', 'rmse')))
                #a nan value will never equal any value, so if the mean outputs nan this will be false
                if nanmean != nanmean:
                    menu.insert(END, str(level) + ': No data found for this level, will not plot')
                else:
                    menu.insert(END, level)

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

        '''Plot time evolution of rmse data for a selected observation type and level

        Keyword arguments:
        event -- an argument passed by any tkinter menu event. Has no influence on output but
        tkinter requires it to be passed to any method called by a menu event.

        '''
        
        level = np.float64(float(self.level_menu.get(self.level_menu.curselection()).split(':')[0]))
        obs_type = self.obs_menu.get(self.obs_menu.curselection())

        forecast = self.reader.filter_single(self.forecast, (self.level_type, level))
        analysis = self.reader.filter_single(self.analysis, (self.level_type, level))

        possible_obs = self.reader.filter_single(forecast, ('copy', 'Nposs'))
        used_obs = self.reader.filter_single(forecast, ('copy', 'Nused'))

        #only need rmse from forecast and analysis

        forecast = self.reader.filter_single(forecast, ('copy', 'rmse'))
        analysis = self.reader.filter_single(analysis, ('copy', 'rmse'))

        #narrow to one region
        
        forecast_region = self.reader.filter_single(forecast, ('region', int(self.region_menu.curselection()[0]) + 1))
        analysis_region = self.reader.filter_single(analysis, ('region', int(self.region_menu.curselection()[0]) + 1))

        possible_obs_region = self.reader.filter_single(possible_obs, ('region', int(self.region_menu.curselection()[0]) + 1))
        used_obs_region = self.reader.filter_single(used_obs, ('region', int(self.region_menu.curselection()[0]) + 1))

        
        #need to change this to tolerate different numbers of subplots (regions)
        fig, ax = plt.subplots(1, 1, figsize = (12, 8))
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
        
        #get rid of nan values by getting masks of only valid values, then indexing into them during plotting

        forecast_mask = np.array(list(filter(lambda v: v == v, forecast_region.values)))
        forecast_mask = ~np.isnan(forecast_region.values)
        forecast_no_nans = forecast_region.values[forecast_mask]
        forecast_times_no_nans = forecast_region.time.values[forecast_mask.flatten()]

        analysis_mask = ~np.isnan(analysis_region.values)
        analysis_no_nans = analysis_region.values[analysis_mask]
        analysis_times_no_nans = analysis_region.time.values[analysis_mask.flatten()]

        #do not plot regions with no values
        if forecast_no_nans.size == 0:
            ax.text(0.5, 0.5, 'No valid rmse data in this region')
            ax.set_title(str(self.reader.bytes_to_string(self.original_data['region_names'].values[index])) + '     ' +
                 'forecast: mean = ' + str(np.nanmean(forecast_region.values.flatten())) + '     ' +
                 'analysis: mean = ' + str(np.nanmean(analysis_region.values.flatten())))
            return

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

        #pad y axis by 20% at top
        y_max = max(np.nanmax(forecast_region.values.flatten()), np.nanmax(analysis_region.values.flatten()))
        y_max = y_max + .20 * y_max
        locs, labels = plt.yticks()
        locs = np.append(locs, 0)
        plt.yticks(locs)
        ax.set_ylim(-.03*y_max, y_max)

        #add horizontal and vertical lines
        for i in range(1, int(y_max)):
            ax.axhline(y = i, ls = ':')

        for time in ax.get_xticks():
            ax.axvline(x = time, ls = ':')


        #subplot title
        ax.set_title(str(obs_type) + ' @ ' + str(level) + '\n' +
                     self.region_menu.get(self.region_menu.curselection()) + '     ' +
                 'forecast: mean = ' + str(round(np.nanmean(forecast_region.values.flatten()), 5)) + '     ' +
                     'analysis: mean = ' + str(round(np.nanmean(analysis_region.values.flatten()), 5)))


        ax.set_ylabel(self.region_menu.get(self.region_menu.curselection()) + '\n' + 'rmse')
        ax.legend(loc = 'upper left', framealpha = 0.25)


        #need to basically plot two plots on top of each other to get 2 y scales
        ax_twin = ax.twinx()
        ax_twin.scatter(x = possible_obs_region.time.values, y = possible_obs_region.values,
                        color = 'blue', marker = 'o', s = 15, facecolors = 'none')
        ax_twin.scatter(x = used_obs_region.time.values, y = used_obs_region.values,
                        color = 'blue', marker = 'x', s = 15)

        y_max = max(max(possible_obs_region.values), max(used_obs_region.values))
        y_max = y_max + .20 * y_max
        
        locs, labels = plt.yticks()
        locs = np.append(locs, 0)
        plt.yticks(locs)
        ax_twin.set_ylim(-.03*y_max, y_max)
        ax_twin.set_ylabel('# of obs: o = poss, x = used', color = 'blue')
        for tick in ax_twin.get_yticklabels():
            tick.set_color('blue')
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
        ax.tick_params(labelsize = 8)

        #need to add plot title and spacing adjustments
        fig.tight_layout(rect=[0.05, 0.1, 0.95, 0.9])


def main(diag):

    '''create a tkinter GUI for plotting time evolution data
    
    Keyword arguments:
    diag -- an obs_diag output netCDF file

    '''
    
    root = Tk()
    root.title("RMSE Time Evolution Plotter")
    widg = GUIObsDiag(root, 0, 0, diag)
    root.mainloop()

if __name__ == '__main__':
    #only cmd line argument is obs diag file name
    main(sys.argv[1])
        
