from tkinter import *
from tkinter import ttk

import xarray as xa

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import numpy as np
from read_diag import ReadDiag
np.set_printoptions(threshold = np.nan) #avoids incomplete StringVar declaration


class GUIVertProf:

    '''
    
    Uses read_diag.py for plotting obs diag info in a GUI
    
    '''
    
    def __init__(self, window, grid_col, grid_row, diag):
        
        '''Initialize GUI with available observation types and regions.
        
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

        #remove obs types that don't actually have vertical profile data
        i = 0
        while i < len(obs_types_sparse):
            try:
                #all obs types with VP data will have a corresponding data variable with this suffix
                self.original_data[obs_types_sparse[i] + '_VPguess']
                i += 1
            except KeyError:
                del obs_types_sparse[i]


        self.obs_type_names = StringVar(value = obs_types_sparse)

        regions = [self.reader.bytes_to_string(name) for name in self.original_data.region_names.values]

        regions_sparse = [name.replace('[', '').replace(']', '').
                            replace(',', '').replace('\'', '')
                            for name in regions]

        self.region_names = StringVar(value = regions_sparse)
        
        self.forecast = None
        self.analysis = None
        
        #GUI config

        #observation selection

        self.obs_frame = ttk.Frame(self.main_frame, padding = "2")
        self.obs_frame.grid(column = 2, row = 1, sticky = "N, S, E, W")
        ttk.Label(self.obs_frame, text = "Observation Type Selection").grid(column = 1, row = 1, sticky = "E, W")
        self.obs_menu = Listbox(self.obs_frame, listvariable = self.obs_type_names,
                                #height = 18, width = 40,
                                exportselection = False)
        self.obs_menu.grid(column = 1, row = 2, rowspan = 1, sticky = "N, S, E, W")
        
        self.obs_menu.bind('<Return>', self.plot)
        
        self.obs_menu.selection_set(0)
        self.obs_menu.event_generate('<<ListboxSelect>>')

        #resizing
        self.obs_frame.grid_columnconfigure(1, weight = 1)
        self.obs_frame.grid_columnconfigure(2, weight = 1)
        self.obs_frame.grid_rowconfigure(1, weight = 1)
        self.obs_frame.grid_rowconfigure(2, weight = 1)
        
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
                                  #height = 18, width = 40,
                                   exportselection = False)
        self.region_menu.grid(column = 1, row = 2, sticky = "N, S, E, W")
        self.region_menu.selection_set(0)
        self.region_menu.bind('<Return>', self.plot)

        #resizing
        self.region_frame.grid_columnconfigure(1, weight = 1)
        self.region_frame.grid_columnconfigure(2, weight = 1)
        self.region_frame.grid_rowconfigure(1, weight = 1)
        self.region_frame.grid_rowconfigure(2, weight = 1)
        
        #region scrollbar
        self.region_bar = ttk.Scrollbar(self.region_frame, orient = VERTICAL, command = self.region_menu.yview)
        self.region_menu.configure(yscrollcommand = self.region_bar.set)
        self.region_bar.grid(column = 2, row = 2, rowspan = 2, sticky = "N, S, E")
        
        #plot button
        self.plot_button = ttk.Button(self.main_frame, text = "Plot", command = self.plot, padding = "2")
        self.plot_button.grid(column = 2, row = 5, sticky = "N, S, E, W")
        
    def plot(self, event = None):

        '''Plot a vertical profile of rmse data given selected obs type and region
        
        Keyword arguments:
        event -- an argument passed by any tkinter menu event. Has no influence on output but
        tkinter requires it to be passed to any method called by a menu event.

        '''

        obs_type = self.obs_menu.get(self.obs_menu.curselection())

        forecast = self.reader.get_variable(self.obs_menu.get(self.obs_menu.curselection()),
                                            'vertical profile', 'forecast', self.original_data)
        analysis = self.reader.get_variable(self.obs_menu.get(self.obs_menu.curselection()),
                                            'vertical profile', 'analysis', self.original_data)

        #narrow to one region
        forecast_region = self.reader.filter_single(forecast, ('region', int(self.region_menu.curselection()[0]) + 1))
        analysis_region = self.reader.filter_single(analysis, ('region', int(self.region_menu.curselection()[0]) + 1))

        possible_obs = self.reader.filter_single(forecast_region, ('copy', 'Nposs'))
        used_obs = self.reader.filter_single(forecast_region, ('copy', 'Nused'))

        #only need rmse from forecast and analysis

        forecast_region = self.reader.filter_single(forecast_region, ('copy', 'rmse'))
        analysis_region = self.reader.filter_single(analysis_region, ('copy', 'rmse'))

        #get level type
        level_type = None
        for coord in forecast_region.coords:
            if coord != 'copy' and coord != 'region':
                level_type = coord

        forecast_levels = forecast_region[level_type]
        analysis_levels = analysis_region[level_type]

        fig, ax = plt.subplots(1, 1, figsize = (8, 8))
        canvas = FigureCanvasTkAgg(fig, master = self.main_frame)
        canvas.get_tk_widget().grid(column = 1, row = 1, rowspan = 4, sticky = "N, S, E, W")
        

        #have to set up a separate toolbar frame because toolbar doesn't like gridding with others
        self.toolbar_frame = ttk.Frame(self.main_frame)
        self.toolbar = NavigationToolbar2TkAgg(canvas, self.toolbar_frame)
        self.toolbar_frame.grid(column = 1, row = 5, sticky = "N, S, E, W")
        self.toolbar.grid(column = 1, row = 1, sticky = "N, S ,E , W")

        #resizing
        self.toolbar_frame.grid_columnconfigure(1, weight = 1)
        self.toolbar_frame.grid_rowconfigure(1, weight = 1)
        
        
        #get rid of nan values by getting masks of only valid values, then indexing into them during plotting
        
        forecast_mask = np.array(list(filter(lambda v: v == v, forecast_region.values.flatten())))
        forecast_mask = ~np.isnan(forecast_region.values)
        forecast_no_nans = forecast_region.values[forecast_mask]
        forecast_levels_no_nans = forecast_levels.values[forecast_mask.flatten()]
        analysis_mask = ~np.isnan(analysis_region.values)
        analysis_no_nans = analysis_region.values[analysis_mask]
        analysis_levels_no_nans = analysis_levels.values[analysis_mask.flatten()]

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
        
        if forecast_levels.values[-1] > forecast_levels.values[0]:
            ax.set_ylim((forecast_levels_no_nans[0] -  pad_y_axis),
                        (forecast_levels_no_nans[-1] +  pad_y_axis))
        else:
            ax.set_ylim((forecast_levels_no_nans[0] + pad_y_axis,
                         forecast_levels_no_nans[-1] - pad_y_axis))
        
        #pad x axis by 20% on right side
        x_max = max(np.nanmax(forecast_region.values.flatten()),
                         np.nanmax(analysis_region.values.flatten()))
        x_min = min(np.nanmin(forecast_region.values.flatten()),
                         np.nanmin(analysis_region.values.flatten()))
        pad_x_axis = .10 * (x_max - x_min)
        ax.set_xlim(0, x_max + pad_x_axis)

        #add horizontal and vertical lines
        for i in range(1, int(x_max)):
            ax.axvline(x = i, ls = ':')

        for lev in ax.get_yticks():
            ax.axhline(y = lev, ls = ':')

        ax.set_xlabel(self.region_menu.get(self.region_menu.curselection()) + '\n' + 'rmse')
        ax.legend(loc = 'upper left', framealpha = 0.25)


        #need to basically plot two plots on top of each other to get 2 x scales
        ax_twin = ax.twiny()
        ax_twin.scatter(y = possible_obs[level_type].values, x = possible_obs.values,
                        color = 'blue', marker = 'o', s = 15, facecolors = 'none')
        ax_twin.scatter(y = used_obs[level_type].values, x = used_obs.values,
                        color = 'blue', marker = 'x', s = 15)

        x_max = max(max(possible_obs.values.flatten()),
                    max(used_obs.values.flatten()))
        x_min = min(min(possible_obs.values.flatten()),
                    min(used_obs.values.flatten()))
        locs, labels = plt.xticks()
        locs = np.append(locs, 0)
        plt.xticks(locs)
        pad_x_axis = .10 * (x_max - x_min)
        ax_twin.set_xlim(-.03 * x_max, x_max + pad_x_axis)
        ax_twin.set_xlabel('# of obs: o = poss, x = used', color = 'blue')
        for tick in ax_twin.get_xticklabels():
            tick.set_color('blue')
        
        fig.suptitle(str(obs_type) + '\n' + self.region_menu.get(self.region_menu.curselection()) + '     ' +
                     'forecast: mean = ' + str(round(np.nanmean(forecast_region.values.flatten()), 5)) + '     ' +
                     'analysis: mean = ' + str(round(np.nanmean(analysis_region.values.flatten()), 5)))
        
        #make sure title doesn't get cutoff
        fig.tight_layout(rect=[0, 0.03, 1, 0.92])


def main(diag):

    '''create a tkinter GUI for plotting vertical profile data
    
    Keyword arguments:
    diag -- an obs_diag output netCDF file
    
    '''
    
    root = Tk()
    root.title("Vertical Profile Plotter")
    widg = GUIVertProf(root, 0, 0, diag)
    root.mainloop()
        
if __name__ == '__main__':
    #only cmd line argument is obs diag file name
    main(sys.argv[1])


        
