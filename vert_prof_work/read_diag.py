#!/Users/wbd1/anaconda3/bin/python3

import xarray as xa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import time
import datetime



class ReadDiag:

    '''
    
    Reading and other non-plotting functionality for obs_diag files.
    
    '''
    
    def __init__(self, file_path):

        '''Initialize reader object with obs diag file opened in xarray

        Keyword arguments:
        file_path -- path of obs_diag output netCDF file
        
        '''
        
        self.full_data = xa.open_dataset(file_path)

    def bytes_to_string(self, bytes):

        '''Convert a byte string variable to a typical Python string format

        Keyword arguments:
        bytes -- byte string to be decoded

        '''


        return ''.join(bytes.decode('UTF-8').split())
                
    def get_variable(self, variable_name, plot_type, data_type, dataset):

        '''Extract DataArray from larger dataset using given parameters

        Keyword arguments:
        variable_name -- name of variable e.g. AIRCRAFT_HORIZONTAL_WIND
        plot_type -- time series or vertical profile
        data_type -- forecast or analysis
        dataset -- xarray Dataset to grab from

        '''
        
        plot_type_dict = {
            'time series' : '_',
            'vertical profile' : '_VP'
        }
        
        data_type_dict = {
            'forecast' : 'guess',
            'analysis' : 'analy'
        }
        

        #concatenate strings and then index into dataset to get matching xarray DataArray
        data_string = variable_name + plot_type_dict[plot_type] + data_type_dict[data_type]
        return dataset[data_string]

    def filter_single(self, data_array, *args):

        '''Filter an xarray DataArray using a single condition

        Keyword arguments:
        data_array -- DataArray to be filtered
        *args -- set of tuples of form (coordinate name, value)
        
        '''
        
        data = data_array

        for (category_name, value) in args:

            #should only be one value in this case
            
            category = getattr(data, category_name)

            if category_name != 'copy':

                if type(category.values[0]) == np.dtype('float64'):
                    data = data.where(abs(category - value) < 1e-4, drop = True)
                else:
                    data = data.where(category == value, drop = True)
                
            else:
                
                value_converted = None

                meta_data = self.full_data['CopyMetaData']

                for i in range(len(meta_data)):
                    if self.bytes_to_string(meta_data.values[i]) == value:
                        value_converted = meta_data['copy'].values[i]
                        break
                    
                #have not tested this line
                data = data.where(data['copy'] == value_converted, drop = True)
                
        #print(data)
        return data
    
    def filter_multiple(self, data_array, *args):

        '''Filter an xarray DataArray using multiple conditions

        Keyword arguments:
        data_array -- DataArray to be filtered
        *args -- set of tuples of form (coordinate name, value)
        
        '''
        
        data = data_array

        for (category_name, values) in args:

            category = getattr(data, category_name)
            print(data.region)
            #below line may not work
            if category_name != 'copy':

                if type(category.values[0]) == np.dtype('float64') or type(category.values[0]) == np.dtype('float32'):
                    #rounding to prevent float comparison mistakes
                    print('FLOAT')
                    mask = np.isin(np.around(category.values, 1), np.around(values, 1))

                else:
                    mask = np.isin(category.values, values)

            else:
                
                values_converted = ()

                meta_data = self.full_data['CopyMetaData']
                for value in values:
                    #convert strings to copy indices
                    for i in range(len(meta_data)):
                        if self.bytes_to_string(meta_data.values[i]) == value:
                            values_converted.append(meta_data['copy'].values[i])

                #this may not work because copy is a reserved word in python
                mask = np.isin(category.values, values_converted)

            print(mask)
            data = data[mask]

        print(data)
        return data
    
    def plot_evolution(self, obs_type, level, dataset):

        '''Replicate plot_evolution.m. Not used in any GUI file
        
        Keyword arguments:
        obs_type -- observation type to be plotted
        level -- observation level to be plotted
        dataset -- xarray Dataset to grab data from

        '''

        #get data of chosen variable
        forecast = self.get_variable(obs_type, 'time series', 'forecast', dataset)
        analysis = self.get_variable(obs_type, 'time series', 'analysis', dataset)

        print('forecast: ', forecast)
        print('analysis: ', analysis)
        #find whether level is plevel, hlevel, or surface

        level_type = None
        
        for coord in forecast.coords:
            #get coordinate type
            if coord.lower() in ('plevel', 'hlevel', 'surface'):
                level_type = coord
                break
            
        print('level_type: ', level_type)
        print(coord)
        #further filter data to proper level
        forecast = self.filter_single(forecast, (coord, level))
        analysis = self.filter_single(analysis, (coord, level))
        print('forecast filtered to 1 level', forecast)
        print('analysis filtered to 1 level', analysis)
            
        possible_obs = self.filter_single(forecast, ('copy', 'Nposs'))
        used_obs = self.filter_single(forecast, ('copy', 'Nused'))
        
        #only need rmse from forecast and analysis

        forecast = self.filter_single(forecast, ('copy', 'rmse'))
        analysis = self.filter_single(analysis, ('copy', 'rmse'))
        print('forecast filtered to rmse', forecast)
        print('analysis filtered to rmse', analysis)
        
        #create 4 subplots
        #need to modularize this and function parameters to allow different numbers of plots
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize = (8,8))
        
        for index, ax in enumerate((ax1, ax2, ax3)):

            forecast_region = forecast.where(forecast.region == index + 1, drop = True)
            analysis_region = analysis.where(analysis.region == index + 1, drop = True)
            possible_obs_region = possible_obs.where(possible_obs.region == index + 1, drop = True)
            used_obs_region = used_obs.where(used_obs.region == index + 1, drop = True)
            #print('forecast region: ', forecast_region)
            #print('analysis region: ', analysis_region)
            #plot both scatter and line for forecast and analysis to achieve connected appearance
            '''
            ax.scatter(x = forecast_region.time.values.astype("M8[ms]").tolist(), y = forecast_region.values, color = 'black')
            ax.plot(forecast_region.time.values.astype("M8[ms]"), forecast_region.values.flatten(), 'ko-')
            ax.scatter(x = analysis_region.time.values.astype("M8[ms]"), y = analysis_region.values, color = 'red')
            ax.plot(analysis_region.time.values.astype("M8[ms]"), analysis_region.values.flatten(), 'ro-')

            ax.set_xlim(forecast_region.time.values[0].astype("M8[ms]"),
                         forecast_region.time.values[forecast_region.time.values.size - 1].astype("M8[ms]"))
            '''

            
            ax.scatter(x = forecast_region.time.values, y = forecast_region.values,
                       edgecolors = 'black', marker = 'x', s = 15)
            ax.plot(forecast_region.time.values, forecast_region.values.flatten(), 'kx-', label = 'forecast')
            ax.scatter(x = analysis_region.time.values, y = analysis_region.values,
                       edgecolors = 'red', marker = 'o', s = 15, facecolors = 'none')
            ax.plot(analysis_region.time.values, analysis_region.values.flatten(), 'ro-', mfc = 'none', label = 'analysis')

            #set min and max x limits to be wider than actual limits for a nicer plot
            #pad x axis by 10% on both sides
            pad_x_axis = .10 * (max(forecast_region.time.values) - min(forecast_region.time.values))
            ax.set_xlim((forecast_region.time.values[0] - pad_x_axis).astype("M8[ms]"),
                        (forecast_region.time.values[-1] + pad_x_axis).astype("M8[ms]"))
            
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
            print('forecast: ', forecast_region.values.flatten())
            print('forecast mean: ', np.nanmean(forecast_region.values.flatten()))
            ax.set_title(str(self.bytes_to_string(dataset['region_names'].values[index])) + '     ' +
                     'forecast: mean = ' + str(np.nanmean(forecast_region.values.flatten())) + '     ' +
                     'analysis: mean = ' + str(np.nanmean(analysis_region.values.flatten())))
            
            
            ax.set_ylabel(str(self.bytes_to_string(dataset['region_names'].values[index])) + '\n' + 'rmse')
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
            
            
            
            
            #possible_obs_region.scatter

        #need to add units
        plt.suptitle(str(obs_type) + ' @ ' + str(level))
        plt.subplots_adjust(hspace = 0.8)
        #fig.tight_layout()
        plt.show()
    
def main():

    '''Initializes a reader and plots from test file'''
    
    reader = ReadDiag('obs_diag_output.nc')
    reader.plot_evolution('GPSRO_REFRACTIVITY', 315.0, reader.full_data)

if __name__ == "__main__":
    main()
            

        
