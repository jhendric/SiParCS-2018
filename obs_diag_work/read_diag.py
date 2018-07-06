#!/Users/wbd1/anaconda3/bin/python3

import xarray as xa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import time

'''

Reading and other non-plotting functionality for obs_diag files.

'''

class ReadDiag:

    def __init__(self, file_path):
        
        self.full_data = xa.open_dataset(file_path)
        
    def get_variable(self, variable_name, plot_type, data_type, dataset):

        '''
        Extract DataArray from larger dataset given parameters
        variable_name: name of variable e.g. AIRCRAFT_HORIZONTAL_WIND
        plot_type: time series or vertical profile
        data_type: forecast or analysis
        dataset: xarray Dataset to grab from
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
        
        data = data_array

        for (category_name, value) in args:

            #should only be one value in this case
            
            category = getattr(data, category_name)

            if category_name != 'copy':

                data = data.where(abs(category - value) < 1e-8, drop = True)
                
            else:
                
                #strings in file need to be decoded into Python strings
                def bytes_to_string(self, bytes):
                    return ''.join(bytes.decode('UTF-8').split())

                
                value_converted = None

                meta_data = self.full_data['CopyMetaData']

                for i in range(len(meta_data)):
                    if bytes_to_string(self, meta_data.values[i]) == value:
                        value_converted = meta_data['copy'].values[i]
                        break
                    
                #have not tested this line
                data = data.where(data['copy'] == value_converted, drop = True)
                
        print(data)
        return data
    
    def filter_multiple(self, data_array, *args):
        
        data = data_array

        for (category_name, values) in args:

            category = getattr(data, category_name)
            print(data.region)
            #below line may not work
            if category_name != 'copy':

                if type(category.values[0]) == np.dtype('float64'):
                    #rounding to prevent flaot comparison mistakes
                    print('FLOAT')
                    mask = np.isin(np.around(category.values, 1), np.around(values, 1))

                else:
                    mask = np.isin(category.values, values)

            else:

                #strings in file need to be decoded into Python strings
                def bytes_to_string(self, bytes):
                    return ''.join(bytes.decode('UTF-8').split())

                
                values_converted = ()

                meta_data = self.full_data['CopyMetaData']
                for value in values:
                    #convert strings to copy indices
                    for i in range(len(meta_data)):
                        if bytes_to_string(self, meta_data.values[i]) == value:
                            values_converted.append(meta_data['copy'].values[i])

                #this may not work because copy is a reserved word in python
                mask = np.isin(category.values, values_converted)

            print(mask)
            data = data[mask]

        print(data)
        return data
    
    def plot_evolution(self, obs_type, level, dataset):
        '''Replicate plot_evolution.m'''

        #get data of chosen variable
        forecast = self.get_variable(obs_type, 'time series', 'forecast', dataset)
        analysis = self.get_variable(obs_type, 'time series', 'analysis', dataset)

        #find whether level is plevel, hlevel, or surface

        level_type = None
        
        for coord in forecast.coords:
            
            if coord.lower() in ('plevel', 'region', 'surface'):

                level_type = coord

        #further filter data to proper level
        forecast = self.filter_single(forecast, (coord, level))
        analysis = self.filter_single(analysis, (coord, level))
        
        possible_obs = self.filter_single(analysis, ('copy', 'Nposs'))
        used_obs = self.filter_single(analysis, ('copy', 'Nused'))
        
        #only need rmse from forecast and analysis

        forecast = self.filter_single(forecast, ('copy', 'rmse'))
        analysis = self.filter_single(analysis, ('copy', 'rmse'))

        #create 4 subplots
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1)
        
        
    
def main():
    
    reader = ReadDiag('../obs_series/obs_diag_output.nc')
    data = reader.get_variable('RADIOSONDE_U_WIND_COMPONENT', 'time series', 'forecast',
                               reader.full_data)
    print(data)
    data = reader.filter_single(data, ('region', 2))
    data = reader.filter_single(data, ('copy', 'rmse'))
    

if __name__ == "__main__":
    main()
            

        
