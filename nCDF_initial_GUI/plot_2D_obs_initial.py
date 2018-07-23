'''#!/Users/wbd1/anaconda3/bin/python3'''

import xarray as xa
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import pandas as pd
import math

class plot_2D_obs:

    '''

    Redo of "3D_plot_initial.py" to have structuring more appropriate for
    GUI work.

    A class for plotting 2D observation data from DART's obs_sequence_tool. 
    Lifted some variables names from read_obs_netcdf.m and plot_obs_netcdf.m.

    Utilizes both pandas and xarray to work around structure of a DART obs sequence file.
    The non-gridded structure is not friendly to xarray. This class could probably be
    modified to never use xarray, but besides a slow GUI boot time everything works fine
    in this format.
    
    '''
    
    def __init__(self, file_path):

        '''Initialize object with relevant information from an obs sequence file opened in xarray

        Keyword arguments:
        file_path -- path of obs sequence file
        
        '''
        
        #mkl.set_num_threads(128)
        
        dataset = xa.open_dataset(file_path)

        #all "strings" from the file actually need to be decoded into Python strings

        def bytes_to_string(self, bytes):
            return ''.join(bytes.decode('UTF-8').split())

        bytes_to_string = np.vectorize(bytes_to_string)

        obs_type_strings = bytes_to_string(self, dataset['ObsTypesMetaData'].values)
        copy_strings = bytes_to_string(self, dataset['CopyMetaData'].values)
        QC_strings = bytes_to_string(self, dataset['QCMetaData'].values)
        obs_types_meta_indexer = dataset['ObsTypes']
        
        times = dataset['time']
        obs_types = dataset['obs_type']
        #not currently using keys but defining it here
        keys = dataset['obs_keys']
        vert_types = dataset['which_vert']
        
        locs = dataset['location']
        obs = dataset['observations']
        qc = dataset['qc']

        my_types = np.unique(obs_types.values)
        time_units = times.dtype
        time_range = times.attrs['valid_range']
        calendar = 'Gregorian' #don't know how to get other calendar types
        #missing some other time stuff

        #only getting things that are the actual obs (not prior and posterior members)

        obs_ind = 0
        obs = obs[:, obs_ind]
        
        lons = locs[:, 0]
        lats = locs[:, 1]
        z = locs[:, 2]

        #first create a pandas DataFrame
        
        pd_array = pd.DataFrame({'obs_types' : obs_types.values,
                                 'times' : times.values,
                                 'lons' : lons.values,
                                 'lats' : lats.values,
                                 'z' : z.values,
                                 'vert_types' : vert_types.values,
                                 'qc_DATA' : qc[:, 0].values,
                                 'qc_DART' : qc[:, 1].values,
                                 'obs' : obs.values
                                 })

        #use set_index to create a MultiIndex so I can access fields that xarray wouldn't normally
        #treat as accessible fields
        
        pd_array.set_index(['obs_types', 'times', 'lons', 'lats', 'z',
                            'qc_DATA', 'qc_DART', 'vert_types'], inplace = True)
        pd_array.sort_index(inplace = True)
        pd_array.reset_index()
        
        ''' Function guideline for creating a proper dimensional xarray DataArray
        if desired at some point in the future.
        

        test_array = np.column_stack((np.array(lons.values, columns = ['lons']),
        np.array(lats.values, columns = ['lats']),
        np.array(obs.values, columns = ['obs'])))
        
        def map(array):
            keys, values = array.sort_values('lons').values.T
            ukeys, index = np.unique(keys, True)
            arrays = np.split(values, index[1:])
            array2 = pd.DataFrame({'lons' : ukeys, 'obs' : [list(a) for a in arrays]})
            return array2

        def map_advanced(array):
            return None
        
        '''
        
        
        #Convert to pd_array to an xarray DataArray
        self.data = xa.DataArray(pd_array)
        self.data = self.data.sortby(self.data.obs_types)


        #map obs type strings to obs type integers that observations have
        obs_type_dict = dict([(type_string, index +1 )
                                   for index, type_string in enumerate(obs_type_strings)])

        #inverse of above dictionary
        self.obs_type_inverse = dict((v, k) for k, v in obs_type_dict.items())

        existing_obs = np.unique(self.data.obs_types.values)

        #dict containing only the types that exist in the file

        self.obs_type_dict = dict([(self.obs_type_inverse[i], i)
                                   for i in existing_obs])

        
        
    def filter_range(self, data_array, conditions):
        '''Takes list of tuples of form ('category_name', min, max) and return lons, lats, and
        observation values satisfying these conditions
        
        Keyword arguments:
        data_array -- xarray DataArray to filter
        conditions -- list of tuples of form ('category_name', min, max) where min and max are the
                      desired limits for the selected category
        
        '''

        data = data_array

        cat_dict = {
            
            'obs_types' : data.obs_types,
            'times' : data.times,
            'lons' : data.lons,
            'lats' : data.lats,
            'z' : data.z,
            'qc_DATA' : data.qc_DATA,
            'qc_DART' : data.qc_DART,
            'vert_types' : data.vert_types
            
        }
        
        for (category_name, min, max) in conditions:
            
            category = cat_dict[category_name]
            
            if min != max:
                if min < max:
                    data = data.where(category >= min)
                    data = data.where(category <= max)
                    data = data.where(np.isnan(data) != True, drop = True)
                else:
                    #for wrapping data (particularly longitudes or times)
                    data = xa.concat([data.where(category >= min, drop = True),
                                      data.where(category <= max, drop = True)], dim = 'dim_0')
            else:
                data = data.where(category == min, drop = True)

        return data
    
        
    def filter(self, data_array, *args):
        '''Takes list of tuples of form ('category_name', (values)) and returns DataArray masked to these values. 
        Uses numpy's isin function, a very fast masking function that accepts lists of values for comparison
        unlike xarray

        Keyword arguments:
        data_array -- xarray DataArray to filter
        *args -- tuples of form ('category_name', (values)). Points which do not have a value in category_name matching
                 any value in values will be removed
        
        '''
        
        data = data_array

        #this can be changed to a getattr command
        cat_dict = {
            
            'obs_types' : data.obs_types,
            'times' : data.times,
            'lons' : data.lons,
            'lats' : data.lats,
            'z' : data.z,
            'qc_DATA' : data.qc_DATA,
            'qc_DART' : data.qc_DART,
            'vert_types' : data.vert_types
            
        }

        data_building = data
        
        for (category_name, values) in args:
            
            category = cat_dict[category_name]
            
            if type(category.values[0]) == np.dtype('float64'):
                #rounding to prevent float comparison mistakes
                mask = np.isin(np.around(category.values, 1), np.around(values, 1))

            else:
                mask = np.isin(category.values, values)
            

            #index data using mask to get needed values
            data = data[mask]

        #sort data for plotting in GUI
        data = data.sortby(data.obs_types)
        return data    

        
