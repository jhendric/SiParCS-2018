#!/Users/wbd1/anaconda3/bin/python3

import mkl
import xarray as xa
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import math
import time
'''

Redo of "3D_plot_initial.py" to have structuring more appropriate for
GUI work.

A class for plotting 2D observation data from DART's obs_sequence_tool. 
Lifted variables names from read_obs_netcdf.m and plot_obs_netcdf.m.


Possible 1-off indexing issues. None obvious right now.

Also string format is weird in nCDF (b's before each entry in large string sets)

Also may have to fiddle with datatypes in initializing args, as well as make sure
inf is float('inf')

'''

class plot_2D_obs:

    def __init__(self, file_path):

        mkl.set_num_threads(128)
        
        dataset = xa.open_dataset(file_path)

        #all "strings" form the file actually need to be decoded into Python strings
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
        
        '''
        MISSING MISSING DATA ACCOUNTABILITY HERE
        '''
        
        locs = dataset['location']
        obs = dataset['observations']
        print(obs.size)
        qc = dataset['qc']

        my_types = np.unique(obs_types.values)
        time_units = times.dtype
        time_range = times.attrs['valid_range']
        calendar = 'Gregorian' #don't know how to get other calendar types
        #missing some other time stuff

        #don't have the verbose section

        #only getting things that are the actual obs (not prior and posterior members)
        obs_ind = 0
        obs = obs[:, obs_ind]
        print(obs.size)
        
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

        
        
        ''' Function guideline for if I ever want to create a proper dimensional xarray DataArray
        
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
        
        
        print(map(pd_array))
        '''
        
        
        #Convert to pd_array to an xarray DataArray
        self.data = xa.DataArray(pd_array)
        self.data = self.data.sortby(self.data.obs_types)
        print(self.data)

        #map obs type strings to obs type integers that observations have
        obs_type_dict = dict([(type_string, index +1 )
                                   for index, type_string in enumerate(obs_type_strings)])

        #inverse of above dictionary
        self.obs_type_inverse = dict((v, k) for k, v in obs_type_dict.items())

        existing_obs = np.unique(self.data.obs_types.values)

        #dict containing only the types that exist in the file
        self.obs_type_dict = dict([(self.obs_type_inverse[i], i)
                                   for i in existing_obs])

        
        
    def filter_range(self, conditions):
        '''Take list of tuples of form ('category_name', min, max) and return lons, lats, and
        observation values satisfying these conditions'''

        data = self.data

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
                           
    def filter_disjoint(self, dataset, *args):
        '''Take list of tuples of form ('category_name', [values]) and return lons, lats, and
        observation values satisfying these conditions
        Discontinued in favor of filter_test which has significantly better performance'''
        

        data = dataset

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

            #this form may not work but maybe it will
            
            #data = data.where(category == values[, drop = True)

            
            #slower less pythonic version
            data_building = data.where(abs(category - values[0]) < 1e-8, drop = True)
            i = 1
            
            while i < len(values):
                
                data_building = xa.concat([data_building, data.where(abs(category-values[i]) < 1e-8, drop = True)], dim = 'dim_0')
                i += 1
                
            data = data_building 

        #sort data for plotting in GUI
        data = data.sortby(data.obs_types)
        print(data)
        return data

    def filter_test(self, dataset, *args):
        '''uses numpy's isin function, a very fast masking function that accepts lists of values for comparison
        unlike xarray'''
        
        data = dataset

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
                print('FLOAT')
                mask = np.isin(np.around(category.values, 1), np.around(values, 1))

            else:
                mask = np.isin(category.values, values)
            
            #index data using mask to get needed values
            data = data[mask]

        #sort data for plotting in GUI
        data = data.sortby(data.obs_types)
        print(data)
        return data

    
    def plot(self, *args):
        '''Each argument represents a range of values to be passed to filter. Any argument given
        should be a tuple ('category_name', min, max) representing the desired coordinate range.
        min and max are inclusive. List of valid arguments: obs_types, times, lons, lats, z,
        qc_DATA, qc_DART, vert_types'''

        print('at plot')
        data = self.filter_disjoint(args)
        print('at plot further')
        ax = plt.axes(projection = ccrs.PlateCarree())
        ax.stock_img()
        ax.gridlines()
        ax.coastlines()
        plt.scatter(data.lons, data.lats, c = data.qc_DART, s = 100,
                    marker = "+", transform = ccrs.PlateCarree())
        plt.tight_layout()
        plt.show()
        
'''
plotter = plot_2D_obs('../obs_series/obs_epoch_001.nc')
#plotter.plot(('obs_types', 10, 20))                
plotter.plot(('obs_types', [3, 4, 5, 6, 7, 8]))                
'''



''' CODE GRAVEYARD DO NOT ENTER 

obs_type_string, region, copy_string, QC_string, max_QC
#def __init__(self, file_path, obs_type_string, region, copy_string, QC_string, max_QC, verbose):   
        self.region = np.array(region) #e.g. [0, 360, -90, 90, -inf, inf] or [0, 360, -90, 90]
        
        'RADIOSONDE_TEMPERATURE',
                  [0, 360, -90, 90, -float('inf'), float('inf')], 'observation', '', 2, 1)

        if copy_string.lower() == 'all':
            self.copy_string = self.dataset['CopyMetaData'].values
        else:
            self.copy_string = copy_string
            
        self.QC_string = QC_string #e.g. 'DART quality control'
        self.max_QC = max_QC #e.g. 2
        self.verbose = verbose #True or False

        self.obs_type_string = obs_type_string #e.g. 'RADIOSONDE_U_WIND_COMPONENT', 'ALL'

        #find observations of the correct type
        if obs_type_string.lower() == 'all':
            self.inds = [i for i in range(0,self.obs.size)]

        else:
            self.my_ind_array = np.where(self.obs_type_strings == self.obs_type_string)[0]
            #print(self.obs_type_strings[0])
            #print(self.obs_type_string)
            #print(self.my_ind_array)
            if self.my_ind_array.size > 0:
                self.my_ind = self.my_ind_array[0] + 1 #first instance of my_ind
                self.inds = np.where(self.obs_type.values == self.my_ind)[0] #all instances of associated type

            else:
                print('FYI - no {} observations \n'.format(self.obs_type_string))
                self.inds = np.array([])


        self.my_obs = self.obs[self.inds, self.obs_ind].values

        self.my_locs = self.location[self.inds].value
        self.my_keys = self.obs_keys[self.inds].values
        self.my_vert_types = self.vert_type[self.inds].values
        self.my_times = self.time[self.inds].values

        
        self.num_copies = self.dataset['copy'].size
        
        if copy_string.lower() == 'all':
            self.my_type_ind = [i for i in range(0,self.num_copies)]
        else:
            #get_copy_index.m
            copy_meta_data = self.dataset['CopyMetaData']
            meta_data_length  = self.dataset['nlines'].size
            copy_string_no_white = ''.join(self.copy_string.split())
            copy_index = -1

            #find matching copy
            for i in range(0, self.num_copies):
                meta_data_no_white = ''.join(bytes_to_string(self, copy_meta_data.values[i].split()))
                if meta_data_no_white == copy_string_no_white:
                    copy_index = i
                    break

            if copy_index < 0:
                raise Exception('No matching copy found')

            self.my_type_ind = copy_index
        
        

        #find QC values for obs

        if self.QC_string != '':
            #get_qc_index.m
            QC_meta_data = self.dataset['QCMetaData']
            #an issue is pointed out in matlab code where order of this array changes
            #based on number of copies
            num_QC_copies = self.dataset['qc_copy'].size
            copy_string_no_white = ''.join(self.QC_string.split())
            copy_index = -1

            #find matching copy
            for i in range(0, self.num_QC_copies):
                meta_data_no_white = ''.join(str(QC_meta_data.values[i].split()))
                if meta_data_no_white == copy_string_no_white:
                    copy_index = i
                    break

            if copy_index < 0:
                raise Exception('No matching copy found')

            self.my_QC_ind = copy_index
            #may be an issue here, confused about QC structure
            self.my_QC = self.QC[self.inds, self.my_QC_ind].values
            
        else:
            self.my_QC = np.array([])



        
        #locations_in_region.m
        if self.region.size == 6:
            z_min = min(self.region[[4,5]])
            z_max = max(self.region[[4,5]])
            
        elif self.region.size < 4:
            raise Exception('Region msut be an array of length 4 or 6')

        else:
            z_min = float(-inf)
            z_max = float(inf)
        
        y_min = min(self.region[[2,3]])
        y_max = max(self.region[[2,3]])
        x_min = self.region[0]%360.0
        x_max = self.region[1]%360.0

        #may not be possible with numpy array, may need some numpy function
        lons = self.my_locs[:, 0]%360.0

        #Need to revise the above code to no longer be a logical (should be a way to do it as one where)
        
        


        #need to redo the following chunk to be similar to the commented chunk below (faster), but not sure how
        #would be more pythonic

        self.loc_inds = []
        if (x_min == x_max):
            for i in range(self.my_locs[:, 0].size):
                #print(i)
                cur_loc = self.my_locs[i]
                #print(i)
                if cur_loc[1] >= y_min and cur_loc[1] <= y_max and cur_loc[2] >= z_min and cur_loc[2] <= z_max:
                    self.loc_inds.append(i)
                    #print(self.loc_inds)
        elif (x_min > x_max):
            for i in range(self.my_locs[:, 0].size):
                cur_loc = self.my_locs[i]
                if ((cur_loc[0] >= x_min or cur_loc[0] <= x_max) and
                    cur_loc[1] >= y_min and cur_loc[1] <= y_max and
                    cur_loc[2] >= z_min and cur_loc[2] <= z_max):

                    self.loc_inds.append(i)
        else:
            for i in range(self.my_locs[:, 0].size):
                cur_loc = self.my_locs[i]
                if ((cur_loc[0] >= x_min and cur_loc[0] <= x_max) and
                    cur_loc[1] >= y_min and cur_loc[1] <= y_max and
                    cur_loc[2] >= z_min and cur_loc[2] <= z_max):

                    self.loc_inds.append(i)
        
        print(np.where(np.bitwise_and(self.my_locs[:, 1] > -40, self.my_locs [:, 1] < 0,
                                      self.my_locs[:, 2] >= z_min, self.my_locs[:, 2] <= z_max))[0].size)
        
        for i in range(self.mylocs[:, 0].size):
            if (x_min == x_max):
                
        if (x_min == x_max):
            self.loc_inds = np.where(self.my_locs[:, 1] >= y_min and self.my_locs[:, 1] <= y_max and
                                     self.my_locs[:, 2] >= z_min and self.my_locs[:, 2] <= z_max)
        else:
            if (x_min > x_max):
                #this line does not follow the sames structure as in matlab, but seems likely to work
                self.loc_inds = np.where((self.my_locs[:, 0] >=x_min or self.my_locs[:, 0] <= x_max) and
                                         self.my_locs[:, 1] >= y_min and self.my_locs[:, 1] <= y_max and
                                         self.my_locs[:, 2] >= z_min and self.my_locs[:, 2] <= z_max)
            else:
                self.loc_inds = np.where((self.my_locs[:, 0] >=x_min and self.my_locs[:, 0] <= x_max) and
                                         self.my_locs[:, 1] >= y_min and self.my_locs[:, 1] <= y_max and
                                         self.my_locs[:, 2] >= z_min and self.my_locs[:, 2] <= z_max)        
        

        print(self.my_locs[:, 0].size)
        print(self.my_locs[:, 0])
        print(np.where((self.my_locs[:, 1] >= y_min)) - (self.my_locs[:, 1] <= y_max).all().size)
        if (x_min == x_max):
            self.loc_inds = np.where(self.my_locs[:, 1] >= y_min and self.my_locs[:, 1] <= y_max and
                                     self.my_locs[:, 2] >= z_min and self.my_locs[:, 2] <= z_max)
        else:
            if (x_min > x_max):
                #this line does not follow the sames structure as in matlab, but seems likely to work
                self.loc_inds = np.where((self.my_locs[:, 0] >=x_min or self.my_locs[:, 0] <= x_max) and
                                         self.my_locs[:, 1] >= y_min and self.my_locs[:, 1] <= y_max and
                                         self.my_locs[:, 2] >= z_min and self.my_locs[:, 2] <= z_max)
            else:
                self.loc_inds = np.where((self.my_locs[:, 0] >=x_min and self.my_locs[:, 0] <= x_max) and
                                         self.my_locs[:, 1] >= y_min and self.my_locs[:, 1] <= y_max and
                                         self.my_locs[:, 2] >= z_min and self.my_locs[:, 2] <= z_max)

        self.loc_inds = np.array(self.loc_inds)
        self.num_obs = self.loc_inds.size
        self.lons = self.my_locs[self.loc_inds, 0]
        self.lats = self.my_locs[self.loc_inds, 1]
        self.z = self.my_locs[self.loc_inds, 2]
        self.obs = self.my_obs[self.loc_inds]
        self.vert_type = self.my_vert_types[self.loc_inds]
        self.keys = self.my_keys[self.loc_inds]
        self.time = self.my_times[self.loc_inds]

        
        if self.my_QC.size > 0:
            self.QC = self.my_QC[self.loc_inds]
        else:
            self.QC = np.array([])

        
        self.vert_types = np.unique(self.vert_type)
        self.num_vert_types = self.vert_types.size

        for vert_type in self.vert_types:
            print(vert_type)
            if vert_type == -2: #Undefined vertical units, but positive is up
                self.vert_pos_dir = 'up'
                self.vert_units = 'unknown'
            elif vert_type == -1: #Vertical units are surface
                self.vert_pos_dir = 'up'
                self.vert_units = 'surface'
            elif vert_type == 1: #Vertical units are level where 1 is uppermost
                self.vert_pos_dir = 'down'
                self.vert_units = 'model level'
            elif vert_type == 2: #Vertical units are pressure
                self.vert_pos_dir = 'down'
                self.vert_units = 'pressure'
            elif vert_type == 3: #Vertical units are heights
                self.vert_pos_dir = 'up'
                self.vert_units = 'height'

                ''' '''A WORLD OCEAN DATABASE CAVEAT THAT I HAVE NOT IMPLEMENTED''' '''

        #print(self.vert_pos_dir, self.vert_units)
        #print(self.z)

'''
        
