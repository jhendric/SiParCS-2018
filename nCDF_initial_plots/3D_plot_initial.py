import xarray as xa
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

'''

A class for plotting 2D observation data from DART's obs_sequence_tool. 
Lifted variables names from read_obs_netcdf.m and plot_obs_netcdf.m.


Possible 1-off indexing issues. None obvious right now.

Also string format is weird in nCDF (b's before each entry in large string sets)

Also may have to fiddle with datatypes in initializing args, as well as make sure
inf is float('inf')

'''

class 2D_plot:

    def __init__(self, file_path, obs_type_string, region, copy_string, QC_string, max_QC, verbose):
        
        self.dataset = xa.open_dataset(file_path)
        self.obs_type_string = obs_type_string #e.g. 'RADIOSONDE_U_WIND_COMPONENT', 'ALL'
        self.region = np.array(region) #e.g. [0, 360, -90, 90, -inf, inf] or [0, 360, -90, 90]
        
        if copy_string.lower() is 'all':
            self.copy_string = self.dataset['CopyMetaData'].values
        else:
            self.copy_string = copy_string
            
        self.QC_string = QC_string #e.g. 'DART quality control'
        self.max_QC = max_QC #e.g. 2
        self.verbose = verbose #True or False

        self.obs_types = self.dataset['ObsTypes']
        self.obs_type_strings = self.dataset['ObsTypesMetaData'].values
        self.copy_strings = self.dataset['CopyMetaData'].values
        self.QC_strings = self.dataset['QCMetaData'].values

        self.time = self.dataset['time'] #may need to change these three lines to .values
        self.obs_type = self.dataset['obs_type']
        self.obs_keys = self.dataset['obs_keys']
        self.vert_type = self.dataset['which_vert']
        
        '''
        MISSING MISSING DATA ACCOUNTABILITY HERE
        '''
        self.location = self.dataset['location']
        self.obs = self.dataset['observations']
        self.QC = self.dataset['qc']

        self.my_types = np.unique(self.obs_type.values)
        self.time_units = self.time.dtype
        self.time_range = self.time.attrs['valid_range']
        self.calendar = 'Gregorian' #don't know how to get other calendar types
        #missing some other time stuff

        #don't have the verbose section

        self.num_copies = self.dataset['copy'].size

        if copy_string.lower() is 'all':
            self.my_type_ind = [i for i in range(0:self.num_copies)]
        else:
            #get_copy_index.m
            copy_meta_data = self.dataset['CopyMetaData']
            meta_data_length  = self.dataset['nlines'].size
            copy_string_no_white = ''.join(self.copy_string.split())
            copy_index = -1

            #find matching copy
            for i in range(0, self.num_copies):
                meta_data_no_white = ''.join(str(copy_meta_data.values[i].split()))
                if meta_data_no_white is copy_string_no_white:
                    copy_index = i
                    break

            if copy_index < 0:
                raise Exception('No matching copy found')

            self.my_type_ind = copy_index

        #find observations of the correct type
        if obs_type_string.lower() is 'all':
            self.inds = [i for i in range(0:self.obs.size)]

        else:
            self.my_ind_array = np.where(self.obs_type_strings == self.obs_type_string)[0]

            if self.my_ind_array.size > 1:
                self.my_ind = self.my_ind_array[0] + 1 #first instance of my_ind
                self.inds = np.where(self.obs_type.values == self.my_ind)[0] #all instances of associated type

            else:
                print('FYI - no {} observations \n'.format(self.obs_type_string))
                inds = np.array([])

        #this may be incorrect. it will currently be an inds by my_type_end array.
        #similar concerns for subsequent assignments
        self.my_obs = self.obs[self.inds, self.my_type_ind].values
        self.my_locs = self.location[self.inds].values
        self.my_keys = self.obs_keys[self.inds].values
        self.my_vert_types = self.which_vert[self.inds].values
        self.my_times = self.time[self.inds].values

        #find QC values for obs

        if self.QC_string is not '':
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
                if meta_data_no_white is copy_string_no_white:
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
        if region.size == 6:
            z_min = min(region[[4,5]])
            z_max = max(region[[4,5]])
            
        elif region.size < 4:
            raise Exception('Region msut be an array of length 4 or 6')

        else:
            z_min = float(-inf)
            z_max = float(inf)
        
        y_min = min(region[[2,3]])
        y_max = max(region[[2,3]])
        x_min = region[0]%360.0
        x_max = region[1]%360.0

        #may not be possible with numpy array, may need some numpy function
        lons = self.my_locs[:, 0]%360.0

        #find correct latitudes and levels
        lat_logical = np.where(self.my_locs[1] >= y_min and self.my_locs[1] <= y_max) #may not work, see command above
        lvl_logical = np.where(self.my_locs[2] >= z_min and self.my_locs[2] <= z_max) #see above

        #Need to revise the above code to no longer be a logical (should be a way to do it as one where)
            

        
                
            
