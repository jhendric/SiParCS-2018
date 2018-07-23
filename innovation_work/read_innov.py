
import xarray as xa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

import cartopy.crs as ccrs


class GenerateInnov:

    '''

    Reading and differencing functionality between two model restart files

    '''
    
    def __init__(self, initial, final):

        '''Initialize object with two model restart files opened in xarray
        
        Keyword arguments:
        initial -- path of pre-DART model restart file
        final -- path of post-DART model restart file

        '''
        
        self.initial = xa.open_dataset(initial, decode_times = False)
        self.final = xa.open_dataset(final, decode_times = False)

    def diff(self, var_name, level, level_name):
        '''get difference between self.initial and self.final for var_name
        at level in category level_name of the dataset
        
        Keyword arguments:
        var_name -- name of data variable to create an xarray DataArray from
        level -- level of data to access
        level_name -- string providing name of category

        '''

        innov_broad = self.final[var_name]-self.initial[var_name]

        #transfer info from original files that xarray will automatically plot (axis labels, title are included here)
        innov_broad.attrs = self.final[var_name].attrs

        if 'long_name' in innov_broad.attrs:
            innov_broad.attrs['long_name'] = 'Difference in ' + innov_broad.attrs['long_name']
        
        print('innov_broad: ', innov_broad)
        
        if type(level) is not np.float64:
            innov_narrow = innov_broad.where(innov_broad[level_name] == level, drop = True)
        else:
            #avoid float equality issues by checking equality to within 3 decimal places
            innov_narrow = innov_broad.where(abs(innov_broad[level_name] - level) < 1e-3, drop = True)

        '''These plots are currently designed to work with files featuring data from only one time
        If a time dimension exists, the time value for the first data point is chosen for making sure all plotted
        data comes from a single time'''
        if 'time' in innov_narrow.dims or 'time' in innov_narrow.coords:
            innov_narrow = innov_narrow.where(innov_narrow.time == innov_narrow.time.values[0], drop = True)
            
        #level_name will definitely be only 1 possible value and should be dropped from the dimensions (not needed?)
        return innov_narrow
        
    def plot(self, innov):
        '''Test plotting function. Not used in any GUI file

        Keyword arguments:
        innov: xarray DataArray that is the difference between two model restart files

        '''
        
        ax = plt.axes(projection = ccrs.PlateCarree())
        innov.plot(ax = ax, transform = ccrs.PlateCarree(), cmap = 'bwr_r')
        ax.gridlines()
        ax.coastlines()
        
def main():

    '''Initializes an innovation generator and plots from test files'''
    
    gen = GenerateInnov('filter_input.0001.nc', 'filter_output.0001.nc')
    innov = gen.diff('T', 1, 'lev')
    gen.plot(innov)
    plt.show()
    
if __name__ == "__main__":
    main()
