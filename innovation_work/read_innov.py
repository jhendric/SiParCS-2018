
import xarray as xa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import time

import cartopy.crs as ccrs


class GenerateInnov:

    def __init__(self, initial, final):

        self.initial = xa.open_dataset(initial, decode_times = False)
        self.final = xa.open_dataset(final, decode_times = False)

    def diff(self, var_name, level, level_name):
        '''get difference between self.initial and self.final for var_name
        at level in category level_name of the datasets'''

        innov_initial = self.final[var_name]-self.initial[var_name]
        innov_final = innov_initial.where(innov_initial[level_name] == level, drop = True)

        return innov_final
        
    def plot(self, innov):
        
        ax = plt.axes(projection = ccrs.PlateCarree())
        innov.plot(ax = ax, transform = ccrs.PlateCarree(), cmap = 'bwr_r')
        ax.gridlines()
        ax.coastlines()
        
def main():
    
    gen = GenerateInnov('filter_input.0001.nc', 'filter_output.0001.nc')
    #print(gen.initial['T'], gen.final['T'])
    #print((gen.final-gen.initial)['T'])
    innov = gen.diff('T', 1, 'lev')
    print(innov)
    gen.plot(innov)
    plt.show()
    
if __name__ == "__main__":
    main()
