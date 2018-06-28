import xarray as xa
import matplotlib
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cartopy.feature
from cartopy.mpl.patch import geos_to_path
import numpy as np
import pandas as pd
import math
from plot_2D_obs_initial import plot_2D_obs
np.set_printoptions(threshold = np.nan)
import time
import itertools


'''

Initial attempt at 3D plotting using cartopy.

'''

plotter = plot_2D_obs('../obs_series/obs_epoch_001.nc')

fig = plt.figure()
a = ccrs.PlateCarree().transform_points(ccrs.PlateCarree(), plotter.data.lons[1:1000].values, plotter.data.lats[1:1000].values)
lons, lats = a[:, 0], a[:, 1]
ax = Axes3D(fig, xlim = [-180, 180], ylim = [-90, 90])

concat = lambda iterable: list(itertools.chain.from_iterable(iterable))

target_projection = ccrs.PlateCarree()

feature = cartopy.feature.NaturalEarthFeature('Physical', 'land', '110m')

geoms = feature.geometries()

geoms = [target_projection.project_geometry(geom, feature.crs) for geom in geoms]

paths = concat(geos_to_path(geom) for geom in geoms)

polys = concat(path.to_polygons() for path in paths)

lc = PolyCollection(polys, edgecolor = 'black', facecolor = 'green', closed = False)

ax.add_collection3d(lc)

ax.scatter(lons, lats,
           plotter.data.z[1:1000], c = plotter.data.qc_DART[1:1000], cmap = plt.get_cmap('gist_ncar', 9))

ax.set_xlabel('lon')
ax.set_ylabel('lat')
ax.set_zlabel('Height')
           
plt.show()