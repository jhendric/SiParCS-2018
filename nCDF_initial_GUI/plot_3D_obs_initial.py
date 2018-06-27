import xarray as xa
import matplotlib
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import math
from plot_2D_obs_initial import plot_2D_obs
np.set_printoptions(threshold = np.nan)
import time
import itertools
import cartopy.feature
from cartopy.mpl.patch import geos_to_path
import cartopy.crs as ccrs
from matplotlib.collections import LineCollection

'''

Initial attempt at 3D plotting using cartopy.

'''

plotter = plot_2D_obs('../obs_series/obs_epoch_001.nc')

fig = plt.figure()
a = ccrs.PlateCarree().transform_points(ccrs.PlateCarree(), plotter.data.lons[1:1000].values, plotter.data.lats[1:1000].values)
lons, lats = a[:, 0], a[:, 1]
ax = Axes3D(fig)
#ax = Axes3D(fig, xlim = [min(plotter.data.lons.values), max(plotter.data.lons.values)], ylim = [-90, 90], zlim = [0, 100000])
#ax = fig.add_subplot(111, projection = 'PlateCarree', xlim = [-180, 180], ylim = [-90, 90])
#ax.set_zlim(bottom = 0)


#print(a)

target_projection = ccrs.PlateCarree()

#ccrs.PlateCarree().transform_points(plotter.data.lons[1:10000], plotter.data.lats[1:10000])

feature = cartopy.feature.NaturalEarthFeature('Physical', 'coastline', '110m')

geoms = feature.geometries()

geoms = [target_projection.project_geometry(geom, feature.crs) for geom in geoms]

paths = list(itertools.chain.from_iterable(geos_to_path(geom) for geom in geoms))

#print(paths)

segments = []
for path in paths:
    vertices = [vertex for vertex, _ in path.iter_segments()]
    vertices = np.asarray(vertices)
    segments.append(vertices)


lc = LineCollection(segments, color = 'black')

ax.add_collection3d(lc)

#X, Y, Z = np.meshgrid(plotter.data.lons[1:10000].values, plotter.data.lats[1:10000].values, plotter.data.z[1:10000].values)

#print(plotter.data.lons[1:10000].values)

ax.scatter(lons, lats,
           plotter.data.z[1:1000], c = "green")
#ax.scatter(X, Y, Z, transform = ccrs.PlateCarree())

ax.set_xlabel('lon')
ax.set_ylabel('lat')
ax.set_zlabel('Height')

plt.show()
