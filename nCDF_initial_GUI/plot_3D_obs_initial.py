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

import cartopy.feature
from cartopy.mpl.patch import geos_to_path
import cartopy.crs as ccrs

'''

Initial attempt at 3D plotting using cartopy.

'''

#plotter = plot_2D_obs('../obs_series/obs_epoch_001.nc')

fig = Figure(figsize = (12, 8))
ax = Axes3D(fig, xlim = [-180, 180], ylim = [-90, 90])
ax.set_zlim(bottom = 0)

target_projection = ccrs.PlateCarree()

feature = cartopy.feature.NaturalEarthFeature('Physical', 'coastline', '110m')

geoms = feature.geometries()

geoms = [target_projection.project_geometry(geom, feature.crs) for geom in geoms]

paths = list(geos_to_path(geom) for geom in geoms)

print(paths)

segments = []
for path in paths:
    vertices = [vertex for vertex, _ in path]
    vertices = np.asarray(vertices)
    segments.append(vertices)


lc = LineCollection(segments, color = 'black')

ax.add_collection3d(lc)

ax.set_xlabel('lon')
ax.set_ylabel('lat')
ax.set_zlabel('Height')

plt.show()
