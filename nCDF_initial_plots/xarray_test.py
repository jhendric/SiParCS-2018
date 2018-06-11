import xarray
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

'''

Work in reading and plotting netCDF files. 
Managed to plot 2D and 3D data from bgrid and cam sample files.

'''

#inno = xarray.open_dataset('Innovation.nc')
#init = xarray.open_dataset('preassim.nc')
#post = xarray.open_dataset('analysis.nc')
#cam_diag = xarray.open_dataset('../obs_series/obs_diag_output.nc')

cam_data = xarray.open_dataset('../obs_series/obs_epoch_001.nc')

#X, Y, Z = np.meshgrid(temp.c
print(cam_data)

quit()


#cam = xarray.open_dataset('caminput.nc')
#print(cam)

combined = xarray.concat([init, post, inno], dim = 'step')
temp = combined['t_mean']
#cam_temp = cam['T']

#print(temp.coords['TmpI'].values)
#print(temp.coords['TmpJ'].values)
#print(temp.coords['lev'].values)

#print(temp.values)


X, Y, Z = np.meshgrid(temp.coords['TmpI'].values, temp.coords['TmpJ'].values,
                         temp.coords['lev'].values)


T = np.swapaxes(np.swapaxes(temp.values[0,0,:,:,:], 0, 2), 0, 1)
#print(X.shape)
#print(Y.shape)
#print(Z.shape)
#print(T.shape)

fig, _ = plt.subplots()

ax = fig.add_subplot(111, projection='3d')
ax.scatter(X, Y, Z, c = T.flatten())

#print(temp.coords['TmpI'].values)

#print(cam_temp)

'''
fig2, ax2 = plt.subplots()

ax2 = plt.axes(projection = ccrs.PlateCarree())

cam_temp[0, 25, :, :].plot(ax = ax2, transform = ccrs.PlateCarree())
                           
ax2.gridlines()
ax2.coastlines()
'''


'''
p3 = cam_temp[0, 21:27, :, :].plot(transform = ccrs.PlateCarree(),
                           subplot_kws = {'projection': ccrs.PlateCarree()},
                           x = 'lon', y = 'lat',
                           col = 'lev')
for ax in p3.axes.flat:
    ax.coastlines()
    ax.gridlines()

'''

'''
p1 = temp[[0, 1], 0, :, :, :].plot(transform = ccrs.PlateCarree(),
                                  subplot_kws = {'projection': ccrs.PlateCarree()},
                                  x = 'TmpI', y = 'TmpJ',
                                 col = 'lev', row = 'step')
for ax in p1.axes.flat:
    ax.coastlines()
    ax.gridlines()
'''

'''
fig1, _  = plt.subplots()
plt.figure(fig1.number)
ax = plt.axes(projection = ccrs.PlateCarree())
temp[2, 0, 0, :, :].plot(ax = ax, transform = ccrs.PlateCarree())
ax.gridlines()
ax.coastlines()
'''

'''
p2 = temp[2, 0, :, :, :].plot(transform  = ccrs.PlateCarree(),
                             subplot_kws = {'projection': ccrs.PlateCarree()},
                             x = 'TmpI', y = 'TmpJ',
                             col = 'lev')

for ax in p2.axes.flat:
    ax.coastlines()
    ax.gridlines()
'''

plt.show()
