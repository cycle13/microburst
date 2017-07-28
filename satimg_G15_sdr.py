#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 10:32:09 2017

@author: ken.pryor
"""
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import urllib

def read_BTD_plot(ncf):
    nc_fid = Dataset(ncf, 'r')
    data = nc_fid.variables["data"][:]  # shape lat, lon as shown above
    lats = nc_fid.variables['lat'][:]  # extract/copy the data
    lons = nc_fid.variables['lon'][:]
    np.savez('btd_3ch.npz', btd=data, lats=lats, lons=lons)
    names = nc_fid.variables.keys()
    nc_fid.close()
    return data, lats, lons, names

"""
WV_IR_BTD_file, headers = urllib.request.urlretrieve('ftp://ftp.star.nesdis.noaa.gov/pub/smcd/opdb/mbimg/GOES0157.nc')
#WV_IR_BTD_file = 'GOES0157.nc'
wvir_btd, lats, lons, names = read_BTD_plot(WV_IR_BTD_file)
wvir_btd = wvir_btd[0,:,:]
print(names)
print('wvir_btd shape', wvir_btd.shape, wvir_btd)
wvir_btd_max = np.amax(wvir_btd)
wvir_btd_min = np.amin(wvir_btd)
print('wvir_btd min, wvir_btd max', wvir_btd_min, wvir_btd_max)
print('lats shape, lons shape', lats.shape, lons.shape)
"""

BTD_file, headers = urllib.request.urlretrieve('ftp://ftp.star.nesdis.noaa.gov/pub/smcd/opdb/mbimg/GOES0150.nc')
#BTD_file = 'GOES0150.nc'
btd, lats, lons, names = read_BTD_plot(BTD_file)
btd = btd[0,:,:]
print(names)
print('btd shape', btd.shape, btd)
btd_max = np.amax(btd)
btd_min = np.amin(btd)
print('btd min, btd max', btd_min, btd_max)
print('lats shape, lons shape', lats.shape, lons.shape)

# CREATE A MAP
"""
fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(projection='merc',llcrnrlat=30,urcrnrlat=45,llcrnrlon=-120,urcrnrlon=-95) # other options include "cyl" projection
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,wvir_btd,shading='gouraud',cmap=plt.cm.jet,latlon=True) # shading can also be "flat"
#m.drawparallels(np.arange(41.,46.,1.),labels=[True,False,False,False])
#m.drawmeridians(np.arange(-118,-108.,2.),labels=[False,False,False,True])
m.drawcoastlines(color="white")
m.drawcountries(color="white")
m.drawstates(color="white")
lon = -112.65173
lat = 43.59413
x, y = m(lon,lat)
print(x,y)
plt.annotate('I', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data', color='white')
lon2 = -117.174
lat2 = 32.714
x2, y2 = m(lon2,lat2)
print(x2,y2)
plt.annotate('SD', xy=(x2, y2), xycoords='data', xytext=(x2, y2), textcoords='data', color='white')
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('GOES-15 SDR WV-IR BTD')
plt.savefig("G15_WVIR_BTD_Map.png",dpi=500,bbox_inches='tight')
plt.show()
"""

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(projection='merc',llcrnrlat=30,urcrnrlat=45,llcrnrlon=-120,urcrnrlon=-95) # other options include "cyl" projection
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,btd,shading='gouraud',cmap=plt.cm.nipy_spectral,latlon=True) # shading can also be "flat"
#m.drawparallels(np.arange(41.,46.,1.),labels=[True,False,False,False])
#m.drawmeridians(np.arange(-118,-108.,2.),labels=[False,False,False,True])
m.drawcoastlines(color="white")
m.drawcountries(color="white")
lon = -112.65173
lat = 43.59413
x, y = m(lon,lat)
print(x,y)
plt.annotate('I', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data', color='white' )
lon2 = -117.174
lat2 = 32.714
x2, y2 = m(lon2,lat2)
print(x2,y2)
plt.annotate('SD', xy=(x2, y2), xycoords='data', xytext=(x2, y2), textcoords='data', color='white')
m.drawstates(color="white")
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('GOES-15 3CH BTD')
plt.savefig("G15_3CH_BTD_Map.png",dpi=500,bbox_inches='tight')
plt.show()
