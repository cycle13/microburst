#!/usr/bin/env python3
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

def read_RAP_500(ncf):
    nc_fid = Dataset(ncf, 'r')

    Z_500 = nc_fid.variables["Geopotential_height_isobaric"][:]  # shape lat, lon as shown above
    RH_500 = nc_fid.variables["Relative_humidity_isobaric"][:]
    T_500 = nc_fid.variables["Temperature_isobaric"][:]
    X = nc_fid.variables["x"][:]
    Y = nc_fid.variables["y"][:]
    lats = nc_fid.variables['lat'][:]  # extract/copy the data
    lons = nc_fid.variables['lon'][:]
    names = nc_fid.variables.keys()
    nc_fid.close()
    return Z_500, RH_500, T_500, X, Y, lats, lons, names

def read_RAP_700(ncf):
    nc_fid = Dataset(ncf, 'r')
    Z_700 = nc_fid.variables["Geopotential_height_isobaric"][:]  # shape lat, lon as shown above
    RH_700 = nc_fid.variables["Relative_humidity_isobaric"][:]
    T_700 = nc_fid.variables["Temperature_isobaric"][:]
    X = nc_fid.variables["x"][:]
    Y = nc_fid.variables["y"][:]
    lats = nc_fid.variables['lat'][:]  # extract/copy the data
    lons = nc_fid.variables['lon'][:]
    names = nc_fid.variables.keys()
    nc_fid.close()
    return Z_700, RH_700, T_700, X, Y, lats, lons, names

def read_RAP_CAPE(ncf):
    nc_fid = Dataset(ncf, 'r')
    CAPE = nc_fid.variables["Convective_available_potential_energy_surface"][:]  # shape lat, lon as shown above
    X = nc_fid.variables["x"][:]
    Y = nc_fid.variables["y"][:]
    lats = nc_fid.variables['lat'][:]  # extract/copy the data
    lons = nc_fid.variables['lon'][:]
    names = nc_fid.variables.keys()
    nc_fid.close()
    return CAPE, X, Y, lats, lons, names

RAP_file_500 = 'RR_CONUS_13km_20170604_2100_500.nc4'
names = read_RAP_500(RAP_file_500)
Z_500, RH_500, T_500, X, Y, lats, lons, names = read_RAP_500(RAP_file_500)
Z_500 = Z_500[0,0,:,:]
Z_500_km = Z_500/1000
RH_500 = RH_500[0,0,:,:]
T_500 = T_500[0,0,:,:]
T_500_C = T_500 - 273.15
#TD_650 = T_650_C - ((100-RH_650)/5)
TD_500 = ((0.198 + (0.0017*T_500_C))*RH_500) + ((0.84*T_500_C)-19.2)
print(names)
print('Z_500 shape', Z_500.shape, Z_500)
print('RH_500 shape', RH_500.shape, RH_500)
print('T_500 shape', T_500.shape, T_500)
print('T_500_C shape', T_500_C.shape, T_500_C)
print('TD_500 shape', TD_500.shape, TD_500)
print('X shape, Y shape', X.shape, Y.shape)
print('lats shape, lons shape', lats.shape, lons.shape)

RAP_file_700 = 'RR_CONUS_13km_20170604_2100_700.nc4'
names = read_RAP_700(RAP_file_700)
Z_700, RH_700, T_700, X, Y, lats, lons, names = read_RAP_700(RAP_file_700)
Z_700 = Z_700[0,0,:,:]
Z_700_km = Z_700/1000
RH_700 = RH_700[0,0,:,:]
T_700 = T_700[0,0,:,:]
T_700_C = T_700 - 273.15
#TD_850 = T_850_C - ((100-RH_850)/5)
TD_700 = ((0.198 + (0.0017*T_700_C))*RH_700) + ((0.84*T_700_C)-19.2)
print(names)
print('Z_700 shape', Z_700.shape, Z_700)
print('RH_700 shape', RH_700.shape, RH_700)
print('T_700 shape', T_700.shape, T_700)
print('T_700_C shape', T_700_C.shape, T_700_C)
print('TD_700 shape', TD_700.shape, TD_700)
print('X shape, Y shape', X.shape, Y.shape)
print('lats shape, lons shape', lats.shape, lons.shape)

RAP_file_CAPE = 'RR_CONUS_13km_20170604_2100_cape.nc4'
names = read_RAP_CAPE(RAP_file_CAPE)
CAPE, X, Y, lats, lons, names = read_RAP_CAPE(RAP_file_CAPE)
CAPE = CAPE[0,:,:]
print(names)
print('CAPE shape', CAPE.shape, CAPE)
print('X shape, Y shape', X.shape, Y.shape)
print('lats shape, lons shape', lats.shape, lons.shape)

def MWPI(Z_500_km, Z_700_km, T_500_C, T_700_C, TD_500, TD_700, CAPE):
        gamma = (T_700_C - T_500_C)/(Z_500_km - Z_700_km)
        DD_500 = T_500_C - TD_500
        DD_700 = T_700_C - TD_700
        DDD = DD_700 - DD_500
        MWPI = (CAPE/100) + gamma + DDD
        MWPIv2 = (CAPE/1000) + (gamma/5) + (DDD/5)
        #WGP = (0.3163 * (MWPI+Tv_c)) + 33.766
        WGP = (0.3163 * MWPI) + 33.766
        WGPw = (5.1532 * MWPIv2) + 27.158
        return gamma, MWPI, MWPIv2, WGP, WGPw
gamma, MWPI, MWPIv2, WGP, WGPw = MWPI(Z_500_km, Z_700_km, T_500_C, T_700_C, TD_500, TD_700, CAPE)
print('gamma shape', gamma.shape, gamma)
print('MWPI shape', MWPI.shape, MWPI)
print('WGP shape', WGP.shape, WGP)
print('MWPIv2 shape', MWPIv2.shape, MWPIv2)
print('WGPwshape', WGPw.shape, WGPw)

# PLOT 2D PLOT
XX, YY = np.meshgrid(X,Y,sparse=True)
print(XX.shape,YY.shape)
XX = XX[0,:]
YY = YY[:,0]

xx, yy = np.meshgrid(lons,lats,sparse=True)
print(xx.shape,yy.shape)
xx = xx[0,:]
yy = yy[:,0]

print(XX.shape,YY.shape)
print(xx.shape,yy.shape)
fig = plt.figure()
fig,ax=plt.subplots()
cf=ax.pcolormesh(XX,YY,gamma,cmap=cm.jet,rasterized=True)
#plt.xlim([-2600,-600])
#plt.ylim([1100,3100])
plt.colorbar(cf)
plt.savefig("gamma_RAP_2100.png",dpi=300)
plt.show()

fig = plt.figure()
fig,ax=plt.subplots()
cf=ax.pcolormesh(XX,YY,MWPI,cmap=cm.jet,rasterized=True)
#plt.xlim([-2600,-600])
#plt.ylim([1100,3100])
plt.colorbar(cf)
plt.savefig("MWPI_RAP_2100.png",dpi=300)
plt.show()

fig = plt.figure()
fig,ax=plt.subplots()
cf=ax.pcolormesh(XX,YY,WGP,cmap=cm.nipy_spectral,rasterized=True)
#plt.xlim([-2600,-600])
#plt.ylim([1100,3100])
plt.colorbar(cf)
plt.savefig("WGP_RAP_2100.png",dpi=300)
plt.show()

fig = plt.figure()
fig,ax=plt.subplots()
cf=ax.pcolormesh(lons,lats,gamma,cmap=cm.jet,rasterized=True)
#plt.xlim([-101,-96])
#plt.ylim([33,36])
plt.colorbar(cf)
plt.savefig("gamma_latlon_2100.png",dpi=300)
plt.show()

fig = plt.figure()
fig,ax=plt.subplots()
cf=ax.pcolormesh(lons,lats,MWPI,cmap=cm.jet,rasterized=True)
#plt.xlim([-101,-96])
#plt.ylim([33,36])
plt.colorbar(cf)
plt.savefig("MWPI_latlon_2100.png",dpi=300)
plt.show()

fig = plt.figure()
fig,ax=plt.subplots()
cf=ax.pcolormesh(lons,lats,WGP,cmap=cm.nipy_spectral,rasterized=True)
#plt.xlim([-101,-96])
#plt.ylim([33,36])
plt.colorbar(cf)
plt.savefig("WGP_latlon_2100.png",dpi=300)
plt.show()


# CREATE A MAP

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(projection='merc',llcrnrlat=35,urcrnrlat=50,llcrnrlon=-123,urcrnrlon=-103) # other options include "cyl" projection
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,gamma,shading='gouraud',cmap=plt.cm.jet,latlon=True) # shading can also be "flat"
m.drawparallels(np.arange(35.,50.,1.),labels=[True,False,False,False])
m.drawmeridians(np.arange(-123.,-103.,5.),labels=[False,False,False,True])
m.drawcoastlines()
m.drawcountries()
m.drawstates()
lon = -112.65173
lat = 43.59413
x, y = m(lon,lat)
print(x,y)
plt.annotate('I', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data')
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('500-700 mb Lapse Rate 2100 UTC 4 June 2017')
plt.savefig("gamma_Map_2100.png",dpi=500)
plt.show()

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(projection='merc',llcrnrlat=35,urcrnrlat=50,llcrnrlon=-123,urcrnrlon=-103) # other options include "cyl" projection
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,CAPE,shading='gouraud',cmap=plt.cm.jet,latlon=True) # shading can also be "flat"
m.drawparallels(np.arange(35.,50.,1.),labels=[True,False,False,False])
m.drawmeridians(np.arange(-123.,-103.,5.),labels=[False,False,False,True])
m.drawcoastlines()
m.drawcountries()
m.drawstates()
lon = -112.65173
lat = 43.59413
x, y = m(lon,lat)
print(x,y)
plt.annotate('I', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data')
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('CAPE 2100 UTC 4 June 2017')
plt.savefig("CAPE_Map_2100.png",dpi=500)
plt.show()

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(projection='merc',llcrnrlat=35,urcrnrlat=50,llcrnrlon=-123,urcrnrlon=-103) # other options include "cyl" projection
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,MWPI,shading='gouraud',cmap=plt.cm.jet,latlon=True) # shading can also be "flat"
m.drawparallels(np.arange(35.,50.,1.),labels=[True,False,False,False])
m.drawmeridians(np.arange(-123.,-103.,5.),labels=[False,False,False,True])
m.drawcoastlines()
m.drawcountries()
m.drawstates()
lon = -112.65173
lat = 43.59413
x, y = m(lon,lat)
print(x,y)
plt.annotate('I', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data')
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('500-700 mb MWPI 2100 UTC 4 June 2017')
plt.savefig("MWPI_Map_2100.png",dpi=500)
plt.show()

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(projection='merc',llcrnrlat=35,urcrnrlat=50,llcrnrlon=-123,urcrnrlon=-103) # other options include "cyl" projection
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,WGP,shading='gouraud',cmap=plt.cm.nipy_spectral,latlon=True) # shading can also be "flat"
m.drawparallels(np.arange(35.,50.,1.),labels=[True,False,False,False])
m.drawmeridians(np.arange(-123.,-103.,5.),labels=[False,False,False,True])
m.drawcoastlines()
m.drawcountries()
m.drawstates()
lon = -112.65173
lat = 43.59413
x, y = m(lon,lat)
print(x,y)
plt.annotate('I', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data')
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('500-700 mb Wind Gust Potential (knots) 2100 UTC 4 June 2017')
plt.savefig("WGP_Map_2100.png",dpi=500)
plt.show()

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(projection='merc',llcrnrlat=35,urcrnrlat=50,llcrnrlon=-123,urcrnrlon=-103) # other options include "cyl" projection
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,MWPIv2,shading='gouraud',cmap=plt.cm.jet,latlon=True) # shading can also be "flat"
m.drawparallels(np.arange(35.,50.,1.),labels=[True,False,False,False])
m.drawmeridians(np.arange(-123.,-103.,5.),labels=[False,False,False,True])
m.drawcoastlines()
m.drawcountries()
m.drawstates()
lon = -112.65173
lat = 43.59413
x, y = m(lon,lat)
print(x,y)
plt.annotate('I', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data')
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('500-700 mb MWPIv2 2100 UTC 4 June 2017')
plt.savefig("MWPIv2_Map_2100.png",dpi=500)
plt.show()

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(projection='merc',llcrnrlat=35,urcrnrlat=50,llcrnrlon=-123,urcrnrlon=-103) # other options include "cyl" projection
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,WGPw,shading='gouraud',cmap=plt.cm.nipy_spectral,latlon=True) # shading can also be "flat"
m.drawparallels(np.arange(35.,50.,1.),labels=[True,False,False,False])
m.drawmeridians(np.arange(-123.,-103.,5.),labels=[False,False,False,True])
m.drawcoastlines()
m.drawcountries()
m.drawstates()
lon = -112.65173
lat = 43.59413
x, y = m(lon,lat)
print(x,y)
plt.annotate('I', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data')
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('WUS Wind Gust Potential (knots) 2100 UTC 4 June 2017')
plt.savefig("WGPw_Map_2100.png",dpi=500)
plt.show()
