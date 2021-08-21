# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 10:29:26 2021
@author: PRJure
"""

#Librerias generales
import os
import sys
from scipy import interpolate
import numpy as np
import pprint
import pandas as pd
import time as t
#%%%%%%%
#PyHDF
import pyhdf
from pyhdf.HDF import *
from pyhdf.VS import *
from pyhdf.SD import SD, SDC  
#%%%
#Librerias geo
import cartopy
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
#%%
#import geopandas as gpd
#%%
from netCDF4 import Dataset
#Librerias para calcular zenith
#%%
import calendar
import logging
from datetime import datetime
#%%
from pyorbital import astronomy
from pyspectral.near_infrared_reflectance import Calculator
#%%
#Recorte
def recorte(data_path,x0 = -555469.8930323641, y0 = 0.0):
  """
  Funcion que recorta la imagen CMI de GOES.
  Parameters
  ----------
  data_path: str. 
  Direccion de los datos GOES.
  x0: float. 
  Coordenada x del limite superior izquierdo en m.
  y0: float. 
  Coordenada y del limite superior izquierdo en m.
  
  Returns
  -------
  im_rec: matriz con los elementos del recorte
  
  Example
  -------
  long-75.194122	lat -0.033069	0	POINT (-21609.527 -3656.536)
  x=-555469.8930323641, y=0.0
  """
  data = Dataset(data_path)
  metadato = data.variables
  banda = metadato['band_id'][:].data[0]
  print('banda =', banda)

  #Parámetros para el recorte
  filas = 1440 # filas del recorte para la de referencia
  columnas = 1440 # filas del recorte para la de referencia
  
  if int(banda)==3:
    psize= 1000
    N = 5424*2 #esc da 1
  else: #bandas 7 y 13
    psize = 2000 #tamaño del pixel en km
    N = 5424 #esc da 1
  
  img_extentr = [x0,x0+columnas*psize,y0-filas*psize,y0] #en proyeccion goes
  print('extent rec en proyeccion goes:', img_extentr)
  esc= int(N/metadato['CMI'][:].shape[0])
  Nx = int(columnas/esc) #numero de puntos del recorte en x
  Ny = int(filas/esc) #numero de puntos del recorte en x
  f0 = int((-y0/psize+N/2-1.5)/esc) #fila del angulo superior izquierdo
  c0 = int((x0/psize+N/2+.5)/esc) #columna del angulo superior izquierdo
  f1 = int(f0+Ny) #fila del angulo inferior derecho
  c1 = int(c0+Nx) #columna del angulo inferior derecho
  print('coordenadas filas, col: ', f0,c0,f1,c1)
  im_rec = metadato['CMI'][:].data[f0:f1,c0:c1]
  return im_rec
#%%
#Hago los recortes
rec03 = recorte('C:/Users/Daniel/OR_ABI-L2-CMIPF-M3C03_G16_s20190021800363_e20190021811129_c20190021811205.nc')
rec07 = recorte("C:/Users/Daniel/OR_ABI-L2-CMIPF-M3C07_G16_s20190021800363_e20190021811141_c20190021811202.nc")
rec13 = recorte("C:/Users/Daniel/OR_ABI-L2-CMIPF-M3C13_G16_s20190021800363_e20190021811141_c20190021811221.nc")
#%%
#Comprobamos tamaños
print(rec03.shape, rec07.shape, rec13.shape)
#%%
#Calculamos zenith para banda 7
def solar_7(ch7,ch13, latlon_extent):
  """"
  Parameters
  ----------
  ch7: matriz (recortada) del canal 7
  ch13: matriz (recortada) del canal 13
  latlon_extent: list
  Lista [x1,y1,x2,y2] de los bordes de la imagen en latitud, longitud donde
      x1=longitud de más al oeste
      y1=latitud de más al sur (punto y inferior)
      x2 = longitud de más al este
      y2=latitud de más al norte (punto y superior)

  Returns
  -------
  data2b: matriz con el cálculo de zenith pixel a pixel
  """
#Calculo del ángulo del sol para banda 7
  lat = np.linspace(latlon_extent[3], latlon_extent[1], ch7.shape[0])
  lon = np.linspace(latlon_extent[0], latlon_extent[2], ch7.shape[1])
  print(lat.shape)
  print(lon.shape)
  zenith = np.zeros((ch7.shape[0], ch7.shape[1]))
  # Calculate the solar zenith angle
  utc_time = datetime(2019, 1, 2, 18, 3) 
  for x in range(len(lat)):
    for y in range(len(lon)):
      zenith[x,y] = astronomy.sun_zenith_angle(utc_time, lon[y], lat[x])
  refl39 = Calculator('GOES-16', 'abi', 'ch7')
  data2b = refl39.reflectance_from_tbs(zenith, ch7, ch13)
  return data2b
#%%
latlon_extent=[-80,-30,-50,0] #lat y longs del recorte
rec07b = solar_7(rec07, rec13, latlon_extent)
rec07b.shape
#%%
# RGB Components
R = rec03 #banda3
G = rec07b #banda7
B = rec13 #banda9
 
# Minimuns and Maximuns
Rmin = 0
Rmax = 1
 
Gmin = 0
Gmax = 0.6
 
Bmin = 203
Bmax = 323
 
#R[R.max] = Rmax
 
#G[G.max] = Gmax
 
#B[B.max] = Bmax
# Choose the gamma -> STANDARIZADAS
gamma_R = 1
gamma_G = 2.5
gamma_B = 1
 
# Normalize the data
R = ((R - Rmin) / (Rmax - Rmin)) ** (1/gamma_R)
G = ((G - Gmin) / (Gmax - Gmin)) ** (1/gamma_G)
B = ((B - Bmin) / (Bmax - Bmin)) ** (1/gamma_B) 
 
# Create the RGB
RGB = np.stack([R, G, B], axis=2)
print(RGB.shape)
#%%
#ploteo
import matplotlib.pyplot as plt
crs=ccrs.Geostationary(central_longitude=-75.0, satellite_height=35786023.0)
fig_dims = (10, 10)
fig, axis = plt.subplots(figsize=fig_dims)
rec_goes_extent =[-555469.8930323641, 2324530.106967636, -2880000.0, 0.0]
axis = plt.axes(projection=crs)
axis.gridlines()
axis.coastlines(resolution='10m',color='green') 
plt.imshow(RGB, origin='upper', extent=rec_goes_extent, vmin=0., vmax=1.)
#%%
ds1 = Dataset("C:/Users/Daniel/OR_ABI-L2-CMIPF-M3C13_G16_s20190021800363_e20190021811141_c20190021811221.nc")
print(ds1.variables.keys())
#%%%
H = ds1.variables['goes_imager_projection'].perspective_point_height
x = ds1.variables['x'][:]
X_goes_proj = x*H
print(X_goes_proj)
