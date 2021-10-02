#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Ago 30 08:33:09 2021

@author: masuelli

Funciones de transformacion de coordenadas Geos a LatLon y viceversa
con verificación basada en Pug3

"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs  # Plot maps
import os
os.chdir('/media/masuelli/Seagate Expansion Drive/Conae/GVT/GOES_R/ImDic2018/GOES_ABI_20183470000355_20183480945355/')


def scan2sat(x, y, lon0=-75., Re=6378000., Rp=6356000., h=3600000.):
    '''
    Transforma coordenadas de scaneo geostacionarias x,y en coordenadas cartesianas con origen en el satelite sx,sy,sz
    En base a 5.2.8.1 de PUG3
    Parameters
    ----------
    x : float, float arr numpy.ma.core.MaskedArray
       coordenada horizontal, en radianes 
    y : float, float arr numpy.ma.core.MaskedArray
       coordenada vertical, en radianes. Paralelo al eje terrestre 
        longitud
    lon0 : float
        longitud del satélite y origen del sistema de coordenadas planas
    Re: float 
        radio ecuatorial, en m
    Rp: float 
        radio polar, en m
    h: float 
        altura del satélite respecto de la superficie, en m    
    Returns
    -------
    sx : float, float arr
        coordenada hacia el centro de la tierra
    sy : float, float arr
        coordenada horizontal
    sz : float, float arr
        coordenada vertical
    '''

    # print (str(type(x))[8:-2])

    if str(type(x))[8:-2] != 'numpy.ma.core.MaskedArray' or str(type(y))[8:-2] != 'numpy.ma.core.MaskedArray':
        x = np.ma.MaskedArray(x)
        y = np.ma.MaskedArray(y)
        # print ("cambia el tipo")
    mask = x.mask

    H = Re + h  # radio orbital del satelite
    a = np.sin(x)**2 + np.cos(x)**2 * (np.cos(y)**2 + (np.sin(y) * Re/Rp)**2)
    b = -2 * H * np.cos(x) * np.cos(y)
    c = H**2 - Re**2

    aux = b**2 - 4 * a * c

    rs = np.zeros(aux.shape)

    sx = np.ma.MaskedArray(np.zeros(aux.shape), mask)
    sy = np.ma.MaskedArray(np.zeros(aux.shape), mask)
    sz = np.ma.MaskedArray(np.zeros(aux.shape), mask)

    rs[aux >= 0] = -(b[aux >= 0] + np.sqrt(aux[aux >= 0])) / (2 * a[aux >= 0])

    sx[aux >= 0] = rs[aux >= 0] * np.cos(x[aux >= 0]) * np.cos(y[aux >= 0])
    sy[aux >= 0] = -rs[aux >= 0] * np.sin(x[aux >= 0])
    sz[aux >= 0] = rs[aux >= 0] * np.cos(x[aux >= 0]) * np.sin(y[aux >= 0])

    return sx, sy, sz


def sat2latlon(sx, sy, sz, lon0=-75., Re=6378000., Rp=6356000., h=3600000.):
    '''
    Transforma coordenadas cartesianas con origen en el satelite sx,sy,sz en coordenadas de latitud/longitud
    En base a 5.2.8.1 de PUG3
    Parameters
    ----------
    sx : float, float arr
        coordenada hacia el centro de la tierra
    sy : float, float arr
        coordenada horizontal
    sz : float, float arr
        coordenada vertical
    lon0 : float
        longitud del satélite y origen del sistema de coordenadas planas
    Re: float 
        radio ecuatorial, en m
    Rp: float 
        radio polar, en m
    h: float 
        altura del satélite respecto de la superficie, en m
    Returns
    -------
    lat : float, float arr
        latitud
    lon : float, float arr
        longitud
    '''
    H = Re + h  # radio orbital del satelite
    gr2rad = np.pi/180

    lat = np.arctan((Re/Rp)**2 * sz / np.sqrt((H-sx)**2 + sy**2)) / gr2rad
    lon = lon0 - np.arctan(sy / (H-sx)) / gr2rad
    return lat, lon


def latlon2scan(lat, lon, lon0=-75., Re=6378000., Rp=6356000., h=3600000.):
    '''
    Transforma coordenadas de latitud/longitud a x/y en proyeccion geoestacionaria
    En base a 5.2.8.2 de PUG3
    Parameters
    ----------
    lat : float, float arr
        latitud
    lon : float, float arr
        longitud
    lon0 : float, float arr
        longitud del satélite y origen del sistema de coordenadas planas
    Re: float 
        radio ecuatorial, en m
    Rp: float 
        radio polar, en m
    h: float 
        altura del satélite respecto de la superficie, en m
    Returns
    -------
    x : float, float arr
       coordenada horizontal, en radianes 
    y : float, float arr
       coordenada vertical, en radianes. Paralelo al eje terrestre 
    '''

    H = Re + h  # radio orbital del satelite
    e = (1-(Rp/Re)**2)**.5  # 0.0818191910435 # excentricidad
    gr2rad = np.pi/180

    latc = np.arctan((Rp/Re)**2 * np.tan(lat*gr2rad))

    rc = Rp / (1 - (e * np.cos(latc))**2)**.5

    sx = H - rc * np.cos(latc) * np.cos((lon-lon0)*gr2rad)
    sy = - rc * np.cos(latc) * np.sin((lon-lon0)*gr2rad)
    sz = rc * np.sin(latc)

    s_norm = np.sqrt(sx**2 + sy**2 + sz**2)

    x = np.arcsin(-sy / s_norm)
    y = np.arctan(sz / sx)

    return x, y


def colfil2scan(col, fil, x0, y0, scale):
    '''
    Parameters
    ----------
    col : int, int arr
        columna
    fil : int, int arr
        fila
    x0 : float
        posición del x[0] en radianes
    y0 : float
        coordenada horizontal del primer punto, en radianes. Paralelo al eje terrestre 
    scale : float
        tamaño del pixel en radianes
    Returns
    -------
    x : float, float arr
       coordenada horizontal, en radianes.
    y : float, float arr
       coordenada vertical, en radianes. Paralelo al eje terrestre
    '''
    x = col * scale + x0
    y = -fil * scale + y0
    return x, y


def scan2colfil(x, y, x0, y0, scale, tipo=0):
    '''
    Parameters
    ----------
    x : float, float arr
       coordenada vertical, en radianes
    x : float
       coordenada vertical del primer punto en radianes
    x0 : float
        posición del x[0] en radianes
    y0 : float
        coordenada horizontal del primer punto, en radianes. Paralelo al eje terrestre 

    scale : float
        tamaño del pixel en radianes
    tipo : TYPE, optional
        tipo de salida. The default is 0 para float, 1 para int.
    Returns
    -------
    col : columna
    fil : fila
    '''
    col = (x - x0) / scale
    fil = - (y - y0) / scale
    if tipo == 0:
        return col, fil
    elif tipo == 1:
        return round(col).astype('int'), round(fil).astype('int')
    else:
        print('error de tipo')

# %% ----------------------------------------------------
# Importacion de una imagen para extracción de parámetros


img_name = [
    'OR_ABI-L1b-RadF-M3C13_G16_s20183480615355_e20183480626133_c20183480626195.nc']
print('Importando la imagen: %s' % img_name[0])

imagenobj = Dataset(img_name[0], 'r')

print('Importando las variables de la imagen: %s' % img_name[0])

metadato = imagenobj.variables

altura = metadato['goes_imager_projection'].perspective_point_height
semieje_may = metadato['goes_imager_projection'].semi_major_axis
semieje_men = metadato['goes_imager_projection'].semi_minor_axis
lat_cen = metadato['goes_imager_projection'].latitude_of_projection_origin
lon_cen = metadato['goes_imager_projection'].longitude_of_projection_origin
scale_factor = metadato['x'].scale_factor
offset = np.array([metadato['x'].add_offset, metadato['y'].add_offset])

lat_lon_extent = metadato['geospatial_lat_lon_extent']

extent_lat_lon = np.array([lat_lon_extent.geospatial_westbound_longitude,
                          lat_lon_extent.geospatial_northbound_latitude,
                          lat_lon_extent.geospatial_eastbound_longitude,
                          lat_lon_extent.geospatial_southbound_latitude])

x_image_bounds = metadato['x_image_bounds'][:].data
y_image_bounds = metadato['y_image_bounds'][:].data


extent_x_y = [altura * x_image_bounds[0], altura * x_image_bounds[1],
              altura * y_image_bounds[1], altura * y_image_bounds[0]]

pol = semieje_may*altura/(semieje_may+altura)
ecu = semieje_men*altura/(semieje_may+altura)
p_size = scale_factor * altura

icanal = int(metadato['band_id'][:])
print('Canal %d' % icanal)

imagen = metadato['Rad'][:].data

print('Calibrando la imagen')
if icanal >= 7:
    # Parámetros de calibracion
    fk1 = metadato['planck_fk1'][0]  # DN -> K
    fk2 = metadato['planck_fk2'][0]
    bc1 = metadato['planck_bc1'][0]
    bc2 = metadato['planck_bc2'][0]

    imag_cal = (fk2 / (np.log((fk1 / imagen) + 1)) - bc1) / \
        bc2-273.15  # DN -> C
    Unit = "Temperatura de Brillo [°C]"
else:
    pendiente = metadato['Rad'].scale_factor
    ordenada = metadato['Rad'].add_offset
    imag_cal = imagen*pendiente+ordenada
    Unit = "Radiancia ["+metadato['Rad'].units+"]"

print("Graficando")

plt.figure()

# proyeccion geoestacionaria para Goes16
crs = ccrs.Geostationary(central_longitude=lon_cen, satellite_height=altura)
ax = plt.axes(projection=crs)
# agrega linea de meridianos y paralelos
ax.gridlines(xlocs=[-135 + i*20 for i in range(8)])
ax.coastlines(resolution='10m', color='blue')  # agrega líneas de costa

img = plt.imshow(imag_cal, extent=extent_x_y,
                 vmin=-100., vmax=30., cmap='Greys')

# Add a colorbar
plt.colorbar(img, label='Brightness Temperatures (°C)',
             extend='both', orientation='vertical', pad=0.05, fraction=0.05)

plt.show()


# %% ---------------------------------------------------------
print('Verifica transformaciones particulares')


print('\nFull disk lat lon 2 fil,col')
lat = extent_lat_lon[1]
lon = -75

x, y = latlon2scan(lat, lon, lon_cen, semieje_may, semieje_men, altura)
print("lat= %.3f, lon=%3f" % (lat, lon))
print("x= %f, y=%f" % (x, y))

col, fil = scan2colfil(x, y, offset[0], offset[1], scale_factor, 1)
print(col, fil)


print('\nFull disk lat lon 2 fil,col')
lat = 0
lon = extent_lat_lon[2]

x, y = latlon2scan(lat, lon, lon_cen, semieje_may, semieje_men, altura)
print("lat= %.3f, lon=%3f" % (lat, lon))
print("x= %f, y=%f" % (x, y))

col, fil = scan2colfil(x, y, offset[0], offset[1], scale_factor, 1)
print(col, fil)

print('\nFull disk fil,col 2 lat,lon')

fil = 2711
col = 5423

x, y = colfil2scan(col, fil, offset[0], offset[1], scale_factor)

print(x, y)

sx, sy, sz = scan2sat(x, y, lon_cen, semieje_may, semieje_men, altura)
print("sx= %.0f ,sy= %.0f, sz= %.0f" % (sx, sy, sz))

llat, llon = sat2latlon(sx, sy, sz, lon_cen, semieje_may, semieje_men, altura)

print(llat, llon)

# %% ---------------------------------------------------------
# Loop de verificación de scan a Latlon a scan (array)

x = np.array([.09, .09])
y = np.array([.06, .06])

print('\nPrueba loop como array')
print(x, y)

sx, sy, sz = scan2sat(x, y, lon_cen, semieje_may, semieje_men, altura)
print(sx, sy, sz)

llat, llon = sat2latlon(sx, sy, sz, lon_cen, semieje_may, semieje_men, altura)
print(llat, llon)

x, y = latlon2scan(llat, llon, lon_cen, semieje_may, semieje_men, altura)
print(x, y)

# %% ---------------------------------------------------------
# Validar con datos de Pug3


print('\nChequeo de Pug3 para Conus')

offset_conus = np.array([-0.101332, 0.128212])  # x, y

x_bounds_conus = np.array([-0.101360, 0.038640])
y_bounds_conus = np.array([0.128240, 0.044240])

fil = 558  # y=0.095340 rad
col = 1539  # x=-0.024052 rad

print(col, fil)

x, y = colfil2scan(col, fil, offset_conus[0], offset_conus[1], scale_factor)

print(x, y)

col, fil = scan2colfil(x, y, offset_conus[0], offset_conus[1], scale_factor, 1)

print(col, fil)

y = 0.095340
x = -0.024052

sx, sy, sz = scan2sat(x, y, lon_cen, semieje_may, semieje_men, altura)
print("sx= %.0f ,sy= %.0f, sz= %.0f" % (sx, sy, sz))

llat, llon = sat2latlon(sx, sy, sz, lon_cen, semieje_may, semieje_men, altura)
print(llat, llon)

x, y = latlon2scan(llat, llon, lon_cen, semieje_may, semieje_men, altura)
print("x= %f, y=%f" % (x, y))


print('\nOrigen Conus según Full Disk')
col0, fil0 = scan2colfil(
    offset_conus[0], offset_conus[1], offset[0], offset[1], scale_factor, 0)

print(col0, fil0)

# %% ---------------------------------------------------------
# Generación de matrices de Lat y de Lon para cada punto scan (con dato)

n = imagen.shape[0]

mask = metadato['Rad'][:].mask

aux = np.tile(np.arange(n), (n, 1))

Col = np.ma.MaskedArray(aux, mask)
Fil = np.ma.MaskedArray(aux.T, mask)

X, Y = colfil2scan(Col, Fil, offset[0], offset[1], scale_factor)

SX, SY, SZ = scan2sat(X, Y, lon_cen, semieje_may, semieje_men, altura)

Lat, Lon = sat2latlon(SX, SY, SZ, lon_cen, semieje_may, semieje_men, altura)
# %%
print(latlon2scan(10., -85.))

X, Y = colfil2scan(Col, Fil, offset[0], offset[1], scale_factor)

SX, SY, SZ = scan2sat(X, Y, lon_cen, semieje_may, semieje_men, altura)

Lat, Lon = sat2latlon(SX, SY, SZ, lon_cen, semieje_may, semieje_men, altura)
# %%
print(latlon2scan(10., -85.))
