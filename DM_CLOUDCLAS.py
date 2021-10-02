# Librerias generales
import sys
from scipy import interpolate
import numpy as np
import pprint
import pandas as pd
import time as t
# %%%%%%%

# PyHDF
import pyhdf
from pyhdf.HDF import *
from pyhdf.VS import *
from pyhdf.SD import SD, SDC
# %%
# Librerias geo
import cartopy
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cf
# %%
from netCDF4 import Dataset
# Librerias para calcular zenith
# %%
import calendar
import logging
from datetime import datetime
# %%
from pyorbital import astronomy
from pyspectral.near_infrared_reflectance import Calculator
# Librerias para graficar
import matplotlib.pyplot as plt
import seaborn as sns

# %%


class GoesClass:
    """
    documentacion
    """

    def __init__(self, file_path):
        self.file_path = file_path
        # guarda desde el nivel L1 o L2
        file_name = self.file_path.split('OR_ABI-')[1]
        start_date = self.file_path.split("s20", 1)[1].split("_", 1)[0]
        self.julian_date = start_date[:5]
        self.gregorian_date = datetime.strptime(
            self.julian_date, '%y%j').date().strftime('%d-%m-%y')
        self.utc_hour = start_date[5:9]

    def __repr__(self):
        return f"GOES obj. Date: {self.gregorian_date}; {self.utc_hour} UTC "

    def recorte(self, filas=1440, columnas=1440, x0=-555469.8930323641, y0=0.0):

        # lat =  0. -> y0
        # lon = -80. -> x0
        """
        Funcion que recorta una imagen tipo CMI de GOES.

        Parameters
        ----------
        data_path: str. 
        Direccion de los datos GOES.
        filas: int.
        Cantidad de filas de pixeles (largo) que tendrá la imagen recortada
        columnas: int.
        Cantidad de columnas de pixeles (ancho) que tendrá la imagen recortada
        x0: float. 
        Coordenada x en sistema geoestacionario GOES del limite superior izquierdo en m.
        y0: float. 
        Coordenada y en sistema geoestacionario GOES del limite superior izquierdo en m.

        Returns
        -------
        im_rec: matriz con los elementos del recorte

        """
        psize = 2000
        N = 5424  # esc da 1
        data = Dataset(self.file_path)  # Abro el archivo netcdf
        metadato = data.variables  # Extraigo todas las variables
        banda = metadato['band_id'][:].data[0]  # Extraigo el nro de banda
        # Extraigo la imagen y la guardo en un array de np
        image = np.array(metadato['CMI'][:].data)

        if int(banda) == 3:
            esc = 0.5
            # escala es 1/2 porque tamaño de pixel de banda 3 = 1 km y tamaño pixel del resto = 2 km
            x = range(0, 10848)
            y = range(0, 10848)
            f = interpolate.interp2d(x, y, image, kind='cubic')
            xnew = np.arange(x[0], x[-1], (x[1]-x[0])/esc)
            ynew = np.arange(y[0], y[-1], (y[1]-y[0])/esc)
            image = f(xnew, ynew)

        # tamaño del recorte en proyeccion goes
        #img_extentr = [x0, x0+columnas*psize, y0 -filas*psize, y0]

        esc = int(N/image.shape[0])
        Nx = int(columnas/esc)  # numero de puntos del recorte en x
        Ny = int(filas/esc)  # numero de puntos del recorte en x
        f0 = int((-y0/psize+N/2-1.5)/esc)  # fila del angulo superior izquierdo
        # columna del angulo superior izquierdo
        c0 = int((x0/psize+N/2+.5)/esc)
        f1 = int(f0+Ny)  # fila del angulo inferior derecho
        c1 = int(c0+Nx)  # columna del angulo inferior derecho
        #print('coordenadas filas, col: ', f0, c0, f1, c1)

        im_rec = image[f0:f1, c0:c1]
        return im_rec

    def reflactance(self, rec07, rec13, latlon_extent=[-80, -30, -50, 0]):
      # Correccion del zenith
        lat = np.linspace(latlon_extent[3], latlon_extent[1], rec07.shape[0])
        lon = np.linspace(latlon_extent[0], latlon_extent[2], rec07.shape[1])
        zenith = np.zeros((rec07.shape[0], rec07.shape[1]))
        # Calculate the solar zenith angle
        utc_time = datetime(2019, 1, 2, int(
            self.utc_hour[:2]), int(self.utc_hour[2:]))
        for x in range(len(lat)):
            for y in range(len(lon)):
                zenith[x, y] = astronomy.sun_zenith_angle(
                    utc_time, lon[y], lat[x])
        refl39 = Calculator(platform_name='GOES-16',
                            instrument='abi', band='ch7')
        data07b = refl39.reflectance_from_tbs(zenith, rec07, rec13)
        return data07b

    def day_microphysicsRGB(self, rec03, rec07, rec13, latlon_extent=[-80, -30, -50, 0]):
        """
        Función que arma una imagen RGB que representa microfísica 
        de día según la guía de la pagina de GOES.

        Parameters
        ----------
        rec03: numpy array
        imágen correctamente procesada de la banda 3
        rec07b: numpy array
        imágen correctamente procesada de la banda 7
        rec13: numpy array
        imágen correctamente procesada de la banda 13
        latlon_extent: list
        lat y longs del recorte
        Returns
        -------
        RGB: numpy array
        Imagen RGB de microfísica de día
        """

        R = rec03  # banda3
        G = rec07  # banda7 con corrección zenith
        B = rec13  # banda13

        # Minimuns and Maximuns
        Rmin = 0
        Rmax = 1

        Gmin = 0
        Gmax = 0.6

        Bmin = 203
        Bmax = 323

        # Choose the gamma -> STANDARIZADAS
        gamma_R = 1
        gamma_G = 2.5
        gamma_B = 1

        # Normalize the data
        R = ((R - Rmin) / (Rmax - Rmin)) ** (1/gamma_R)
        G = ((G - Gmin) / (Gmax - Gmin)) ** (1/gamma_G)
        B = ((B - Bmin) / (Bmax - Bmin)) ** (1/gamma_B)

        # Normalizamos (matplotlib lo normaliza de todas formas)
        RR = np.copy(R)
        BB = np.copy(B)
        GG = np.copy(G)

        RR[RR < 0] = 0.
        RR[RR > 1] = 1.
        BB[BB < 0] = 0.
        BB[BB > 1] = 1.
        GG[GG < 0] = 0.
        GG[GG > 1] = 1.

    # Create the RGB
        RGB = np.stack([R, G, B], axis=2)
        # el axis está para que el shape sea fil col dim y no dim col fil
        RRGB = np.stack([RR, GG, BB], axis=2)
        print(RGB.shape)
        return RRGB
# %%


goes3 = GoesClass(
    "C:\\Users\\Paula\\Downloads\\GOES_L2\\2019-003\\OR_ABI-L2-CMIPF-M3C03_G16_s20190031700362_e20190031711129_c20190031711199.nc")

goes3

banda3rec = goes3.recorte()

goes7 = GoesClass(
    "C:\\Users\\Paula\\Downloads\\GOES_L2\\2019-003\\OR_ABI-L2-CMIPF-M3C07_G16_s20190031700362_e20190031711140_c20190031711200.nc")


banda7rec = goes7.recorte()

goes13 = GoesClass(
    "C:\\Users\\Paula\\Downloads\\GOES_L2\\2019-003\\OR_ABI-L2-CMIPF-M3C13_G16_s20190031700362_e20190031711140_c20190031711224.nc")


banda13rec = goes13.recorte()

banda7cor = goes7.reflactance(banda7rec, banda13rec)

rgb = goes3.day_microphysicsRGB(banda3rec, banda7cor, banda13rec)
# %%
plt.figure(1)
plt.imshow(rgb)

# %%


def mask(rgb):
    """This function returns a labeled-by-color image according to the interpretation
    of the product Day Microphysics 
    (https://weather.msfc.nasa.gov/sport/training/quickGuides/rgb/QuickGuide_DtMicroRGB_NASA_SPoRT.pdf)

    Parameters:
    -----------
    rgb: numpy array
    Numpy Array object containig the Day Microphysics product

    Returns:
    -------
    img_mask: numpy array

    """

    img_mask = np.zeros(rgb.shape)

    # Large drops, Low clouds-> pink/magenta
    lc_rfilter = rgb[:, :, 0] > 0.7  # R>0.8
    lc_gfilter = rgb[:, :, 1] < 0.4  # G
    lc_bfilter = rgb[:, :, 2] > 0.6  # B
    lc_filter = lc_rfilter*lc_gfilter*lc_bfilter
    #Mask= magenta
    img_mask[lc_rfilter, 0] = 1.
    img_mask[lc_rfilter, 1] = 0.
    img_mask[lc_rfilter, 2] = 1.

    # Stratus/Stratoculumus (small drops, low clouds) -> bright green/blue
    st_rfilter = (rgb[:, :, 0] > 0.3)*(rgb[:, :, 0] < 0.45)  # R
    st_gfilter = (rgb[:, :, 1] > 0.5)*(rgb[:, :, 1] < 0.8)  # G
    st_bfilter = (rgb[:, :, 2] < 0.7)
    st_filter = st_rfilter*st_gfilter*st_bfilter
    # Mask=Light blue
    img_mask[st_filter, 0] = 0.
    img_mask[st_filter, 1] = 1.
    img_mask[st_filter, 2] = 1.

    # CumuloNimbis (high clouds) -> red, dark orange
    cb_rfilter = rgb[:, :, 0] > 0.7  # R>0.8
    cb_gfilter = rgb[:, :, 1] < 0.3  # G
    cb_bfilter = rgb[:, :, 2] < 0.3  # B
    cb_filter = cb_rfilter*cb_gfilter*cb_bfilter
    # Mask=Red
    img_mask[cb_filter, 0] = 1.
    img_mask[cb_filter, 1] = 0.
    img_mask[cb_filter, 2] = 0.

    # Cirrus (high cliuds)-> green, dark green
    cr_rfilter = rgb[:, :, 0] < 0.3  # R
    cr_gfilter = rgb[:, :, 1] > 0.7  # G
    cr_bfilter = rgb[:, :, 2] < 0.3  # B
    cr_filter = cr_rfilter*cr_gfilter*cr_bfilter
    #Mask= Green
    img_mask[cr_filter, 0] = 0.
    img_mask[cr_filter, 1] = 1.
    img_mask[cr_filter, 2] = 0.

    # supercooled clouds Thick, small drops, medium clouds-> yellow
    super_rfilter = rgb[:, :, 0] > 0.8
    super_gfilter = rgb[:, :, 1] > 0.8
    super_bfilter = rgb[:, :, 2] < 0.2  # amarillo
    super_filter = super_rfilter*super_gfilter*super_bfilter
    # Mask=Yellow
    img_mask[super_filter, 0] = 1.
    img_mask[super_filter, 1] = 1.
    img_mask[super_filter, 2] = 0.

    return img_mask[:, :, [0, 1, 2]]


# %%
mascara = mask(rgb)
crs = ccrs.Geostationary(central_longitude=-75.0, satellite_height=35786023.0)
img_extentr = [-555469.8930323641, 2324530.106967636, -2880000.0, 0.0]

fig, axis = plt.subplots()
axis = plt.axes(projection=crs)
axis.gridlines()
axis.coastlines(resolution='10m', color='white')
axis.add_feature(cf.BORDERS)
axis.add_feature(cf.LAKES)
axis.add_feature(cf.RIVERS)
plt.imshow(mascara, origin='upper', extent=img_extentr)

# %%


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


# %%
print(latlon2scan(10., -85.))


# %%
print(latlon2scan(10., -85.))
