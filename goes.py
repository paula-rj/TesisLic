# Librerias generales
import sys

from scipy import interpolate
import numpy as np
import pprint
import time as t

# Librerias geo
import cartopy
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cf

from netCDF4 import Dataset

# Librerias para calcular zenith
import calendar
import logging
from datetime import datetime

from pyorbital import astronomy
from pyspectral.near_infrared_reflectance import Calculator

import matplotlib.pyplot as plt

import latlon2geos8

# %%


class GoesClass:
    """
    documentacion
    """

    def __init__(self, file_path):
        self.file_path = file_path
        # guarda desde el nivel L1 o L2
        file_name = self.file_path.split("OR_ABI-")[1]
        start_date = self.file_path.split("s20", 1)[1].split("_", 1)[0]
        self.julian_date = start_date[:5]
        self.gregorian_date = (
            datetime.strptime(self.julian_date, "%y%j")
            .date()
            .strftime("%d-%m-%y")
        )
        self.utc_hour = start_date[5:9]

    def __repr__(self):
        return f"GOES obj. Date: {self.gregorian_date}; {self.utc_hour} UTC "

    def recorte(
        self, rows=2500, cols=2500, lat_sup=10., lon_west=-80.
    ):

        # lat =  0. -> y0
        # lon = -80. -> x0
        """
        Funcion que recorta una imagen tipo CMI de GOES.

        Parameters
        ----------
        data_path: str.
        Direccion de los datos GOES.
        rows: int.
        Cantidad de filas de pixeles (largo) que tendrá la imagen recortada
        cols: int.
        Cantidad de columnas de pixeles (ancho) que tendrá la imagen recortada
        lon_west: float.
        Longitud maxima al oeste del recorte
        lat_sup: float.
        Latitud superior del recorte.

        Returns
        -------
        im_rec: matriz con los elementos del recorte

        """
        psize = 2000  # Tamaño del pixel en m
        N = 5424  # Tamaño de imagen con psize=2000 m

        data = Dataset(self.file_path)  # Abro el archivo netcdf
        metadato = data.variables  # Extraigo todas las variables
        banda = metadato["band_id"][:].data[0]  # Extraigo el nro de banda
        altura = metadato['goes_imager_projection'].perspective_point_height
        semieje_may = metadato['goes_imager_projection'].semi_major_axis
        semieje_men = metadato['goes_imager_projection'].semi_minor_axis
        lat_cen = metadato['goes_imager_projection'].latitude_of_projection_origin
        lon_cen = metadato['goes_imager_projection'].longitude_of_projection_origin
        scale_factor = metadato['x'].scale_factor
        offset = np.array([metadato['x'].add_offset, metadato['y'].add_offset])

        pto_sup_izq = latlon2geos8.latlon2scan(lat_sup, lon_west, lon_cen, Re=semieje_may,
                                               Rp=semieje_men, h=altura)
        x0 = pto_sup_izq[1]*3600000.0
        y0 = pto_sup_izq[0]*3600000.0

        # Extraigo la imagen y la guardo en un array de np
        image = np.array(metadato["CMI"][:].data)

        if int(banda) == 3:
            esc = 0.5
            # escala es 1/2 porque tamaño de pixel de banda 3 = 1 km y tamaño pixel del resto = 2 km
            x = range(0, 10848)
            y = range(0, 10848)
            f = interpolate.interp2d(x, y, image, kind="cubic")
            xnew = np.arange(x[0], x[-1], (x[1] - x[0]) / esc)
            ynew = np.arange(y[0], y[-1], (y[1] - y[0]) / esc)
            image = f(xnew, ynew)

        # tamaño del recorte en proyeccion goes
        # img_extentr = [x0, x0+columnas*psize, y0 -filas*psize, y0]

        esc = int(N / image.shape[0])
        Nx = int(cols / esc)  # numero de puntos del recorte en x
        Ny = int(rows / esc)  # numero de puntos del recorte en x
        f0 = int(
            (-y0 / psize + N / 2 - 1.5) / esc
        )  # fila del angulo superior izquierdo
        # columna del angulo superior izquierdo
        c0 = int((x0 / psize + N / 2 + 0.5) / esc)
        f1 = int(f0 + Ny)  # fila del angulo inferior derecho
        c1 = int(c0 + Nx)  # columna del angulo inferior derecho
        # print('coordenadas filas, col: ', f0, c0, f1, c1)

        im_rec = image[f0:f1, c0:c1]
        return im_rec

    def reflactance(self, rec07, rec13, latlon_extent=[-80, -55, -34, 10]):
        # Correccion del zenith
        lat = np.linspace(latlon_extent[3], latlon_extent[1], rec07.shape[0])
        lon = np.linspace(latlon_extent[0], latlon_extent[2], rec07.shape[1])
        zenith = np.zeros((rec07.shape[0], rec07.shape[1]))

        # Calculate the solar zenith angle
        utc_time = datetime(
            2019, 1, 2, int(self.utc_hour[:2]), int(self.utc_hour[2:])
        )
        for x in range(len(lat)):
            for y in range(len(lon)):
                zenith[x, y] = astronomy.sun_zenith_angle(
                    utc_time, lon[y], lat[x]
                )
        refl39 = Calculator(
            platform_name="GOES-16", instrument="abi", band="ch7"
        )
        data07b = refl39.reflectance_from_tbs(zenith, rec07, rec13)

        return data07b

    def day_microphysicsRGB(
        self, rec03, rec07, rec13, latlon_extent=[-80, -55, -34, 10]
    ):
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
        R = ((R - Rmin) / (Rmax - Rmin)) ** (1 / gamma_R)
        G = ((G - Gmin) / (Gmax - Gmin)) ** (1 / gamma_G)
        B = ((B - Bmin) / (Bmax - Bmin)) ** (1 / gamma_B)

        # Normalizamos (matplotlib lo normaliza de todas formas)
        RR = np.copy(R)
        BB = np.copy(B)
        GG = np.copy(G)

        RR[RR < 0] = 0.0
        RR[RR > 1] = 1.0
        BB[BB < 0] = 0.0
        BB[BB > 1] = 1.0
        GG[GG < 0] = 0.0
        GG[GG > 1] = 1.0

        # Create the RGB
        RGB = np.stack([R, G, B], axis=2)
        # el axis está para que el shape sea fil col dim y no dim col fil
        RRGB = np.stack([RR, GG, BB], axis=2)
        print(RGB.shape)
        return RRGB


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
    lc_filter = lc_rfilter * lc_gfilter * lc_bfilter
    # Mask= magenta
    img_mask[lc_rfilter, 0] = 1.0
    img_mask[lc_rfilter, 1] = 0.0
    img_mask[lc_rfilter, 2] = 1.0

    # Stratus/Stratoculumus (small drops, low clouds) -> bright green/blue
    st_rfilter = (rgb[:, :, 0] > 0.3) * (rgb[:, :, 0] < 0.45)  # R
    st_gfilter = (rgb[:, :, 1] > 0.5) * (rgb[:, :, 1] < 0.8)  # G
    st_bfilter = rgb[:, :, 2] < 0.7
    st_filter = st_rfilter * st_gfilter * st_bfilter
    # Mask=Light blue
    img_mask[st_filter, 0] = 0.0
    img_mask[st_filter, 1] = 1.0
    img_mask[st_filter, 2] = 1.0

    # CumuloNimbis (high clouds) -> red, dark orange
    cb_rfilter = rgb[:, :, 0] > 0.7  # R
    cb_gfilter = rgb[:, :, 1] < 0.3  # G
    cb_bfilter = rgb[:, :, 2] < 0.3  # B
    cb_filter = cb_rfilter * cb_gfilter * cb_bfilter
    # Mask=Red
    img_mask[cb_filter, 0] = 1.0
    img_mask[cb_filter, 1] = 0.0
    img_mask[cb_filter, 2] = 0.0

    # Cirrus (high clouds)-> green, dark green
    cr_rfilter = rgb[:, :, 0] < 0.3  # R
    cr_gfilter = rgb[:, :, 1] > 0.7  # G
    cr_bfilter = rgb[:, :, 2] < 0.3  # B
    cr_filter = cr_rfilter * cr_gfilter * cr_bfilter
    # Mask= Green
    img_mask[cr_filter, 0] = 0.0
    img_mask[cr_filter, 1] = 1.0
    img_mask[cr_filter, 2] = 0.0

    # supercooled clouds Thick, small drops, medium clouds-> yellow
    super_rfilter = rgb[:, :, 0] > 0.8
    super_gfilter = rgb[:, :, 1] > 0.8
    super_bfilter = rgb[:, :, 2] < 0.2  # amarillo
    super_filter = super_rfilter * super_gfilter * super_bfilter
    # Mask=Yellow
    img_mask[super_filter, 0] = 1.0
    img_mask[super_filter, 1] = 1.0
    img_mask[super_filter, 2] = 0.0

    return img_mask[:, :, [0, 1, 2]]
# %%


# Main
main_path = "C:\\Users\\Paula\\Downloads\\GOES_L2\\2019-004\\"

goes3 = GoesClass(main_path +
                  "OR_ABI-L2-CMIPF-M3C03_G16_s20190041600363_e20190041611130_c20190041611201.nc"
                  )

goes3

banda3rec = goes3.recorte()

goes7 = GoesClass(main_path +
                  "OR_ABI-L2-CMIPF-M3C07_G16_s20190041600363_e20190041611142_c20190041611202.nc"
                  )


banda7rec = goes7.recorte()

goes13 = GoesClass(main_path +
                   "OR_ABI-L2-CMIPF-M3C13_G16_s20190041600363_e20190041611142_c20190041611224.nc"
                   )


banda13rec = goes13.recorte()

banda7cor = goes7.reflactance(banda7rec, banda13rec)  # banda 7 corregida

rgb = goes3.day_microphysicsRGB(banda3rec, banda7cor, banda13rec)

plt.figure(1)
plt.imshow(rgb)


mascara = mask(rgb)

crs = ccrs.Geostationary(central_longitude=-75.0, satellite_height=35786023.0)
img_extentr = [
    -5017302.234120312,  # left
    9456653.737754004,  # right
    -22360446.075856283,  # bottom
    10288448.967580065,  # top
]

fig, axis = plt.subplots()
axis = plt.axes(projection=crs)
axis.gridlines()
axis.coastlines(resolution="10m", color="white")
axis.add_feature(cf.BORDERS)
axis.add_feature(cf.LAKES)
axis.add_feature(cf.RIVERS)
plt.imshow(mascara, origin="upper", extent=img_extentr)

plt.figure(2)
plt.imshow(mascara)
plt.show()
