# CON CLOUDSAT
import os

import numpy as np

import matplotlib.pyplot as plt

import geopandas as gpd

import pandas as pd

import seaborn as sns

# PyHDF
import pyhdf
from pyhdf.HDF import *
from pyhdf.VS import *
from pyhdf.SD import SD, SDC


class CldClass:
    def __init__(self, path_total):
        self.path_total = path_total
        self.file_name = os.path.split(self.path_total)[-1]
        date = self.file_name.split("_")[0]
        self.year = date[:4]
        self.julian_day = date[4:7]
        self.hour_utc = date[7:9]
        self.light = ""
        if int(self.hour_utc) > 10:
            self.light = "day"
        else:
            self.light = "night"

    def open_cldclass(self, sur=True):
        """
          Función que toma un path donde debe estar el archivo CLDCLASS de Cloudsat y devuelve un dataframe de Pandas
        con los datos del archivo cldclass guardados en columnas del dataframe: latitud, longitud, capa0, capa1,..., capa9;
        con cada tipo de nube segun cada capa.
        Parameters:
        -----------
        sur : Bool. Si se quiere abrir solo la parte de sudamerica. Default=True.

        Returns:
        -----------
        layers_df: Pandas Dataframe lat, lon, tipo de nube en cada una de las capas (son 10)
        """
        # Read v data
        hdf_file = HDF(self.path_total, HC.READ)
        vs = hdf_file.vstart()
        vdata = (
            vs.vdatainfo()
        )  # es una lista de tuplas. acá estan lat y long y cloud layers

        vd_lat = vs.attach("Latitude", write=0)
        lat = vd_lat[:]
        vd_lat.detach

        vd_lon = vs.attach("Longitude", write=0)
        lon = vd_lon[:]
        vd_lon.detach
        vs.end()
        # hdf_file.close()

        latitud = np.array(lat).flatten()
        longitud = np.array(lon).flatten()

        # Read sd data
        file = SD(self.path_total, SDC.READ)
        start_point = 0
        end_point = 36951
        if sur == True:
            if self.light == "night":
                end_point = 6000
            else:
                end_point = 20000
                if self.hour_utc == (15):
                    start_point = 6000
                else:  # 16,17,18 utc
                    start_point = 10000
            latitud = latitud[start_point:end_point]
            longitud = longitud[start_point:end_point]
            cld_layertype = file.select("CloudLayerType")[
                start_point:end_point
            ]
        else:  # grafica toda la orbita
            cld_layertype = file.select("CloudLayerType")[:]

        layers_df = pd.DataFrame(
            {
                "Longitude": longitud,
                "Latitude": latitud,
                "capa0": cld_layertype[:, 0],
                "capa1": cld_layertype[:, 1],
                "capa2": cld_layertype[:, 2],
                "capa3": cld_layertype[:, 3],
                "capa4": cld_layertype[:, 4],
                "capa5": cld_layertype[:, 5],
                "capa6": cld_layertype[:, 6],
                "capa7": cld_layertype[:, 7],
                "capa8": cld_layertype[:, 8],
                "capa9": cld_layertype[:, 9],
            }
        )

        return layers_df

    def plot_cldclass_geos(self, layers_data, capa_n):
        """
          Función que grafica los tipos de nubes en la capa asignada en la pasada del Cloudsat.
        Parameters:
        -----------
        layers_data : Pandas Dataframe.
          Dataframe de Pandas que contiene Latitud,Longitud,y la capa que se quiere graficar.
        capa_n : int
          Numero de capaque se quiere graficar. Van de 0 a 9.

        Returns:
        -----------
        Plot. Tipo de nube en la capa 9 con órbita proyectada en proyeccion geoestacionaria de GOES16,
        con mapa de costas de fondo.
        """
        capa_str = "capa" + str(capa_n)

        # Generamos geodataframe a partir del pd dataframe de entrada
        geo_df = gpd.GeoDataFrame(
            layers_data.loc[:, ["Longitude", "Latitude", capa_str]],
            geometry=gpd.points_from_xy(
                layers_data.Longitude, layers_data.Latitude
            ),
        )
        geo_df.crs = {
            "init": "epsg:4326"
        }  # EPSG 4326 corresponds to coordinates in latitude and longitude
        # Reprojecting into GOES16 geostationary projection
        geodf_GOESproj = geo_df.to_crs("+proj=geos +h=35786023.0 +lon_0=-75.0")
        crs = ccrs.Geostationary(
            central_longitude=-75.0, satellite_height=35786023.0
        )  # proyeccion geoestacionaria para Goes16
        fig_dims = (10, 10)
        fig, axis = plt.subplots(figsize=fig_dims)
        axis = plt.axes(projection=crs)
        axis.gridlines
        axis.coastlines(resolution="10m", color="blue")
        sns.scatterplot(
            x="Longitude",
            y="Latitude",
            data=geodf_GOESproj,
            hue=capa_str,
            palette="bright",
            s=2,
            transform=ccrs.PlateCarree(),
        )
        axis.set_title(
            """year {}; day {}; hour {}; {}""".format(
                self.year, self.julian_day, self.hour_utc, self.light
            )
        )
        plt.show()

    def plotlatlon_cld(self, capa_n, layers_data_df):
        """
        Grafica la capa capa_n con latitud en el eje y y longitud en el eje x, sin ninguna proyeccion
         Parameters
         ----------
         capa_n : int
             El numero de capa que quiero dibujar. Entre 0 y 9.
         layers_data_df : pandas DataFrame
             DataFrame de Pandas o Geopandas que incluya latitud, longitud y la capa que vamos a dibujar.
         Returns
         ----------
         Plot."""

        capa_str = "capa" + str(capa_n)
        fig_dims = (6, 6)
        fig, ax = plt.subplots(figsize=fig_dims)
        sns.scatterplot(
            x="Longitude",
            y="Latitude",
            data=layers_data_df,
            hue="capa2",
            palette="bright",
            marker="o",
            s=1,
        )
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        plt.show()
