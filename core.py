import goes
import cldsat
import latlon2geos8


def superpix(cs_fila, cs_col, npix, lat_lon_extent):
    """
    Define un superpixel, es decir una sección de la imagen
    con cierto tamaño, pequeño, de forma tal que se aprecien
    las coincidencias de colores de GOES con las clasificaciones
    de Cloudsat para esos pixeles.

    Parameters
    ----------
    cs_fila: int
    Fila en la que se encuentra el pixel de Cloudsat
    cs_col: int
    Columna en la que se encuentra el pixel de Cloudsat
    npix= int
    cantidad de pixeles que tendra
    lat_lon_extent

    Returns:
    -------
    superpix: np array
    El punto de cloudsat + superpixel, pero todo de GOES ya en RRGB
    """
    banda3 = Dataset(GOES_path_list[0])
    banda7 = Dataset(GOES_path_list[1])
    banda13 = Dataset(GOES_path_list[2])
    metadato1 = banda3.variables
    metadato2 = banda7.variables
    metadato3 = banda13.variables
    superpix1 = metadato1["CMI"][:].data[
        cs_fila - npix : cs_fila + npix + 1, cs_col - npix : cs_col + npix + 1
    ]
    superpix2 = metadato2["CMI"][:].data[
        cs_fila - npix : cs_fila + npix + 1, cs_col - npix : cs_col + npix + 1
    ]
    superpix3 = metadato3["CMI"][:].data[
        cs_fila - npix : cs_fila + npix + 1, cs_col - npix : cs_col + npix + 1
    ]

    superpix2b = solar_7(superpix2, superpix3, lat_lon_extent)

    RGBsuperpix = dayRGB(superpix1, superpix2b, superpix3)
    return RGBsuperpix

    # main


lista_paths = [
    "/content/gdrive/MyDrive/Kaggle/OR_ABI-L2-CMIPF-M3C03_G16_s20190021800363_e20190021811129_c20190021811205.nc",
    "/content/gdrive/MyDrive/Kaggle/OR_ABI-L2-CMIPF-M3C07_G16_s20190021800363_e20190021811141_c20190021811202.nc",
    "/content/gdrive/MyDrive/Kaggle/OR_ABI-L2-CMIPF-M3C13_G16_s20190021800363_e20190021811141_c20190021811221.nc",
]
lat_lon_extent = [-69.890762, -24.111683, -69.879143, -23.976511]
a = superpix(3979, 2965, lat_lon_extent, lista_paths)

fig_dims = (6, 6)
cmap = mpl.colors.ListedColormap(["navy", "crimson", "limegreen", "gold"])
fig, axis = plt.subplots(figsize=fig_dims)
# sns.scatterplot(x='col', y='fil', data=day2df.iloc[5903:6119,:], hue='capa0', palette="bright",s=2)
plt.imshow(
    a,
    vmin=0.0,
    vmax=1.0,
)
scatter = axis.scatter(
    x=day2df.iloc[5903:6119, :]["col"].to_numpy() - 2965 + 50,
    y=day2df.iloc[5903:6119, :]["fil"].to_numpy() - 3979 + 50,
    c=day2df.iloc[5903:6119, :]["capa0"].to_numpy(),
    cmap=cmap,
)

fig.colorbar(scatter, ticks=np.linspace(0, 3, 4))
plt.show()
