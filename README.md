# Procesamiento-de-datos-Cloudsat-y-GOES16
Este repositorio incluye codigo en Python para abrir y graficar archivos Cldclass de Cloudsat, archivos de GOES16 y cómo graficar archivos correspondientes de ambos satelites.
Se utilizan los productos de nivel 2 CLDCLASS de CloudSat (http://www.cloudsat.cira.colostate.edu/data-products/level-2b/2b-cldclass?term=88) y los productos de nivel 2 CMIPF (Cloud and moisture: https://www.goes-r.gov/products/baseline-cloud-moisture-imagery.html) de GOES16 para generar un producto RGB de nivel 3 denominado microfísica de día. 

# Páginas útiles

## Relacionadas con Cloudsat:
* Para *descargar archivos de Cloudsat*, crear un usuario y:
* * Opción 1: Descargarlos desde la página oficial: http://www.cloudsat.cira.colostate.edu/order-data
* * Opcion 2: Util en Windows y si no funciona la página, utilizar un programa estilo FileZilla o similar para descargar ftp
* Notacion de los nombres de los archivos de Cloudsat: http://www.cloudsat.cira.colostate.edu/data-products
* Código para graficar tipos de nubes en cada capa segun la altura vs la latitud: https://moonbooks.org/Codes/Plot-cldclass-lidar-granule-vertical-profile-using-python-3/
* Código para  graficar tipo de nube por altura vs latitud de archivos cloudclass: https://docserver.gesdisc.eosdis.nasa.gov/public/project/MEaSUREs/Fetzer/README.AIRS_CloudSat.pdf (es un pdf, e código está al final)

## Relacionadas con GOES
* Para *descargar archivos de GOES16*: si son pocas y de alguna fecha específica recomiendo esta página, hay que llenar el formulario y descargar la im{agen de la banda y hora que se desee: https://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/goes16_download.cgi
* Páginas oficiales de GOES16 con todos sus productos: https://www.goes-r.gov/products/overview.html , https://www.star.nesdis.noaa.gov/goesr/product_cp_cloud.php
* Guía rápida del producto *Daytime Microphysics*: https://weather.msfc.nasa.gov/sport/training/quickGuides/rgb/QuickGuide_DtMicroRGB_NASA_SPoRT.pdf . Nota: este producto no viene "listo para descargar", hay que generarlo como dice en la guía. 
* Códigos para generar productos RGB mediante la libreria GDAL: https://geonetcast.wordpress.com/2019/07/03/python-script-examples-to-generate-goes-16-rgbs/ (al final está el RGB de microfísica de día)

