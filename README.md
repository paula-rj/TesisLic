# Procesamiento-de-datos-Cloudsat-y-GOES16
Este repositorio incluye codigo en Python para abrir y graficar archivos Cldclass de Cloudsat, archivos de GOES16 y cómo graficar archivos correspondientes de ambos satelites.
Se utilizan los productos de nivel 2 CLDCLASS de CloudSat (http://www.cloudsat.cira.colostate.edu/data-products/level-2b/2b-cldclass?term=88) y los productos de nivel 2 CMIPF (Cloud and moisture: https://www.goes-r.gov/products/baseline-cloud-moisture-imagery.html) de GOES16 para generar un producto RGB de nivel 3 denominado microfísica de día. 

# Motivación
Las nubes juegan un rol fundamental en los procesos meteorológicos de nuestro planeta, por lo tanto entenderlos y caracterizarlos es de vital importancia para poder predecir fenómenos meteorológicos. 
Los nubes se clasifican, basicamente, según su desarrollo vertical y la fase en la que se enuentran las gotas que las componen (líquida o sólida), entre otras variables.
Cuando la radiación electromagnética interactúa con las gotas, producen resultados distintos dependientes de la longitud de onda, en el marco de la teoría de la transferencia radiativa para la atmósfera (Rees 2012). Observando las nubes en distintas bandas espectrales, se pueden idear algoritmos para diferenciar tipos de nubes a partir de imágenes espectrales (Kidder and Vonder Haar 1995, Marshak and Davis 2005).
CloudSat es un satélite artificial lanzado por la NASA en 2006, cuyo instrumento principal es un radar que mide las microondas reflejadas por las nubes en función de su distancia. Su órbita es heliosíncrona y actualmente sólo puede tomar perfiles durante el día debido a una falla en su batería.
En cambio, GOES16 O GOES-R es un satélite geocentrico que toma imágenes en 16 bandas de todo el continente americano. Combinando algunas de estas bandas se pueden generar, por ejemplo, imágenes RGB que permitan distinguir tipos de nubes según los colores. 
Combinando productos de ambos satélites, se puede obtener una clasificación de las nubes que se observan desde el espacio exterior cercano al planeta. En particular, nos vamos a enfocar en la región de Sudamérica.


# Páginas útiles

## Relacionadas con Cloudsat:
* Para *descargar archivos de Cloudsat*, crear un usuario y:
  * Opción 1: Descargarlos desde la página oficial: http://www.cloudsat.cira.colostate.edu/order-data
  * Opcion 2: Util en Windows y si no funciona la página, utilizar un programa estilo FileZilla o similar para descargar ftp
* Notacion de los nombres de los archivos de Cloudsat: http://www.cloudsat.cira.colostate.edu/data-products
* Código para graficar tipos de nubes en cada capa segun la altura vs la latitud: https://moonbooks.org/Codes/Plot-cldclass-lidar-granule-vertical-profile-using-python-3/
* Código para  graficar tipo de nube por altura vs latitud de archivos cloudclass: https://docserver.gesdisc.eosdis.nasa.gov/public/project/MEaSUREs/Fetzer/README.AIRS_CloudSat.pdf (es un pdf, e código está al final)

## Relacionadas con GOES
* Para *descargar archivos de GOES16*: si son pocas y de alguna fecha específica recomiendo esta página, hay que llenar el formulario y descargar la im{agen de la banda y hora que se desee: https://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/goes16_download.cgi
* Páginas oficiales de GOES16 con todos sus productos: https://www.goes-r.gov/products/overview.html , https://www.star.nesdis.noaa.gov/goesr/product_cp_cloud.php
* Guía rápida del producto *Daytime Microphysics*: https://weather.msfc.nasa.gov/sport/training/quickGuides/rgb/QuickGuide_DtMicroRGB_NASA_SPoRT.pdf . Nota: este producto no viene "listo para descargar", hay que generarlo como dice en la guía. 
* Códigos para generar productos RGB mediante la libreria GDAL: https://geonetcast.wordpress.com/2019/07/03/python-script-examples-to-generate-goes-16-rgbs/ (al final está el RGB de microfísica de día)

