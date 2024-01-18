from .package_installer import install_packages
try:
    from qgis.core import QgsProject,QgsVectorLayer,QgsFeature,QgsGeometry,QgsPoint
    from qgis.PyQt.QtCore import QDate
    import json
    import os
    import io
    import ee
    import glob
    import shutil
    import requests
    import geemap
    import numpy as np
    import folium
    from .codes import download
    # from .parameters import aoi, date
    import rasterio
    from rasterio.merge import merge
    from rasterio.mask import mask
    from natsort import natsorted
    from pyproj import CRS
except ImportError:
    install_packages()
    
import warnings
warnings.filterwarnings("ignore")
def initialize():
    json_file="ee-brianchelloti-f3bf0b9227be.json"
    script_directory = os.path.dirname(os.path.abspath(__file__))
    json_file_path = os.path.join(script_directory, json_file)
    service_account = 'brian-python@ee-briansimiyu.iam.gserviceaccount.com'
    credentials = ee.ServiceAccountCredentials(service_account, json_file_path)
    print(credentials)
    ee.Initialize(credentials)

def download_image(dlg):
    # Create temporary folder
    # if not os.path.exists('./output/satellite-image/pre-processed-images'):
    #     os.makedirs('./output/satellite-image/pre-processed-images')
    initialize()
    # ------ Main Code Execution ------ #
    start_date = dlg.DownloadStartDateEdit.date().toString("yyyy-MM-dd")
    end_date = dlg.DownloadEndDateEdit.date().toString("yyyy-MM-dd")
    output_path=dlg.DownloadOutputLineEdit.text()

    # Create filename
    start_yyyy = start_date[:4]
    start_mm = start_date[5:7]
    end_yyyy = end_date[:4]
    end_mm = end_date[5:7]
    file_name = start_yyyy+start_mm+'_'+end_yyyy+end_mm
    print(file_name)

    # Get the selected shapefile
    shp_name=dlg.BoundComboBox.currentText()
    bound=QgsProject.instance().mapLayersByName(shp_name)[0]

    if not bound.isValid():
        print("layer is invalid")
    else:

        polygon=next(bound.getFeatures())
        polygon_geometry=polygon.geometry()
        if not polygon_geometry.isEmpty():
            polygon_coordinates = []
            # print(polygon_geometry)
            if polygon_geometry.isMultipart():
                for part in polygon_geometry.asMultiPolygon():
                    for point in part[0]:
                        polygon_coordinates.append([point.x(),point.y()])
            else:
                for point in polygon_geometry.asPolygon()[0]:
                    polygon_coordinates.append([point.x(),point.y()])
            polygon_json={
                "type":"Polygon",
                "coordinates":[polygon_coordinates]
            }
            print(polygon_json)
    polygon_json=ee.Geometry.Polygon(polygon_json['coordinates'])
    # ------ Calculate EPSG ------ #
    # epsg_code = '4326' # WGS 84
    lon = polygon_json.geometries().getInfo()[0]['coordinates'][0][0][0]
    lat = polygon_json.geometries().getInfo()[0]['coordinates'][0][0][1]
    epsg_code = "3857"
    dlg.progressBar.setValue(10)
    # Split the geometry to small grid 
    grid = geemap.fishnet(polygon_json, h_interval=0.1, v_interval=0.1, delta=1)
    gridList = grid.toList(grid.size())
    grid_num = grid.toList(grid.size()).length()

    # Create list of each grid feature
    ls_feature = []
    for i in range(grid_num.getInfo()):
        feature = ee.Feature(gridList.get(i)).geometry().bounds()
        ls_feature.append(feature)

    # Create temporary folder
    if not os.path.exists('./content/temp/grid'):
        os.makedirs('./content/temp/grid')

    # Download image by grid
    for i in range(grid_num.getInfo()):
        # Extract image from GEE and apply pre-process function
        image = download.Sentinel_no_clouds(ls_feature[i], start_date, end_date)

        # Set band IDs based on landsat collection
        BandIDs = ['B11', 'B8', 'B4', 'B3', 'B2']

        # Download pre-processed image
        download_id = ee.data.getDownloadId({
            'image': image,
            'bands': BandIDs,
            'region': ls_feature[i],
            'scale': 10,
            'format': 'GEO_TIFF',
            'crs' : 'EPSG:'+str(epsg_code),
        })

        response = requests.get(ee.data.makeDownloadUrl(download_id))
        with open('./content/temp/grid/image_grid_'+str(i)+'.tif', 'wb') as fd:
            fd.write(response.content)
    dlg.progressBar.setValue(30)
    # ------ Merge all grid image to one image ------ #
    # Make a search criteria to select the ndvi files
    q = os.path.join("./content/temp/grid/image*.tif") 

    # sorted files by name
    fp = natsorted(glob.glob(q)) 

    # List for storing the raster image
    src_files = []

    # Open each raster files by iterating and then append to our list
    for raster in fp:
        # open raster file
        files = rasterio.open(raster)

        # add each file to our list
        src_files.append(files)

    # Merge function returns a single mosaic array and the transformation info
    mosaic, out_trans = merge(src_files)

    # Set metadata
    out_meta = src_files[0].meta.copy()
    out_meta.update({"driver": "GTiff",
                    "dtype": "float32",
                    "nodata": None,
                    "height": mosaic.shape[1],
                    "width": mosaic.shape[2],
                    "transform": out_trans,
                    "count": 5,
                    "crs": CRS.from_epsg(int(epsg_code))
                    }
                    )
                    
    # Write the mosaic raster
    # Create temporary folder
    if not os.path.exists('./content/temp/main'):
        os.makedirs('./content/temp/main')
    output = os.path.join('./content/temp/main/image_snrgb.tif')
    with rasterio.open(output, "w", **out_meta) as dest:
        dest.write(mosaic.astype(np.float32))
    dlg.progressBar.setValue(50)
    # Clip image to aoi geometry
    img_grid = rasterio.open('./content/temp/main/image_snrgb.tif')
    aoi_epsg = polygon_json.transform(ee.Projection('EPSG:'+str(epsg_code)), 1)

    # Mask raster
    clip, clip_transform = mask(img_grid, aoi_epsg.geometries().getInfo(), crop=True)
    dlg.progressBar.setValue(60)
    # Set metadata
    out_meta = img_grid.meta.copy()
    out_meta.update({"driver": "GTiff",
                    "dtype": "float32",
                    "nodata": 0 and None,
                    "height": clip.shape[1],
                    "width": clip.shape[2],
                    "transform": clip_transform,
                    "count": 5,
                    "crs": img_grid.crs
                    }
                    )
    dlg.progressBar.setValue(80)
    # Write the clip raster
    output = os.path.join(output_path+'/'+file_name+'.tif')
    with rasterio.open(output, "w",tiled=True, compress='lzw', **out_meta) as dest:
        dest.write(clip.astype(np.float32))

    # Shut down temporary directory
    # shutil.rmtree('./content/temp')
    
    #---------------------------------------------------------------------
    if img_grid.read().any() == 0:
        print(' Error: Images are not available for this area within the given date.')
    else:
        print(' Download completed!')
    dlg.progressBar.setValue(100)
