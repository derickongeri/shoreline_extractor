from .package_installer import install_packages
try:
    import os
    import ee
    import glob
    import requests
    import geemap
    import shutil
    import numpy as np
    from .Plugin_download import waterMask_download, Sentinel_no_clouds
    import rasterio
    from rasterio.merge import merge
    from rasterio.features import shapes
    from shapely.ops import unary_union
    from .codes import shoreline
    # from rasterio.mask import mask
    from natsort import natsorted
    from pyproj import CRS
    import geopandas as gpd
    # import shapely
    from shapely.geometry import Polygon, MultiPolygon
    from shapely.geometry.polygon import LinearRing
    from qgis.core import QgsProject,QgsVectorLayer,QgsRasterLayer,QgsFeature,QgsGeometry,QgsPoint
except ImportError:
    install_packages()
import warnings
warnings.filterwarnings("ignore")

#Create a function to plot vector maps on the QGIS project
def addMapLayer(path):
    project = QgsProject.instance()
    filename=path.split('/')[-1][0:-5]
    layer=QgsVectorLayer(path,filename,'ogr')
    if layer.isValid():
        project.instance().addMapLayer(layer)
    else:
        print("Failed to load layer")

#create a function to plot the raster maps on the QGIS project
def addRasterLayer(path):
    project=QgsProject.instance()
    filename=path.split('/')[-1][0:-4]
    layer=QgsRasterLayer(path,filename,'gdal')
    if layer.isValid():
        project.instance().addMapLayer(layer)
    else:
        print("Failed to load layer")

def initialize():
    json_file="ee-brianchelloti-f3bf0b9227be.json"
    script_directory = os.path.dirname(os.path.abspath(__file__))
    json_file_path = os.path.join(script_directory, json_file)
    service_account = 'brian-python@ee-briansimiyu.iam.gserviceaccount.com'
    credentials = ee.ServiceAccountCredentials(service_account, json_file_path)
    print(credentials)
    ee.Initialize(credentials)

def extract_SAR_Shoreline(dlg):
    initialize()
    #Get input field values
    start_date = dlg.DownloadStartDateEditSAR.date().toString("yyyy-MM-dd")
    end_date = dlg.DownloadEndDateEditSAR.date().toString("yyyy-MM-dd")
    output_folder=dlg.DownloadOutputLineEditSAR.text()
    
    #Get bounds
    shp_name=dlg.BoundComboBoxSAR.currentText()
    bound=QgsProject.instance().mapLayersByName(shp_name)[0]
    # Define the working CRS
    epsg_code = "3857"

    #Requires PyCRS
    dlg.progressBar.setValue(10)
    # Load the roi shapefile
    roi = geemap.shp_to_ee(bound.source())
    print(bound.source())
    # Create temporary folder
    if not os.path.exists(output_folder+'/output/WaterMask/pre-processed-images'):
        os.makedirs(output_folder+'/output/WaterMask/pre-processed-images')

    # ------ Main Code Execution ------ #
    # Split the geometry to small grid 
    grid = geemap.fishnet(roi, h_interval=0.05, v_interval=0.05, delta=1)
    gridList = grid.toList(grid.size())
    grid_num = grid.toList(grid.size()).length()

    # Create list of each grid feature
    ls_feature = []
    for i in range(grid_num.getInfo()):
        feature = ee.Feature(gridList.get(i)).geometry()
        ls_feature.append(feature)

    dlg.progressBar.setValue(20)
    # Create temporary folder
    if not os.path.exists(output_folder+'/output/Tiles/watermask'): 
        os.makedirs(output_folder+'/output/Tiles/watermask')

    if not os.path.exists(output_folder+'/output/Tiles/rgb/'): 
        os.makedirs(output_folder+'/output/Tiles/rgb/')
    dlg.progressBar.setValue(25)
    # 2. Define start date and end date
    date = [start_date, end_date]

    # ------ Main Code Execution ------ #
    start_date = date[0]
    end_date = date[1]

    # Create filename
    start_yyyy = start_date[:4]
    start_mm = start_date[5:7]
    end_yyyy = end_date[:4]
    end_mm = end_date[5:7]
    file_name = start_yyyy+start_mm+'_'+end_yyyy+end_mm

    # Download SAR data
    #Download image by grid
    #Download image by grid
    BandIDs_watermask = ['waterMask']
    BandIDs_rgb = ['VV', 'VH', 'VV-VH']
    for i in range(grid_num.getInfo()):
  # Extract image from GEE and apply pre-process function
        image = waterMask_download(ls_feature[i], start_date, end_date)

        # Download pre-processed image
        download_watermask = ee.data.getDownloadId({
            'image': image,
            'bands': BandIDs_watermask,
            'region': ls_feature[i],
            'scale': 10,
            'format': 'GEO_TIFF',
            'crs' : 'EPSG:'+str(epsg_code),
        })
        
        response_watermask = requests.get(ee.data.makeDownloadUrl(download_watermask))
        with open(output_folder+'/output/Tiles/watermask/image_grid_'+str(i)+'.tif', 'wb') as fd:
            fd.write(response_watermask.content)
        
        # Download pre-processed image
        download_rgb = ee.data.getDownloadId({
            'image': image,
            'bands': BandIDs_rgb,
            'region': ls_feature[i],
            'scale': 10,
            'format': 'GEO_TIFF',
            'crs' : 'EPSG:'+str(epsg_code),
        })

        response_rgb = requests.get(ee.data.makeDownloadUrl(download_rgb))
        with open(output_folder+'/output/Tiles/rgb/image_grid_'+str(i)+'.tif', 'wb') as fd:
            fd.write(response_rgb.content)
    dlg.progressBar.setValue(30)   
            # ------ Merge all grid image to one image ------ #
    # Make a search criteria to select the image files
    watermask_infolder = output_folder+'/output/Tiles/watermask/image_grid_*.tif'
    rgb_infolder = output_folder+'/output/Tiles/rgb/image_grid_*.tif'

    watermask_outfolder = os.path.join(output_folder+'/output/WaterMask/pre-processed-images'+'/WaterMask_'+file_name+'.tif')
    rgb_outfolder = os.path.join(output_folder+'/s1rgb_'+file_name+'.tif')

    def MergeGrids(in_folder,out_folder, bands):
        
        q = in_folder
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
                        "nodata": 0,
                        "height": mosaic.shape[1],
                        "width": mosaic.shape[2],
                        "transform": out_trans,
                        "count": bands,
                        "crs": CRS.from_epsg(int(epsg_code))})
                        
        # Write the mosaic raster
        # output = os.path.join('./output/'+folder+'/WaterMask/pre-processed-images'+'/WaterMask_'+file_name+'.tif')
        output = out_folder
        with rasterio.open(output, "w", tiled=True, compress='lzw', **out_meta) as dest:
            dest.write(mosaic.astype(np.float32))
        
    MergeGrids(watermask_infolder, watermask_outfolder, 1)
    dlg.progressBar.setValue(40)   
    MergeGrids(rgb_infolder, rgb_outfolder, 3)
    addRasterLayer(output_folder+'/s1rgb_'+file_name+'.tif')
    dlg.progressBar.setValue(50)   

        # Create folder
    if not os.path.exists(output_folder+'/output/WaterMask/post-processed-images'):
        os.makedirs(output_folder+'/output/WaterMask/post-processed-images')
    
    # ------ Shoreline extraction ------ #
    # Read input image data
    Image = rasterio.open(output_folder+'/output/WaterMask/pre-processed-images'+'/WaterMask_'+file_name+'.tif')

    print(Image.crs)
    if Image.read().any() == 0:
        print('Warning: The image is empty, so shoreline cannot be extracted.')

    else:
        print('Extracting shoreline')
    
    rescale_image, transform = shoreline.resampling(image=Image, scale_factor=5)

    # Georeference
    horizontal_step = 0     # (+ Positive) Move to right // (- Negative) Move to left
    vertical_step = 0       # (+ Positive) Move to top // (- Negative) Move to bottom

    ncol = -horizontal_step
    nrow = vertical_step
    buffer_rate = 0

        # Slice column number
    if ncol < 0:
        # Negative
        rescale_image = np.delete(rescale_image, np.s_[rescale_image.shape[2]-abs(ncol):rescale_image.shape[2]], axis=2)
        new_col = np.empty((rescale_image.shape[0],rescale_image.shape[1], abs(ncol)))
        new_col.fill(np.nan)
        rescale_image = np.array([np.c_[new_col[i], rescale_image[i]] for i in range(rescale_image.shape[0])])
    else:
        # Positive
        rescale_image = np.delete(rescale_image, slice(0,ncol), axis=2)
        new_col = np.empty((rescale_image.shape[0],rescale_image.shape[1], ncol))
        new_col.fill(np.nan)
        rescale_image = np.c_[rescale_image, new_col]

    # Slice row number
    if nrow < 0:
        # Negative
        rescale_image = np.delete(rescale_image, np.s_[rescale_image.shape[1]-abs(nrow):rescale_image.shape[1]], axis=1)
        new_row = np.empty((rescale_image.shape[0], abs(nrow), rescale_image.shape[2]))
        new_row.fill(np.nan)
        rescale_image = np.array([np.r_[new_row[i], rescale_image[i]] for i in range(rescale_image.shape[0])])
    else:
        # Positive
        rescale_image = np.delete(rescale_image, slice(0,abs(nrow)), axis=1)
        new_row = np.empty((rescale_image.shape[0], abs(nrow), rescale_image.shape[2]))
        new_row.fill(np.nan)
        rescale_image = np.array([np.r_[rescale_image[i], new_row[i]] for i in range(rescale_image.shape[0])])


    crs ='EPSG:3857'
    # Set metadata
    out_meta = Image.meta.copy()
    out_meta.update({"driver": "GTiff",
                    "dtype": "float32",
                    "nodata": 0,
                    "height": rescale_image.shape[1],
                    "width": rescale_image.shape[2],
                    "transform": transform,
                    "count": 1,
                    "crs": crs#Image.crs
                    }
                    )
    dlg.progressBar.setValue(60)   
    # Write the clip raster
    output = os.path.join(output_folder+'/output/WaterMask/post-processed-images/WaterMask_'+file_name+'.tif')
    with rasterio.open(output, "w",tiled=True, compress='lzw', **out_meta) as dest:
        dest.write(rescale_image.astype(np.float32))
        
    rescale_image = rasterio.open(output_folder+'/output/WaterMask/post-processed-images/WaterMask_'+file_name+'.tif')

    mask = rescale_image.read(1)

        # ------ Save result as GeoJSON file ------ #
    # Export result
    polygons = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) 
            in enumerate(
                shapes(mask.astype('uint8'), mask=None, transform=transform)))
    geometry = list(polygons)

    # Create new geodataframe
    geom_dataframe  = gpd.GeoDataFrame.from_features(geometry)
    dlg.progressBar.setValue(70)   
    # Set projection of dataframe
    geom = geom_dataframe.set_crs(Image.crs)

    # Remove no-data geometries
    geom = geom[geom.raster_val != 0]

    # Extract  boundaries
    # list_poly = []
    # for p in geom['geometry']:
    #     for interior in p.interiors:
    #         list_poly.append(Polygon(interior))

    list_poly = []
    for p in geom['geometry']:
        list_poly.append(Polygon(p.exterior))
    dlg.progressBar.setValue(80)   
    smooth_poly = []
    for i in range(len(list_poly)):
        poly_line = list_poly[i]
        outbuffer = poly_line.buffer(10, resolution=5, cap_style=1, join_style=1, mitre_limit=2, single_sided=True)
        inbuffer = outbuffer.buffer(-10.5, resolution=5, cap_style=1, join_style=1, mitre_limit=2, single_sided=True)
        simplified = inbuffer.simplify(1, preserve_topology=False) # False: Use Douglas-Peucker algorithm
        inbuffer2 = simplified.buffer(buffer_rate, resolution=5, cap_style=1, join_style=1, mitre_limit=2, single_sided=True)
        
        if type(inbuffer2) == MultiPolygon:
            ring = [LinearRing(inbuffer2.geoms[k].exterior) for k in range(len(inbuffer2.geoms))]
            ring = unary_union(ring)
            smooth_poly.append(ring)
        else:
            ring = LinearRing(inbuffer2.exterior)
            smooth_poly.append(ring)
    dlg.progressBar.setValue(90)   
    # Create new geodataframe for exterior boundaries
    geo_shoreline = gpd.GeoDataFrame({'geometry':smooth_poly}, crs=Image.crs)
    geo_shoreline = geo_shoreline.dropna().reset_index(drop=True)
    geo_shoreline['id'] = geo_shoreline.index

    # Save to geojson file
    outfp = output_folder+'/shoreline_'+file_name+'.json'
    geo_shoreline.to_file(outfp, driver='GeoJSON')
    addMapLayer(outfp)
    dlg.progressBar.setValue(100)   

    # Disconnect file processing and Shut down temporary directory
    # Close files
    Image.close()
    rescale_image.close()
    del Image
    del rescale_image
    files = None
    src_files = None
    img_grid = None
    shutil.rmtree(output_folder+'/output')