import os
import glob
import natsort
import rasterio
from rasterio.plot import show
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.patches import Patch
from shapely.ops import unary_union
from .codes.shoreline import Create_points, ExtrapolateOut, ExtrapolateIn
from .codes.shoreline import create_union_polygon, create_shoreline_change_points
from .codes.shoreline import merge_shoreline_change, linearring_to_polygon
from shapely.geometry import MultiPolygon, Polygon
import warnings
warnings.filterwarnings("ignore")
from qgis.core import QgsProject,QgsVectorLayer
def addMapLayer(path):
    project = QgsProject.instance()
    filename=path.split('/')[-1][0:-5]
    layer=QgsVectorLayer(path,filename,'ogr')
    if layer.isValid():
        project.instance().addMapLayer(layer)
    else:
        print("Failed to load layer")
def shoreline_analysis(dlg):
    #read combobox text
    shoreline_fp_1 = dlg.baselineShorelineComboBox.currentText()
    shoreline_fp_2=dlg.comparisonShorelineComboBox.currentText()
    # changeType=dlg.changeTypeComboBox.currentText()

    year1=int(dlg.baselineYearLineEdit.text())
    year2=int(dlg.comparisonYearLineEdit.text())
    total_years=year2-year1
    outputpath=dlg.outputSClineEdit.text()
    print(outputpath)
    #search for layernames marching the combobox text
    baseShoreline=QgsProject.instance().mapLayersByName(shoreline_fp_1)[0]
    comparisonShoreline=QgsProject.instance().mapLayersByName(shoreline_fp_2)[0]

    #obtain the directory and read each geojson with geopandas
    baseShoreline=baseShoreline.source().split('|')[0]
    comparisonShoreline=comparisonShoreline.source().split('|')[0]
    shl_past = gpd.read_file(baseShoreline).dropna().reset_index(drop=True)
    shl_present = gpd.read_file(comparisonShoreline).dropna().reset_index(drop=True)
    
    # Convert LinearRing to Polygon
    shl_past = linearring_to_polygon(shl_past)
    shl_present = linearring_to_polygon(shl_present)
    dlg.progressBar.setValue(20)
    
   
    retreat = gpd.overlay(shl_past, shl_present, how='difference', keep_geom_type=False)
    retreat.to_file(outputpath+'/retreat.json', driver='GeoJSON')
    addMapLayer(outputpath+'/retreat.json')
    
    growth = gpd.overlay(shl_present, shl_past, how='difference', keep_geom_type=False)
    growth.to_file(outputpath+'/growth.json', driver='GeoJSON')
    addMapLayer(outputpath+'/growth.json')
    # Export growth and retreat geometry to GeoJSON
    
    

    # # Create union polygon from geodata of growth area 
    growth_poly = create_union_polygon(growth)

    # # Create union polygon from geodata of retreat area 
    retreat_poly = create_union_polygon(retreat)

    # Create shoreline change as points along shoreline
    growth_shoreline_change = create_shoreline_change_points(shl_present, growth_poly)
    retreat_shoreline_change = create_shoreline_change_points(shl_present, retreat_poly)
    dlg.progressBar.setValue(40)
    # Export shoreline growth and retreat to GeoJSON
    growth_shoreline_change.to_file(outputpath+'/growth_points.json', driver='GeoJSON')
    retreat_shoreline_change.to_file(outputpath+'/retreat_points.json', driver='GeoJSON')
    addMapLayer(outputpath+'/growth_points.json')
    addMapLayer(outputpath+'/retreat_points.json')
    # Calculate total year
    # total_year = int(shoreline_fp[-1][-21:-17]) - int(shoreline_fp[0][-21:-17]) + 1

    # Create shoreline change
    change_distance = merge_shoreline_change(growth_shoreline_change, retreat_shoreline_change)
    dlg.progressBar.setValue(80)
    shoreline_change = retreat_shoreline_change.drop(columns=['change_m'])
    shoreline_change['total change_m'] = change_distance
    shoreline_change['rate per year_m'] = (shoreline_change['total change_m']/total_years).round(2)
    # print(outputpath)
    # # Export shoreline change to GeoJSON

    shoreline_change.to_file(outputpath+'/shorelineChange.json', driver='GeoJSON')
    addMapLayer(outputpath+'/shorelineChange.json')
    dlg.progressBar.setValue(100)