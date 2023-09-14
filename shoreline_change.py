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

def shoreline_analysis(dlg):
    shoreline_fp_1 = 
    shl_past = gpd.read_file(shoreline_fp[0]).dropna().reset_index(drop=True)
    shl_present = gpd.read_file(shoreline_fp[-1]).dropna().reset_index(drop=True)