import sys
import subprocess

packages_to_install = ['folium', 'rasterio', 'rtree', 'pandas', 'geopandas', 'geemap', 'mapclassify', 'contextily', 'matplotlib_scalebar', 'opencv-python', 'natsort', 'scikit-learn']

for package in packages_to_install:
    try:
        subprocess.run(["pip","install", package],check=True)
        print(f"{package} installed Successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error installing the dependencies: {e}")