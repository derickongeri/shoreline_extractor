import pip

packages_to_install = ['folium', 'rasterio', 'rtree', 'pandas', 'geopandas', 'geemap', 'mapclassify', 'contextily', 'matplotlib_scalebar', 'opencv-python', 'natsort', 'scikit-learn']

for package in packages_to_install:
    try:
        pip.main(["install", package])
        print(f"{package} installed Successfully")
    except:
        print(f"Error installing {package}")