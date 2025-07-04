import geopandas as gpd
import rasterio
from rasterio.features import shapes
from shapely.geometry import Polygon
import numpy as np
from shapely.ops import unary_union
import pandas as pd
import os


# Open the raster and read data into memory
def process_raster_and_points(raster_path, points, ID, ScenarioID):
    print(raster_path)
    with rasterio.open(raster_path) as src:
        asc_crs = points.crs
        asc_transform = src.transform
        asc_data = src.read(1)  # Reading the first band (raster data into memory)

    # Step 1: Convert the non-zero raster values to polygons
    mask = asc_data != 0  # Mask non-zero values
    shapes_list = list(shapes(asc_data, mask=mask, transform=asc_transform))

    # Step 2: Convert raster shapes into polygons
    polygons = [Polygon(shape['coordinates'][0]) for shape, value in shapes_list if value != 0]

    # Step 3: Merge all polygons into a single polygon (or multipolygon)
    merged_polygon = unary_union(polygons)

    # Create a GeoDataFrame for the polygon
    polygon_gdf = gpd.GeoDataFrame(geometry=[merged_polygon], crs=asc_crs)

    # Step 4: Perform a spatial join to get points within the polygon
    points_in_polygon = points[points.geometry.within(merged_polygon)]

    # Count the number of points inside the polygon
    number_of_points = len(points_in_polygon)

    # Step 5: Add the point count as a new attribute in the polygon GeoDataFrame
    polygon_gdf['NUMPOINTS'] = number_of_points
    polygon_gdf['USTABILEFJELLID'] = ID
    polygon_gdf['SCENARIOID'] = ScenarioID

    return polygon_gdf


def run_out_consequences(out_folder_SLBL, output_shapefile_path, Housing_shp_file, IDs, ScenarioIDs):

    print(out_folder_SLBL)
    out_folder_SLBL = r'C:\\Users\\Fiolleau_Sylvain\\Desktop'
    points = gpd.read_file(Housing_shp_file)


    combined_polygons = gpd.GeoDataFrame(columns=['geometry', 'USTABILEFJELLID', 'SCENARIOID', 'NUMPOINTS'],
                                             crs=points.crs)  # Initialize with a default CRS
    # List of raster files to process
    with os.scandir(out_folder_SLBL) as it:
        for entry in it:
            if entry.name.startswith('.') or entry.is_file():

                continue  # skips if its not a site
            if IDs: ## check if choose only some site with IDS
                if entry.name not in IDs:

                    continue
            with os.scandir(os.path.join(out_folder_SLBL, entry.name)) as it2:
                for entry2 in it2:
                    print(entry2.name)
                    if entry2.name.startswith('.') or entry.is_file():

                        continue  # skips if its not a scenario
                    if ScenarioIDs: ## check if choose only some scenarios
                        if entry2.name not in ScenarioIDs:

                            continue
                    if not os.path.exists(
                            os.path.join(out_folder_SLBL, entry.name, entry2.name, 'Outputs', 'com1DFA', 'peakFiles')):

                        continue  # skips if results are missing

                    print('Evaluating consequences, site: ' + entry.name + ' - ' + entry2.name, flush=True)

                    # find name of simulation output
                    with os.scandir(
                            os.path.join(out_folder_SLBL, entry.name, entry2.name, 'Outputs', 'com1DFA', 'peakFiles')) as it3:
                        for entry3 in it3:
                            if "pft.asc" in entry3.name and 'aux' not in entry3.name:
                                break  # stops the loop when the simulation output is found
                    # Step 6: Process each raster and append polygons to the combined GeoDataFrame
                    print(entry3.path, entry.name, entry2.name)
                    raster_polygon_gdf = process_raster_and_points(entry3.path, points, entry.name, entry2.name)
                    combined_polygons = pd.concat([combined_polygons, raster_polygon_gdf], ignore_index=True)


    # Step 7: Save the combined polygons (with point count) to a shapefile or GeoPackage
    print(combined_polygons)
    combined_polygons.to_file(output_shapefile_path, driver='GPKG')
    print(f"Polygon boundary with point count saved to {output_shapefile_path}")
