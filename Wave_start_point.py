import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import reproject, Resampling

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from scipy.ndimage import label, center_of_mass
import pandas as pd
import os



def compute_wave_start_point(runout_path, water_path, ID, scenarioID):
    # Open the second raster (extent provider)
    # print(runout_path, water_path)
    WaveId = 1
    with rasterio.open(runout_path) as src_volume:
        bounds = src_volume.bounds

        # print(bounds)
        volume_data = src_volume.read(1)
        volume_transform = src_volume.transform
        volume_crs = src_volume.crs if src_volume.crs else None  # Get CRS from volume raster, if it exists
        print(volume_crs)
        # Open the first raster (to be cropped)
        with rasterio.open(water_path) as src_water:
            # Create a window based on the bounds of the second raster
            window = from_bounds(bounds.left, bounds.bottom, bounds.right, bounds.top, src_water.transform)
            # print('window ', window)
            # Read the data within the window
            water_mask_init = src_water.read(window=window)

            water_mask_init = water_mask_init[0]
            src_water.window_bounds(window)
            water_transform_init = src_water.transform
            water_transform = water_transform_init * rasterio.Affine.translation(window.col_off, window.row_off)
            water_crs = src_water.crs
            res_factor = np.sqrt(water_transform[0] - volume_transform[0])

            print(water_mask_init.shape)
            # resampling water raster so it matches volume resolution
            print(f'resolution of {ID} volume is {volume_transform[0]}')
            resample_height = int(water_mask_init.shape[0] * (water_transform[0] / volume_transform[0]))
            resample_width = int(water_mask_init.shape[1] * (water_transform[0] / volume_transform[0]))

            new_transform = water_transform * water_transform.scale(
                (water_mask_init.shape[1] / resample_width), (water_mask_init.shape[0] / resample_height)
            )

            # Create an empty array for the resampled data
            water_mask = np.empty((resample_height, resample_width), dtype=water_mask_init.dtype)

            # Reproject and resample the data
            reproject(
                source=water_mask_init,
                destination=water_mask,
                src_transform=water_transform,
                dst_transform=new_transform,
                src_crs=src_water.crs,
                dst_crs=src_water.crs,
                resampling=Resampling.nearest  # or use nearest, cubic, etc.
            )
            window_bounds = src_water.window_bounds(window)
            window_bounds = tuple(round(val, 2) for val in window_bounds)
            # Ensure the two rasters have the same spatial extent (same bounding box)
            rounded_bounds = rasterio.coords.BoundingBox(*map(lambda x: round(x, 2), bounds))
            # windows_bounds = rasterio.coords.BoundingBox(*map(lambda x: round(x, 1), window_bounds))
            if window_bounds != rounded_bounds:
                print(window_bounds, src_volume.bounds)
                raise ValueError("The bounding boxes of the two rasters do not match.")

            water_mask = (water_mask == 0)  # Convert to boolean mask for water (0 = water)
            min_pix = 15 * res_factor
            filtered_water_mask = np.zeros_like(water_mask, dtype=bool)
            Ini_labeled_mask, Ininum_features = label(water_mask)
            for feature_num in range(1, Ininum_features + 1):
                feature_mask = (Ini_labeled_mask == feature_num)
                if np.sum(feature_mask) >= min_pix:
                    filtered_water_mask[feature_mask] = True

            labeled_mask, num_features = label(filtered_water_mask)
            gdf = gpd.GeoDataFrame(columns=['Id', 'ScenarioID', 'WaveID', 'Volume', 'X', 'Y', 'geometry'])

            for region_id in range(1, num_features + 1):
                # Initialize total volume for the entire water body
                total_volume = 0
                weighted_sum_x = 0
                weighted_sum_y = 0
                total_weight = 0
                totvol1 = 0
                # Loop through all pixels in the water mask to compute the total volume and center of mass
                # rows, cols = np.where(water_mask_t)  # Get row and column indices of all water pixels
                rows, cols = np.where(labeled_mask == region_id)
                for row, col in zip(rows, cols):
                    # Convert pixel coordinates (row, col) to world coordinates using the water transform
                    x, y = new_transform * (col, row)

                    # Now find the corresponding pixel in the volume raster
                    # We use the inverse of the volume transform to convert from world coordinates to raster indices
                    col_vol, row_vol = ~volume_transform * (x, y)

                    # Check if the column and row are within the bounds of the volume raster
                    if 0 <= col_vol < src_volume.width and 0 <= row_vol < src_volume.height:
                        vol1 = volume_data[int(row_vol), int(col_vol)]
                        totvol1 += vol1
                        volume_value = volume_data[int(row_vol), int(col_vol)] * (
                                    volume_transform[0] * volume_transform[0])

                        total_volume += volume_value  # Add the volume value

                        # Weighted sums for center of mass
                        weighted_sum_x += volume_value * x
                        weighted_sum_y += volume_value * y
                        total_weight += volume_value

                        # Calculate center of mass for the entire region
                if total_weight > 0:  # Avoid division by zero
                    com_x = weighted_sum_x / total_weight
                    com_y = weighted_sum_y / total_weight
                else:
                    com_x, com_y = None, None  # In case there's no volume
                if com_x != None and com_y != None:
                    # Create a GeoDataFrame with a single point representing the center of mass for the entire water region
                    new_point = gpd.GeoDataFrame({
                        'Id': ID,  # Single ID for all contiguous water areas
                        'ScenarioID': scenarioID,  # Scenario ID, modify if needed
                        'WaveID': WaveId,  # Wave ID, modify if needed
                        'Volume': [total_volume],  # Total volume for the entire water body
                        'X': com_x,
                        'Y': com_y,
                        'geometry': [Point(com_x, com_y)]  # Center of mass as a point geometry
                    })
                    WaveId += 1
                    gdf = pd.concat([gdf, new_point], ignore_index=True)
                    print(totvol1)

            if not gdf.empty:
                # Set the coordinate reference system (CRS) to match the volume raster
                gdf.set_crs(water_crs, inplace=True)
                print(total_volume)
                print(f"Center of mass for {ID}, scenario {scenarioID} computed")

            else:
                print(f"{ID}, scenario {scenarioID} not in water")

        return gdf


def Launch_wave_start_points(fkb_water_path, input_RunOuts_path, output_shapefile, IDs, ScenarioIDs):
    ## initialize geodataframe for starting points
    starting_points = gpd.GeoDataFrame(columns=['Id', 'ScenarioID', 'WaveID', 'Volume', 'X', 'Y', 'geometry'])
    gdf = gpd.GeoDataFrame(columns=['Id', 'ScenarioID', 'WaveID', 'Volume', 'X', 'Y', 'geometry'])
    # List of raster files to process
    with os.scandir(input_RunOuts_path) as it:
        for entry in it:
            if entry.name.startswith('.') or entry.is_file():
                continue  # skips if its not a site
            if IDs:  ## check if choose only some site with IDS
                if entry.name not in IDs:
                    continue
            with os.scandir(os.path.join(input_RunOuts_path, entry.name)) as it2:
                for entry2 in it2:
                    if entry2.name.startswith('.') or entry.is_file():
                        continue  # skips if its not a scenario
                    if ScenarioIDs:  ## check if choose only some scenarios
                        if entry2.name not in ScenarioIDs:
                            continue
                    if not os.path.exists(
                            os.path.join(input_RunOuts_path, entry.name, entry2.name, 'Outputs', 'com1DFA',
                                         'peakFiles')):
                        continue  # skips if results are missing

                    print('Evaluating consequences, site: ' + entry.name + ' - ' + entry2.name, flush=True)

                    # find name of simulation output
                    with os.scandir(
                            os.path.join(input_RunOuts_path, entry.name, entry2.name, 'Outputs', 'com1DFA',
                                         'peakFiles')) as it3:
                        for entry3 in it3:
                            if "FT.asc" in entry3.name and 'aux' not in entry3.name:
                                print('simulation found')
                                break  # stops the loop when the simulation output is found
                    # Step 6: Process each raster and append polygons to the combined GeoDataFrame

                    new_starting_point = compute_wave_start_point(entry3.path, fkb_water_path, entry.name, entry2.name)
                    ## append new points:
                    if not new_starting_point.empty:
                        gdf = pd.concat([gdf, new_starting_point], ignore_index=True)

    starting_points = gpd.GeoDataFrame(gdf, geometry='geometry')
    # Save Shapefile
    starting_points.to_file(output_shapefile)
    print('---done computing wave start points---')

