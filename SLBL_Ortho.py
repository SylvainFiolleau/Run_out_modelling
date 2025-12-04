import numpy as np
import rasterio
import geopandas as gpd
from copy import deepcopy
from rasterio.transform import from_origin
import math
import os
import sys
sys.path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))


## reprojecting part
from rasterio import warp
from rasterio.enums import Resampling
from rasterio.crs import CRS
from rasterio.features import geometry_mask
from rasterio.windows import from_bounds
import Correct_Thickness as CT



def read_vrt_window_resampled(vrt_path, bbox, target_resolution):
    with rasterio.open(vrt_path) as dataset:
        # Convert the bounding box to a window in the VRT's coordinate system
        window = from_bounds(*bbox, transform=dataset.transform)

        # Read the portion of the VRT specified by the window
        data_window = dataset.read(1, window=window)  # Assuming we're reading the first band

        # Calculate the original resolution (pixel size)
        original_resolution_x = dataset.transform[0]  # Pixel size in the x direction
        original_resolution_y = -dataset.transform[4]  # Pixel size in the y direction (negative)

        # Calculate the scaling factor based on the target resolution
        scale_factor_x = target_resolution / original_resolution_x
        scale_factor_y = target_resolution / original_resolution_y

        # Calculate the new shape (width and height) for the resampled data
        new_width = int(data_window.shape[1] / scale_factor_x)
        new_height = int(data_window.shape[0] / scale_factor_y)

        # Resample the windowed data to the new resolution
        data_resampled = dataset.read(
            1,
            window=window,
            out_shape=(new_height, new_width),
            resampling=Resampling.nearest  # You can choose other resampling methods (e.g., nearest)
        )
        new_transform = dataset.window_transform(window) * dataset.transform.scale(
            (window.width / new_width), (window.height / new_height)
        )
        print(new_transform)

        new_meta = dataset.meta.copy()
        new_meta.update({
            'height': new_height,
            'width': new_width,
            'transform': new_transform
        })
        # Calculate the new transform for the resampled data

        return data_resampled, new_transform, new_meta



def extend_bbox(source_area, factor=10):
    # Load the shapefile


    # Get the bounding box as a tuple (minx, miny, maxx, maxy)
    minx, miny, maxx, maxy = source_area.bounds

    # Calculate the center of the original bounding box
    center_x = (minx + maxx) / 2
    center_y = (miny + maxy) / 2

    # Extend the bounding box by the factor (e.g., 10 times)
    width = (maxx - minx) * factor
    height = (maxy - miny) * factor

    # Calculate new bounding box by extending equally around the center
    new_minx = center_x - width / 2
    new_maxx = center_x + width / 2
    new_miny = center_y - height / 2
    new_maxy = center_y + height / 2

    # Return the extended bounding box
    return new_minx, new_miny, new_maxx, new_maxy

def read_raster_with_mask(dtm_path, sources_areas, resample_resolution = None):

    if resample_resolution is not None:
        # Extract and extend the bounding box by 10 times
        extended_bbox = extend_bbox(sources_areas, factor=10)
        data, dem_new_transform, dem_meta = read_vrt_window_resampled(dtm_path, extended_bbox, resample_resolution)
        # data_from_vrt, dem_meta = FRR.read_vrt_window(dtm_path, extended_bbox)
        print('done with reading')
        # data, dem_meta = FRR.resample_raster_to_resolution(data_from_vrt,dem_meta, target_resolution)
        dtm_res = resample_resolution
        buffer_distance = 3 * dtm_res  # Buffer distance, scaled by the raster's pixel size
        buffered_areas = sources_areas.buffer(buffer_distance, resolution=2)  # Create buffer around polygons
        print('done with reading2')

        # Prepare output shape (round to nearest integer) based on the resampled data
        out_shape = (dem_meta['height'], dem_meta['width'])  # Now it directly reflects the new shape after resampling

        # Create a mask for the windowed data based on the polygon, using the updated transform
        mask = geometry_mask([sources_areas],
                             transform=dem_new_transform,  # Use the new transform after resampling
                             invert=True,
                             out_shape=out_shape)  # The mask shape must match the resampled shape
        # # Define the metadata for the output raster (same as input raster)
        with rasterio.open(dtm_path) as src:
            metadata = src.meta.copy()
        metadata.update({
            'driver': 'GTiff',
            'count': 2,  # Two bands (data and mask)
            'dtype': 'float32',  # Adjust depending on the data type
            'width': out_shape[1],
            'height': out_shape[0],
            'crs': src.crs,
            'transform': dem_new_transform
        })
        print('done with reading3')

    else:

        with rasterio.open(dtm_path) as src:

            dtm_res = src.res
            buffer_distance = 3 * dtm_res[0]  # Buffer distance, scaled by the raster's pixel size
            buffered_areas = sources_areas.buffer(buffer_distance, resolution=2)  # Create buffer around polygons

            # Step 2: Calculate bounding box of the buffered geometries
            minx_b, miny_b, maxx_b, maxy_b = buffered_areas.bounds  # Bounding box of the buffered geometries
            # minx_b, miny_b, maxx_b, maxy_b = np.round([minx_b, miny_b, maxx_b, maxy_b])
            # Create a window for the region of interest
            window = from_bounds(minx_b, miny_b, maxx_b, maxy_b, transform=src.transform)

            # Read only the window (subset of the VRT) from the raster
            data = src.read(1, window=window)  # Assuming we're reading the first band (change as needed)

            out_shape = (int(window.height.round()), int(window.width.round()))
            # Create a mask for the windowed data based on the polygon
            # mask = geometry_mask([geom for geom in sources_areas.geometry],
            #                      transform=src.window_transform(window),
            #                      invert=True,
            #                      out_shape=(window.height, window.width))
            mask = geometry_mask([sources_areas],
                                 transform=src.window_transform(window),
                                 invert=True,
                                 out_shape=out_shape)

            # # Define the metadata for the output raster (same as input raster)
            metadata = src.meta.copy()
            metadata.update({
                'driver': 'GTiff',
                'count': 2,  # Two bands (data and mask)
                'dtype': 'float32',  # Adjust depending on the data type
                'width': out_shape[1],
                'height': out_shape[0],
                'crs': src.crs,
                'transform': src.window_transform(window)
            })
        # output_raster_path = r"C:\Users\Fiolleau_Sylvain\Documents\Avaframe_simulations\test_WQGIS\slbl_15011_1501101_base__2.tif"
        # # Save the data and mask to the new raster
        # with rasterio.open(output_raster_path, 'w', **metadata) as dest:
        #     # Write the data (first band)
        #     dest.write(data, 1)
        #     # Write the mask (second band, inverted to mask the area outside)
        #     dest.write(mask.astype('float32'), 2)

    return data, mask, buffered_areas

# Step 2: Extract Points along Lines (every 1 meter)
# Create a function to extract points along lines at regular distances
def extract_points_along_line(line, distance=1.0):
    points = []
    length = line.length
    num_points = int(length // distance)

    for i in range(num_points):
        # Calculate point along the line at regular intervals
        point = line.interpolate(i * distance)
        points.append(point)

    return points
def reproject(input_file, target_crs, data_type='shapefile'):
    """
    Reproject either a shapefile or a raster file to a new CRS and return the reprojected data.

    Parameters:
    - input_file: Path to the input file (shapefile or raster).
    - target_crs: The target CRS (could be a string, EPSG code, or CRS object).
    - data_type: Specify 'shapefile' for shapefiles or 'raster' for raster files.

    Returns:
    - For shapefile: Returns a GeoDataFrame with reprojected data.
    - For raster: Returns a tuple (reprojected_raster_array, new_transform, new_crs).
    """

    if data_type.lower() == 'shapefile':
        # Reproject a shapefile
        gdf = gpd.read_file(input_file)

        # Reproject to the target CRS
        gdf = gdf.to_crs(target_crs)

        # Return the reprojected GeoDataFrame
        return gdf

    elif data_type.lower() == 'raster':
        # Reproject a raster file
        with rasterio.open(input_file) as src:
            # Get the source CRS and transform
            src_crs = src.crs
            transform = src.transform
            width = src.width
            height = src.height
            dtype = src.dtypes[0]

            # Prepare for the reprojected output file
            target_crs_obj = CRS.from_user_input(target_crs)
            output_transform, output_width, output_height = warp.calculate_default_transform(
                src.crs, target_crs_obj, width, height, *src.bounds)

            # Create an empty array to hold the reprojected raster
            reprojected_raster = np.empty((output_height, output_width), dtype=dtype)

            # Perform the reproject operation
            warp.reproject(
                source=rasterio.band(src, 1), destination=reprojected_raster,
                src_transform=transform, src_crs=src_crs,
                dst_transform=output_transform, dst_crs=target_crs_obj,
                resampling=Resampling.nearest)

            # Return the reprojected raster array, transform, and CRS
            return reprojected_raster, output_transform, target_crs_obj
    else:
        raise ValueError("Invalid data_type. Please specify 'shapefile' or 'raster'.")

def SLBL(grid_dem, grid_mask, tol, maxt, maxv, z_min, cellSize, nb_neigh, stop, criteria, inverse):
    # initilizes thickness and difference grids
    s = grid_dem.shape
    grid_thickn = np.zeros(s)
    grid_diff = np.ones(s)

    # Creates a matrice to store the values of the neighbouring cells in the previous iteration
    mat_neigh = np.zeros(s)
    mat_neigh = np.expand_dims(mat_neigh, axis=2)
    if nb_neigh == 4:
        mat_neigh = np.tile(mat_neigh, (1, 1, 4))
    else:
        mat_neigh = np.tile(mat_neigh, (1, 1, 8))

    # Creates a matrice where the proposed value and previous value are stored for comparison
    mat_comp = np.zeros(s)
    mat_comp = np.expand_dims(mat_comp, axis=2)
    mat_comp = np.tile(mat_comp, (1, 1, 2))

    # Initiate the slbl grid (altitude)
    grid_slbl = deepcopy(grid_dem)

    nb_iter = 0
    volume = 0.

    # The SLBL strarts here
    while np.nanmax(grid_diff) > stop and np.nanmax(grid_thickn) < maxt and volume < maxv:
        nb_iter = nb_iter + 1
        grid_thickn_prev = deepcopy(grid_thickn)
        grid_slbl_prev = deepcopy(grid_slbl)

        # writes the values of the neighbourings cells in the 3rd dimension of the matrix
        mat_neigh[:-1, :, 0] = grid_slbl_prev[1:, :]
        mat_neigh[1:, :, 1] = grid_slbl_prev[:-1, :]
        mat_neigh[:, :-1, 2] = grid_slbl_prev[:, 1:]
        mat_neigh[:, 1:, 3] = grid_slbl_prev[:, :-1]

        # diagonals
        if nb_neigh == 8:
            mat_neigh[:-1, :-1, 4] = grid_slbl_prev[1:, 1:]
            mat_neigh[:-1, 1:, 5] = grid_slbl_prev[1:, :-1]
            mat_neigh[1:, 1:, 6] = grid_slbl_prev[:-1, :-1]
            mat_neigh[1:, :-1, 7] = grid_slbl_prev[:-1, 1:]

        if criteria == 'minmax':
            mat_max = np.nanmax(mat_neigh, axis=2)
            mat_min = np.nanmin(mat_neigh, axis=2)
            mat_mean = (mat_max + mat_min) / 2
        elif criteria == 'average':
            mat_mean = np.mean(mat_neigh, axis=2)
        mat_mean = mat_mean + tol
        if np.isfinite(z_min):
            mat_mean = np.maximum(mat_mean, z_min)

        mat_comp[:, :, 0] = mat_mean
        mat_comp[:, :, 1] = grid_slbl

        # Check if the new value should be kept
        if inverse == 'true':
            grid_slbl = np.nanmax(mat_comp, axis=2)
        else:
            grid_slbl = np.nanmin(mat_comp, axis=2)

        # Replaces the values of the SLBL by the original values outside the masked area
        grid_slbl[~grid_mask] = grid_dem[~grid_mask]

        grid_thickn = np.absolute(grid_dem - grid_slbl)
        grid_diff = np.absolute(grid_thickn - grid_thickn_prev)

        volume = (np.sum(grid_thickn) * cellSize ** 2)


    return grid_slbl, grid_thickn, nb_iter

def prepare_and_compute_SLBL(sources_areas, dtm_path, out_folder, tolerance, IDs, SCENARIOID, resample_resolution = None):
    with rasterio.open(dtm_path) as src_dtm:
        dtm_transform = src_dtm.transform
        dtm_crs = src_dtm.crs
        dtm_res = src_dtm.res
    print(dtm_path, dtm_crs, dtm_res)
    if dtm_crs != sources_areas.crs:
        print('reproject sources areas to match dtm crs')
        return

    sources_areas['SourceDip'] = None
    sources_areas['SourceDir'] = None
    sources_areas['Volume'] = None
    sources_areas['Tol'] = None

    i = 0
    print(sources_areas.shape[0])
    for k in range(sources_areas.shape[0]):
        print(dtm_path, dtm_crs, dtm_res)

        sources_area = sources_areas.iloc[k]
        sources_area_geometry = sources_areas.geometry.iloc[k]  # Replace 0 with the desired index
        ust_id = sources_area['USTABILEFJELLID']
        scen_id = sources_area['SCENARIOID']
        print(ust_id)

        if (math.isnan(ust_id) or math.isnan(scen_id)):
            continue


        if IDs and SCENARIOID:
            if str(int(ust_id)) in IDs and str(int(scen_id)) in SCENARIOID:
                    print('Start configuring : ' + str(int(ust_id)) + ' / ' + str(int(scen_id)))
            else:
                continue
        elif IDs and not SCENARIOID:
            if str(int(ust_id)) in IDs:
                    print('Start configuring : ' + str(int(ust_id)) + ' / ' + str(int(scen_id)))
            else:
                continue
        ust_id = int(ust_id)
        scen_id = int(scen_id)
        if os.path.exists(os.path.join(out_folder, f"slbl_{ust_id}_{scen_id}_ortho.tif")):
            print(f"{ust_id}-{scen_id} ortho exist, skip")
            continue

        try:
            # write progress in a text file (usefull for debuging / finding a problematic site for which SLBL struggle)
            with open(os.path.join(out_folder, 'log.txt'), 'a') as f:
                f.write(str(ust_id) + ' - ' + str(scen_id) + '\n')
        except:
            print("processing source polygon without USTABILEFJELLID attributes")
        print(f'computation slbl for {ust_id}-{scen_id} ')
        # Convert polygons to lines
        lines_layer = sources_areas.boundary.iloc[k]
        # Extract Points along each line in the lines_layer

        all_points = []
        points = extract_points_along_line(lines_layer, distance=1.0)  # Adjust distance as needed
        all_points.extend(points)


        # Convert points to a GeoDataFrame
        points_gdf = gpd.GeoDataFrame(geometry=all_points)


        print(sources_area_geometry.area)
        # Step 3: Add Geometry Columns (X, Y Coordinates)
        points_gdf['xcoord'] = points_gdf.geometry.x
        points_gdf['ycoord'] = points_gdf.geometry.y
        # Read the DEM file
        with rasterio.open(dtm_path) as src_dtm:
            # Extract raster values at point locations
            raster_values = list(src_dtm.sample([(point.x, point.y) for point in points_gdf.geometry]))
            dtm_res = src_dtm.res
            # Add elevation values to the points DataFrame
            points_gdf['zcoord'] = [val[0] for val in raster_values]  # Assuming single-band raster


        # Step 5: Filter Points Based on Z-Values (e.g., Quantiles for Foot and Crest)
        # Calculate quantiles of the 'zcoord' field
        first_quartile = points_gdf['zcoord'].quantile(0.25)
        third_quartile = points_gdf['zcoord'].quantile(0.75)


        # Filter points for foot (zcoord <= first_quartile) and crest (zcoord >= third_quartile)
        points_foot = points_gdf[points_gdf['zcoord'] <= first_quartile]
        points_crest = points_gdf[points_gdf['zcoord'] >= third_quartile]



        # Step 6: Compute Statistics for Foot and Crest Points
        # Example of calculating mean x, y, z for foot and crest
        foot_stats = {
            'mean_x': points_foot['xcoord'].mean(),
            'mean_y': points_foot['ycoord'].mean(),
            'mean_z': points_foot['zcoord'].mean()
        }
        foot_stats_min = {
            'min_x': points_foot['xcoord'].min(),
            'min_y': points_foot['ycoord'].min(),
            'min_z': points_foot['zcoord'].min()
        }

        crest_stats = {
            'mean_x': points_crest['xcoord'].mean(),
            'mean_y': points_crest['ycoord'].mean(),
            'mean_z': points_crest['zcoord'].mean()
        }

        # Print statistics (example)
        print("Foot Stats:", foot_stats)
        print("Crest Stats:", crest_stats)

        # Step 7: Calculate Dip and Direction (based on foot and crest statistics)
        dx = foot_stats['mean_x'] - crest_stats['mean_x']
        dy = foot_stats['mean_y'] - crest_stats['mean_y']
        dz = crest_stats['mean_z'] - foot_stats['mean_z']

        # Calculate horizontal distance (planimetric)
        dxy = np.sqrt(dx ** 2 + dy ** 2)

        # H/L ratio and dip
        HoL = dz / dxy
        dip = np.degrees(np.arctan(HoL))
        print("Dip (in degrees):", dip)

        # Direction (azimuth)
        dir = np.degrees(np.arccos(np.dot([dx,dy],[0,1])/dxy))
        if dx < 0:
            dir = 360 - dir

        print("Direction (in degrees):", dir)

        # Step 8: Rasterize Polygon (mask creation from the polygon)
        # For creating a mask raster based on the polygon geometry
        from rasterio.features import geometry_mask
        if resample_resolution is not None:
            #extract mask and clip dtm
            clipped_raster, mask, buffered_areas = read_raster_with_mask(dtm_path, sources_area_geometry, resample_resolution)
            # clipped_raster = clipped_raster[0,:,:]

        else:
            clipped_raster, mask, buffered_areas = read_raster_with_mask(dtm_path, sources_area_geometry)
        #calculate SLBL based on Pierrick's script
        #Load dtm via gdal
        grid_dem  = clipped_raster#numpy array
        grid_mask = mask
        grid_mask = (grid_mask[:] > 0)#boolean
        if tolerance == -99:
            tol       = 2*(1-2**0.5)*dz*dtm_res[0]**2/dxy**2 #intermediate tolerance, as in slbl script, simply half of max
            print('tolerance tune in respect to source area')
        else:
            tol = tolerance

        maxt      = np.inf
        maxv      = np.inf
        z_min     = foot_stats_min['min_z']
        nb_neigh  = 4
        stop      = 0.00001
        criteria  = 'average'
        inverse   = 'false'
        if resample_resolution == None:
            cellSize  = dtm_res[0]
        else:
            cellSize = resample_resolution
        print(grid_dem.shape, grid_mask.shape, cellSize)
        print('StartSLBL')
        grid_slbl, grid_thickn, nb_iter = SLBL(grid_dem, grid_mask, tol, maxt, maxv, z_min, cellSize,
                                               nb_neigh, stop, criteria, inverse)

        volume = (np.sum(grid_thickn) * cellSize ** 2)

        print('Volume:' + str(volume) + ', Tol:' + str(tol))


        # Extracting necessary data
        Lx = dtm_res[0]  # DTM X resolution in meters
        Ly = dtm_res[1]  # DTM Y resolution in meters
        ncols = grid_thickn.shape[1]  # number of columns in the raster
        nrows = grid_thickn.shape[0]  # number of rows in the raster
        xll, ymin, xur, yur = buffered_areas.bounds
        # Output file path
        fn = os.path.join(out_folder, f"slbl_{ust_id}_{scen_id}")

        orthogonal_thickness = CT.Calc_Norm(grid_slbl, grid_thickn, cellSize)
        volume_ortho = (np.sum(orthogonal_thickness) * cellSize ** 2)
        print('Volume Ortho:' + str(volume_ortho))
        # Creating an empty raster for 'thickn'
        with rasterio.open(
            fn + ".tif", 'w', driver='GTiff', width=ncols, height=nrows, count=1, dtype='float32',
            crs=dtm_crs, transform=from_origin(xll, yur, Lx, Ly)
        ) as dst:
            dst.write(grid_thickn, 1)  # Writing the grid_thickn data to the raster

        with rasterio.open(
            fn + "_ortho.tif", 'w', driver='GTiff', width=ncols, height=nrows, count=1, dtype='float32',
            crs=dtm_crs, transform=from_origin(xll, yur, Lx, Ly)
        ) as dst:
            dst.write(orthogonal_thickness, 1)  # Writing the grid_thickn data to the raster


        # Creating an empty raster for 'slbl'
        with rasterio.open(
            fn + "_base.tif", 'w', driver='GTiff', width=ncols, height=nrows, count=1, dtype='float32',
            crs=dtm_crs, transform=from_origin(xll, yur, Lx, Ly)
        ) as dst:
            dst.write(grid_slbl, 1)  # Writing the grid_slbl data to the raster
