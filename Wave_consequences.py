import rasterio
from rasterio.merge import merge
from rasterio.plot import show
import numpy as np
import os
from rasterio.features import shapes
import geopandas as gpd
# Example: Open raster with rasterio
import glob
from scipy.ndimage import binary_erosion, generate_binary_structure, binary_dilation
from skimage.morphology import disk
import math
from rasterio.mask import mask
from shapely.geometry import box, mapping, Polygon, GeometryCollection
from scipy.ndimage import binary_dilation, distance_transform_edt
from rasterio.transform import from_origin
import pandas as pd

def contract_runup(input_raster_path, output_raster_path, radius_pixels=2):
    # Open the input raster
    with rasterio.open(input_raster_path) as src:
        # Read the raster data into an array
        data = src.read(1)

        # Convert the data into a binary mask (assuming non-zero values are water areas)
        binary_mask = data > 0

        # Create a structuring element for erosion (disk shape, with radius equal to radius_pixels)
        structuring_element = disk(radius_pixels)

        # Perform erosion (contract the water areas by 2 pixels)
        contracted_mask = binary_erosion(binary_mask, structure=structuring_element)

        # Replace the non-zero areas with the contracted mask
        contracted_data = np.where(contracted_mask, data, None)

        # Copy metadata and update to reflect new shape
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": contracted_data.shape[0],
            "width": contracted_data.shape[1],
            "transform": src.transform
        })
    count_value = np.count_nonzero(contracted_data)
    # Write the contracted raster to the output path
    if count_value > 0:
        with rasterio.open(output_raster_path, 'w', **out_meta) as dst:
            dst.write(contracted_data, 1)
        return 1
    else:
        return 0

    print(f"Contracted raster saved to {output_raster_path}")

def merge_rasters(raster_paths, out_path, method='max', range_vals=None, weights=None):
    src_files_to_mosaic = []

    # Open the rasters and append them to the list for mosaic
    for path in raster_paths:
        src = rasterio.open(path)
        src_files_to_mosaic.append(src)

    # Merge rasters
    mosaic, out_trans = merge(src_files_to_mosaic)
    for src in src_files_to_mosaic:
        no_data_value = src.nodata
        if no_data_value is not None:
            mosaic[mosaic == no_data_value] = None  # Replace NoData with 0


    # Apply method to combine raster values, e.g., maximum, average, etc.
    if method == 'max':
        mosaic = np.max(mosaic, axis=0)
    elif method == 'min':
        mosaic = np.min(mosaic, axis=0)
    elif method == 'mean':
        mosaic = np.mean(mosaic, axis=0)
    elif method == 'quantile' and weights is not None:
        mosaic = np.quantile(mosaic, q=weights, axis=0)
    else:
        raise ValueError(f"Method {method} not supported or weights missing for quantile")

    # Apply range restriction if specified
    if range_vals is not None:
        min_val, max_val = range_vals
        mosaic = np.clip(mosaic, min_val, max_val)

    # Prepare metadata for output
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[0],
        "width": mosaic.shape[1],
        "transform": out_trans
    })

    # Write the result to the output path
    with rasterio.open(out_path, "w", **out_meta) as dest:
        dest.write(mosaic, 1)  # Write the final mosaic to band 1


def raster_to_polygon(raster_path, out_path):
    with rasterio.open(raster_path) as src:
        # Read the first band of the raster
        image = src.read(1)
        # Mask out NoData values
        mask = image != src.nodata
        # Set NoData pixels to a constant value (e.g., NaN or any suitable placeholder)
        image = np.where(mask, image, np.nan)

        # Extract shapes (polygons) from the raster
        results = (
            {'properties': {'thickness': v}, 'geometry': s}
            for s, v in shapes(image, mask=mask, transform=src.transform)
        )
        # Create a GeoDataFrame from the results
        gdf = gpd.GeoDataFrame.from_features(list(results), crs=src.crs)
        gdf_filtered = gdf[gdf['thickness'] > 0]

        # Save the polygons to a GeoPackage (or other formats like shapefile)
        gdf_filtered.to_file(out_path, driver="GPKG")  # You can change the driver to 'ESRI Shapefile' if needed

        print(f"Polygonized output saved to: {out_path}")


def count_points_in_polygon(polygon_gdf, pointsI, indexes=[]):
    # Perform spatial join to count points within each polygon

    points = pointsI.drop(indexes)
    joined = gpd.sjoin(points, polygon_gdf, how="left", predicate="within")
    joined = joined[~joined['thickness'].isna()]
    joined = joined.drop_duplicates(subset=['oest_koo', 'nord_koo'])
    numPoints = joined['bosatte'][~joined['thickness'].isna()].sum()
    numPoints = joined['bosatte'].sum()

    Index = joined.index

    return numPoints, Index


def adjust_raster_to_bbox(src_raster_path, bbox, out_raster_path, nodata_value=None):
    with rasterio.open(src_raster_path) as src:
        # Get the original raster properties
        raster_bounds = src.bounds
        raster_width = src.width
        raster_height = src.height
        raster_transform = src.transform
        raster_nodata = src.nodata

        # Bounding box coordinates
        xmin, ymin, xmax, ymax = bbox

        # Check if the bounding box is larger than the raster
        if (xmax - xmin > raster_bounds[2] - raster_bounds[0]) or (ymax - ymin > raster_bounds[3] - raster_bounds[1]):
            # Padding case: Bounding box is larger than raster
            print("Padding raster to fit the bounding box.")

            # Calculate new raster dimensions
            pixel_width = raster_transform[0]  # The pixel width (resolution)
            pixel_height = abs(raster_transform[4])  # The pixel height (resolution)

            # Calculate the new number of pixels for the output raster
            new_width = int((xmax - xmin) / pixel_width)
            new_height = int((ymax - ymin) / pixel_height)

            # Create the new transform for the padded raster
            new_transform = from_origin(xmin, ymax, pixel_width, pixel_height)

            # Create an empty array for the new raster (padded with nodata)
            padded_image = np.full((new_height, new_width), nodata_value, dtype=np.float32)

            # Read the original raster data
            original_window = src.window(raster_bounds[0], raster_bounds[1], raster_bounds[2], raster_bounds[3])
            original_data = src.read(1, window=original_window)

            # Calculate the offset to place the original data inside the padded array
            offset_x = int((raster_bounds[0] - xmin) / pixel_width)
            offset_y = int((ymax - raster_bounds[3]) / pixel_height)

            # Place the original data inside the padded array
            padded_image[offset_y:offset_y + raster_height, offset_x:offset_x + raster_width] = original_data

            # Write the padded raster to the output file
            with rasterio.open(out_raster_path, 'w', driver='GTiff', height=new_height, width=new_width,
                               count=1, dtype=np.float32, crs=src.crs, transform=new_transform,
                               nodata=-9999) as dest:
                dest.write(padded_image, 1)
            print(f"Raster padded and saved to: {out_raster_path}")

        else:
            # Clipping case: Bounding box is smaller than raster
            print("Clipping raster to the bounding box.")

            # Ensure the bounding box is within the raster bounds
            xmin = max(xmin, raster_bounds[0])
            ymin = max(ymin, raster_bounds[1])
            xmax = min(xmax, raster_bounds[2])
            ymax = min(ymax, raster_bounds[3])

            # Get the window to clip the raster
            window = src.window(xmin, ymin, xmax, ymax)
            clipped_data = src.read(1, window=window)

            # Create a new transform for the clipped raster
            clipped_transform = src.window_transform(window)

            # Write the clipped raster to the output file
            with rasterio.open(out_raster_path, 'w', driver='GTiff', height=clipped_data.shape[0],
                               width=clipped_data.shape[1], count=1, dtype=np.float32, crs=src.crs,
                               transform=clipped_transform, nodata=None) as dest:
                dest.write(clipped_data, 1)
            print(f"Raster clipped and saved to: {out_raster_path}")

def clip_raster_by_extent(runup, raster, output_path, dtm_grid, max_ingressT = 500):
    # Round max_ingress to the nearest multiple of 50
    max_ingress = math.ceil(max_ingressT / 50) * 50

    with rasterio.open(dtm_grid) as src_dtm:
        # Get the bounds (xmin, ymin, xmax, ymax)
        xmin_dtm, ymin_dtm, xmax_dtm, ymax_dtm = src_dtm.bounds

    # Open the raster file to get its bounds
    with rasterio.open(runup) as src1:
        # Get the bounds (xmin, ymin, xmax, ymax)
        xmin, ymin, xmax, ymax = src1.bounds

    # Adjust the extent based on the max_ingress and buffer (2 * 50)
    xmin = xmin - (max_ingress + 2 * 50)
    xmax = xmax + (max_ingress + 2 * 50)
    ymin = ymin - (max_ingress + 2 * 50)
    ymax = ymax + (max_ingress + 2 * 50)

    if xmin < xmin_dtm:
        xmin = xmin_dtm
    if xmax > xmax_dtm:
        xmax = xmax_dtm
    if ymin < ymin_dtm:
        ymin = ymin_dtm
    if ymax > ymax_dtm:
        ymax = ymax_dtm

    # Clip raster by the extent
    # Define the bounding box for clipping
    adjust_raster_to_bbox(raster, [xmin, ymin, xmax, ymax], output_path, nodata_value=None)


def Add_Height(runup_pathIn, dtm_local, runup_pathOut):
    with rasterio.open(runup_pathIn) as src_runup, rasterio.open(dtm_local) as src_dtm_local:
        runup = src_runup.read(1)
        out_meta = src_runup.meta.copy()

        dtm = src_dtm_local.read(1)
        # print(runup.shape, dtm.shape)
        if runup.shape != dtm.shape:
            raise ValueError("Runup and DTM rasters must have the same shape.")

    Runup_Updated = runup + dtm
    out_meta.update({
        "nodata": np.nan,  # If you want to set a NoData value, you can specify it
    })
    with rasterio.open(runup_pathOut, 'w', driver='GTiff',
                       height=src_runup.height, width=src_runup.width,
                       count=1, dtype=np.float32, crs=src_runup.crs,
                       transform=src_runup.transform, nodata=np.nan) as dst:
        dst.write(Runup_Updated, 1)

def aggregate_polygons(polygon_path, out_path, crs=None):
    # Read polygon layer
    gdf = gpd.read_file(polygon_path)

    # Aggregate by dissolving all polygons into one and summing "thickness"
    gdf_agg = gdf.dissolve(by=None, aggfunc={"thickness": "sum"}).reset_index()
    if crs:
        gdf_agg = gdf_agg.set_crs(crs, allow_override=True)
    # Save output
    gdf_agg.to_file(out_path, driver="GPKG")
   # print(f"Aggregation done. Output saved to: {out_path}")

    return gdf_agg


def OnLand(runup_raster_path, dtm_raster_path, output_raster_path):
    with rasterio.open(runup_raster_path) as src_runup, rasterio.open(dtm_raster_path) as src_dtm:
        runup_data = src_runup.read(1)  # Read first band
        dtm_data = src_dtm.read(1)  # Read first band
        if runup_data.shape != dtm_data.shape:
            raise ValueError("Runup and DTM rasters must have the same shape.")

        result_data = np.where(runup_data >= dtm_data, 1.0, np.nan).astype(np.float32)

        with rasterio.open(output_raster_path, 'w', driver='GTiff',
                           height=src_runup.height, width=src_runup.width,
                           count=1, dtype=np.float32, crs=src_runup.crs,
                           transform=src_runup.transform, nodata=np.nan) as dst:
            dst.write(result_data, 1)


def grow_raster(runup_path, max_ingressT, out_path):
    with rasterio.open(runup_path) as src:
        runup = src.read(1)
        profile = src.profile
        transform = src.transform
        res = transform.a  # cell size

    max_ingress = math.ceil(max_ingressT / res) + 2
    print('max_ingress:  ', max_ingress)

    # Identify valid runup data
    mask_valid = np.isfinite(runup) & (runup != 0)

    # Invert mask for distance transform (True = background)
    distance, indices = distance_transform_edt(~mask_valid, return_indices=True)

    # Max ingression in number of pixels
    max_pixels = max_ingress

    # Create a mask where distance is within the radius
    grow_mask = distance <= max_pixels

    # Extrapolate by copying values from nearest valid pixel
    extrapolated = np.zeros_like(runup)
    extrapolated[grow_mask] = runup[tuple(idx[grow_mask] for idx in indices)]

    # Set everything else to nodata
    extrapolated[~grow_mask] = profile.get('nodata', -9999)

    profile.update(dtype='float32', nodata=profile.get('nodata', -9999))

    with rasterio.open(out_path, 'w', **profile) as dst:
        dst.write(extrapolated.astype('float32'), 1)


# Function to extrapolate runup elevation
def extrapolate_runup(runup_raster_path, max_ingress, output_raster_path):
    # Open the raster file
    with rasterio.open(runup_raster_path) as src:
        runup_data = src.read(1)  # Read the first band (assuming single-band raster)
        transform = src.transform  # Get the affine transform
        nodata_value = src.nodata  # Get the nodata value
        res = src.res[0]

        width, height = src.width, src.height

    # Calculate the number of cells for max_ingress
    max_ingress_cells = math.ceil(max_ingress / res) + 2

    # Create a binary mask (1 where values are not nodata)
    print('nodata_value:', nodata_value)
    mask = runup_data != nodata_value

    # Apply binary dilation to simulate the onland ingression
    dilated_mask = binary_dilation(mask, structure=np.ones((max_ingress_cells, max_ingress_cells)))

    # Create a new runup data array with the dilated mask
    extrapolated_runup = np.copy(runup_data)
    extrapolated_runup[dilated_mask] = runup_data.max()  # Assign max value for extrapolated cells

    # Save the extrapolated runup as a new raster
    with rasterio.open(output_raster_path, 'w', driver='GTiff',
                       count=1, dtype=extrapolated_runup.dtype,
                       crs=src.crs, transform=transform,
                       width=width, height=height,
                       nodata=nodata_value) as dst:
        dst.write(extrapolated_runup, 1)

    print("done extrapolate")



def Wave_consequence(raster_paths, people, dtm50, out_folder):
    # out_folder = r'C:\Users\Fiolleau_Sylvain\Documents\Avaframe_simulations\Wave_Outputs'
    # out_folder = r'/lustre02/sandbox/sylvain/'
    runup1 = os.path.join(out_folder, 'temp_runup1.tif')
    #merge_rasters(raster_paths, runup1, method='max', range_vals=None, weights=None)
    runup2 = os.path.join(out_folder, 'temp_runup2.tif')
    print (runup1)
    i = 0
    print(raster_paths)
    for runup1 in raster_paths:
        ##2 contract
        IsData = contract_runup(runup1, runup2, radius_pixels=2)
        if IsData == 0:
            return

        ## 3 adjust
        ## adjust runup
        runup3 = os.path.join(out_folder, 'temp_runup3.tif')

        clip_raster_by_extent(runup2, runup2, runup3, dtm50, max_ingressT = 500)
        ## adjust dtm
        dtm_local =os.path.join(out_folder, 'temp_dtm.tif')
        clip_raster_by_extent(runup2, dtm50, dtm_local, dtm50)

        runup4 = os.path.join(out_folder, 'temp_runup4.tif')

        # add height of lake or sea to the runup height and increase the extent to prepare for the next step
        Add_Height(runup3, dtm_local, runup4)


        runup5 = os.path.join(out_folder, 'temp_runup5.tif')


        max_ingress = 500  # Max ingression distance (in meters)

        # extrapolate_runup(runup4, max_ingress, runup5)
        grow_raster(runup4, max_ingress, runup5)#, metric=0, old=None, new=None)


        runup6 = os.path.join(out_folder, 'temp_runup6.tif')
        OnLand(runup5, dtm_local, runup6)


        out_poly = os.path.join(out_folder, 'temp_runPoly.gpkg')
        raster_to_polygon(runup6, out_poly)
        out_poly_agg = os.path.join(out_folder, 'temp_runPoly_agg.gpkg')
        # Step 1: Aggregate polygons
        with rasterio.open(runup6) as src:
            crs = src.crs
        aggregated_gdf = aggregate_polygons(out_poly, out_poly_agg, crs)

        # Step 2: Count points in polygon (with geometry fix if needed)
        points = gpd.read_file(people)
        if i == 0:
            try:
                result_gdf, index = count_points_in_polygon(aggregated_gdf, points)
            except Exception as e:
                print(f"Error encountered: {e}. Fixing geometries and retrying...")

                # Fix geometries using buffer(0) trick
                aggregated_gdf["geometry"] = aggregated_gdf["geometry"].buffer(0)

                # Retry point counting
                result_gdf, index = count_points_in_polygon(aggregated_gdf, points)

            geometry = GeometryCollection(aggregated_gdf["geometry"].values)  # Convert to GeometryCollection

            numpoint = result_gdf
            previous_index = index
        else:
            try:
                result_gdf, index = count_points_in_polygon(aggregated_gdf, points, previous_index)
            except Exception as e:
                print(f"Error encountered: {e}. Fixing geometries and retrying...")

                # Fix geometries using buffer(0) trick
                aggregated_gdf["geometry"] = aggregated_gdf["geometry"].buffer(0)

                # Retry point counting
                result_gdf, index = count_points_in_polygon(aggregated_gdf, points, previous_index)

            geometry = GeometryCollection(list(geometry.geoms) + list(aggregated_gdf["geometry"].values))


            numpoint = numpoint + result_gdf

        i += 1
    print(type(geometry))
    gdf = gpd.GeoDataFrame(geometry=[geometry], columns=["geometry"])
    gdf["USTABILEFJELLID"] = raster_paths[0].split('\\')[-3]
    gdf["SCENARIOID"] = raster_paths[0].split('\\')[-2]
    gdf["Consequences"] = numpoint

    #
    ## remove temp files
    files_to_remove = [runup2, runup3, runup4, runup5, runup6, dtm_local, out_poly, out_poly_agg]

    for file in files_to_remove:
        if os.path.exists(file):  # Check if file exists before deleting
            os.remove(file)
        else:
            print(f"File not found: {file}")  # Optional: Print missing files
    return gdf
        #        return result_gdf
# # Save final result
# final_output = "path_to_output_folder/final_result.gpkg"
# result_gdf.to_file(final_output, driver="GPKG")
# if __name__==__main__:
#     columns_to_keep = ["USTABILEFJELLID", "SCENARIOID", "Consequences", "geometry"]
#     layer_output = gpd.GeoDataFrame(columns=columns_to_keep, geometry="geometry", crs="EPSG:25833")
#     dtm50_path = r'N:\Prosjekter\391400_Nasjonal_oversikt_over_ustabile_fjellområder\3_Konsekvensanalyse\Waves\2023\dtm_50m.tif'
#     AllScenarios = r'C:\Users\Fiolleau_Sylvain\Documents\Avaframe_simulations\Wave_Outputs_python'
#     Ids = os.listdir(AllScenarios)
#     IdsFold = [item for item in Ids if os.path.isdir(os.path.join(AllScenarios, item))]
#
#     people = r"C:\Users\Fiolleau_Sylvain\Documents\Avaframe_simulations\Bosatte2019_25833.shp"
#     for Id in IdsFold:
#         scenarioIds = os.listdir(os.path.join(AllScenarios,Id))
#         # print(scenarioIds)
#         for scenarioId in scenarioIds:
#             files = '*.tif'
#             raster_paths = glob.glob(os.path.join(AllScenarios,Id, scenarioId, files))
#             Consequences = Wave_consequence(raster_paths, people, dtm50_path)
#             print(Consequences)
#
#             # print('Consequences ', Consequences.columns)
#             if Consequences is not None and not Consequences.empty:
#                 layer_output = pd.concat([layer_output, Consequences[columns_to_keep]], ignore_index=True)
#     layer_output['geometry'] = layer_output['geometry'].buffer(0)
#     layer_output.to_file(r'C:\Users\Fiolleau_Sylvain\Documents\Avaframe_simulations\Wave_Outputs_python\WaveConsequences.gpkg', driver='GPKG')
#


if __name__=='__main__':
    runUp_Gp = r"R:\Arbeidsomr\Temp\Vanja\Piggtind\Wave_Consequences_ArcGis.gpkg"
    aggregated_gdf = gpd.read_file(runUp_Gp)
    people = r"C:\Users\Fiolleau_Sylvain\Documents\Avaframe_simulations\Bosatte2019_25833.shp"
    points = gpd.read_file(people)

    try:
        result_gdf, index = count_points_in_polygon(aggregated_gdf, points)
    except Exception as e:
        print(f"Error encountered: {e}. Fixing geometries and retrying...")

        # Fix geometries using buffer(0) trick
        aggregated_gdf["geometry"] = aggregated_gdf["geometry"].buffer(0)

        # Retry point counting
        result_gdf, index = count_points_in_polygon(aggregated_gdf, points)

    geometry = GeometryCollection(aggregated_gdf["geometry"].values)  # Convert to GeometryCollection

    numpoint = result_gdf
    previous_index = index

    print(type(geometry))
    gdf = gpd.GeoDataFrame(geometry=[geometry], columns=["geometry"])
    gdf["USTABILEFJELLID"] = raster_paths[0].split('\\')[-3]
    gdf["SCENARIOID"] = raster_paths[0].split('\\')[-2]
    gdf["Consequences"] = numpoint


