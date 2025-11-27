# this script prepare the dtms for the on-land runouts. It requires those inputs:
# input sites from the database to be selected in the layers of the project
# slbl folder location where the tickness outputs have the name of the "USTABILEFJELLID_t.tif"
# the dtm mosaic to be clipped

# If it crashes at running a GRASS tool, run r.null tool from the processing toolbox to activate GRASS prior to run this script

# Author: Francois Noel, 2024. Open source CC BY 4.0 (free to use, modify, sell if citing original author) modified by S. Fiolleau
import numpy as np
import rasterio
import fiona
import geopandas as gpd
from shapely.geometry import Polygon, box, Point
from rasterio.enums import Resampling
from rasterio.transform import Affine, from_bounds
from rasterio.mask import mask
import sys
import os

from osgeo import gdal

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import rasterstats as rs
import FuncRunAva as FRA
import tempfile
import math
import pandas as pd


# Define the volume classification
def classify_by_volume(row, mu_voell):
    vol = row["pre-vol"]
    res = 5
    # voellmy = 'VoellmyMinShear'
    voellmy = 'voellmy'

    if vol < (0.250 * 10 ** 6):
        return res, 'coulomb', np.tan(np.radians(31)), 2000
    elif vol < (0.750 * 10 ** 6):
        return res, voellmy, mu_voell, int(10000 * (5 / res))
    elif vol < (1 * 10 ** 6):
        return res, voellmy, mu_voell, int(20000 * (5 / res))
    elif vol < (10 * 10 ** 6):
        return res, voellmy, mu_voell, int(10000 * (5 / res))
    elif vol < (50 * 10 ** 6):
        return res, voellmy, mu_voell, int(50000 * (5 / res))
    elif vol < (100 * 10 ** 6):
        return res, voellmy, mu_voell, int(100000 * (5 / res))
    else:
        return res, voellmy, mu_voell, int(300000 * (5 / res))


def resample_raster(input_raster_path, output_raster_path, out_res, method="bilinear"):
    # Open the input raster
    src = gdal.Open(input_raster_path)

    # Define resampling method (nearest, bilinear, cubic, lanczos)
    resampling_methods = {
        'nearest': gdal.GRA_NearestNeighbour,
        'bilinear': gdal.GRA_Bilinear,
        'cubic': gdal.GRA_Cubic,
        'lanczos': gdal.GRA_Lanczos
    }

    # Resample to desired output resolution
    options = gdal.WarpOptions(
        xRes=out_res, yRes=out_res, resampleAlg=resampling_methods[method]
    )

    gdal.Warp(output_raster_path, src, options=options)

    # Close the raster source
    src = None


def fill_nodata(input_raster_path, output_raster_path, fill_value=0):
    # Open the raster
    with rasterio.open(input_raster_path) as src:
        # Read the raster data
        raster_data = src.read(1)

        # Get metadata from the raster
        meta = src.meta.copy()

        # Get NoData value
        nodata_value = src.nodata

        # Replace NoData values with the fill value (e.g., 0)
        filled_data = np.where(raster_data == nodata_value, fill_value, raster_data)

        # Update metadata
        meta.update({'nodata': fill_value})

        # Write the filled raster to a new file
        with rasterio.open(output_raster_path, 'w', **meta) as dst:
            dst.write(filled_data, 1)


def raster_to_ascii(input_raster_path, output_ascii_path, nodata_value=-9999):
    # Open the filled raster with GDAL
    src_raster = gdal.Open(input_raster_path)

    # Translate the raster to ESRI ASCII grid
    gdal.Translate(
        output_ascii_path,
        src_raster,
        format='AAIGrid',
        noData=nodata_value
    )


def save_polygon_to_shp(input_polygon, output_shp_path, crs, existing_properties=None):
    # Check if input_polygon has geometry (this assumes input_polygon is a row from a GeoDataFrame)
    if not hasattr(input_polygon, 'geometry') or input_polygon.geometry is None:
        raise ValueError("Input polygon is not a valid geometry.")

    # Dynamically create schema by inspecting data types of each property in the input_polygon
    properties_types = {key: type(input_polygon[key]).__name__ for key in input_polygon.index if key != 'geometry'}

    # Map Python types to Fiona types
    type_mapping = {
        'int32': 'int32',
        'float64': 'float64',
        'str': 'str',
        'bool': 'bool',
        'NoneType': 'str',
        'Timestamp': 'str',
        'NaTType': 'str',
        'Point': 'str',
    }

    # Create schema with correct types for each property
    schema = {'geometry': 'Polygon',
              'properties': {key: type_mapping.get(properties_types[key], 'str') for key in properties_types}}

    # Construct the properties dictionary (excluding the geometry field)
    properties = {key: str(input_polygon[key]) for key in input_polygon.index if key != 'geometry'}

    # Debugging: Check the properties being written and schema
    print(f"Schema being used: {schema}")
    print(f"Properties being written: {properties}")

    # Open the shapefile for writing
    with fiona.open(output_shp_path, 'w', driver='ESRI Shapefile', crs=crs, schema=schema) as dst:
        # Use Shapely's mapping function to convert geometry to dict format
        # Fiona expects geometry to be in dict format (GeoJSON-like)
        dst.write({
            'geometry': input_polygon.geometry,  # Convert geometry to dict using Shapely's mapping function
            'properties': properties,  # The properties are now the updated ones
        })


def Prep_Ava_Simul(sources_areas, dtm_path, ava_folder, out_folder, slbl_folder, IDs, SCENARIOID, autoparam=True,
                   ava_time_step=0.1, fric_model='Voellmy', ava_mu_user=0.06, ava_t_coef=400, nb_part=10000, out_res=5,
                   ava_density=2700, run_simul=0):
    Skiped_Ids = {}

    oversize_area = 3  # how many times the dtm area should expand beyound the one covered by the cone method

    # time step to save if want multiple frames # saving time step, i.e.  time in seconds (first and last time step are always saved)
    # option 1: give an interval with start:interval in seconds (tStep = 0:5 - this will save desired results every 5 seconds for the full simulation)
    # option 2: explicitly list all desired time steps (closest to actual computational time step) separated by | (example tSteps = 1|50.2|100)
    tSteps = 1
    sphOption = 3
    sphOptionIni = 2  # initially was set to 2

    # print('Outres: ', out_res)
    # MPPDH= mass per particles through release thickness. massPerPart = rho x cellarea x deltaTh
    # MPPDIR= mass per particle direct. The massPerPart value is taken from the configuration and is the same for all cells.
    # MPPKR= mass per particles through number of particles per kernel radius.
    massPerParticleDeterminationMethod = 'MPPDIR'

    # is computed with: nPPK = nPPK0 * (sphKR/sphKR0)^aPPK
    # where sphKR is the sphKernelRadius specified further up
    # reference kernel radius [m]
    sphKR0 = out_res
    # reference number of particles per kernel radius
    nPPK0 = 35  # 15

    sources_areas['pre-vol'] = np.nan

    #################### loop over all sites #################
    for k in range(sources_areas.shape[0]):
        sources_area = sources_areas.iloc[k]
        print(sources_area)
        # sources_area_geometry = sources_areas.loc[k, 'geometry']  # Replace 0 with the desired index
        sources_area_geometry = sources_areas.iloc[k]['geometry']
        ust_id = sources_area['USTABILEFJELLID']
        scen_id = sources_area['SCENARIOID']
        # SCENARIOID = ['1800102', '1800103', '1800107', '1800108']
        # skips if not a USTABILEFJELLID
        if (math.isnan(ust_id) or math.isnan(scen_id)):
            Skiped_Ids[str(ust_id)] = ['No UstabilID or scenarioID']
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

        path_slbl = slbl_folder + "/slbl_" + str(int(ust_id)) + "_" + str(int(scen_id)) + "_ortho.tif"

        if not os.path.isfile(path_slbl):
            Skiped_Ids[str(int(ust_id))] = [str(int(scen_id)), 'no sblb file, wrong path']
            continue

        USTABILEFJELLID = str(int(ust_id))
        print("Processing site: " + str(scen_id))
        # Add slbl layer to project

        # Open the DTM raster
        with rasterio.open(dtm_path) as src_dtm:
            # dtm_array = src_dtm.read(1)  # Read the first band of the raster
            dtm_transform = src_dtm.transform
            dtm_xmin, dtm_ymin, dtm_xmax, dtm_ymax = src_dtm.bounds
            dtm_crs = src_dtm.crs
            dtm_res = src_dtm.res
        with rasterio.open(path_slbl) as src_slbl:
            slbl_res = src_slbl.res
        zonal_stats_result_dtm = rs.zonal_stats(sources_area_geometry, dtm_path, stats=["max"])
        zonal_stats_result_slbl = rs.zonal_stats(sources_area_geometry, path_slbl, stats=["sum"])

        sources_areas.loc[k, "pre-max"] = zonal_stats_result_dtm[0]["max"]
        sources_areas.loc[k, "pre-vol"] = zonal_stats_result_slbl[0]["sum"] * abs(slbl_res[0]) * abs(slbl_res[1])
        sources_areas.loc[k, "centroid"] = sources_area_geometry.centroid
        try:
            if autoparam == True:
                sources_areas.loc[k, ['out_res', 'ava_model', 'ava_mu', 'ava_nb_part']] = classify_by_volume(
                    sources_areas.iloc[k], ava_mu_user)
            else:
                sources_areas.loc[
                    k, ['out_res', 'ava_model', 'ava_mu', 'ava_nb_part']] = out_res, fric_model, ava_mu_user, nb_part
        except:
            continue
        out_res = sources_areas.loc[k, 'out_res']
        ava_model = sources_areas.loc[k, 'ava_model']
        ava_mu = sources_areas.loc[k, 'ava_mu']
        ava_nb_part = sources_areas.loc[k, 'ava_nb_part']
        ## make folders if doesn't exist

        layer_ustabilefjellid = str(int(ust_id))
        layer_scenarioid = str(int(scen_id))

        # Create directories
        base_dir = os.path.join(out_folder, str(layer_ustabilefjellid), str(layer_scenarioid))
        directories = ['Inputs', 'Outputs', 'Work']
        subdirectories = ['ENT', 'LINES', 'POINTS', 'REL', 'RELTH', 'RES', 'SECREL']

        for dir_name in directories:
            dir_path = os.path.join(base_dir, dir_name)
            os.makedirs(dir_path, exist_ok=True)
            if dir_name == 'Inputs':
                for sub_dir in subdirectories:
                    os.makedirs(os.path.join(dir_path, sub_dir), exist_ok=True)

        # Assuming 'layer' is a GeoDataFrame with the geometries of features
        bounding_box = sources_area_geometry.bounds  # (minx, miny, maxx, maxy)

        # Calculate width and height
        width = bounding_box[2] - bounding_box[0]
        height = bounding_box[3] - bounding_box[1]
        source_size = [width, height]

        #################### produce clip rectangle of site #################
        V = sources_areas.loc[k, "pre-vol"]
        phi = min(31, np.degrees(np.arctan(10 ** (0.623419 - 0.15666 * np.log10(V)))))  # truncated Scheidegger
        H = sources_areas.loc[k, 'pre-max']  # assuming H max being the max elevation of the source
        # print(sources_areas.loc[k, "pre-vol"], sources_areas.loc[k, "pre-max"])
        L = H / (np.tan(np.radians(phi)))  # runout length

        ## calculate the dtm surface in square km if smaller than 20 clipping box is the same size of dtm
        # Extract the coordinates of the extent

        dtm_surface = (dtm_xmax - dtm_xmin) * (dtm_ymax - dtm_ymin) / 10 ** 6
        if dtm_surface > 20:

            box_width = oversize_area * (H / np.tan(np.radians(phi)) + source_size[
                0])  # go slighly oversized to limit out of bounds simulations
            box_height = oversize_area * (H / np.tan(np.radians(phi)) + source_size[
                1])  # go slighly oversized to limit out of bounds simulations

            # Create a rectangular clipping box around the centroid
            bbox_dtm = (
            sources_area_geometry.centroid.x - box_width / 2, sources_area_geometry.centroid.y - box_height / 2,
            sources_area_geometry.centroid.x + box_width / 2, sources_area_geometry.centroid.y + box_height / 2)

            with rasterio.open(dtm_path) as src_dtm:
                # Clip the raster using the mask
                window_dtm = src_dtm.window(*bbox_dtm)
                clip_dtm_transform = from_bounds(*bbox_dtm, window_dtm.width, window_dtm.height)
                clip_dtm = src_dtm.read(1, window=window_dtm)

                clip_slbl = clip_dtm * 0
                clip_slbl_transform = clip_dtm_transform
                height_clip_slbl, width_clip_slbl = clip_slbl.shape
                height_clip_dtm, width_clip_dtm = clip_dtm.shape
                clip_dtm_meta = src_dtm.meta.copy()
                clip_dtm_meta.update({
                    "driver": "GTiff",
                    "count": 1,
                    "dtype": "float32",
                    "crs": src_dtm.crs,
                    "transform": clip_dtm_transform,
                    "width": width_clip_dtm,  # Update the width based on the clipped raster
                    "height": height_clip_dtm  # Update the height based on the clipped raster
                })

        else:

            dtm_center_x = (dtm_xmin + dtm_xmax) / 2
            dtm_center_y = (dtm_ymin + dtm_ymax) / 2
            box_width = dtm_xmax - dtm_xmin
            box_height = dtm_ymax - dtm_ymin
            # Define the box around the center point
            bbox_dtm = (dtm_center_x - box_width / 2, dtm_center_y - box_height / 2,
                        dtm_center_x + box_width / 2, dtm_center_y + box_height / 2)

            with rasterio.open(dtm_path) as src_dtm:
                # Clip the raster using the mask
                window_dtm = src_dtm.window(*bbox_dtm)
                clip_dtm_transform = from_bounds(*bbox_dtm, window_dtm.width, window_dtm.height)
                clip_dtm = src_dtm.read(1, window=window_dtm)

                clip_slbl = clip_dtm * 0
                clip_slbl_transform = clip_dtm_transform
                height_clip_slbl, width_clip_slbl = clip_slbl.shape

                height_clip_dtm, width_clip_dtm = clip_dtm.shape
                clip_dtm_meta = src_dtm.meta.copy()
                clip_dtm_meta.update({
                    "driver": "GTiff",
                    "count": 1,
                    "dtype": "float32",
                    "crs": src_dtm.crs,
                    "transform": clip_dtm_transform,
                    "width": width_clip_dtm,  # Update the width based on the clipped raster
                    "height": height_clip_dtm  # Update the height based on the clipped raster
                })

        with rasterio.open(path_slbl) as src_slbl:

            # Clip the raster using the mask
            window_slbl = src_slbl.window(*bounding_box)
            clip_slbl1_transform = from_bounds(*bounding_box, window_slbl.width, window_slbl.height)
            clip_slbl1 = src_slbl.read(1, window=window_slbl)

            # Update the metadata after the first clipping
            clip_slbl1_meta = src_slbl.meta.copy()
            clip_slbl1_meta.update({
                "driver": "GTiff",
                "count": 1,
                "dtype": "float32",
                "crs": src_slbl.crs,
                "transform": clip_slbl1_transform
            })
        # Define the offset where slbl1 should be placed in clip_slbl
        x_offset = int((clip_slbl1_transform[2] - clip_slbl_transform[2]) / clip_slbl_transform[
            0])  # Example offset for the y-direction
        y_offset = int((clip_slbl_transform[5] - clip_slbl1_transform[5]) / clip_slbl_transform[
            0])  # Example offset for the x-direction

        clip_slbl[y_offset:y_offset + clip_slbl1.shape[0], x_offset:x_offset + clip_slbl1.shape[1]] = clip_slbl1
        # Update the metadata for the second clipping
        clip_slbl_meta = clip_slbl1_meta.copy()
        clip_slbl_meta.update({
            "transform": clip_slbl_transform,
            'height': clip_dtm.shape[0],
            'width': clip_dtm.shape[1],
            'count': 1  # Assuming it's a single band raster
        })

        # Calculate the extent of the clipped raster using the transform
        xmin_clip_slbl, ymax_clip_slbl = clip_slbl_transform * (0, 0)  # Top-left corner (0,0) in pixel coordinates
        xmax_clip_slbl, ymin_clip_slbl = clip_slbl_transform * (width_clip_slbl, height_clip_slbl)

        # Calculate the extent of the clipped raster using the transform
        xmin_clip_dtm, ymax_clip_dtm = clip_dtm_transform * (0, 0)  # Top-left corner (0,0) in pixel coordinates
        xmax_clip_dtm, ymin_clip_dtm = clip_dtm_transform * (width_clip_dtm, height_clip_dtm)

        # Create extents as tuples
        extent_slbl = (xmin_clip_slbl, ymin_clip_slbl, xmax_clip_slbl, ymax_clip_slbl)
        extent_dtm = (xmin_clip_dtm, ymin_clip_dtm, xmax_clip_dtm, ymax_clip_dtm)

        if extent_dtm != extent_slbl:
            print('clipping extents do not match slbl and dtm')

        clip_dtm2 = clip_dtm - clip_slbl
        with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as temp_clip_dtm:
            # Open the temporary file to write clip_dtm
            with rasterio.open(temp_clip_dtm.name, 'w', **clip_dtm_meta) as dst:
                dst.write(clip_dtm2, 1)

        # Create temporary files
        with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as temp_clip_slbl:

            # Open the temporary file to write clip_slbl
            with rasterio.open(temp_clip_slbl.name, 'w', **clip_slbl_meta) as dst:
                dst.write(clip_slbl, 1)

        ## resample to out_res and store in ascii format the dtm
        temp_resample_dtm = out_folder + 'resampled_clip_dtm.tif'
        resample_raster(temp_clip_dtm.name, temp_resample_dtm, out_res, 'bilinear')
        #
        temp_resample_filled_dtm = out_folder + 'filled_clip_dtm.tif'
        fill_nodata(temp_resample_dtm, temp_resample_filled_dtm, fill_value=0)

        output_ascii_path_dtm = os.path.join(base_dir, 'Inputs', 'clip_dtm.asc')
        raster_to_ascii(temp_resample_filled_dtm, output_ascii_path_dtm, nodata_value=-9999)

        #################### produce source thickness of site #################
        temp_resample_slbl = out_folder + 'resampled_clip_slbl.tif'
        resample_raster(temp_clip_slbl.name, temp_resample_slbl, out_res=out_res, method='bilinear')

        temp_resample_filled_slbl = out_folder + 'filled_clip_slbl.tif'
        fill_nodata(temp_resample_slbl, temp_resample_filled_slbl, fill_value=0)

        output_ascii_path_slbl = os.path.join(base_dir, 'Inputs', 'RELTH', 'thickness.asc')
        raster_to_ascii(temp_resample_filled_slbl, output_ascii_path_slbl, nodata_value=-9999)

        #################### produce source polygon of site #################

        output_layer_path = os.path.join(base_dir, 'Inputs', 'REL', 'ava_source.shp')
        crs = sources_areas.crs
        # Assuming you have the existing properties (maybe from an existing shapefile)
        existing_properties = {key: value for key, value in sources_areas.iloc[k].items()}
        save_polygon_to_shp(sources_areas.iloc[k], output_layer_path, crs, existing_properties)

        #################### produce AvaFrame parameter file #################

        # opens reference parameter file to modify
        cellsurface = out_res ** 2  # cellsize^2
        with open(os.path.join(ava_folder, 'com1DFA', 'com1DFACfg.ini'), 'r') as fp:
            data = fp.readlines()  # read a list of lines into data
            fp.seek(0, 0)  # goes back to the start of the file
            for l_no, line in enumerate(fp):
                if '#' not in line[0]:  # skip line if comment
                    if 'resType =' in line:  # search string
                        data[l_no] = 'resType = pft|pfv|FT\n'  # set to the desired outputs
                    elif 'tSteps =' in line:  # time steps saved
                        data[l_no] = 'tSteps = ' + str(
                            tSteps) + '\n'  # set the time step to save # option 1: give an interval with start:interval in seconds (tStep = 0:5 - this will save desired results every 5 seconds for the full simulation)
                    # option 2: explicitly list all desired time steps (closest to actual computational time step) separated by | (example tSteps = 1|50.2|100)
                    elif 'sphOptionIni =' in line:  # time steps saved
                        data[l_no] = 'sphOptionIni = ' + str(sphOptionIni) + '\n'
                    elif 'rho =' in line:  # search string
                        data[l_no] = 'rho = ' + str(ava_density) + '\n'  # set to the desired density
                    elif 'relThFromShp =' in line:  # search string
                        data[l_no] = 'relThFromShp = False\n'  # set to read the thickness from a raster
                    elif 'relThFromFile =' in line:  # search string
                        data[l_no] = 'relThFromFile = True\n'  # set to read the thickness from a raster
                    elif 'dt =' in line:  # search string
                        data[l_no] = 'dt = ' + str(
                            ava_time_step) + '\n'  # time step for the SPH (more particles requires smaller steps)
                    elif 'sphOption =' in line:  # search string
                        data[l_no] = 'sphOption = ' + str(sphOption) + '\n'  # set the SPH merthod
                    elif 'massPerParticleDeterminationMethod =' in line:  # search string
                        data[
                            l_no] = 'massPerParticleDeterminationMethod = ' + massPerParticleDeterminationMethod + '\n'  # set the SPH merthod
                    elif 'massPerPart =' in line:  # search string
                        data[l_no] = 'massPerPart = ' + str(ava_density * cellsurface * V / (
                                    cellsurface * ava_nb_part)) + '\n'  # set the corresponding mass per particle
                    elif 'deltaTh =' in line:  # search string
                        data[l_no] = 'deltaTh = ' + str(V / (
                                    cellsurface * ava_nb_part)) + '\n'  # set to the desired turbulence coefficient for the Voellmy model
                    elif 'sphKR0 =' in line:  # search string
                        data[l_no] = 'sphKR0 =' + str(sphKR0) + '\n'  # reference kernel radius [m]
                    elif 'nPPK0 =' in line:  # search string
                        data[l_no] = 'nPPK0 =' + str(nPPK0) + '\n'  # reference number of particles per kernel radius
                    elif 'splitOption =' in line:  # search string
                        data[l_no] = 'splitOption = ' + str(0) + '\n'  # set to adjust the nb of particle based on mass
                        # data[l_no] = 'splitOption = ' + str(1) + '\n'#set to adjust the nb of particle for a proper number given the kernel size (1=constant nb of particle per kernel)
                    elif 'meshCellSize =' in line:  # search string
                        data[l_no] = 'meshCellSize = ' + str(out_res) + '\n'  # adjust the size of the mesh
                    elif 'frictModel =' in line:  # search string
                        data[
                            l_no] = 'frictModel = ' + ava_model + '\n'  # set to the desired friction type #model to use (samosAT, samosATSmall, samosATMedium, samosATAuto, Coulomb, Voellmy, VoellmyMinShear, CoulombMinShear, wetsnow)
                    elif 'muvoellmy =' in line:  # search string
                        data[l_no] = 'muvoellmy = ' + str(
                            ava_mu) + '\n'  # set to the desired friction (5.7 deg gives a mu of 0.1)
                    elif 'xsivoellmy =' in line:  # search string
                        data[l_no] = 'xsivoellmy = ' + str(
                            ava_t_coef) + '\n'  # set to the desired turbulence coefficient for the Voellmy model
                    # elif 'muvoellmyminshear =' in line:  # search string
                    #     data[l_no] = 'muvoellmyminshear = ' + str(ava_mu) + '\n'  # set to the desired friction (5.7 deg gives a mu of 0.1)
                    # elif 'xsivoellmyminshear =' in line:  # search string
                    #     data[l_no] = 'xsivoellmyminshear = ' + str(ava_t_coef) + '\n'  # set to the desired turbulence coefficient for the Voellmy model
                    # elif 'tau0voellmyminshear =' in line:  # search string
                    #     data[l_no] = 'tau0voellmyminshear = ' + str(tau0voellmyminshear) + '\n'  # set to the desired turbulence coefficient for the Voellmy model
                    elif 'mucoulomb =' in line:  # search string
                        data[l_no] = 'mucoulomb = ' + str(
                            ava_mu) + '\n'  # set to the desired friction (5.7 deg gives a mu of 0.1)
                    elif 'dam =' in line:  # search string
                        data[l_no] = 'dam = False\n'  # do not use a dam input
                    elif 'plotFields =' in line:  # search string
                        data[l_no] = 'plotFields = pft|pfv|FT\n'  # do not use a dam input


        # write new modified parameter file to the site folder
        with open(os.path.join(base_dir, "local_com1DFACfg.ini"), 'w') as fp:
            fp.writelines(data)

        #################### produce AvaFrame settings file #################

        # opens reference parameter file to modify
        with open(os.path.join(ava_folder, 'avaframeCfg.ini'), 'r') as fp:
            data = fp.readlines()  # read a list of lines into data
            fp.seek(0, 0)  # goes back to the start of the file
            for l_no, line in enumerate(fp):
                if '#' not in line[0]:  # skip line if comment
                    if 'avalancheDir =' in line:  # search string
                        data[l_no] = 'avalancheDir = ' + base_dir + '\n'  # set to the desired outputs

        # write new modified parameter file to the site folder
        with open(os.path.join(base_dir, 'local_avaframeCfg.ini'), 'w') as fp:
            fp.writelines(data)
        try:
            os.remove(temp_resample_dtm)
            os.remove(temp_resample_filled_dtm)
            os.remove(temp_resample_slbl)
            os.remove(temp_resample_filled_slbl)
            os.remove(temp_clip_dtm.name)
            os.remove(temp_clip_slbl.name)
            os.remove(os.path.join(out_folder, 'temp_clip_slbl3.tif'))
            os.remove(os.path.join(out_folder, 'temp_clip_slbl4.tif'))
            os.remove(os.path.join(out_folder, 'temp_clip_slbl5.tif'))
            os.remove(os.path.join(out_folder, 'temp_clip_slbl6.tif'))
        except:
            print('did not manage to delete temp files')

        ## run simulation
        if run_simul == 1:
            FRA.run_avaframe(base_dir, ava_folder)
        # remove temp tif files (due to issues saving temp files on macos)

        print('Done generating the files for AvaFrame\n')


