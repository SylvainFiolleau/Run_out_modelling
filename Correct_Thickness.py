# This python script converts depths (DoDs) into thickness values (along normals)

import numpy as np
import os
from osgeo import gdal
import re


def process_files(folder_path):
    for filename in os.listdir(folder_path):
        match = re.match(pattern, filename)
        if match:
            stars = match.group(1)
            hashes = match.group(2)
            yield stars, hashes
    # Dictionary to store **** and corresponding ####
    file_info = {}

    # Loop over all files in the folder
    for filename in os.listdir(folder_path):
        match = re.match(pattern, filename)
        if match:
            # Extract **** and ####
            stars = match.group(1)
            hashes = match.group(2)
            # Store in the dictionary
            file_info[stars] = hashes
    In_folder = folder_path

def resample(src_ds, Output, res):

    resampling_options = gdal.WarpOptions(
        xRes=res,
        yRes=res,
        resampleAlg='bilinear'  # Options: 'nearest', 'bilinear', 'cubic', 'cubicspline', etc.
    )
    gdal.Warp(Output, src_ds, options=resampling_options)

def cut_raster(RasterToCut, RefRaster, output):
    ref_geotrans = RefRaster.GetGeoTransform()
    xmin = ref_geotrans[0]
    ymax = ref_geotrans[3]
    xmax = xmin + ref_geotrans[1] * RefRaster.RasterXSize
    ymin = ymax + ref_geotrans[5] * RefRaster.RasterYSize
    # Define the gdal.Warp options with the extent from the reference raster
    warp_options = gdal.WarpOptions(
        outputBounds=(xmin, ymin, xmax, ymax),  # Set clipping extent
        resampleAlg='bilinear',                  # Resampling method
        cropToCutline=True                      # Ensures it clips exactly to the boundary
    )
    gdal.Warp(output, RasterToCut, options=warp_options)
    print('HI')
def write_Output_Raster(raster_data, output_asc, geo_transform, projection):
    # Specify paths for output files
    temp_tiff = 'temp_output.tif'

    # Create a new raster (GeoTIFF) as an intermediate step
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(temp_tiff, raster_data.shape[1], raster_data.shape[0], 1,
                                   gdal.GDT_Float32)

    if output_dataset is None:
        raise Exception("Failed to create the temporary GeoTIFF dataset.")

    # Set geo-transform and projection
    output_dataset.SetGeoTransform(geo_transform)
    output_dataset.SetProjection(projection)

    # Write the array to the temporary GeoTIFF
    output_band = output_dataset.GetRasterBand(1)
    output_band.WriteArray(raster_data)
    output_band.FlushCache()
    output_dataset = None

    # Now translate the GeoTIFF to .asc format using gdal.Translate()
    gdal.Translate(output_asc, temp_tiff, format="AAIGrid")
    ds = gdal.Open(output_asc)
    ds.Close()
    ds = None
    # output_asc.GetRasterBand(1)

    print(f"Raster data successfully written to {output_asc}")

    output_asc = None
    # Optional: Delete the temporary file if no longer needed
    if os.path.exists(temp_tiff):
        os.remove(temp_tiff)
        print(f"Temporary file {temp_tiff} deleted.")


def write_Output_Raster_tif(raster_data, output_tiff, geo_transform, projection):
    # Create a new raster (GeoTIFF)
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_tiff, raster_data.shape[1], raster_data.shape[0], 1,
                                   gdal.GDT_Float32)

    if output_dataset is None:
        raise Exception("Failed to create the GeoTIFF dataset.")

    # Set geo-transform and projection
    output_dataset.SetGeoTransform(geo_transform)
    output_dataset.SetProjection(projection)

    # Write the array to the GeoTIFF
    output_band = output_dataset.GetRasterBand(1)
    output_band.WriteArray(raster_data)
    output_band.FlushCache()
    output_dataset = None

    print(f"Raster data successfully written to {output_tiff}")

def Calc_Norm(source_dtm, source_dod, cellsize):
    ###### Normals calculation ######
    # calculates normals of the dtm based on the neighboors (E-W and N-S cells: b-d-f-h, ignoring the D8 corners: a-c-g-i)
    # inputs: dtm grid and its header information
    #    ___ ___ ___
    #   | a | b | c |
    #   |___|___|___|
    #   | d | e | f |
    #   |___|___|___|
    #   | g | h | i |
    #   |___|___|___|
    b = np.array(source_dtm)
    d = np.array(source_dtm)
    f = np.array(source_dtm)
    h = np.array(source_dtm)
    h = np.array(source_dtm)
    b[1:, :] = b[0:-1, :]  # shift cells down
    d[:, 1:] = d[:, 0:-1]  # shift cells right
    f[:, 0:-1] = f[:, 1:]  # shift cells left
    h[0:-1, :] = h[1:, :]  # shift cells up

    # calculate eigenvectors (vi = null(A-lambda_i*Identity_matrix))
    dx = b - h
    dy = f - d
    sz = (b - source_dtm) ** 2 + (d - source_dtm) ** 2 + (source_dtm - f) ** 2 + (source_dtm - h) ** 2
    # clearvars b d f h

    # v2 = N':
    ny = -(2 * dx) / ((4 * dx ** 2 + 4 * dy ** 2 + sz ** 2 - 4 * sz + 4) ** 0.5 - sz + 2)
    nx = -(2 * dy) / ((4 * dx ** 2 + 4 * dy ** 2 + sz ** 2 - 4 * sz + 4) ** 0.5 - sz + 2)
    # clearvars dx dy sz
    nz = np.ones((np.size(nx, 0), np.size(nx, 1)))

    # scale based on cellsize
    h_ratio = 1 / cellsize
    nx *= h_ratio
    ny *= h_ratio

    # normalise
    norm_n = (nx ** 2 + ny ** 2 + nz ** 2) ** 0.5
    nx /= norm_n
    ny /= norm_n
    nz /= norm_n
    # clearvars norm_n

    # flip if not pointing up
    nx[nz < 0] *= -1
    ny[nz < 0] *= -1
    nz[nz < 0] *= -1

    orthogonal_thickness = source_dod / norm_n

    # end of normal function with outputs: nx, ny and nz
    return orthogonal_thickness
