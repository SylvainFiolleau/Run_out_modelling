import os
import numpy as np
from numpy.ma.extras import hstack
from osgeo import gdal
import rasterio
from rasterio.transform import from_origin
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
import sys
import Func_Waves_Angles as CorAngle

# Global vars for multiprocessing
global_elev_model = None
global_is_nodata = None
global_OpenSeaData = None


# Function to load the digital terrain model (DTM)
def dtm_import():
    # Load the water raster
    # Open file dialog to select file
    root = tk.Tk()
    root.withdraw()  # Hide Tkinter root window
    filepath = filedialog.askopenfilename(
        filetypes=[("Raster Files", "*.tif *.tiff *.asc *.txt"), ("All Files", "*.*")])
    with rasterio.open(filepath) as dataset:
        elev_model = {
            'dtm': dataset.read(1),
            'nrows': dataset.height,
            'ncols': dataset.width,
            'cellsize': dataset.res[0],  # Assuming square pixels
            'xllcorner': dataset.bounds.left,
            'yllcorner': dataset.bounds.bottom,
        }
    return elev_model


# Define functions for handling GeoTIFF writing
def geotiffwrite(filepath, data, dtype, cellsize, xllcorner, yurcorner):
    transform = from_origin(xllcorner, yurcorner, cellsize, cellsize)
    with rasterio.open(filepath, 'w', driver='GTiff', height=data.shape[0], width=data.shape[1], count=1, dtype=dtype,
                       transform=transform) as dst:
        dst.write(data, 1)


def gdal_write_geotiff(filename, data, cellsize, xll, yur):
    driver = gdal.GetDriverByName('GTiff')
    rows, cols = data.shape
    out_tif = driver.Create(filename, cols, rows, 1, gdal.GDT_Float32)
    out_tif.SetGeoTransform((xll, cellsize, 0, yur, 0, -cellsize))
    out_tif.GetRasterBand(1).WriteArray(data)
    out_tif.FlushCache()


def init_worker(elev, nodata):
    """Initializer to pass global data to worker processes."""
    global global_elev_model, global_is_nodata
    global_elev_model = elev
    global_is_nodata = nodata


def process_site(k, sites, result_dir, max_vol, min_wave, max_wave, IDs, ScenarioIDs, lim_azimuth=0, save_vid=False):
    global global_elev_model, global_is_nodata

    # use global_elev_model instead of passing as arg
    elev_model = global_elev_model
    is_nodata = global_is_nodata

    # print(f'k{k + 1}: site id {sites.Id[k]} of {sites.shape[0]}')
    # Create directories for the sites
    site_dir = f'{result_dir}/{sites.Id[k]}/{sites.ScenarioID[k]}'

    # Check if the result exists
    if os.path.isfile(f'{site_dir}/runup_id{int(sites.WaveID[k])}.tif'):
        return
    else:
        print(f'k{k + 1}: site id {sites.Id[k]} of {sites.shape[0]}')
    os.makedirs(site_dir, exist_ok=True)
    # Initialize distance, azimuth, and other arrays

    distance = np.ones_like(elev_model['dtm'], dtype=np.float32) * 9999999

    azimuth = np.zeros_like(distance)
    # Assuming that elev_model and sites are already defined as arrays
    p_row = elev_model['nrows'] - np.floor((sites.Y[k] - elev_model['yllcorner']) / elev_model['cellsize']).astype(int)
    p_col = np.floor((sites.X[k] - elev_model['xllcorner']) / elev_model['cellsize']).astype(int) + 1

    particles = np.array([[p_row, p_col, 0, 0]])
    w_ext_subs = np.array([[p_row, p_col], [p_row, p_col]])
    distance[p_row, p_col] = 0
    elev_source = elev_model['dtm'][p_row, p_col]
    if save_vid:
        import matplotlib.pyplot as plt
        import cv2
        vid_file = f'{site_dir}/wave_id{int(sites.WaveID[k])}.mp4'
        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        v = cv2.VideoWriter(vid_file, fourcc, 30, (640, 480))
    # While loop to propagate
    is_change = True
    nbvl = elev_model['ncols'] * elev_model['nrows']
    subs_24 = np.array([
        [-3, -3, -3, -3, -3, -3, -3, -2, -2, -1, -1, 0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, -2, -2, -2, -2, -2, -1,
         -1, 0, 0, 1, 1, 2, 2, 2, 2, 2, -1, -1, -1, 0, 0, 1, 1, 1],
        [-3, -2, -1, 0, 1, 2, 3, -3, 3, -3, 3, -3, 3, -3, 3, -3, 3, -3, -2, -1, 0, 1, 2, 3, -2, -1, 0, 1, 2, -2, 2,
         -2, 2, -2, 2, -2, -1, 0, 1, 2, -1, 0,
         1, -1, 1, -1, 0, 1]
    ])
    subs_8 = subs_24[:, 40:48]

    dist_opt = np.linalg.norm((subs_24 * elev_model['cellsize']).T, axis=1)
    m = 1
    while is_change:
        particles = np.repeat(particles, subs_8.shape[1], axis=0) + np.tile(
            np.hstack((subs_8.T, np.zeros([subs_8.shape[1], 2]))), (particles.shape[0], 1))

        ## check if particles in water
        is_inbound = (particles[:, 0] >= 0) & (particles[:, 0] < distance.shape[0]) & (particles[:, 1] >= 0) & (
                    particles[:, 1] < distance.shape[1])
        particles = particles[is_inbound]
        particles = particles[
            np.logical_and(~is_nodata[(particles[:, 0] - 1).astype(int), (particles[:, 1] - 1).astype(int)],
                           particles[:, 2] <= distance[
                               (particles[:, 0] - 1).astype(int), (particles[:, 1] - 1).astype(int)])]
        # Sort and update particles
        #######################################
        repmat_part = np.tile(np.hstack([particles[:, :3] * np.array([1, 1, 0]),
                                         np.arange(1, particles.shape[0] + 1).reshape(-1, 1)])
                              , (1, subs_24.shape[1]))
        reshape_repmat = repmat_part.reshape([particles.shape[0] * subs_24.shape[1], particles.shape[1]])
        repmat_part2 = np.tile(np.hstack([subs_24.T, dist_opt.reshape(-1, 1), np.zeros([subs_24.shape[1], 1])]),
                               (particles.shape[0], 1))
        particles_opt = reshape_repmat + repmat_part2
        is_inbound = (particles_opt[:, 0] >= 1) & (particles_opt[:, 0] <= distance.shape[0]) & (
                    particles_opt[:, 1] >= 1) & (particles_opt[:, 1] <= distance.shape[1])
        particles_opt = particles_opt[is_inbound]
        particles_opt[:, 2] = distance[
                                  particles_opt[:, 0].astype(int), particles_opt[:, 1].astype(int)] + particles_opt[:,
                                                                                                      2]
        is_opt = ~is_nodata[(particles_opt[:, 0] - 1).astype(int), (particles_opt[:, 1] - 1).astype(int)]
        particles_opt = particles_opt[is_opt]
        particles_opt = particles_opt[np.argsort(particles_opt[:, 2])]
        unique_indices = np.unique(particles_opt[:, 3], axis=0, return_index=True)[1]
        particles_opt = particles_opt[unique_indices]
        particles_opt = particles_opt[np.argsort(particles_opt[:, 3])]

        particles = particles[(particles_opt[:, 3] - 1).astype(int), :]
        particles[:, 2] = particles_opt[:, 2]

        # Calculate azimuth and update distances
        p_azimuth = particles_opt[:, :2] - particles[:, :2]

        particles[:, 3] = np.degrees(np.arccos(
            np.sum(np.column_stack((p_azimuth[:, 1], - p_azimuth[:, 0])) * np.tile([0, -1], [p_azimuth.shape[0], 1]),
                   axis=1)
            / np.linalg.norm(p_azimuth, axis=1))) * (1 - (p_azimuth[:, 1] > 0) * 2) + (p_azimuth[:, 1] > 0) * 360
        particles[particles[:, 2] < lim_azimuth, 3] = 0
        particles = particles[np.argsort(particles[:, 2])]
        unique_indices = np.unique([particles[:, 0], particles[:, 1]], axis=1, return_index=True)[1]
        particles = particles[np.sort(unique_indices)]

        runup_opt = 18.093 * min(max_vol, sites['Volume'][k] / 10 ** 6) ** 0.57110 * (
                    particles[:, 2] / 1000) ** -0.74189
        ind_opt = particles[:, 0].astype(int), particles[:, 1].astype(int)
        is_opt = (distance[ind_opt] > particles[:, 2]) & ((elev_model['dtm'][ind_opt] - runup_opt) < elev_source) & (
                    runup_opt >= min_wave)
        particles = particles[is_opt]
        distance[particles[:, 0].astype(int), particles[:, 1].astype(int)] = particles[:, 2]
        azimuth[particles[:, 0].astype(int), particles[:, 1].astype(int)] = particles[:, 3]  # update azimuth

        # Update is_change and azimuth, clip to extent
        is_change = len(particles) > 0
        if np.nanmax(distance[distance < 9999990]) > 100000:
            print("One or more particles are farther than 100km. Stopping wave propagation.")
            is_change = False

        if is_change:
            w_ext_subs = np.array(
                [[max(w_ext_subs[0, 0], np.max(particles[:, 0])), min(w_ext_subs[0, 1], np.min(particles[:, 1]))],
                 [min(w_ext_subs[1, 0], np.min(particles[:, 0])), max(w_ext_subs[1, 1], np.max(particles[:, 1]))]])

        if m == 1:
            subs_24 = subs_24[:, :-8]
            dist_opt = dist_opt[:-8]
        elif m == 2:
            subs_24 = subs_24[:, :-16]
            dist_opt = dist_opt[:-16]
        m += 1
        # Save video frame
        if save_vid:
            if m == 10:
                plt.scatter(particles[:, 1], -particles[:, 0], 30, np.log10(np.minimum(50, runup_opt[is_opt])),
                            marker='.')
                plt.gca().set_aspect('equal', adjustable='box')
                plt.draw()
                frame = np.frombuffer(plt.gcf().canvas.tostring_rgb(), dtype=np.uint8).reshape(
                    plt.gcf().canvas.get_width_height()[::-1] + (3,))
                v.write(cv2.cvtColor(frame, cv2.COLOR_RGB2BGR))
                m = 3

    if save_vid:
        v.release()
    distance = distance[w_ext_subs[1, 0].astype(int):w_ext_subs[0, 0].astype(int) + 1,
               w_ext_subs[0, 1].astype(int):w_ext_subs[1, 1].astype(int) + 1]
    azimuth = azimuth[w_ext_subs[1, 0].astype(int):w_ext_subs[0, 0].astype(int) + 1,
              w_ext_subs[0, 1].astype(int):w_ext_subs[1, 1].astype(int) + 1]

    # Calculate runup and export results
    runup = 18.093 * np.minimum(max_vol, sites.Volume[k] / 1e6) ** 0.57110 * (distance / 1000) ** -0.74189
    runup = np.minimum(max_wave, runup)
    runup = np.where(runup < min_wave, 0, runup)

    Lim = 100000

    runup[distance > Lim] = 0
    matrix = runup.copy()
    max_index = np.nanargmax(matrix)
    starting_point = np.unravel_index(max_index, matrix.shape)
    threshold = 100

    lake = np.where(matrix > 0.1, 1, 0)
    print(f'Done runup now correction k{k + 1}: site id {sites.Id[k]} of {sites.shape[0]}')
    try:
        final_raster_stack = CorAngle.compute_lake_angles(lake, starting_point, threshold)
        for i in range(len(final_raster_stack[:, 0, 0])):
            alpha = final_raster_stack[i, :, :]
            if i == 0:
                # Assuming matrix and alpha are already defined
                FinalRunup = np.where(
                    alpha <= 90,
                    matrix * (0.7 + (0.3 * np.cos(np.radians(alpha)) ** 2)),  # For alpha <= 90 degrees
                    matrix * (0.4 + (0.3 * np.cos(np.radians(alpha - 90)) ** 2))  # For alpha > 90 degrees
                )
            else:
                FinalRunup = np.where(
                    alpha <= 90,
                    FinalRunup * (0.7 + (0.3 * np.cos(np.radians(alpha)) ** 2)),  # For alpha <= 90 degrees
                    FinalRunup * (0.4 + (0.3 * np.cos(np.radians(alpha - 90)) ** 2))  # For alpha > 90 degrees
                )
    except Exception as e:
        error_log_path = "/lustre02/home/sylvain/error_log_angles.txt"
        with open(error_log_path, "a") as f:
            f.write(f"Error in process_sites k{k + 1}: site id {sites.Id[k]} of {sites.shape[0]}: {str(e)}\n")
    xll = elev_model['xllcorner'] + elev_model['cellsize'] * (w_ext_subs[0, 1] - 1)
    yll = elev_model['yllcorner'] + elev_model['cellsize'] * (elev_model['nrows'] - w_ext_subs[0, 0])
    yur = yll + (w_ext_subs[0, 0] - w_ext_subs[1, 0] + 1) * elev_model['cellsize']

    # Write to geotiff
    output_file = f'{site_dir}/runup_id{int(sites.WaveID[k])}.tif'
    print(output_file, np.sum(runup))

    geotiffwrite(output_file, FinalRunup, 'float32', elev_model['cellsize'], xll, yur)


def process_sites(sites, elev_model, is_nodata, result_dir, max_vol, min_wave, max_wave, IDs, ScenarioIDs,
                  lim_azimuth=0, save_vid=False):
    with ThreadPoolExecutor(max_workers=10, initializer=init_worker, initargs=(elev_model, is_nodata)) as executor:
        futures = [executor.submit(process_site, k, sites, result_dir, max_vol, min_wave, max_wave, IDs, ScenarioIDs,
                                   lim_azimuth=0, save_vid=False)
                   for k in range(sites.shape[0])]
        for future in as_completed(futures):
            result = future.result()


