
import Wave_Propa_par as WPp

import numpy as np
import time
import geopandas as gpd
import rasterio


# from processing import init_worker, process_site

def Wave_simulation(output_dir, wave_start_points_path, max_vol, min_wave, max_wave, save_vid, IDs, ScenarioIDs,
                    lim_azimuth, water_raster_path='', par=1):
    # Get the current time
    current_time1 = time.localtime()

    if not water_raster_path:
        interactive = 'ok'
    else:
        interactive = ''

    # load list of sites
    #    print(wave_start_points_path)
    sites_gpd = gpd.read_file(wave_start_points_path)
    sites = sites_gpd.drop(columns='geometry')
    # load raster water
    if 'elev_model' not in locals():
        print('load water tiff')
        if interactive == 'ok':
            elev_model = WP.dtm_import()
        else:
            with rasterio.open(water_raster_path) as dataset:
                elev_model = {
                    'dtm': dataset.read(1),
                    'nrows': dataset.height,
                    'ncols': dataset.width,
                    'cellsize': dataset.res[0],  # Assuming square pixels
                    'xllcorner': dataset.bounds.left,
                    'yllcorner': dataset.bounds.bottom,
                }

        is_nodata = np.logical_or(elev_model['dtm'] == 9999, elev_model['dtm'] <= -9999)
        elev_model['dtm'] = np.maximum(0, elev_model['dtm'])

    print(f"Start at: {current_time1.tm_year}-{current_time1.tm_mon:02d}-{current_time1.tm_mday:02d} "
          f"{current_time1.tm_hour:0.0f}h{current_time1.tm_min:0.0f}min{current_time1.tm_sec:0.0f}s")

 #   if par == 0:
 #       WPp.process_sites(sites, elev_model, is_nodata, output_dir, max_vol, min_wave, max_wave,
 #                        IDs, ScenarioIDs, lim_azimuth, save_vid=save_vid)
 #   elif par == 1:
    WPp.process_sites(sites, elev_model, is_nodata, output_dir, max_vol, min_wave, max_wave,
                          IDs, ScenarioIDs, lim_azimuth, save_vid=save_vid)

    current_time2 = time.localtime()
    Proccessing_time = time.mktime(current_time2) - time.mktime(current_time1)
    print(f"End wave simulation process took: {Proccessing_time}s")
