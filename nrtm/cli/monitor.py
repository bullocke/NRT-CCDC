""" Module for creating probability of change maps based on previous time series results
"""
from datetime import datetime
import logging
import os
import re
import click
import csv
import numpy as np
from osgeo import gdal, gdal_array
import patsy
from ..ccdc_monitor import *
import options
from monitor_map import make_map, write_shapefile
from ..config_parser import parse_config_file
from ..utils import get_mon_outputs, get_output_names, find_results, iter_records, write_output
gdal.AllRegister()
gdal.UseExceptions()

logger = logging.getLogger('yatsm')
logger_algo = logging.getLogger('yatsm_algo')

pattern = '*npz'
_result_record = 'yatsm_r*'

@click.command(short_help='Monitor for changes give up to 1 year of new images')
@options.arg_config_file
@options.arg_mon_csv
@options.opt_date_format
@options.opt_nodata
@options.opt_format
@click.pass_context


def monitor(ctx, config, mon_csv, gdal_frmt, date_frmt, ndv=0):
    """Command line interface to handle monitoring of new imagery. This program will not
     pre-process the data, which is done in yatsm.process_modis. This program will calculate
     the change probabilities in time-sequential order for all images in input monitoring log.
     Currently, the output is written to shapefiles for tileing on Mapbox. """

    logger_algo.setLevel(logging.DEBUG)

    #Parse config and open csv with images previously processed
    cfg = parse_config_file(config)
    done_csv = cfg['dataset']['input_file']
    read=csv.reader(open(done_csv,"rb"),delimiter=',')
    done_array = list(read)

    try:
	begin_monitor = cfg['NRT']['begin_monitor']
    except:
	begin_monitor = 2016001

    #Get first and last dates
    last = int(done_array[-1][0])
    veryfirst = int(done_array[1][0])

    #Read monitor csv with images to use in monitoring
    read_mon=csv.reader(open(mon_csv,"rb"),delimiter=',')
    monitor_array = list(read_mon)

    if monitor_array is None:
        logger.error('Incorrect path to monitor csv')
        raise click.Abort()

    if len(monitor_array) == 0:
        logger.error('Not new images to monitor')
        raise click.Abort()
    first = int(monitor_array[0][0])

    #Loop over each date in monitor list. Check again if the date is in input list
    num_monitor=len(monitor_array)
    for i in range(num_monitor):
        cur_image = monitor_array[i]
	date = int(cur_image[0])
	image_path = cur_image[1]
	if date <= last:
            logger.error('Previous results processed past image date. Skipping.')
            continue

        #Read the image as an array.
        try:
            image_ds = gdal.Open(image_path, gdal.GA_ReadOnly)
        except:
            logger.error('Could not open new image for reading')
            raise click.Abort()


        #Do monitor
        logger.info('Doing image %s' % image_path)
        out = ccdc_monitor(cfg, date, image_ds)

        #Get output file names
        output_lp_today, output_hp_today, output_lp, output_hp, output_conf, output_conf_today, master = get_mon_outputs(cfg, date)

	#Write out the shapefiles. Currently a copy is saved in the daily folders in addition to a master version.
        if np.any(out['lowprob'] > begin_monitor):
            write_shapefile(out['lowprob'], output_lp_today,image_ds, gdal_frmt,
     	    	            ndv, band_names=None)
	    if os.path.isfile(output_lp):
	        os.remove(output_lp)
            write_shapefile(out['lowprob'], output_lp,image_ds, gdal_frmt,
     	    	            ndv, band_names=None)
        if np.any(out['highprob'] > begin_monitor):
            write_shapefile(out['highprob'], output_hp_today,image_ds, gdal_frmt,
   	    	            ndv, band_names=None)
	    if os.path.isfile(output_hp):
	        os.remove(output_hp)
            write_shapefile(out['highprob'], output_hp,image_ds, gdal_frmt,
     	    	            ndv, band_names=None)
        if np.any(out['confirmed_today'] > begin_monitor):
            write_shapefile(out['confirmed_today'], output_conf_today,image_ds, gdal_frmt,
   	    	            ndv, band_names=None)
        if np.any(out['confirmed'] > begin_monitor):
	    if os.path.isfile(master):
	        os.remove(master)
            write_shapefile(out['confirmed'], master,image_ds, gdal_frmt,
   	    	            ndv, band_names=None)

	#update processed image csv
	out_log = [str(date),'Com',image_path]
	done_array.append(out_log)
	with open(done_csv, 'wb') as f:
	    writer = csv.writer(f)
	    writer.writerows(done_array)

	output_rast = None

