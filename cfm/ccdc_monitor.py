""" Module for creating probability of change maps based on previous time series results
"""
import logging
import os
import re
from datetime import datetime
import click
import numpy as np
from osgeo import gdal
import patsy
from config_parser import parse_config_file
from utils import design_to_indices, find_results, iter_records

gdal.AllRegister()
gdal.UseExceptions()

logger = logging.getLogger('yatsm')
logger_algo = logging.getLogger('yatsm_algo')

pattern = '*npz'
_result_record = 'yatsm_r*'


def ccdc_monitor(cfg, date, image_ds):
    """Update previous CCDC results based on newly available imagery
       This implementation transfers the work-flow from line-based to
       2-d array based. When no regression is needed, indexing with
       arrays instead of iterating over lines improves performance
       exponentially.
    """
    logger_algo.setLevel(logging.DEBUG)

    #Open new mage as array
    image_ar = image_ds.ReadAsArray()

    #Find previous result location containing numpy-saved results sorted by row.
    result = cfg['dataset']['output']

    #Convert the date to proleptic Gregorian ordinal date
    date = datetime.strptime(str(date), '%Y%j').toordinal()

    #Do the monitoring
    consec = do_monitor(
            date, result, image_ds, image_ar, cfg
        )

    return consec



def get_mon_changes(scores, consec_new, consec, ndvi, mask, threshold, date, confirmed, previous_ds):
    """Return pixels that exceed threshold of deforestation on given date
    """

    #Flag change for areas above threshold. Note: Threshold is negative because we only want large decreases in NDVI.
    consec_new[scores[:,:,ndvi] < (-threshold)] = 1

    #Mask clouds
    consec_new[mask > 0] = 0

    #Set previously confirmed pixels to 0
    consec_new[confirmed > 0] = 0
    consec_new[previous_ds > 10] = 0

    #Add newly flagged pixels to old 'consecutive' array
    consec += consec_new

    #Turn back to 0 if cloud free and not increasing
    consec[(np.logical_and(mask == 0, consec_new == 0))] = 0

    #Create arrays to hold probabilities and confirmed changes
    lowprob = np.zeros_like(consec)
    highprob = np.zeros_like(consec)
    confirmed_today = np.zeros_like(consec)
    d = datetime.fromordinal(date)
    date = int(d.strftime("%Y%j"))

    #Assign probabilities
    lowprob[consec == 3] = date
    highprob[consec == 4] = date
    confirmed_today[np.logical_and(consec_new == 1, consec == 5)] = date

    return consec, lowprob, highprob, confirmed_today



def mask_clouds(cloud_mask,scores, image_ar, mask_values, cloudthreshold, green, swir, shadowthreshold, mask_band, vza_band, vza_threshold, test_indices):
    """Convert masked pixels to 1"""

    #Mask anything that is not 1 in QA band
    cloud_mask[image_ar[mask_band,:,:] != 1] = 1

    #Multiple-temporal cloud and shadow screening based on change scores
    cloud_mask[scores[:,:, green] > cloudthreshold] = 1
    cloud_mask[scores[:,:, swir] < shadowthreshold] = 1

    #Mask pixels above view-zenith angle threshold
    cloud_mask[image_ar[vza_band,:,:] > vza_threshold] = 1

    return cloud_mask

def mask_nonforest(f_mask, image_ar, coef, ndvi, rmse, ndvi_threshold, slope_threshold, rmse_threshold, test_indices):
    """Convert nonforest pixels to 1"""

    #Mask pixels whose starting model is below NDVI threshold
    f_mask[coef[:,:,0] < ndvi_threshold] = 1

    #Mask pixels whose starting model is above slope threshold. It is unlikely a stable forest will have a high slope.
    f_mask[coef[:,:,1] > slope_threshold] = 1

    #Mask pixels whose starting model is above RMSE threshold.
    f_mask[rmse[:,:,ndvi] > rmse_threshold] = 1

    return f_mask

def get_mon_scores(raster, image_ar, n_bands, i_bands, rmse, scores):

    """Return change scores based on model RMSE. Also returns VZA weight."""

    scores=np.zeros((np.shape(raster)[0],np.shape(raster)[1], n_bands))
    for _band in i_bands:
        scores[:,:,_band]=(image_ar[_band,:,:].astype(float) - raster[:,:,_band])/rmse[:,:,_band]
    return scores



def find_mon_result_attributes(results, bands, coefs, config, prefix=''):
    """ Returns attributes about the dataset from result files"""

    _coef = 'coef'
    _rmse = 'rmse'

    n_bands, n_coefs = None, None
    design = None
    for r in results:
        try:
            _result = np.load(r)
            rec = _result['record']
            if 'metadata' in _result.files:
                logger.debug('Finding X design info for version>=v0.5.4')
                md = _result['metadata'].item()
                design = md['YATSM']['design']
                design_str = md['YATSM']['design_matrix']
            else:
                logger.debug('Finding X design info for version<0.5.4')
                design = _result['design_matrix'].item()
                design_str = _result['design'].item()
        except:
            continue

        if not rec.dtype.names:
            continue

        if _coef not in rec.dtype.names or _rmse not in rec.dtype.names:
            logger.error('Could not find coefficients ({0}) and RMSE ({1}) '
                         'in record'.format(_coef, _rmse))
            if prefix:
                logger.error('Coefficients and RMSE not found with prefix %s. '
                             'Did you calculate them?' % prefix)
            raise click.Abort()

        try:
            n_coefs, n_bands = rec[_coef][0].shape
        except:
            continue
        else:
            break

    if n_coefs is None or n_bands is None:
        logger.error('Could not determine the number of coefficients or bands')
        raise click.Abort()
    if design is None:
        logger.error('Design matrix specification not found in results.')
        raise click.Abort()

    # How many bands does the user want?
    if bands == 'all':
        i_bands = range(0, n_bands)
    else:
        # NumPy index on 0; GDAL on 1 -- so subtract 1
        i_bands = [b - 1 for b in bands]
        if any([b > n_bands for b in i_bands]):
            logger.error('Bands specified exceed size of bands in results')
            raise click.Abort()

    # How many coefficients did the user want?
    use_rmse = False
    if coefs:
        if 'rmse' in coefs or 'all' in coefs:
            use_rmse = True
        i_coefs, coef_names = design_to_indices(design, coefs)
    else:
        i_coefs, coef_names = None, None

    logger.debug('Bands: {0}'.format(i_bands))
    if coefs:
        logger.debug('Coefficients: {0}'.format(i_coefs))

    return (i_bands, i_coefs, use_rmse, coef_names, design_str, design)

def find_mon_indices(record, date, after=False, before=False):
    """ Yield indices matching time segments for a given date"""

    #Find most recent model that starts before the monitoring date
    index = np.where(record['start'] <= date)[0]
    index = index[::-1]
    _, _index = np.unique(record['px'][index], return_index=True)
    index = index[_index]
    yield index


def get_prediction(index, rec, i_coef, n_i_bands, X, i_bands, raster):
    """ Get prediction for date of input based on model coeffiecients. """

    _coef = rec['coef'].take(index, axis=0).\
            take(i_coef, axis=1).take(i_bands, axis=2)
    return np.tensordot(_coef, X, axes=(1, 0))


def do_monitor(date, result_location, image_ds, image_ar, cfg,
		   bands='all', prefix='', ndv=-9999, pattern=_result_record):

    """ Get change prediction for date of input image.
	Outputs:
		Consec: 2-d array indicating consecutive days above threshold
		Confirmed: 2-d array indicating date of confirmed change
		lowprob: 2-d array indicating date of low confidence change
		highprob: 2-d array indicating date of high confidence change
	Probabilities:
		lowprob: 3 consecutive days above threshold
		highprob: 4 consecutive days above threshold
		Confirmed: 5 conscutive days above threshold
"""

    #Define parameters. This should be from the parameter file.
    threshold = cfg['CCDCesque']['threshold']
    test_ind = cfg['CCDCesque']['test_indices']
    cloudthreshold = cfg['CCDCesque']['screening_crit']
    green = cfg['dataset']['green_band'] - 1
    swir = cfg['dataset']['swir1_band'] - 1
    shadowthreshold = 0 - (cfg['CCDCesque']['screening_crit'])
    mask_band = cfg['dataset']['mask_band'] - 1
    mask_values = cfg['dataset']['mask_values']
    vza_band = cfg['NRT']['view_angle_band'] - 1
    vza_threshold = cfg['NRT']['view_angle_threshold'] * 100
    try:
        ndvi = cfg['NRT']['ndvi']
    except:
	ndvi = cfg['test_indices'][0] #Can only use one band at moment
    rmse_threshold = cfg['NRT']['rmse_threshold']
    ndvi_threshold = cfg['NRT']['ndvi_threshold']
    slope_threshold = cfg['NRT']['slope_threshold']
    consec_file = cfg['NRT']['out_file']
    shapefile_dir = cfg['NRT']['master_shapefile_dir']
    previous = cfg['NRT']['previousresult']
    try:
	begin_monitor = cfg['NRT']['begin_monitor']
    except:
	begin_monitor = 2016001
    # Find results
    records = find_results(result_location, pattern)

    # Find result attributes to extract
    i_bands, _, _, _, design, design_info = find_mon_result_attributes(
        records, bands, None, cfg, prefix=prefix)
    n_bands = len(i_bands)
    n_i_bands = len(i_bands)
    im_Y = image_ds.RasterYSize
    im_X = image_ds.RasterXSize

    #Initiate arrays for prediction, change scores, rmse, coefficients, and consecutive days above threshold
    raster = np.zeros((im_Y, im_X, n_bands), dtype=np.float16) * int(ndv)
    cloud_mask = np.zeros((im_Y, im_X), dtype=np.int32) * int(ndv)
    f_mask = np.zeros((im_Y, im_X), dtype=np.int32) * int(ndv)
    rmse = np.zeros((im_Y, im_X, n_bands), dtype=np.float16) * int(ndv)
    scores = np.zeros((im_Y, im_X, n_bands), dtype=np.float16) * int(ndv)
    consec_new = np.zeros((im_Y, im_X), dtype=np.int32) * int(ndv)

    # Create X matrix from date -- ignoring categorical variables

    if re.match(r'.*C\(.*\).*', design):
        logger.warning('Categorical variable found in design matrix not used'
                       ' in predicted image estimate')
    design = re.sub(r'[\+\-][\ ]+C\(.*\)', '', design)
    X = patsy.dmatrix(design, {'x': date}).squeeze()

    i_coef = []
    for k, v in design_info.iteritems():
        if not re.match('C\(.*\)', k):
            i_coef.append(v)
    i_coef = np.asarray(i_coef)

    coef = np.zeros((im_Y, im_X, len(i_coef)), dtype=np.float16) * int(ndv)

    logger.info('Creating prediction')
    for z, rec in iter_records(records, warn_on_empty=True):

       for index in find_mon_indices(rec, date):
           if index.shape[0] == 0:
                continue

           #Get x and y location from records
	   _px = rec['px'][index]
	   _py = rec['py'][index]

	   #First, get prediction
           raster[_py ,_px, :n_i_bands] = get_prediction(index, rec, i_coef,
					        	 n_i_bands, X, i_bands, raster)
           #Then RMSE
           rmse[_py, _px, :n_i_bands] = rec['rmse'][index][:, i_bands]

          #Finally, the coefficients
           rec['coef'][:, 0, :] += (
           (rec['start'] + rec['end'])
            / 2.0)[:, np.newaxis] * \
            rec['coef'][:, 1, :]

           coef[_py ,_px, :len(i_coef)] = rec['coef'][:,:,ndvi][index]

    #Calculate change score for each band
    logger.info('Calculating change scores')
    scores = get_mon_scores(raster, image_ar, n_bands, i_bands, rmse, scores)

    # Mask out everything based on thresholds
    logger.info('Creating cloud mask')
    cloud_mask = mask_clouds(cloud_mask, scores, image_ar, mask_values, cloudthreshold, green, swir, shadowthreshold, mask_band, vza_band, vza_threshold, test_ind)

    #We're only interested in places that are deforestation
    logger.info('Creating forest mask')
    f_mask = mask_nonforest(f_mask, image_ar, coef, ndvi, rmse, ndvi_threshold, slope_threshold, rmse_threshold, test_ind)

    #Combine the masks
    mask = cloud_mask + f_mask
    #load monitor file
    if os.path.isfile(consec_file):
        consec = np.load(consec_file)['consec']
        confirmed = np.load(consec_file)['confirmed']
    else:
        consec = np.zeros((im_Y, im_X), dtype=np.int32) * int(ndv)
        confirmed = np.zeros((im_Y, im_X), dtype=np.int32) * int(ndv)
    #Reset previous changes
    consec[consec[:,:] >= 5] = 0

    #TODO Hack: don't overlap results from 2015
    previous_im = gdal.Open(previous, gdal.GA_ReadOnly)
    previous_ds = previous_im.ReadAsArray()
    confirmed[confirmed < begin_monitor] = 0

    #Calculate change probabilities for the day
    consec, lowprob, highprob, confirmed_today = get_mon_changes(scores, consec_new, consec, ndvi, mask, threshold, date, confirmed, previous_ds)

    #Add today's confirmed changes to previous confirmed changes
    confirmed[:,:] += confirmed_today[:,:]

    #Save the results
    out = {}
    for k, v in z.iteritems(): # Get meta data items from z
        out[k] = v
    out['consec'] = consec
    out['lowprob'] = lowprob
    out['highprob'] = highprob
    out['confirmed_today'] = confirmed_today
    out['confirmed'] = confirmed
    np.savez(consec_file, **out)
    return out

