#!/usr/bin/env python2
""" Utils for pre-processing and compositing MODIS images
Modified version of scripts by Chris Holden
https://github.com/ceholden/misc/tree/master/composites
"""
from __future__ import division, print_function

import logging
import math

import click
import numpy as np
try:
    import progressbar
except:
    _has_progressbar = False
else:
    _has_progressbar = True
import rasterio
from rasterio.rio.options import _cb_key_val, creation_options
import snuggs
import fnmatch
import multiprocessing
import os
import sys
from docopt import docopt
import numexpr as ne
from osgeo import gdal
import scipy.ndimage


#Original author of these functions
__author__ = 'Chris Holden (ceholden@gmail.com)'

logging.basicConfig(format='%(asctime)s.%(levelname)s: %(message)s',
                    level=logging.INFO,
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

# Predefined algorithms -- name: snuggs expression
_ALGO = {
    'maxNDVI': '(max (/ (- nir red) (+ nir red)))',
    'medianNDVI': '(median (/ (- nir red) (+ nir red)))',
    'ZheZhu': '(max (/ nir blue))',
    'minBlue': '(min blue)',
    'minVZA': '(min vza)',
    'maxNIR': '(max nir)'
}

_context = dict(
    token_normalize_func=lambda x: x.lower(),
    help_option_names=['--help', '-h']
)


gdal.AllRegister()
gdal.UseExceptions()



compress_algos = ['LZW', 'PACKBITS', 'DEFLATE', 'None']

def _valid_band(ctx, param, value):
    if value is None:
        return None
    try:
        band = int(value)
        assert band >= 1
    except:
        raise click.BadParameter('Band must be integer above 1')
    return band

def image_composite(inputs, algo, output, oformat, vza, mask_band, mask_val):

    """ Create image composites based on some criteria
    Output image composites retain original values from input images that meet
    a certain criteria. For example, in a maximum NDVI composite with 10 input
    images, all bands for a given pixel will contain the band values from the
    input raster that had the highest NDVI value.
    Users can choose from a set of predefined compositing algorithms or may
    specify an Snuggs S-expression that defines the compositing criteria.
    Normalized Differenced indexes can be computed using "(normdiff a b)" for
    the Normalized Difference between "a" and "b" (or "nir" and "red").
    See https://github.com/mapbox/snuggs for more information on Snuggs
    expressions.
    The indexes for common optical bands (e.g., red, nir, blue) within the
    input rasters are included as optional arguments and are indexed in
    wavelength sequential order. You may need to overwrite the default indexes
    of bands used in a given S-expression with the correct band index.
    Additional bands may be identified and indexed using the
    '--band NAME=INDEX' option.
    Currently, input images must be "stacked", meaning that they contain the
    same bands and are the same shape and extent.
    Example:
    1. Create a composite based on maximum NDVI
        Use the built-in maxNDVI algorithm:
        \b
        $ image_composite.py --algo maxNDVI image1.gtif image2.gtif image3.gtif
            composite_maxNDVI.gtif
        or with S-expression:
        \b
        $ image_composite.py --expr '(max (/ (- nir red) (+ nir red)))'
            image1.gtif image2.gtif image3.gtif composite_maxNDVI.gtif
        or with S-expressions using the normdiff shortcut:
        \b
        $ image_composite.py --expr '(max (normdiff nir red))'
            image1.gtif image2.gtif image3.gtif composite_maxNDVI.gtif
    2. Create a composite based on median EVI (not recommended)
        With S-expression:
        \b
        $ evi='(median (/ (- nir red) (+ (- (+ nir (* 6 red)) (* 7.5 blue)) 1)))'
        $ image_composite.py --expr "$evi"  image1.gtif image2.gtif image3.gtif
            composite_medianEVI.gtif
    3. Create a composite based on median NBR
        With S-expression:
        \b
        $ image_composite.py --expr '(median (normdiff nir sswir))'
            image1.gtif image2.gtif image3.gtif composite_maxNBR.gtif
    """
    verbose = True
    if verbose:
        logger.setLevel(logging.DEBUG)
    elif quiet:
        logger.setLevel(logging.ERROR)

        expr = _ALGO[algo]
    if algo is not None:
        logger.debug('Using predefined algorithm: {}'.format(algo))
        expr = _ALGO[algo]


    # Setup band keywords
    _bands = {'vza': vza}

    # Find only the band names and indexes required for the composite criteria
    crit_indices = {k: v - 1 for k, v in _bands.iteritems() if k in expr}

    # Enhance snuggs expressions to return index of value matching function
    snuggs.func_map['max'] = lambda a: np.argmax(a, axis=0)
    snuggs.func_map['min'] = lambda a: np.argmin(a, axis=0)
    snuggs.func_map['median'] = lambda a: np.argmin(
        np.abs(a - np.median(a, axis=0)), axis=0)
    snuggs.func_map['normdiff'] = lambda a, b: snuggs.eval(
        '(/ (- a b) (+ a b))', **{'a':a, 'b':b})

    with rasterio.drivers():

        # Read in the first image to fetch metadata
        with rasterio.open(inputs[0]) as first:
            meta = first.meta
            if 'transform' in meta:
                meta.pop('transform')  # remove transform since deprecated
            meta.update(driver=oformat)
            if len(set(first.block_shapes)) != 1:
                click.echo('Cannot process input files - '
                           'All bands must have same block shapes')
                raise click.Abort()
            block_nrow, block_ncol = first.block_shapes[0]
            windows = first.block_windows(1)
            n_windows = math.ceil(meta['height'] / block_nrow *
                                  meta['width'] / block_ncol)

            # Ensure mask_band exists, if specified
            if mask_band:
                if mask_band <= meta['count'] and mask_band > 0:
                    mask_band -= 1
                else:
                    click.echo('Mask band does not exist in INPUT images')
                    raise click.Abort()

        # Initialize output data and create composite
        with rasterio.open(output, 'w', **meta) as dst:
            # Process by block
            dat = np.ma.empty((len(inputs), meta['count'],
                               block_nrow, block_ncol),
                              dtype=np.dtype(meta['dtype']))
            mi, mj = np.meshgrid(np.arange(block_nrow), np.arange(block_ncol),
                                 indexing='ij')
            # Open all source files one time
            srcs = [rasterio.open(fname) for fname in inputs]

            logger.debug('Processing blocks')
            if _has_progressbar:
                widgets = [
                    progressbar.Percentage(),
                    progressbar.BouncingBar(
                        marker=progressbar.RotatingMarker())
                ]
                pbar = progressbar.ProgressBar(widgets=widgets).start()

            for i, (idx, window) in enumerate(windows):
                # Update dat and mi, mj only if window changes
                nrow = window[0][1] - window[0][0]
                ncol = window[1][1] - window[1][0]
                if dat.shape[-2] != nrow or dat.shape[-1] != ncol:
                    dat = np.ma.empty((len(inputs), meta['count'],
                                       nrow, ncol),
                                      dtype=np.dtype(meta['dtype']))
                    mi, mj = np.meshgrid(np.arange(nrow), np.arange(ncol),
                                         indexing='ij')
                for j, src in enumerate(srcs):
                    dat[j, ...] = src.read(masked=True, window=window)
                    # Mask values matching mask_vals if mask_band
                    if mask_band and mask_val:
                        dat[j, ...].mask = np.logical_or(
                            dat[j, ...].mask,
                            np.in1d(dat[j, mask_band, ...], mask_val,).reshape(
                                dat.shape[-2], dat.shape[-1])
                        )

                # Find indices of files for composite
                crit = {k: dat[:, v, ...] for k, v in crit_indices.iteritems()}
                crit_idx = snuggs.eval(expr, **crit)

                # Create output composite
                # Use np.rollaxis to get (nimage, nrow, ncol, nband) shape
                composite = np.rollaxis(dat, 1, 4)[crit_idx, mi, mj]

                # Write out
                for i_b in range(composite.shape[-1]):
                    dst.write(composite[:, :, i_b], indexes=i_b + 1,
                              window=window)
                if _has_progressbar:
                    pbar.update(int((i + 1) / n_windows * 100))

#if __name__ == '__main__':
 #   image_composite()

def find_MODIS_pairs(location, pattern='M*09*hdf'):
    """ Finds matching sets of M[OY]D09GQ and M[OY]D09GA within location

    Args:
      location (str): directory of stored data
      pattern (str, optional): glob pattern to limit search

    Returns:
      pairs (list): list of tuples containing M[OY]D09GQ and M[OY]D09GA

    """
    files = [os.path.join(location, f) for f in
             fnmatch.filter(os.listdir(location), pattern)]

    if len(files) < 2:
        raise IOError('Could not find any MODIS image pairs')

    # Parse out product and acquisition date
    products = []
    dates = []
    for f in files:
        s = os.path.basename(f).split('.')
        products.append(s[0])
        dates.append(s[1])

    products = np.array(products)
    dates = np.array(dates)

    # Retain dates if there are matching MOD09GA/MOD09GQ or MYD09GA/MYD09GQ
    pairs = []
    for d in np.unique(dates):
        i = np.where(dates == d)[0]
        prods = products[i]

        i_mod09ga = np.core.defchararray.startswith(prods, 'MOD09GA')
        i_mod09gq = np.core.defchararray.startswith(prods, 'MOD09GQ')
        i_myd09ga = np.core.defchararray.startswith(prods, 'MYD09GA')
        i_myd09gq = np.core.defchararray.startswith(prods, 'MYD09GQ')

        if i_mod09gq.sum() == 1 and i_mod09ga.sum() == 1:
            pairs.append((files[i[i_mod09gq]], files[i[i_mod09ga]]))
        if i_myd09gq.sum() == 1 and i_myd09ga.sum() == 1:
            pairs.append((files[i[i_myd09gq]], files[i[i_myd09ga]]))

    return pairs


def get_output_names(pairs, outdir):
    """ Get output filenames for pairs of MODIS images

    Output filename pattern will be:

        M[OY]D_A[YYYY][DOY]_stack.gtif

    Args:
      pairs (list): list of tuples containing M[OY]D09GQ and M[OY]D09GA
      outdir (str): output directory

    Returns:
      output (list): list of output filenames

    """
    output_names = []
    for pair in pairs:
        _temp = os.path.basename(pair[0]).split('.')
        out_name = _temp[0][0:3] + '_' + _temp[1] + '_stack.gtif'
        output_names.append(os.path.join(outdir, out_name))

    return output_names


def check_resume(pairs, output):
    """ Returns pairs and output with pre-existing results removed

    Args:
      pairs (list): list of tuples containing M[OY]D09GQ and M[OY]D09GA
      output (list): list of output filenames for each entry in pairs

    Returns:
      pairs, output (tuple): pairs and outputs with pre-existing files removed

    """
    out_pairs = []
    out_output = []
    for p, o in zip(pairs, output):
        if os.path.isfile(o):
            try:
                ds = gdal.Open(o, gdal.GA_ReadOnly)
                ds = None
            except:
                logger.warning('File {f} already exists but cannot be opened. \
                    Overwriting'.format(f=o))
            else:
                logger.debug('Skipping output {f}'.format(f=o))
                continue
        out_pairs.append(p)
        out_output.append(o)

    return (out_pairs, out_output)


def create_stack(pair, output, ndv=-28672, compression='None',
                 tiled=False, blockxsize=None, blockysize=None):
    """ Create output stack from MODIS image pairs (M[OY]D09GQ & M[OY]D09GA)

    Args:
      pair (tuple): pairs of images (M[OY]D09GQ & M[OY]D09GA)
      output (str): location to output stack
      ndv (int, optional): NoDataValue for output
      compression (str, optional): compression algorithm to use
      tiled (bool, optional): use tiles instead of stripes
      blockxsize (int, optional): tile width
      blockysize (int, optional): tile or strip height

    Stack will be formatted as:
        Band        Definition
        ----------------------
        1           Normalized Difference Vegetation Index (NDVI) * 10000
        2           500m blue from M[OY]D09GA
        3           500m SWIR1 from M[OY]D09GA
        4           1000m Mask band from M[OY]D09GA
        5           1000m VZA * 100 from M[OY]D09GA

    Mask values     Definition
    --------------------------
        0           Not-clear, or not-land surface

    """
    # Outut band descriptions
    bands = ['NDVI * 10000',
             '500m blue * 10000', '500m SWIR1 * 10000',
             '1km QAQC Mask', '1km VZA * 100']

    ga_state = 1  # 1km state flags
    ga_vza = 2  # 1km view zenith angle
    ga_blue = 12 # 500m blue band
    ga_green = 13  # 500m green band
    ga_swir1 = 15  # 500m swir1 band
    ga_qc = 17  # 500m QC band
    gq_red = 1  # 250m red band
    gq_nir = 2  # 250m NIR band
    gq_qc = 3  # 250m QC band

    compress_algos = ['LZW', 'PACKBITS', 'DEFLATE', 'None']

    # Open and find subdatasets
    gq_ds = gdal.Open(pair[0], gdal.GA_ReadOnly)
    ga_ds = gdal.Open(pair[1], gdal.GA_ReadOnly)

    # Read in datasets
    ds_red = gdal.Open(gq_ds.GetSubDatasets()[gq_red][0], gdal.GA_ReadOnly)
    ds_nir = gdal.Open(gq_ds.GetSubDatasets()[gq_nir][0], gdal.GA_ReadOnly)

    ds_blue = gdal.Open(ga_ds.GetSubDatasets()[ga_blue][0], gdal.GA_ReadOnly)
    ds_state = gdal.Open(ga_ds.GetSubDatasets()[ga_state][0], gdal.GA_ReadOnly)
    ds_vza = gdal.Open(ga_ds.GetSubDatasets()[ga_vza][0], gdal.GA_ReadOnly)
    ds_swir1 = gdal.Open(ga_ds.GetSubDatasets()[ga_swir1][0], gdal.GA_ReadOnly)

    # Setup options
    opts = []
    if tiled:
        opts.append('TILED=YES')

    if blockxsize and tiled:
        opts.append('BLOCKXSIZE=%s' % blockxsize)
    if blockysize:
        opts.append('BLOCKYSIZE=%s' % blockysize)

    if compression:
        opts.append('COMPRESS=%s' % compression)
        if compression == 'LZW':
            opts.append('PREDICTOR=2')

    # Create file and setup metadata
    driver = gdal.GetDriverByName('GTiff')

    out_ds = driver.Create(output,
                           ds_red.RasterYSize, ds_red.RasterXSize, 5,
                           gdal.GDT_Int16,
                           options=opts)

    out_ds.SetProjection(ds_red.GetProjection())
    out_ds.SetGeoTransform(ds_red.GetGeoTransform())

    out_ds.SetMetadata(ga_ds.GetMetadata())
    out_ds.GetRasterBand(1).SetNoDataValue(-28672)

    # Write out image
    red = ds_red.GetRasterBand(1).ReadAsArray().astype(np.int16)
    nir = ds_nir.GetRasterBand(1).ReadAsArray().astype(np.int16)
    blue = enlarge(ds_blue.GetRasterBand(1).ReadAsArray().astype(np.int16),
                    2)
    swir1 = enlarge(ds_swir1.GetRasterBand(1).ReadAsArray().astype(np.int16),
                    2)

    # Calculate SR and NDVI, masking invalid values
    eqn = '(red <= 0) | (red >= 10000) | (nir <= 0) | (nir >= 10000) | ' \
        ' (swir1 <= 0) | (swir1 >= 10000) | ' \
        '(blue <= 0) | (blue >= 10000)'
    invalid = ne.evaluate(eqn)

    ndvi = ne.evaluate('(nir - red) / (nir + red) * 10000').astype(np.int16)


    # Apply valid range mask to data
    blue[invalid] = ndv
    red[invalid] = ndv
    nir[invalid] = ndv
    swir1[invalid] = ndv

    # Write out
    out_ds.GetRasterBand(1).WriteArray(ndvi)
    out_ds.GetRasterBand(2).WriteArray(blue)
    out_ds.GetRasterBand(3).WriteArray(swir1)

    # Perform masking -- 1 for land and VZA in separate band
    mask = get_mask(ds_state.GetRasterBand(1).ReadAsArray()).astype(np.int16)
    out_ds.GetRasterBand(4).WriteArray(enlarge(mask, 4))

    vza = ds_vza.GetRasterBand(1).ReadAsArray().astype(np.int16)
    out_ds.GetRasterBand(5).WriteArray(enlarge(vza, 4))

    # Write data
    for i, desc in enumerate(bands):
        out_ds.GetRasterBand(i + 1).SetDescription(desc)
        out_ds.GetRasterBand(i + 1).SetMetadataItem('Description', desc)

    # Close
    ga_ds = None
    gq_ds = None

    ds_red = None
    ds_nir = None
    ds_vza = None
    ds_state = None
    ds_swir1 = None
    ds_blue = None

    out_ds = None


def get_mask(modQA, dilate=1):
    """ Return a mask image from an input QA band from M[OY]D09G[AQ]

    Args:
      modQA (ndarray): input QA/QC band image
      dilate (int, optional): pixels around aerosols and clouds to buffer

    Returns:
      mask (ndarray): output mask image with only good observations masked

    Porting note:
        MATLAB's 'bitshift' shifts to the left

    """
    # Identify land from water
    land = (np.mod(np.right_shift(modQA, 3) + 6, 8) / 7).astype(np.uint8)
    # Identify cloud
    cloud = (np.mod(modQA, 8) |  # unsure!
             np.mod(np.right_shift(modQA, 8), 4) |  # cirrus == '00' (none)
             np.mod(np.right_shift(modQA, 10), 2) |  # cloud mask == '0'
             np.mod(np.right_shift(modQA, 13), 2)) > 0  # adjacent to cloud

    cloud_buffer = scipy.ndimage.morphology.binary_dilation(
        cloud, structure=np.ones((dilate, dilate)))

    return ((cloud_buffer == 0) * land)


def enlarge(array, scaling):
    """ Enlarge an array by a scaling factor

    Args:
      array (ndarray): array to be scaled
      scaling (int): amount of scaling

    Returns:
      scaled (ndarray): scaled array

    """
    return np.kron(array, np.ones((scaling, scaling)))


if __name__ == '__main__':
    args = docopt(__doc__, version=_version)

    logger.setLevel(logging.INFO)
    if args['--quiet']:
        logger.setLevel(logging.WARNING)
    if args['--verbose']:
        logger.setLevel(logging.DEBUG)

    location = args['<input_location>']
    outdir = args['<output_location>']
    if not os.path.isdir(outdir):
        # Make output directory
        try:
            os.makedirs(outdir)
        except:
            if os.path.isdir(outdir):
                pass
            else:
                logger.error('Output directory does not exist and could '
                             'not create output directory')
                raise

    # Optional arguments
    pattern = args['--pattern']

    ndv = int(args['--nodata'])

    resume = args['--resume']

    skip_error = args['--skip_error']

    compression = args['--compression']
    if compression.lower() == 'none':
        compression = None
    else:
        if compression not in compress_algos:
            logger.error('Unknown compression algorithm. Available are:')
            for _algo in compress_algos:
                logger.error('    %s' % _algo)
            sys.exit(1)

    tiled = args['--tiled']

    blockxsize = args['--blockxsize']
    if blockxsize.lower() == 'none':
        blockxsize = None
    else:
        blockxsize = int(blockxsize)

    blockysize = args['--blockysize']
    if blockysize.lower() == 'none':
        blockysize = None
    else:
        blockysize = int(blockysize)

    ncpu = int(args['--ncpu'])
    if ncpu > 1:
        logger.warning('NCPU only supported for `numexpr` so far...')
    ne.set_num_threads(ncpu)

    logger.debug('Finding pairs of MODIS data')

    pairs = find_MODIS_pairs(location, pattern)
    output_names = get_output_names(pairs, outdir)

    logger.info('Found {n} pairs of M[OY]D09GQ and M[OY]D09GA'.format(
        n=len(pairs)))

    if resume:
        pairs, output_names = check_resume(pairs, output_names)
        logger.info('Resuming calculation for {n} files'.format(n=len(pairs)))

    failed = 0
    for i, (p, o) in enumerate(zip(pairs, output_names)):
        logger.info('Stacking {i} / {n}: {p}'.format(
            i=i, n=len(pairs), p=os.path.basename(p[0])))
        try:
            create_stack(p, o,
                         ndv=ndv, compression=compression,
                         tiled=tiled,
                         blockxsize=blockxsize, blockysize=blockysize)
        except RuntimeError:
            if skip_error:
                logger.error('Could not preprocess pair: {p1} {p2}'.format(
                    p1=p[0], p2=p[1]))
                failed += 1
            else:
                raise

    if failed != 0:
        logger.error('Could not process {n} pairs'.format(n=failed))
    logger.info('Complete')

def get_mon_outputs(cfg, date):
    "Get output shapefile names"
    out_dir = cfg['NRT']['master_shapefile_dir']
    if not os.path.isdir(out_dir):
	os.mkdir(out_dir)
    date_path = '%s/%s' % (out_dir, date)
    if not os.path.isdir(date_path):
	os.mkdir(date_path)
    output_lp_today = '%s/lowprob.shp' % (date_path)
    output_hp_today = '%s/highprob.shp' % (date_path)
    output_lp = '%s/lowprob.shp' % (out_dir)
    output_hp = '%s/highprob.shp' % (out_dir)
    output_conf = '%s/confirmed.shp' % (date_path)
    output_conf_today = '%s/confirmed.shp' % (date_path)
    master = cfg['NRT']['master_shapefile']
    return output_lp_today, output_hp_today, output_lp, output_hp, output_conf, output_conf_today, master


def design_to_indices(design_matrix, features):
    """ Return indices of coefficients for features in design matrix. Taken
        from Chris Holden's YATSM package.

    Args:
      design_matrix (OrderedDict): OrderedDict containing design features keys
        and indices of coefficient matrix as values
      features (list): list of feature coefficients to extract

    Return:
      tuple: list of indices and names for each feature specified in `features`

    """
    if 'all' in features:
        features = design_coefs[1:]

    i_coefs = []
    coef_names = []
    for c in features:
        if c == 'intercept':
            k = _key_lookup_ignorecase(design_matrix, 'intercept')
            i_coefs.append(design_matrix.get(k))
            coef_names.append(k)
        elif c == 'slope':
            k = _key_lookup_ignorecase(design_matrix, 'x')
            i_coefs.append(design_matrix.get(
                _key_lookup_ignorecase(design_matrix, 'x')))
            coef_names.append(k)
        elif c == 'seasonality':
            i = [k for k in design_matrix.keys() if 'harm' in k]
            i_coefs.extend([design_matrix[_i] for _i in i])
            coef_names.extend(i)
        elif c == 'categorical':
            i = [k for k in design_matrix.keys() if 'C' in k]
            i_coefs.extend([design_matrix[_i] for _i in i])
            coef_names.extend(i)

    i_coefs = [i for i in i_coefs if i is not None]
    coef_names = [n for n in coef_names if n is not None]

    return i_coefs, coef_names

def find_results(location, pattern):
    """ Create list of result files and return sorted

    Args:
      location (str): directory location to search
      pattern (str): glob style search pattern for results

    Returns:
      results (list): list of file paths for results found

    """
    # Note: already checked for location existence in main()
    records = []
    for root, dirnames, filenames in os.walk(location):
        for filename in fnmatch.filter(filenames, pattern):
            records.append(os.path.join(root, filename))

    if len(records) == 0:
        logger.error('Error: could not find results in: {0}'.format(location))
        sys.exit(1)

    records.sort()

    return records


def iter_records(records, warn_on_empty=False, yield_filename=False):
    """ Iterates over records, returning result NumPy array

    Args:
      records (list): List containing filenames of results
      warn_on_empty (bool, optional): Log warning if result contained no
        result records (default: False)
      yield_filename (bool, optional): Yield the filename and the record

    Yields:
      np.ndarray or tuple: Result saved in record and the filename, if desired

    """
    n_records = len(records)

    for _i, r in enumerate(records):
        # Verbose progress
        if np.mod(_i, 100) == 0:
            logger.debug('{0:.1f}%'.format(_i / n_records * 100))
        z = np.load(r)
        rec = z['record']
        if rec.shape[0] == 0:
            # No values in this file
            if warn_on_empty:
                logger.warning('Could not find results in {f}'.format(f=r))
            continue

        if yield_filename:
            yield rec, r
        else:
            yield z, rec

def write_output(raster, output, image_ds, gdal_frmt, ndv, band_names=None):
    """ Write raster to output file """
    from osgeo import gdal, gdal_array

    logger.debug('Writing output to disk')

    driver = gdal.GetDriverByName(str(gdal_frmt))

    if len(raster.shape) > 2:
        nband = raster.shape[2]
    else:
        nband = 1

    ds = driver.Create(
        output,
        image_ds.RasterXSize, image_ds.RasterYSize, nband,
        gdal_array.NumericTypeCodeToGDALTypeCode(raster.dtype.type)
    )

    if band_names is not None:
        if len(band_names) != nband:
            logger.error('Did not get enough names for all bands')
            sys.exit(1)

    if raster.ndim > 2:
        for b in range(nband):
            logger.debug('    writing band {b}'.format(b=b + 1))
            ds.GetRasterBand(b + 1).WriteArray(raster[:, :, b])
            ds.GetRasterBand(b + 1).SetNoDataValue(ndv)

            if band_names is not None:
                ds.GetRasterBand(b + 1).SetDescription(band_names[b])
                ds.GetRasterBand(b + 1).SetMetadata({
                    'band_{i}'.format(i=b + 1): band_names[b]
                })
    else:
        logger.debug('    writing band')
        ds.GetRasterBand(1).WriteArray(raster)
        ds.GetRasterBand(1).SetNoDataValue(ndv)

        if band_names is not None:
            ds.GetRasterBand(1).SetDescription(band_names[0])
            ds.GetRasterBand(1).SetMetadata({'band_1': band_names[0]})

    ds.SetProjection(image_ds.GetProjection())
    ds.SetGeoTransform(image_ds.GetGeoTransform())

    ds = None


