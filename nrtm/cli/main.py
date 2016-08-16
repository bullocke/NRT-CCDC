""" Loads all commands for nrtm command line interface

Modeled after very nice `click` interface for YATSM, which was modeled after `rasterio`:
https://github.com/mapbox/rasterio/blob/master/rasterio/rio/main.py

"""
import logging
from pkg_resources import iter_entry_points

import click
import click_plugins
#import pdb; pdb.set_trace()
#from click_plugins import with_plugins

import nrtm


# nrtm CLI group
_context = dict(
    token_normalize_func=lambda x: x.lower(),
    help_option_names=['--help', '-h']
)


@click_plugins.with_plugins(ep for ep in
                            iter_entry_points('nrtm.nrtm_commands'))
@click.group(help='nrtm command line interface', context_settings=_context)
@click.option('--verbose', '-v', is_flag=True, help='Be verbose')
@click.option('--quiet', '-q', is_flag=True, help='Be quiet')
@click.pass_context
def cli(ctx, verbose, quiet):
    # Logging config
    logger = logging.getLogger('yatsm')
    if verbose:
        logger.setLevel(logging.DEBUG)
    if quiet:
        logger.setLevel(logging.WARNING)
