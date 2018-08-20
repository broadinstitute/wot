# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

import os
import sys

sys.path.insert(0, os.path.join(os.path.abspath('..'), 'wot'))


# -- Project information -----------------------------------------------------

project = 'wot'
copyright = '2018, Broad Institute'
author = 'Broad Institute'

# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = ''

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinxcontrib.napoleon'
]

source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

html_theme = 'sphinx_rtd_theme'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
