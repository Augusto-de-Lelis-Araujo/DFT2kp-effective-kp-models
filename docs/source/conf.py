import sys
import os

from pydft2kp.__version import __version__ as version

project = 'dft2kp'
copyright = '2023, João V. V. Cassiano, Gerson J. Ferreira'
author = 'João V. V. Cassiano, Gerson J. Ferreira'
release = version

extensions = ["myst_parser",
              "sphinx.ext.duration",
              "sphinx_gallery.load_style",
              "sphinx.ext.autosectionlabel",
              "sphinx_rtd_theme",
              'sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.napoleon',
              'sphinx.ext.mathjax',
              'sphinx.ext.viewcode',
              'sphinx_multiversion'
              ]

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_sidebars = {'**': ['versions.html',]}


numpydoc_show_class_members = True 
numpydoc_class_members_toctree = True
numpydoc_show_inherited_class_members = True

autodoc_member_order = 'bysource'
