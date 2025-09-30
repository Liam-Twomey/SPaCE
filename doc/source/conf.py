# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- imports for autofunctions
import sys
from pathlib import Path
#from importlib.metadata import version as vn
#import mock

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
#html_css_files = ["eprsim.css"]
html_show_sourcelink = False 
#html_favicon = "_static/eprsim_favicon.png"
#html_logo = "_static/eprsim_logo.png"
#html_theme_options = {
#	"logo": {
#		"image_light": "_static/eprsim_logo_light.svg",
#		"image_dark": "_static/eprsim_logo_dark.svg",
#	}
#}

# For web document autobuild with readthedocs.io
html_context = {
	"github_user": "zhasanbasri-hash",
	"github_repo" : "SPaCE",
	"github_version" : "main",
	"doc_path":"doc",
}

# -- Set path for autofunctions and autoset version info
#
pkg_path  = str(Path('../..').resolve())
print('package location:', pkg_path)
sys.path.insert(0,pkg_path) 
#ver = vn('eprsim')

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SPaCE'
copyright = '2025, Zikri Hasanbasri' 
author = 'Zikri Hasanbasri'
#release = '.'.join(ver.split('.'))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
	'sphinx.ext.napoleon', # for numpy in autodoc
	'sphinx.ext.mathjax',  # for latex support in html 
	'sphinx.ext.autodoc',  # for automatic documentation of docstrings
#	'autoapi',  # for latex support in html 
	'sphinx_copybutton',   # allow code copy button in examples
	'sphinx_design'		   # enable design elements like grids
]
#autoapi_dirs = ['../..','.']
language = "en"
templates_path = ['_templates']

#toMock = ['Validate_input_parameter','Tools']
#for modName in toMock:
#	sys.modules[modName] = mock.Mock()

