# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- imports for autofunctions
import sys
from pathlib import Path
from importlib.metadata import version as vn

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_show_sourcelink = False 

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
ver = vn('SPaCE')

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SPaCE'
copyright = '2025, Zikri Hasanbasri' 
authors = 'Zikri Hasanbasri'
release = '.'.join(ver.split('.'))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
	'sphinx.ext.napoleon', # for numpy in autodoc
	'sphinx.ext.mathjax',  # for latex support in html 
	'sphinx.ext.autodoc',  # for automatic documentation of docstrings
	'sphinx_copybutton',   # allow code copy button in examples
	'sphinx_design'		   # enable design elements like grids
]
language = "en"
templates_path = ['_templates']
