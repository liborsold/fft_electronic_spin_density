# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# INCLUDE package in path
import os
import sys
sys.path.insert(0, os.path.abspath('../../src/'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'fft_electronic_spin_density'
copyright = '2025, Libor Vojáček'
author = 'Libor Vojáček'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx_copybutton'
]

autosummary_generate = True  # This triggers autosummary to auto-create stub files

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

templates_path = ['_templates']
html_static_path = ['_static']
html_logo = "logo_with_text.png"
exclude_patterns = []

html_theme = 'sphinx_book_theme'
html_theme_options = {
    "repository_url": "https://github.com/liborsold/fft_electronic_spin_density",
    "use_source_button": False,
    "use_repository_button": True,
    'logo_only': True,
    "repository_branch": "master",
    "path_to_docs": "docs/source/",
}