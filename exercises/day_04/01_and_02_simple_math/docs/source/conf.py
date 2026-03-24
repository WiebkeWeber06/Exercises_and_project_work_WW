# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# Configuration file for the Sphinx documentation builder.

from pathlib import Path
import sys

HERE = Path(__file__).resolve()
PROJECT_ROOT = HERE.parents[2]   # .../04_day_4
sys.path.insert(0, str(PROJECT_ROOT))

print("Sphinx conf.py loaded")
print("PROJECT_ROOT =", PROJECT_ROOT)
print("sys.path[0] =", sys.path[0])

# -- Project information -----------------------------------------------------

project = 'simple_math'
copyright = '2026, Wiebke Weber'
author = 'Wiebke Weber'
release = '1.0'

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "numpydoc",
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

html_theme = 'alabaster'
html_static_path = ['_static']
