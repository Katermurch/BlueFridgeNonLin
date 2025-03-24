# Configuration file for the Sphinx documentation builder.
#

import os
import sys

sys.path.insert(0, os.path.abspath(".."))


# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "BlueFridgeNonLin"
copyright = "2025, Murch Lab"
author = "Murch Lab"
release = "0.1"

autodoc_mock_imports = [
    "hardware_control.atsapi",
    "hardware_control.tewx",
    "hardware_control.wx_programs",
    "usbtmc",
    "dg535_control",
    "AlazarTech",
    "ats",
    "pyvisa",
    "usbtmc",
    "atsapi",
    "dg535_control",
    # If you're using CDLL("libATSApi.so") under this name
]


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon"]


templates_path = ["_templates"]
exclude_patterns = []

language = "English"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
html_static_path = ["_static"]
