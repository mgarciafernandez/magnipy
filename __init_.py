"""
Package:
	magnipy

Modules:
	catutils:
		Functions to handle catalogs and masks.
		Allows to do txt-fits conversion and define masks as well as mask data.
	magnipy:
		Statistical functions to compare with theory.
	jsonw:
		Json interface to store 2pacf.
"""

__version__ = "0.1"

def version():
	return __version__

import catutils
import magnipy
import jsonw
import plotutils
