# MAGNIPY

Cosmic Magnification Tools on Python.
Author: Manuel Garcia-Fernandez <manuel.garcia-fernandez@ciemat.es>
Date: 02/14/2016

## File Description

The file `catutils.py` contains some function utilities to handle catalogs and masks.
The file `jsonw.py` contains a class description of the 2pacf and implements a json file way of save.

Inside each function there should be enough explanation.

## Usefull Command-Line Instructions.

The function has been designed to be used as API. Nevertheless, their generality allow to perform simple catalog-handling tasks in a single line command. Here we provide some examples (assuming the library is located at '$PYTHONPATH')

* Convert database csv file to fits: `python -c "import catutils; catutils.TxtToFits('name.csv','name.fits')"`
* Append mask values to an existing fits file : `python -c "import catutils; catutils.Mask('filename.fits',['healpix#1.fits','healpix#2'],['namemask#1','namemask#2'])"`
