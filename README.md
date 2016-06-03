# MAGNIPY

Cosmic Magnification Tools on Python 2.x.

Author: Manuel Garcia-Fernandez <manuel.garcia-fernandez@ciemat.es>

### Installation Procedure
No need to compile, nor perfom any installation.

Just copy the `/magnipy` directory to anywhere and add its path to `$PYHTONPATH`. At Bash:
```
export PYTHONPATH=${PYTHONPATH}:/somepath/magnipy
```
Besides the Python standard libraries the dependencies are:
* numpy
* scipy
* matplotlib
* json
* healpy
* astropy
* scikit-learn
* ROOT a.k.a. PyROOT, the wrap from CERN's libraries (https://root.cern.ch).

Currently trying to avoid as much as possible ROOT dependencies since this is not installed at all computers.

### File Description

The file `./catutils.py` contains some function utilities to handle catalogs and masks.
The file `./jsonw.py` contains a class description of the 2pacf and implements a json file way of save.

Inside each function there should be enough explanation.

### Useful Command-Line Instructions.

The function has been designed to be used as API. Nevertheless, their generality allow to perform simple catalog-handling tasks in a single line command. Here we provide some examples.

* Convert database csv file to fits:
```
python -c "import catutils; catutils.TxtToFits('name.csv','name.fits')"
```
* Append mask values to an existing fits file :
```
python -c "import catutils; catutils.Mask('filename.fits',['healpix#1.fits'],['namemask#1'])"
```
