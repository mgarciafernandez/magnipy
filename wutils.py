from astropy.io import fits
import healpy, numpy, math, multiprocessing, sklearn, collections
import catutils
from ROOT import *

__NCPU__ = multiprocessing.cpu_count()

def DoRandomFlat(N=1,ra=[],dec=[],masks=[],masknames=[]):
	"""
	Creates an BinTable with uniform distributed points trough space with mask values appended.
	-Input:
		N (int): The number of points.
		ra (list): The ra edges of the box where do draw the points.
		dec (list): The dec edges of the box where do draw the points.
		masks (list): The list of the masks to append.
	-Output:
		rand (bintable HDU): The catalog.
	"""

	rnd = TRandom3(0)

	ra_r  = []
	dec_r = []
	while len(ra_r) < N:
		ra_tmp  = rnd.Uniform(ra[0],ra[1])
		cth_tmp = rnd.Uniform(math.cos((90.-dec[0])*numpy.pi/180.),math.cos((90.-dec[1])*numpy.pi/180.))
		dec_tmp = 90.-math.acos(cth_tmp)*180./numpy.pi

		ra_r.append( ra_tmp )
		dec_r.append( dec_tmp )

	mask_array = []
	for filemask_ in masks:
		mask_array.append( catutils.GetMaskArray(filemask_,ra=ra_r,dec=dec_r) )	

	colnames = ['ra','dec']+masknames
	typles   = ['E']*len(colnames)
	arg      = [ra_r,dec_r]+mask_array

	columnlist = map(lambda name_,format_,array_: fits.Column( name=name_,format=format_,array=array_ ),colnames,types,arg)
	cols  = fits.ColDefs(columnlist)
	tbhdu = fits.BinTableHDU.from_columns(cols)
	return tbhdu
