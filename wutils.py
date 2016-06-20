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

def ReweightKNN(array_to_match,array_to_reweight,keys,nn=100):
	"""
	Computes the weights of an given a table with N objects on an k-dim space with the KNN approach.
	The weights are computed such that the k-dim space of M objects of another reference table.
	-Input:
		array_to_match (structured array): The reference table.
		array_to_reweight (structured array): The array for which the weights are going to be computed.
		keys (list of str): the names of the fields of the array that will compose the k-dim space.
		nn (int): the number of nearest neighbors.
	-Output:
		w_norm (float array): the weights normalized at the range (0,1).
	"""

	array_to_match    = numpy.vstack([array_to_match[key_] for key_ in keys]).T
	array_to_reweight = numpy.vstack([array_to_reweight[key_] for key_ in keys]).T

	tree_match = sklearn.neighbors.NearestNeighbors(n_neighbors=nn+1).fit(array_to_match)
	tree_torew = sklearn.neighbors.NearestNeighbors(n_neighbors=nn+1).fit(array_to_reweight)

	dist_torew,indx_torew = tree_torew.kneighbors(array_to_reweight,n_neighbors=nn+1)

	which_nn = numpy.array( [nn]*len(dist_torew) )
	distance = dist_torew[numpy.arange(len(dist_torew)),which_nn]
	dist_match,indx_match = tree_match.radius_neighbors(array_to_reweight)

	ww = [float(len(d_)) for d_ in dist_match]
	w_norm = [ w_/sum(ww) for w_ in ww ]

	return w_norm

