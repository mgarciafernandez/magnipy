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

def ReweightKNN(array_to_match,array_to_reweight,keys,nn=100,ncpu=None):
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

	if ncpu is None:
		ncpu = __NCPU__-1
	elif not ncpu < __NCPU__:
		ncpu = __NCPU__-1

	array_to_match    = numpy.vstack([array_to_match[key_] for key_ in keys]).T
	array_to_reweight = numpy.vstack([array_to_reweight[key_] for key_ in keys]).T

	tree_torew = sklearn.neighbors.NearestNeighbors(n_neighbors=nn+1,n_jobs=ncpu).fit(array_to_reweight)
	dist_torew,indx_torew = tree_torew.kneighbors(array_to_reweight,n_neighbors=nn+1)

	which_nn = numpy.array( [nn]*len(dist_torew) )
	distance = dist_torew[numpy.arange(len(dist_torew)),which_nn]

	tree_match = sklearn.neighbors.NearestNeighbors(radius=distance,n_jobs=ncpu).fit(array_to_match)
	dist_match,indx_match = tree_match.radius_neighbors(array_to_reweight,radius=distance)

	ww = [float(len(d_)) for d_ in dist_match]
	w_norm = [ w_/sum(ww) for w_ in ww ]

	return w_norm

def Get2pacf(data_lens=[],data_sour=None,rand_lens=[],rand_sour=None):
	"""
	Computes the 2pacf density-density cross/auto-correlation.
	-Input:
		data_lens (list): a list containing [ra,dec,weight] for the lens. If weight is not provided, ti will be assumed as 1.
		data_sour (list): a list containing [ra,dec,weight] for the source. If weight is not provided, ti will be assumed as 1.
				  if None, the auto-correlation will be computed instead of the cross.
		rand_lens (list): a list containing [ra,dec,weight] for the random sample associated with the lens.
		rand_sour (list): a list containing [ra,dec,weight] for the random sample associated with the source.
	-Output:
		xi (list): a list (xi,sigma) such that xi contains the value of the 2pacf and sigma the Poisson error.
		th (list): a list containing the angle values of the bins where the 2pacf is computed.
		
	"""

	if len(data_lens) == 2:
		lens_cat = treecorr.Catalog(ra=data_lens[0],dec=data_lens[1],ra_units='degrees',dec_units='degrees')
	elif len(data_lens) == 3:
		lens_cat = treecorr.Catalog(ra=data_lens[0],dec=data_lens[1],w=data_lens[2],ra_units='degrees',dec_units='degrees')
	else:
		raise ValueError('Too many lens parameters.')
	if len(rand_lens) == 2:
		lens_cat = treecorr.Catalog(ra=rand_lens[0],dec=rand_lens[1],ra_units='degrees',dec_units='degrees')
	elif len(rand_lens) == 3:
		lens_rnd = treecorr.Catalog(ra=rand_lens[0],dec=rand_lens[1],w=rand_lens[2],ra_units='degrees',dec_units='degrees')
	else:
		raise ValueError('Too many lens parameters.')

	if data_sour is None and not rand_sour is None:
		raise Exception('Not data of sources provided')
	if rand_sour is None and not data_sour is None:
		raise Exception('Not rand of sources provided')

	if data_sour is None and rand_sour is None:
		data_sour = data_lens
		rand_sour = rand_lens

	if len(data_sour) == 2:
		sour_cat = treecorr.Catalog(ra=data_sour[0],dec=data_sour[1],ra_units='degrees',dec_units='degrees')
	elif len(data_lens) == 3:
		sour_cat = treecorr.Catalog(ra=data_sour[0],dec=data_sour[1],w=data_sour[2],ra_units='degrees',dec_units='degrees')
	else:
		raise ValueError('Too many lens parameters.')
	if len(rand_sour) == 2:
		sour_rnd = treecorr.Catalog(ra=rand_sour[0],dec=rand_sour[1],ra_units='degrees',dec_units='degrees')
	elif len(rand_sour) == 3:
		sour_rnd = treecorr.Catalog(ra=rand_sour[0],dec=rand_sour[1],w=rand_sour[2],ra_units='degrees',dec_units='degrees')
	else:
		raise ValueError('Too many lens parameters.')


	dd = treecorr.NNCorrelation(nbins=6,min_sep=0.01,max_sep=1.0,sep_units='degrees')
	dr = treecorr.NNCorrelation(nbins=6,min_sep=0.01,max_sep=1.0,sep_units='degrees')
	rd = treecorr.NNCorrelation(nbins=6,min_sep=0.01,max_sep=1.0,sep_units='degrees')
	rr = treecorr.NNCorrelation(nbins=6,min_sep=0.01,max_sep=1.0,sep_units='degrees')

	dd.process(lens_cat,sour_cat)
	dr.process(lens_cat,sour_rnd)
	rd.process(lens_rnd,sour_cat)
	rr.process(lens_rnd,sour_rnd)

	xi = dd.calculateXi(rr=rr,dr=dr,rd=rd)
	th = dd.meanr

	return xi,th
