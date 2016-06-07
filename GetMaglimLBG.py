#!/usr/bin/env python

from astropy.io import fits
import healpy, numpy, argparse

if __name__ == "__main__":	

	parser = argparse.ArgumentParser()
	parser.add_argument("--sim","-s",help="Path to Balrog simulated file.",type=str)
	parser.add_argument("--out","-o",help="Name of the mask to write.",type=str)
	parser.add_argument("--band","-b",help="Band to get the magnitude limit",type=str)
	args = parser.parse_args()

	balrog = fits.open(args.sim)[1].data

	for band in ['g','r','i','z']:
		balrog = balrog[ balrog['flux_auto_'+band]>2*balrog['fluxerr_auto_'+band] ]
	balrog = balrog[ balrog['spread_model_i']+2*balrog['spreaderr_model_i']<0.003 ]

	mask   = numpy.zeros(healpy.nside2npix(4096))

	for i_ in xrange(len(balrog)):
		ra  = balrog['ra'][i_]
		dec = balrog['dec'][i_]

		if ra > 180.:
			ra -= 360.
		tht = (90.-dec)*numpy.pi/180.
		phi = ra*numpy.pi/180.

		pix = healpy.ang2pix(4096,tht,phi,nest=False)

		if balrog['mag_auto_'+args.band][i_] > mask[pix]:
			mask[pix] = balrog['mag_auto_'+args.band][i_]

	healpy.write_map(args.out,mask)
