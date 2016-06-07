#!/usr/bin/env python

from astropy.io import fits
import healpy, numpy, sys


balrog = fits.open(sys.argv[1])[1].data
mask   = numpy.zeros(healpy.nside2npix(4096))

for i_ in xrange(len(balrog)):
	ra  = balrog['ra'][i_]
	dec = balrog['dec'][i_]

	mag     = [ balrog['mag_auto_g'][i_], balrog['mag_auto_r'][i_], balrog['mag_auto_i'][i_], balrog['mag_auto_z'][i_] ]
	flux    = [ balrog['flux_auto_g'][i_], balrog['flux_auto_r'][i_], balrog['flux_auto_i'][i_], balrog['flux_auto_z'][i_] ]
	fluxerr = [ balrog['fluxerr_auto_g'][i_], balrog['fluxerr_auto_r'][i_], balrog['fluxerr_auto_i'][i_], balrog['fluxerr_auto_z'][i_] ]

	if not balrog['spread_model_i'] + 3*balrog['spread_model_i'] < 0.003:
		continue


	if all( map(lambda flux_,fluxerr_:flux_>2*fluxerr_,flux,fluxerr) ):
		if ra > 180.:
			ra -= 360.
		tht = (90.-dec)*numpy.pi/180.
		phi = ra*numpy.pi/180.

		pix = healpy.ang2pix(4096,tht,phi,nest=False)

		if mag[0] > mask[pix]:
			mask[pix] = mag[0]

healpy.write_map(sys.argv[2],mask)
