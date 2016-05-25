from astropy.io import fits
import healpy, numpy, math


def WriteFitsTable(filename=None,colnames=[],types=[],*arg):

	"""
	Given a set of arrays/lists creates a fits file containing it.

	-Input:
		filename (str): The name of the file to write.
		colnames (list): An ordered list containing the names of the columns.
		types (list): An ordered list containg the type of each column. The types must be FITSIO tipes.
		*arg: The next arguments are each one a list containing the data.
	-Output:
		tbhdu (BinTableHDU): The BinTable resulting from this construction.

	"""

	if len(colnames) == 0:
		raise ValueError('No colnames given.')
	if len(types) == 0:
		raise ValueError('No types given.')
	if not len(colnames) == len(arg):
		raise Exception('The number of data-arrays is different than the names.')

	for data_ in arg:
		if not data_ == arg[0]:
			raise Exception('The size of the data-arrays is different.')

	columnlist = map(lambda name_,format_,array_: fits.Column( name=name_,format=format_,array=array_ ),colnames,types,arg)

	cols  = fits.ColDefs(columnlist)
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.writeto(filename)

	return tbhdu

def GetMaskArray(filename=None,ra=[],dec=[],units='degrees'):

	"""
	Given a a set of positons at the sky, returns an ordered array with the value of the mask at that position.
	A recentering is performed to ensure the coordinates to be -180 deg < ra < 180 deg.

	-Input:
		filename (str): The path to the HEALPix mask.
		ra (array): An array of floats containing the ra of each point of the sky.
		dec (array): An array of floats containing the dec of each point of the sky.
		units (str): Units of the input ra and dec. Degrees, Arcmin or Radians.

	-Output:
		maskvalues (list): a list containing the value of the HEALPIx mask on each point. 
	"""

	if len(ra) == 0:
		raise ValueError('No ra given.')
	if len(dec) == 0:
		raise ValueError('No dec given.')
	if not len(ra) == len(dec):
		raise Exception('The size of the ra and dec arrays is different.')
	if not units in ['degrees','radians','arcmin']:
		raise ValueError('No recognized angular units.')

	if units == 'degrees':
		ra  = map(lambda ra_ : ra_ *numpy.pi/180.,ra )
		dec = map(lambda dec_: dec_*numpy.pi/180.,dec)

	elif units == 'arcsec':
		ra  = map(lambda ra_ : ra_ /60.,ra )
		dec = map(lambda dec_: dec_/60.,dec)
	
	ordering = fits.open(filename)[1].header['ordering']
	nside    = fits.open(filename)[1].header['naxis1']
	
	if ordering == 'NESTED':
		isnest = True:
	elif ordering == 'RING':
		isnest = False
	else:
		raise Exception('Wrong ordering value '+ordering)

	mask = healpy.read_map(filename,nest=isnest)

	maskvalues = []
	for ra_,dec_ in zip(ra,dec):
		tht = numpy.pi/2.-dec_
		phi = ra_
		pix = healpy.ang2pix( nside,tht,phi,nest=isnest )

		maskvalues.append( mask[pix] )

	return maskvalues	

