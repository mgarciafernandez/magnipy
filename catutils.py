from astropy.io import fits
import healpy
import numpy
import math
import multiprocessing
import sklearn.neighbors
import collections


__NCPU__ = multiprocessing.cpu_count()

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
		if not len(data_) == len(arg[0]):
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
		units (str): Units of the input ra and dec. Degrees, arcmin or rad

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
	nside    = fits.open(filename)[1].header['nside']
	
	if ordering == 'NESTED':
		isnest = True
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

def Mask(filecat=None,filemask=[],maskname=[]):
	"""
	Given a Fits file, appends the value of the mask to the table fits on a new file 
	-Input:
		filecat (str): The name of the file with the catalog.
		filemask (list): List of the names of the file with the mask.
		maskname (list): List of the names of the new field
	"""

	if len(filemask) == 0:
		raise ValueError('No mask given.')
	if not len(filemask) == len(maskname):
		raise ValueError('The number of files and headers does not match.')

	catalog       = fits.open(filecat)[1].data
	columnnames   = fits.open(filecat)[1].columns.names
	columnformats = fits.open(filecat)[1].columns.formats

	masklist = []
	for file_ in filemask:
		masklist.append( GetMaskArray(file_,catalog['ra'],catalog['dec']) )


	columns  = [ catalog[col_] for col_ in columnnames ]
	columns       += masklist
	columnnames   += maskname
	columnformats += [ 'E' for name_ in maskname ]

	columnlist = map(lambda name_,format_,array_: fits.Column( name=name_,format=format_,array=array_ ),columnnames,columnformats,columns)

	cols  = fits.ColDefs(columnlist)
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.writeto(filecat+'_'.join(maskname))


def GetMasterMask(logic=None,*arg):
	"""
	Builds a binary mask such that 0 means not in the mask and 1 means in the mask.
	The mask is build as a HEALPix map form other HEALPix masks.

	-Input:
		logic (str): A string cotaining the logic to apply to the subsequent masks.
			The logic is as follows: xmin1,xmax1;xmin2,xmax2
			*The conditions of different masks must be separated by ;
			*The condition xmin,xmax is equivalent to xmin < x < xmax
			*The condition xmin~,xmax~ is equivalent to xmin <= x <= xmax
			*If no bound on a side that side must put inf, that is xmin,inf is equivalent to xmin < x
			*To demand to be equal to a quantity xmin~,xmin~ is the same as xmin = x
		arg (array): a variable number of entries, each one a HEALPix mask
	-Output:
		mask (array): the master mask
	"""

	conditions = logic.split(';')
	for condition_ in conditions:
		if not len(condition_.split(',')) == 2:
			raise SyntaxError('Condition syntax incorrect.')
	if not len(conditions) == len(arg):
		raise ValueError('Not given the same number of conditions as masks.')
	for mask_ in arg:
		if not len(mask_) == len(arg[0]):
			raise Exception('Lengths of the mask does not match.')

	mask = numpy.ones(len(arg[0]))

	for condition_,mask_ in zip(conditions,arg):
		xmin,xmax = map(lambda eq:eq.split('~'),condition_.split(','))
		if len(xmin) == 2 and len(xmin) == 2:
			for pix_ in xrange(len(mask_)):
				if not float(xmin[0]) <= mask_[pix_] <= float(xmax[0]):
					mask[pix_] *= 0

		elif len(xmin) == 2 and len(xmax) == 1:
			for pix_ in xrange(len(mask_)):
				if not float(xmin[0]) <= mask_[pix_] < float(xmin[0]):
					mask[pix_] *= 0
			
		elif len(xmin) == 1 and len(xmax) == 2:
			for pix_ in xrange(len(mask_)):
				if not float(xmin[0]) < mask_[pix_] <= float(xmin[0]):
					mask[pix_] *= 0

		elif len(xmin) == 1 and len(xmax) == 1:
			for pix_ in xrange(len(mask_)):
				if not float(xmin[0]) < mask_[pix_] < float(xmin[0]):
					mask[pix_] *= 0

	return mask

def TxtToFits(filein=None,fileout=None):
	"""
	Converts .csv (as read from desdb with E. Sheldon library) to a FITS table.

	-Input:
		filein (str): The name of the file to read
		fileout(str): The name of the file to write
	"""

	with open(filein) as csv:
		line = csv.readline()
	colnames = [ line_.translate(None,'\n') for line_ in line.split(',') ]

	duplicated = [item for item, count in collections.Counter(colnames).items() if count > 1]
	if not len(duplicated) is 0:
		raise Exception('There are columns repeated: '+','.join(duplicated) )

	table = [ [] for _ in colnames ]
	types = ['E' for _ in colnames ]

	with open(filein) as csv:
		line = csv.readline()
		line = csv.readline()

		while not line == '':
			cols = line.split(',')
			line = csv.readline()

			for col_ in range(len(cols)):
				try:
					if cols[col_] is 'ra' and float(cols[col_]) > 180.:
						table[col_].append( float(cols[col_])-360. )
					else:
						table[col_].append( float(cols[col_]) )
				except ValueError:
					table[col_].append( cols[col_] )
					types[col_] = '100A'
		csv.close()


	columnlist = map(lambda name_,format_,array_: fits.Column( name=name_,format=format_,array=array_ ),colnames,types,table)
	
	cols  = fits.ColDefs(columnlist)
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.writeto(fileout)

	return tbhdu

def FitsToTxT(filein=None,fileout=None):
	"""
	Converts to .csv (as read from desdb with E. Sheldon library) from a FITS table.

	-Input:
		fileout (str): The name of the file to read
		filein(str): The name of the file to write
	"""

	catalog  = fits.open(filein)[1].data
	colnames = fits.open(filein)[1].columns.names

	with open(fileout,'w') as write_file:
		write_file.write( ','.join(colnames)+'\n' )

		for i_ in xrange(len(catalog)):
			line = []
			for col_ in colnames:
				line.append( str(catalog[col_][i_]) )
			write_file.write( ','.join( line )+'\n' )

		write_file.close()

def Recenter(filename=None,colname=[]):
	"""
	Given a catalog, recenters a set of columns from the 0 < ra < 360 to -180 < ra < 180.
	-Input:
		filename (str): the file to recenter.
		colname (list): a list containing the colnames to recenter.
	"""

	if len(colname) is 0:
		raise Exception('WARNING: no colname given to recenter')

	catalog  = fits.open(filename)[1].data
	colnames = fits.open(filename)[1].columns.names
	types    = fits.open(filename)[1].columns.formats


	newcol = {}
	for col_ in colname:
		newcol[col_] = []
		for ra_ in catalog[col_]:
			if ra_ > 180.:
				newcol[col_].append( ra_ - 360. )
			else:
				newcol[col_].append( ra_ )

	table = []
	for i_ in range(len(colnames)):
		if colnames[i_] in colname:
			table.append( newcol[colnames[i_]] )
		else:
			table.append( catalog[colnames[i_]] )

	columnlist = map(lambda name_,format_,array_: fits.Column( name=name_,format=format_,array=array_ ),colnames,types,table)
	
	cols  = fits.ColDefs(columnlist)
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.writeto(filename+'_centered')

def RenameColumn(filein,keysin=[],keysout=[]):
	"""
	Renames the columns of a bintable on a fits.
	-Input:
		filein (str): The name of the file containing the bintable.
		keysin (str list): A list containing the names to rename.
		keysoyt (str list): A list containing the new names on the same order.
	
	"""
	data          = fits.open(filein)[1].data
	columnnames   = fits.open(filein)[1].columns.names
	columnformats = fits.open(filein)[1].columns.formats

	if not len(keysin) == len(keysout):
		raise Exception('The length of the keys lists are different.')
	for key_ in keysin:
		if not key_ in columnnames:
			raise ValueError('The key '+key_+' does not belong to the bintable at '+filein)


	newcolnames = []
	newdata     = []
	for col_ in columnnames:
		if col_ in keysin:
			newcolnames.append( keysout[keysin.index(col_)] )
		else:
			newcolnames.append( col_ )
		newdata.append( data[col_] )

	columnlist = map(lambda name_,format_,array_: fits.Column( name=name_,format=format_,array=array_ ),newcolnames,columnformats,newdata)

	cols  = fits.ColDefs(columnlist)
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.writeto(filein+'_renamed')

	return tbhdu
