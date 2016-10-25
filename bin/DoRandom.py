#!/usr/bin/python

import sys
import random
import numpy
import math

if __name__ == "__main__":
	_filename_ = sys.argv[1]

	wf = open(_filename_,'w')
	wf.write('ra,dec\n')
	ii_= 0
	while ii_ < 1000000000:
		ra_  = random.uniform(0,90)
		cth_ = random.uniform(0,1)
		dec_ = 90.-math.acos(cth_)*180./numpy.pi

		wf.write( str(ra_)+','+str(dec_)+'\n')

		ii_ += 1
	wf.close()
