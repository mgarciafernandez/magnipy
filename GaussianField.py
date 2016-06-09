import numpy, math, cmath, scipy.interpolate, scipy.integrate, argparse
from ROOT import *

class MatterDensityContrast(object):
	"""
	Object defining a matter cosmological density contrast.
	The field is just defined according to a gaussian random field. Not precise, but fast for testing.
	The algorithm is as follows:
		1.- Define a grid of NxN harmonic points.
		2.- Draw on each point a normal distributed random number,
		3.- Do a 2D-FFT to these random numbers.
		4.- Apply the power-spectrum as a filter function:
 		That is for each mode with k =  ||kx**2+ky**2|| multiply by the amplitude of Pk(k)
		5.- Do the 2D-FFT inverse transform.
	"""

	def __init__(self,cambfile='',N=500,L=1.,seed=626):

		kPk   = numpy.loadtxt(cambfile)
		k,Pk  = zip(*kPk)
		
		kPkf = scipy.interpolate.interp1d(k,Pk,bounds_error=False,fill_value=0.)
		kx   = map(abs,numpy.fft.fftfreq(N))
		ky   = map(abs,numpy,fft.fftfreq(N))
	
		rnd = TRandom3(seed)
		delta_R = [ [ rnd.Gaus(0,1) for kx_ in kx ] for ky_ in ky ]
		delta_F = numpy.fft2(delta_R)
	
		A_F = [ [ kPkf(math.sqrt(kx_**2+ky_**2)) for kx_ in kx ] for ky_ in ky ]
		
		for i_ in xrange(N):
			for j_ in xrange(N):
				delta_F[i_][j_] *= A_F[i_][j_]/sum(A_F[0])
	
		delta_M = numpy.fft.ifft2(delta_F)

		self.N  = N
		self.L  = L
		self.k  = k
		self.Pk = Pk
		self.delta_M = delta_M
	
	def GetGalaxies(self,Ngal=100000,seed=626):

		galaxies = []
		density_map = TH2F('density_map','',self.N,0,self.L,self.N,0,self.L)
		for i_ in xrange(self.N):
			for j_ in xrange(self.N):
				density_map.SetBinContent(i_+1,j_+1,self.delta_M.real)

		rnd = TRandom3(seed)
		niter = 0
		while niter < Ngal :
			ra  = rnd.Uniform(0,self.L)
			dec = rnd.Uniform(0,self.L)
			
			ii = density_map.GetXaxis().FindBin(ra)
			jj = density_map.GetYaxis().FindBin(dec)

			P = density_map.GetBinContent(ii,jj)
			p = rnd.Uniform(0,1)
			if p < P :
				galaxies.append( (ra,dec) )
				niter += 1
		return galaxies

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("-sm","Set seed of the random number generator of the gaussian random field.",type=int)
	parser.add_argument("-sg","Set seed of the random number generator of the galaxy coordinates.",type=int)
	parser.add_argument("-N","Set the number of modes of the harmonic space to test.",type=int)
	parser.add_argument("-L","Set the length of the squared box.",type=float)
	parser.add_argument("-pk","File path to CAMB.",type=str)
	parser.add_argument("--Ngal","Number of galaxies",type=int)
	parser.add_argument("-o","Outfile",type=str)
	args = parser.parse_args()

	matter_field = MatterDensityContrast(cambfile=args.pk,N=args.N,L=args.L,seed=args.sm)
	galaxies = matter_field.GetGalaxies(Ngal=args.Ngal,seed=args.sg)

	with open(args.o) as f:
		for ra_,dec_ in galaxies:
			f.write( str(ra)+' '+str(dec)+'\n' )
		f.close()

