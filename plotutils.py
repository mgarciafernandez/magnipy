from astropy.io import fits
from matplotlib import pyplot, use
import numpy


class TH1F(object):

	def __init__(self,name='',title='',nbins=10,xlow=None,xup=None):
		self.ndata = 0
		self.nbins = nbins
		self.xlow = xlow
		self.xup  = xup
		self.BinContent = []
		self.BinEdges   = []

	def Fill(self,data):
		self.BinContent,self.BinEdges = numpy.histogram(data,bins=self.nbins,range=(self.xlow,self.xup))
		self.ndata = len(data)

	def Integral(self,xlow=None,xup=None,width=True):
		if xlow is None:
			xlow = self.xlow
		if xup is None:
			xup = self.xup
		if width is not True:
			integral = sum(self.BinContent)
		else:
			integral = numpy.trapz(self.BinContent,self.BinEdges)
		return integral

	def Scale(self,c=1.0):
		for content_ in self.BinContent:
			content_ /= c

	def Draw(self,filename=None):
		if filename is None:
			pyplot.plot(self.BinEdges,self.BinContent)
		else:
			use('agg')
			pyplot.savefig(filename,filename.split('.')[-1])
		
