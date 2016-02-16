import numpy, math, json, os, ROOT

__GRAPHIX__ = 'ROOT'


class CorrelationFunction(object):

	def __init__(self,name=''):
		self.angle_ = numpy.zeros([0])
		self.w_     = numpy.zeros([0])
		self.error_ = numpy.zeros([0])
		self.Nth_   = len(self.w_)
		self.name_  = name
		self.plot_  = None

		if __GRAPHIX__ == 'ROOT':
			self.plot_ = ROOT.TGraphErrors()

	def __str__(self):
		return str(zip(self.angle_,self.w_,self.error_))

	def __repr__(self):
		return str(zip(self.angle_,self.w_,self.error_))

	def __len__(self):
		return self.Nth_

	def __getitem__(self,ii):
		return (self.angle_[ii],self.w_[ii],self.error_[ii])

	def __hash__(self):

		data = {}
		data['angle'] = self.angle_
		data['w']     = self.w_
		data['error'] = self.error_
		data['name']  = self.name_
		data['Nth']   = self.Nth_

		return 0

	def Plot(self,opts=''):
		if self.plot_ is None:
			print 'Graphic library not defined.'
			raise Exception

		elif self.Nth_ == 0:
			print 'Empty object. No points to draw.'
			raise Exception

		elif not (self.Nth_ == len(self.w_) == len(self.angle_)):
			print 'Dimension does not agree. Chech angle and w.'
			raise Exception

		elif __GRAPHIX__ == 'ROOT' :
			for th in xrange(self.Nth_):
				self.plot_.SetPoint(th,self.angle_[th],self.w_[th])
				self.plot_.SetPointError(th,0.,self.error_[th])

class DataW(CorrelationFunction):

	def __init__(self,name=''):
		CorrelationFunction.__init__(self,name)
		self.covariance_     = numpy.zeros([0,0])
		self.pathWtheta_     = ''
		self.pathCovariance_ = ''

	def ReadAthenaFunction(self,path):

		tmp = numpy.loadtxt(path,usecols=[0])
		if self.pathCovariance_ != '' and tmp.size != self.Nth_ :
			print 'Number of points does not agree with previous covariance file: ',self.pathCovariance_
			raise Exception
		else:
			self.angle_  = numpy.loadtxt(path,usecols=[0])
			self.w_      = numpy.loadtxt(path,usecols=[1])
			self.error_  = numpy.loadtxt(path,usecols=[2])
			self.Nth_    = len(self.w_)

			if path[0] != '/':
				if path[0] == '.':
					self.pathW_ = os.getcwd()+'/'+path[1:]
				else:
					self.pathW_ = os.getcwd()+'/'+path[:]
			else:
				self.pathW_ = path[:]

	def SetDiagonalCovariance(self):

		self.covariance_ = numpy.zeros([self.Nth_,self.Nth_])
		for th in xrange(self.Nth_):
			self.covariance_[th][th] = self.error_[th]**2
		self.pathCovariance_ = ''


	def ReadAthenaCovariance(self,path):
		tmp = numpy.loadtxt(path)
		if self.pathWtheta_ != '' and tmp.size != self.Nth_**2 :
			print 'Number of points does not agree with previous wtheta file: ',self.pathWtheta_
			raise Exception
		else:
			self.covariance_ = numpy.loadtxt(path)
			self.Nth_        = math.sqrt(self.covariance_.size)

			if path[0] != '/':
				if path[0] == '.':
					self.pathCovariance_ = os.getcwd()+'/'+path[1:]
				else:
					self.pathCovariance_ = os.getcwd()+'/'+path[:]
			else:
				self.pathCovariance_ = path[:]


	def GetChi(self,w_theory):

		if self.Nth_ == 0:
			print 'Empty data object! Can not do fit.'
			raise Exception
		elif len(w_theory) != self.Nth_:
			print 'Length of input array ',len(w_theory),' does not agree with those of data ',self.Nth_,'.'
			raise Exception

		errorM = numpy.linalg.inv(self.covariance_)
		chisq = 0.
		for th1 in xrange(self.Nth_):
			for th2 in xrange(self.Nth_):
				chisq += ( self.w_[th1]-w_theory[th1] )*( self.w_[th2]-w_theory[th2] )*errorM[th1][th2]
		return chisq

class TheoMagW(CorrelationFunction):

	def __init__(self,name='',bias=1.,alpha=1.):
		CorrelationFunction.__init__(self,name)
		self.w0_         = numpy.zeros([0])
		self.alpha_      = alpha
		self.bias_       = bias
		self.pathWtheta_ = ''

	def ReadFunction(self,path):
		self.angle_ = numpy.loadtxt(path,usecols=[0])
		self.w0_    = numpy.loadtxt(path,usecols=[1])
		self.Nth_   = len(self.w0_)
		self.w_     = map(lambda x: x*self.alpha_*self.bias_,self.w0_)
		self.error_ = numpy.zeros([self.Nth_])
		
		if path[0] != '/':
			if path[0] == '.':
				self.pathWtheta_ = os.getcwd()+'/'+path[1:]
			else:
				self.pathWtheta_ = os.getcwd()+'/'+path[:]
		else:
			self.pathWtheta_ = path[:]

	def Update(self):
		self.w_ = map(lambda x: x*self.alpha_*self.bias_,self.w0_)
