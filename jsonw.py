import json

class MeasuredWTheta(object):

	def __init__(self,theta=[],DD=[],DR=[],RD=[],RR=[],Nd1=[],Nd2=[],Nr1=[],Nr2=[]):

		if not len(theta) == len(DD):
			raise ValueError('Lenth of inputs does not agree.')
		if not len(DR) == len(DD):
			raise ValueError('Lenth of inputs does not agree.')
		if not len(RD) == len(DD):
			raise ValueError('Lenth of inputs does not agree.')
		if not len(RR) == len(DD):
			raise ValueError('Lenth of inputs does not agree.')

		self.theta = theta
		self.DD    = DD
		self.DR    = DR
		self.RD    = RD
		self.RR    = RR
		self.Nd1   = Nd1
		self.Nd2   = Nd2
		self.Nr1   = Nr1
		self.Nr2   = Nr2
		self.w     = []

		for th_ in theta:
			w_  = 1.
			w_ += (DD/RR) * ((Nr1*Nr2)/(Nd1*Nd2)))
			w_ -= (DR/RR) * (Nr1/Nd1)
			w_ -= (RD/RR) * (Nr2/Nd2)

	def __len_(self):

		return len(self.w)

	def __str__(self):

		return str(zip(self.theta,self.w))

	def __str__(self):

		return str(zip(self.theta,self.w))

	def __getitem__(self,ii):

		return zip(self.theta,self.w)[ii]

	def __hash__(self):

		tohash = {}
		tohash['theta'] = list(self.theta)
		tohash['w']     = list(self.w)
		tohash['DD']    = list(self.DD)
		tohash['DR']    = list(self.DD)
		tohash['RD']    = list(self.DD)
		tohash['RR']    = list(self.DD)
		tohash['Nd1']    = list(self.Nd1)
		tohash['Nd2']    = list(self.Nd2)
		tohash['Nr1']    = list(self.Nr1)
		tohash['Nr2']    = list(self.Nr2)

		return tohash

	def SaveJson(self,filename):
		
		tohash = {}
		tohash['theta'] = list(self.theta)
		tohash['w']     = list(self.w)
		tohash['DD']    = list(self.DD)
		tohash['DR']    = list(self.DD)
		tohash['RD']    = list(self.DD)
		tohash['RR']    = list(self.DD)
		tohash['Nd1']    = list(self.Nd1)
		tohash['Nd2']    = list(self.Nd2)
		tohash['Nr1']    = list(self.Nr1)
		tohash['Nr2']    = list(self.Nr2)

		with open(filename) as data_file:
			json.dump(tohash,data_file)

