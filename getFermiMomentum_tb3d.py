#!/usr/bin/python

import argparse
import os
import scipy
import numpy
import pylab
from scipy import special
from scipy import integrate
from scipy import optimize
import math

def main():
	parser = argparse.ArgumentParser(description='Get Friedel wavelength for given fillings of a 3d-tight-binding-band')
	parser.add_argument('-c', '--filling', default=0.01, help='Band filling', type=float)
	parser.add_argument('-p', '--plot', action='store_true', help='Band filling')
	parser.add_argument('-o', '--output', default="output/", help='Output folder')
	args = parser.parse_args()
	n=args.filling


	# integrand for dos (analytical formula with Besselfunction)
	def integrand(x, w, eta, t):
		return 2*scipy.imag(complex(0,1)*scipy.exp(complex(0,1)*x*complex(w, eta))*pow(special.jv(0,2*t*x), 3))/math.pi

	# get dos
	def occ(epsilon, eta, t):
		result, err=integrate.dblquad(integrand, -1, epsilon, lambda x:0, lambda x:1E2, args=(eta, t))
		return result

	# small imaginary part
	eta=1E-7
	# hopping
	t=1.0/6.0

	if (args.plot):
		if not os.path.exists(args.output):
			os.makedirs(args.output)
		e_list=numpy.arange(-1,1,0.1)
		n_list=[]
		for e in e_list:
			o=occ(e, eta, t)
			print e, o
			n_list.append(o)
	
		# save occupation number vs. energy to file
		e_list_shifted=e_list+1.0
		data=zip(e_list_shifted, n_list)
		occ_filename=args.output + "/occVsEnergy.dat"
		numpy.savetxt(occ_filename, data, fmt="%0.17e")
	
		pylab.plot(e_list, n_list, 'b')
		pylab.show()

	# find energy to given filling
	def root_func(e, n, eta, t):
		return occ(e, eta, t) - n

	epsilon0=optimize.brentq(root_func, -1, 1, args=(n,eta,t))
	# shift to positive region
	epsilon0=epsilon0+1.0
	print "Filling = ", n
	print "Energy= ", epsilon0 
	kf=math.acos((1.0-epsilon0)/(2.0/6.0)-2)
	print "Momentum kf (1/a)= ", kf
	l=2.0*math.pi/(2.0*kf)
	print "Friedel wavelength lambda=2*pi/(2*kf) (a)= ", l
	
if __name__=="__main__":
	main()
