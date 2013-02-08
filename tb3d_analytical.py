#!/usr/bin/python

import argparse
import scipy
import numpy
import pylab
from scipy import special
from scipy import integrate
import math

def main():
	parser = argparse.ArgumentParser(description='Get 3d-tight-binding density of states by analytical formula')
	parser.add_argument('-o', '--output', default="output/", help='Output folder')
	args = parser.parse_args()

	# small imaginary part
	eta=1E-7
	# hopping
	t=1.0/6.0

	def integrand(x, w, eta, t):
		return scipy.imag(complex(0,1)*scipy.exp(complex(0,1)*x*complex(w, eta))*pow(special.jv(0,2*t*x), 3))/math.pi
	w_list=numpy.arange(-1.0, 1.0, 0.01)
	dos_list=[]
	for w in w_list:
		print w
		dos, err=integrate.quadrature(integrand, 0, 1E2, args=(w,eta, t), maxiter=100)
		dos_list.append(dos)

		# print integrand for
		#xlist=numpy.arange(0.0,100,1)
		#ylist=[]
		#for x in xlist:
		#	ylist.append(integrand(x, w, eta, t))
		#pylab.plot(xlist, ylist)
	#pylab.show()


	# save dos to file
	if not os.path.exists(args.output):
		os.makedirs(args.output)
	w_list_shifted=w_list+1.0
	data=zip(w_list_shifted, dos_list)
	dos_filename=args.output + "/dos.dat"
	numpy.savetxt(dos_filename, data, fmt="%0.17e")

	# print dos
	pylab.xlabel("omega")
	pylab.ylabel("DoS")
	pylab.plot(w_list, dos_list, 'b')
	pylab.show()
	
if __name__=="__main__":
	main()
