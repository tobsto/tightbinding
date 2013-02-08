#!/usr/bin/python

import argparse
from math import *
def main():
	parser = argparse.ArgumentParser(description='Get Friedel wavelength for given fillings of a 1d-tight-binding-band')
	parser.add_argument('-c', '--filling', default=0.01, help='Band filling', type=float)
	args = parser.parse_args()
	n=args.filling
	k=pi*n
	l=2*pi/(2*k)
	print "Filling = ", n
	print "Fermi momentum (1/a)= ", k
	print "Friedel wavelength (a)= ", l

	
if __name__=="__main__":
	main()
