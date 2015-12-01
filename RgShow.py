#! /usr/bin/python
"""
Reads radius of gyration data from crysol file and plots histograms
"""
__author__ = "Wojtek Potrzebowski"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "Wojciech.Potrzebowski@biochemistry.lu.se"

import optparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import string
import pylab as P


def read_rg_data(filename):
	f = open(filename)
	lines = f.readlines()
	Rg_Model = {}
	for line in lines:
		line_arr=string.split(line)
		model_name = line_arr[1]
		rg = float(line_arr[3])
		Rg_Model[model_name] = rg
	return Rg_Model


def show_rg(Rg_all, Rg_selected):
	"""
	Makes plots from layer lines and store them in image file
	"""
	rg_all = Rg_all.values()
	rg_selected = Rg_selected.values()

	# the histogram of the data with histtype='step'
	n1, bins1, patches1 = P.hist(rg_all, 50, normed=1, histtype='stepfilled')
	P.setp(patches1, 'facecolor', 'r', 'alpha', 0.75)
	n, bins, patches = P.hist(rg_selected, 50, normed=1, histtype='stepfilled')
	P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

	#l = P.plot(bins, 'k--', linewidth=1.5)
	P.show()



if __name__=="__main__":
	doc = """
	    Reads radius of gyration data from crysol file and plots histograms
	    Usage: python RgShow.py --help
	"""
	print doc
	usage = "usage: %prog [options] args"
	option_parser_class = optparse.OptionParser
	parser = option_parser_class( usage = usage, version='0.1' )

	parser.add_option("-a", "--all", dest="all",
                      help="Rgs for entire pool [OBLIGATORY]")
	parser.add_option("-s", "--selectec", dest="selected",
                      help="Rgs for selected models [OBLIGATORY]")
	options, args = parser.parse_args()
	all = options.all
	selected = options.selected

	All = read_rg_data(all)
	Selected = read_rg_data(selected)

	show_rg(All, Selected)
