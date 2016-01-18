#! /usr/bin/python
"""

Usage: runVBW.py -k -s layer_lines(simulated)

"""
__author__ = "Wojtek Potrzebowski"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "Wojciech.Potrzebowski@biochemistry.lu.se"

import optparse
import os
import sys
import vbw

if __name__=="__main__":
	doc = """
	    Python interface to Variational Bayesian algorithm
	    Usage: python runVBW.py --help
	"""
	print doc
	usage = "usage: %prog [options] args"
	option_parser_class = optparse.OptionParser
	parser = option_parser_class( usage = usage, version='0.1' )

	parser.add_option("-r", "--restart", dest="restart",default = 0,
                      type = 'int',
                      help="Restart or not")
	parser.add_option("-n", "--number_of_strcutures", dest="nstruct",default = None,
                      type = 'int',
                      help="Number of strcutural model [OBLIGATORY]")
	parser.add_option("-p", "--priors", dest="priors",
                      help="Prior weights [OBLIGATORY]")
        parser.add_option("-k", "--number_of_measurements", dest="measures",default = None,
                      type = 'int',
                      help="Number of exp measurements [OBLIGATORY]")
        parser.add_option("-s", "--simulated", dest="simulated",
                      help="Simulated SAXS curves [OBLIGATORY]")
        parser.add_option("-e", "--experimental", dest="experimental",
                      help="Experimental SAXS curves [OBLIGATORY]")
        parser.add_option("-d", "--delta", dest="errors",
                      help="Experimental SAXS errors [OBLIGATORY]")
        parser.add_option("-o", "--output", dest="output",
                      help="Output file [OBLIGATORY]")
        parser.add_option("-c", "--cores", dest="nprocs",default = None,
                      type = 'int',
                      help="Number of proccessors [OBLIGATORY]")
	parser.add_option("-w", "--weights", dest="weight_cut",default = None,
                      type = 'float',
                      help="Weight cutoff [OBLIGATORY]")
	options, args = parser.parse_args()
	vbw.run_vbw(options.restart, options.nstruct, options.priors,\
		options.measures, options.simulated, options.experimental,\
		options.errors, options.output, options.nprocs, options.weight_cut)
