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
import vbwSC

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
    parser.add_option("-m", "--measures", dest="measures",default = None,
                      type = 'int',
                      help="Number of measure on which algorithm is evaluated [OBLIGATORY]")
    parser.add_option("-v", "--number_of_curves", dest="ncurves",default = 1,
                      type = 'int',
                      help="Number of strcutural model [OBLIGATORY]")
    parser.add_option("-p", "--priors", dest="priors",
                      help="Prior weights [OBLIGATORY]")
    	parser.add_option("-s", "--simulated", dest="simulated",
                      help="Simulated SAXS curves [OBLIGATORY]")
    	parser.add_option("-e", "--experimental", dest="experimental",
                      help="Experimental SAXS curves [OBLIGATORY]")
    	parser.add_option("-o", "--output", dest="output",
                      help="Output file [OBLIGATORY]")
    	parser.add_option("-c", "--cores", dest="nprocs",default = None,
                      type = 'int',
                      help="Number of proccessors [OBLIGATORY]")
    parser.add_option("-w", "--weights", dest="weight_cut",default = None,
                      type = 'float',
                      help="Weight cutoff [OBLIGATORY]")
    	parser.add_option("-k", "--skip_vbw", dest="skip_vbw",default = 0,
                      type = 'int',
                      help="Skipping VBW step goes to model evidence directly")
    options, args = parser.parse_args()
    vbwSC.run_vbw(options.restart, options.nstruct, options.priors,\
		options.measures, options.simulated, options.ncurves, options.experimental,\
		options.output, options.nprocs, options.weight_cut, options.skip_vbw)
