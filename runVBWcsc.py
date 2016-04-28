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
import vbwCSC

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
    parser.add_option("-m", "--saxs_measures", dest="saxs_measures",default = None,
                      type = 'int',
                      help="Number of CS points [OBLIGATORY]")
    parser.add_option("-k", "--cs_measures", dest="cs_measures",default = None,
                      type = 'int',
                      help="Number of CS points [OBLIGATORY]")
    parser.add_option("-v", "--number_of_curves", dest="ncurves",default = 1,
                      type = 'int',
                      help="Number of strcutural model [OBLIGATORY]")
    parser.add_option("-p", "--priors", dest="priors",
                      help="Prior weights [OBLIGATORY]")
    parser.add_option("-s", "--saxs_simulated", dest="saxs_simulated",
                      help="Simulated SAXS curves [OBLIGATORY]")
    parser.add_option("-e", "--saxs_experimental", dest="saxs_experimental",
                      help="Experimental SAXS curves [OBLIGATORY]")
    parser.add_option("-t", "--cs_simulated", dest="cs_simulated",
                      help="Simulated CS curves [OBLIGATORY]")
    parser.add_option("-u", "--cs_rms", dest="cs_rms",
                      help="Simulated CS rms [OBLIGATORY]")
    parser.add_option("-f", "--cs_experimental", dest="cs_experimental",
                      help="Experimental CS curves [OBLIGATORY]")
    parser.add_option("-o", "--output", dest="output",
                      help="Output file [OBLIGATORY]")
    parser.add_option("-c", "--cores", dest="nprocs",default = None,
                      type = 'int',
                      help="Number of proccessors [OBLIGATORY]")
    parser.add_option("-w", "--weights", dest="weight_cut",default = None,
                      type = 'float',
                      help="Weight cutoff [OBLIGATORY]")
    options, args = parser.parse_args()
    vbwCSC.run_vbw(options.restart, options.nstruct, options.priors,\
		options.saxs_measures, options.saxs_measures, options.ncurves,
        options.saxs_simulated,  options.saxs_experimental,\
        options.cs_simulated, options.cs_rms,  options.cs_experimental,\
		options.output, options.nprocs, options.weight_cut)
