import argparse
import os
import subprocess
import sys
import tempfile

from svsim.reads.metasim import MetaSimSimulator
from svsim.reads.dwgsim import DwgsimSimulator
from svsim.reads.lognsim import LogNSimulator

##
# Returns a simulator that can be used to simulate reads.
#
# @param simulator A string representing the simulator to use,
#                  either 'metasim' or 'dwgsim'.
#
# @return A matching simulator, or raises ValueError if no simulator
#         was found.
#
def get_simulator(simulator):
    if simulator == "metasim":
        return MetaSimSimulator( )
    elif simulator == "dwgsim":
        return DwgsimSimulator( )
    elif simulator == "lognsim":
        return LogNSimulator( )
    else:
        raise ValueError( "No such simulator found" )

if __name__ == '__main__':
    parser = argparse.ArgumentParser( description="Simulates illumina reads with metasim or simulater from the more general logNormal distribution." )
    parser.add_argument( 'genome_file', type=str, help='Path to the genome' )
    parser.add_argument( 'output_prefix', type=str, help='Output prefix for the paired end files, will apped _pe1.fa and pe2.fa to this.' )
    parser.add_argument( '-c', type=float, help='Coverage.', default=10.0 )
    parser.add_argument( '-m', type=float, help='Mean of the library distribution.', default=550.0 )
    parser.add_argument( '-s', type=float, help='Standard deviation of the library distribution.', default=50.0 )
    parser.add_argument( '-t', type=str, choices=[ "metasim", "dwgsim", "lognsim" ], help="Type of simulator 'metasim' or 'dwgsim'.", required=True )

    args = parser.parse_args( )

    simulator = get_simulator( args.t )
    simulator.coverage = args.c
    simulator.mean = args.m
    simulator.std = args.s
    simulator.read_length = 100

    simulator.simulate( args.genome_file, args.output_prefix )
