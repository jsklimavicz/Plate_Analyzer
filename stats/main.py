# stats/main.py
# This file is a component of Plate_Analyzer, which can count Drosophila
# L1 larvae, classify them as alive or dead, and determine dose-repsonse
# information based on live/dead count data. 

# Copyright (C) 2021 James Klimavicz

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from stats.merlinanalyzer import MerlinAnalyzer
import sys
import argparse
import os
import errno
from stats.utils import parse_config_file as parse_config

parser = argparse.ArgumentParser(description="""Produces pdf of figures.""")
parser.add_argument("-i", "--input", metavar="/path/to/input.csv", 
	default = None, required = True,
	help="Input live/count file for Merlin in .csv format.")
parser.add_argument("-k", "--key", metavar="/path/to/key.csv", 
	default = None, required = True,
	help="Key file that contains the compound codes and names.")
parser.add_argument("-p", "--path", metavar="/path/for/saving", 
	default = None, required = True,
	help="Directory in which the output will be saved.")
parser.add_argument("-o", "--output_csv", metavar="output_filename.csv", default = "output.csv",
	help="Name of the output .csv file.")
parser.add_argument("-a", "--archive_path", metavar="/path/to/pickle/archive", default = '.',
	help="Name of the output .csv file.")
parser.add_argument("-g", "--output_pdf", metavar="graphs_filename", default = "graphs",
	help="Name of the output file containing the graphs WITHOUT file extension.")


#For use by the LArval counter AI; keywords are specified via the config file as 
#opposed to the command line
def analyze_data(config_path = os.path.abspath('./stats/analysis_config.txt'), **kwargs):
	return MerlinAnalyzer(config_path = config_path, **kwargs)
	
#For command line running. 
if __name__ == "__main__":
	args = parser.parse_args()
	args.archive_path = args.path
	MA = MerlinAnalyzer()
	MA.full_process(new_datafile = args.input, 
		key_file = args.key, 
		out_path=args.path, 
		archive_path = args.archive_path,
		csv_outfile = args.output_csv,
		pdf_outfile = args.output_pdf)
