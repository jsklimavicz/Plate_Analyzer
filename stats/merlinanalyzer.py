# stats/merlinanalyzer.py
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

import numpy as np
import os
import sys
import pandas as pd
from stats.curvell import CI_finder
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from stats.compound import Compound
import math
from copy import deepcopy
import pickle
import hashlib
import hmac
from stats.merlin_grapher import MerlinGrapher
import stats.utils as utils
from time import time
from stats.latex_writer import LatexWriter


class MerlinAnalyzer:
	'''
	Class for processing all bioassay data.
	'''
	archivefilename = "stats.pickle"
	picklesha1hash = ".statshash"
	sha_key = b"merlin-data"
	cache_path = os.path.abspath('./stats/cache')

	def __init__(self, 
					*args, 
					config_path = os.path.abspath('./stats/analysis_config.txt'), 
					**kwargs):
		self.cmpd_data = {}
		self.options = utils.parse_config_file(config_path = config_path, **kwargs)
		self.message = ""
		self.error_message = None
		self.progress = 0.0
		# print(self.options)

	def column_name_modifier(self, filename):
		'''
		Method to ensure all column names are consistent between files in case anyone makes minor modifications.
		'''
		new_data = pd.read_csv(filename, header = 0)
		column_names = list(new_data.columns)
		new_col_names = {}
		for col_name in column_names:
			#live counts
			if col_name.lower() in ["alive", 'live', 'living']: new_col_names[col_name] = "Live"
			#dead counts
			elif col_name.lower() in ["dead"]: new_col_names[col_name] = "Dead"
			#total count
			elif col_name.lower() in ["total", 'count', 'sum']: new_col_names[col_name] = "Count"
			#ID for a compound
			elif col_name.lower() in ["ref.id", 'ref id', 'ref', 'id']: new_col_names[col_name] = "ID"
			elif 'col' in col_name.lower(): new_col_names[col_name] = "Column"
			elif 'row' in col_name.lower(): new_col_names[col_name] = "Row"
			#Concentration of compound
			elif col_name.lower() in ["ppm", "conc", "concentration"]: new_col_names[col_name] = "Conc"
			elif col_name.lower() in ["plate"]: new_col_names[col_name] = "Plate"
			elif col_name.lower() in ["rep"]: new_col_names[col_name] = "Rep"
			elif col_name.lower() in ["class"]: new_col_names[col_name] = "Class"
			elif col_name.lower() in ["date", 'day']: new_col_names[col_name] = "Date"
			else: new_col_names[col_name] = col_name
		new_data.rename(columns=new_col_names, inplace=True)
		return new_data

	def read_new_data(self, filename, key_file):
		'''
		Driver to read in new data based on a .csv filename and a key.csv file
		'''

		self.message = "Reading new data."
		#If read directly from the plate reader, we must go through all the permitted folders. 
		if self.options["MASK_RCNN_SAVE"]:
			dirs = [x[0] for x in os.walk(filename)]
			readable_files = []
			for path in dirs:
				csvpath = os.path.join(path, self.options["MASK_RCNN_SAVE_NAME"])
				if os.path.exists(csvpath): readable_files.append(csvpath)
			cmpd_data = self.new_file_reader_helper(readable_files[0])
			for csv in readable_files[1:]:
				new_data = self.new_file_reader_helper(csv)
				cmpd_data = cmpd_data.append(new_data, ignore_index=True)

		else:
			cmpd_data = self.new_file_reader_helper(filename)
		#read key and merge in compound names
		if key_file is not None: key = self.read_key(key_file)
		cmpd_data = cmpd_data.merge(key, how='left', on='ID')
		cmpd_data['Compound'] =  cmpd_data['Compound'].str.lower()
		#Remove any rows with missing count or compound data.
		cmpd_data.dropna(subset = ['Live', 'Dead', 'Compound'], inplace=True)
		cmpd_data = cmpd_data[cmpd_data.Count != 0]
		return self.process_compounds(cmpd_data)


	def new_file_reader_helper(self, filename):
		new_data = self.column_name_modifier(filename) #check/change filenames
		cmpd_data = deepcopy(new_data[new_data["Conc"] != 0])
		#Find control mortalities
		ctrl_data = deepcopy(new_data[new_data["Conc"] == 0])
		ctrl_live_sum = sum(np.array(ctrl_data["Live"].tolist()))
		ctrl_dead_sum = sum(np.array(ctrl_data["Dead"].tolist()))
		ctrl_ave_mort = (ctrl_dead_sum *1.)/(ctrl_dead_sum + ctrl_live_sum *1.)
		cmpd_data["ctrl_mort"] = ctrl_ave_mort

		return cmpd_data


	def read_key(self, filename):
		#Select important columns for merging key file
		return self.column_name_modifier(filename)[["Compound", "ID", "Class"]]
		
	def read_archive(self):
		'''
		Reads in old pickle file after making sure that the file is not corrupted/modified
		by means of using a sha1 hash
		'''
		self.message = "Reading archived data."
		with open(os.path.join(self.cache_path, self.picklesha1hash), 'r') as file:
			pickle_hash = file.read().strip()
		with open(os.path.join(self.cache_path, self.archivefilename), 'rb') as file:
			pickled_data = file.read()
		digest =  hmac.new(self.sha_key, pickled_data, hashlib.sha1).hexdigest()
		self.progress = 2.0 #update progress for each compound

		if pickle_hash == digest:
			unpickled_data = pickle.loads(pickled_data)
			return unpickled_data
		else:
			print('Pickled data as been compromised. Old data cannot be loaded.')

	def merge_old_new(self, new_datapath, key_file):
		if os.path.exists(os.path.join(self.cache_path, self.archivefilename)): 
			self.cmpd_data = self.read_archive()
			# print(self.options)
			for cmpd in self.cmpd_data.keys():
				# print(self.cmpd_data[cmpd].options)
				dict_change = utils.check_library_change(self.cmpd_data[cmpd].options, self.options)
				# print(self.cmpd_data[cmpd].__dict__)
				if dict_change: 
					self.cmpd_data[cmpd] = self.cmpd_data[cmpd].reset_curves()
				#update dictionaries
				for k, v in self.options.items():
					self.cmpd_data[cmpd].options[k] = v
					# self.cmpd_data[cmpd].curve_data.options[k] = v


		if new_datapath is not None:
			new_cmpd_dict = self.read_new_data(new_datapath, key_file)
			#merge
			for k, v, in new_cmpd_dict.items():
				if k in self.cmpd_data: self.cmpd_data[k] = self.cmpd_data[k] + new_cmpd_dict[k]
				else: self.cmpd_data[k] = v
				# self.cmpd_data[k].test_print()

		self.progress = 4.0 #update progress for each compound

	def save_archive(self):
		saveable_lib = {}
		for k, v in self.cmpd_data.items():
			saveable_lib[k] = self.cmpd_data[k].saveable_cmpd()
		pickle_data = pickle.dumps(saveable_lib)
		digest =  hmac.new(self.sha_key, pickle_data, hashlib.sha1).hexdigest()
		header = '%s' % (digest)
		if not os.path.exists(self.cache_path): os.makedirs(self.cache_path)
		with open(os.path.join(self.cache_path, self.picklesha1hash), 'w') as file:
			file.write(header)
		with open(os.path.join(self.cache_path, self.archivefilename), 'wb') as file:
			file.write(pickle_data)
		
	def process_compounds(self, new_data, *args, **kwargs):
		self.message = "Processing compounds with new data."
		new_compound_dict = {}
		self.number_of_compounds = len(new_data["Compound"].unique())
		for cmpd_id in new_data["Compound"].unique():

			cmpd_data = new_data[new_data["Compound"] == cmpd_id].copy()
			#UID = name_date_plate_row_ID
			unique_ids = ["_".join([cmpd_id,x,str(y),str(z),w]) for x,y,z,w in zip(cmpd_data["Date"].tolist(), 
																		cmpd_data["Plate"].tolist(), 
																		cmpd_data["Row"].tolist(),
																		cmpd_data["ID"].tolist())]
			new_compound_dict[cmpd_id] = Compound(name = cmpd_id, 
							ids = cmpd_data["ID"].unique(), 
							test_dates = cmpd_data["Date"].unique(),
							max_conc = max(cmpd_data["Conc"].tolist()),
							min_conc = min(cmpd_data["Conc"].tolist()),
							n_trials = len(set(unique_ids)),
							column_IDs = cmpd_data["Column"].tolist(),
							row_IDs = cmpd_data["Row"].tolist(),
							conc = np.log(np.array(cmpd_data["Conc"].tolist()))/math.log(2),
							live_count = np.array(cmpd_data["Live"].tolist()),
							dead_count = np.array(cmpd_data["Dead"].tolist()),
							plate_ids = cmpd_data["Plate"].tolist(),
							reps = cmpd_data["Rep"].tolist(),
							ctrl_mort = np.array(cmpd_data["ctrl_mort"].tolist()),
							unique_plate_ids = unique_ids,
							include_now = [True] * len(unique_ids),
							*args, **kwargs)

			self.progress += 2.0/self.number_of_compounds #update progress for each compound
		return new_compound_dict

	def set_diallowed(self, disallowed_uids): 
		for cmpd, key in self.cmpd_data.items():
			for ind, uid in enumerate(key.data['unique_plate_ids']):
				key.data['include_now'][ind] = False if uid in disallowed_uids else True

	def get_uid_list(self):
		uid_list = []
		for k, v in self.cmpd_data.items(): 
			uid_list += list(v.data['unique_plate_ids'])
		return uid_list

	def save_csv(self, filename, header, body):
		with open(filename, 'w') as file:
			file.write(",".join(header) + "\n")
			file.write("\n".join(body))
	
	def generate_csv_data_lines(self, header):
		output = []
		for cmpd_name, cmpd in self.cmpd_data.items():
			#check to see if the data is included at all. 
			if not any(cmpd.data['include_now']):
				continue

			line = []
			comment = " "
			good_curve = True
			#Check to make sure that the slope is positive enough
			slope_info = cmpd.curve_data.get_slope_CI(CI_val = 0.95)

			if slope_info[1] < 1.e-2:
				comment += "Fitted slope is too shallow. "
				good_curve = False

			#Check to make sure that the LC50 value is close enough to the 
			LC50_info = np.power(2.,cmpd.curve_data.get_LC50_CI(CI_val=0.95, log = True))
			# print("The outer get_LC_CIs function", LC50_info)

			lc_vals = utils.format_LC_to_CSV(cmpd.get_LC_CIs()) 
			# print(cmpd_name, LC50_info[1], 2**(1 + cmpd.data["max_conc"]) ,  2**(cmpd.data["min_conc"] - 1))

			if LC50_info[1] > self.options['EXTRAPOLATION_FACTOR']**(1 + cmpd.data["max_conc"]) or \
					LC50_info[1]< self.options['EXTRAPOLATION_FACTOR']**(cmpd.data["min_conc"] - 1) : 
				# print(LC50_info[1], {lc_vals[0]}, cmpd.data["min_conc"]/2., 2*cmpd.data["max_conc"])
				comment += f"Calculated LC50 ({utils.format_lc_val(LC50_info[1])}) out of bounds. "
				lc_vals = ['NA'] * len(lc_vals)
				good_curve = False
				
			if good_curve:
				rel_pot = self.compare_LC(cmpd = cmpd_name, n_bs = 100000)
				rel_pot = utils.format_LC_to_CSV(rel_pot)
			else: 
				lc_vals = ['NA'] * len(lc_vals)
				rel_pot = ["NA"] * 2*len(self.options['LC_VALUES'])

			rel_done = False
			LC_done = False

			biol, tech = cmpd.get_rep_counts()

			for item in header:
				if item.lower() in 'compound': line.append(cmpd.data["name"])
				elif 'rows' in item.lower(): line.append(str(tech))
				elif item.lower() in 'slope': 
					line.append(utils.format_lc_val(slope_info[1]))
				elif 'slope' in item.lower() and 'ci' in item.lower(): 
					line.append(utils.CI_to_string(slope_info[0], slope_info[2]))
				elif self.options['REFERENCE_COMPOUND'].lower() in item.lower():
					if rel_done: continue
					else:
						line = [*line, *rel_pot]
						rel_done = True
				elif self.options['REFERENCE_COMPOUND'].lower() not in item.lower() and "lc" in item.lower():
					if LC_done: continue
					else:
						line = [*line, *lc_vals]
						LC_done = True
				elif "codes" in item.lower(): line.append('"' + ", ".join(list(set([i for i in cmpd.data["ids"]]))) + '"')
				elif "date" in item.lower(): line.append('"' + ", ".join(list(set([i for i in cmpd.data["test_dates"]]))) + '"')
				elif "bio" in item.lower(): line.append(f"{biol}")
				elif item == "R2": line.append(utils.format_lc_val(cmpd.curve_data.r2))
				elif "comment" in item.lower():
					comment = comment if len(comment) == 1 else comment[1:len(comment)] 
					line.append(comment)
				else: line.append(" ")
			output.append(",".join(line))
		return output

	def generate_csv_header(self):
		LC_title_names = []
		for idx, LC in enumerate(self.options['LC_VALUES']): 
			LC_title_names.append("LC" + str(round(self.options['LC_VALUES'][idx] * 100)))
			LC_title_names.append("LC" + str(round(self.options['LC_VALUES'][idx]* 100)) + 
						" " + str(round(self.options['LC_CI']* 100)) + "%CI")
		rel_pot_col = []
		if self.options['REL_POT_TO_REF']:
			init_name = self.options['REFERENCE_COMPOUND'] + ' Pot. Rel. to Cmpd at LC'
		else:
			init_name = "Pot. of Cmpd Rel. to " + self.options['REFERENCE_COMPOUND'] + " at LC"
		for idx, LC in enumerate(self.options['LC_VALUES']): 
			rel_pot_col.append(init_name + str(round(self.options['LC_VALUES'][idx] * 100)))
			rel_pot_col.append(init_name + str(round(self.options['LC_VALUES'][idx] * 100)) + 
						" " + str(round(self.options['REL_POT_CI']* 100)) + "%CI")
		header = ['Compound', 
					'Biological Reps',
					'Total Rows', 
					*LC_title_names,
					*rel_pot_col,
		 			'R2',
		 			'Slope', 
		 			'Slope CI',  
		 			'Tested Codes', 
		 			'Test Dates',
		 			'Comments']
		return header

	@utils.surpress_warnings
	def compare_LC(self, cmpd, n_bs = 100000):
		'''
		Calculates relative potency using random samples from the kernel density of the LCx distribution for
		cmpd1 and cmpd2.
		'''
		lcs = self.options['LC_VALUES']
		vals = np.zeros((len(lcs),3))
		for ind, lc in enumerate(lcs):
			cmpd_kernel = self.cmpd_data[cmpd].curve_data.LC_kernel(
				LC_val = 1. - lc).sample(n_samples=n_bs)
			ref_kernel= self.cmpd_data[self.options['REFERENCE_COMPOUND']].curve_data.LC_kernel(
				LC_val = 1. - lc).sample(n_samples=n_bs)
			diff = cmpd_kernel - ref_kernel if self.options['REL_POT_TO_REF'] else kernel2 - kernel1
			func = utils.calc_ET_CI if self.options["CI_METHOD"].lower() in utils.ET_VARS else utils.calc_HPDI_CI
			val = func(diff, CI_level = self.options['REL_POT_CI'])
			vals[ind,:] = np.power(2., val).T
		return vals

	def collect_data(self,
						new_datapath = None, 
						key_file = None,
						archive_path = None):
		
		new_datapath =  new_datapath if new_datapath is not None else self.options["INPUT_PATH"]
		archive_path = self.cache_path if archive_path is not None else self.options["ARCHIVE_PATH"]
		key_file = key_file if key_file is not None else self.options["KEY_PATH"]

		self.merge_old_new(new_datapath = new_datapath, key_file = key_file)


	def full_process(self, 
						new_datapath = None, 
						key_file = None,
						csv_outfile = None, 
						out_path = None,
						pdf_outfile = None,
						archive_path = None,
						*args, 
						**kwargs):

		
		archive_path = self.cache_path if archive_path is not None else self.options["ARCHIVE_PATH"]
		key_file = key_file if key_file is not None else self.options["KEY_PATH"]
		out_path = out_path if out_path is not None else self.options["SAVE_PATH"]
		pdf_outfile = pdf_outfile if pdf_outfile is not None else self.options["OUTPUT_PDF_NAME"]
		csv_outfile = csv_outfile if csv_outfile is not None else self.options["OUTPUT_CSV_NAME"]

		#Set progress to zero; this is needed if running the stats more than once to prevent 100% completions.
		self.progress = 0.

		image_dir = os.path.abspath(os.path.join(out_path, 'images'))
		pdf_dir = os.path.abspath(os.path.join(image_dir, 'pdf'))

		if not os.path.exists(image_dir):
			os.makedirs(image_dir)
		if not os.path.exists(pdf_dir):
			os.makedirs(pdf_dir)

		LW = LatexWriter(img_folder = pdf_dir)

		for cmpd_name, cmpd in self.cmpd_data.items():
			self.progress += 90.0/self.number_of_compounds #update progress for each compound

			#check to see if the data is included at all. 
			if not any(cmpd.data['include_now']):
				continue

			self.message = f"Calculating LC data and dose-response curves for {cmpd_name:s}."
			cmpd.fit_data(options = self.options)
			cmpd.make_plot()
			pdf_path = cmpd.plot.save_plot(cmpd_name, image_dir, pdf_dir = pdf_dir)
			lc50lb, lc50med, lc50ub = cmpd.curve_data.get_LC50_CI(CI_val=0.95).squeeze()

			if lc50med > self.options['EXTRAPOLATION_FACTOR']**(1 + cmpd.data["max_conc"]) or \
				lc50med< self.options['EXTRAPOLATION_FACTOR']**(cmpd.data["min_conc"] - 1): lcTrue = False 
			else: lcTrue = True 
				
			lc50med = utils.format_lc_val(lc50med)
			lc50CI = '[' + utils.format_lc_val(lc50lb) + ', ' + utils.format_lc_val(lc50ub) + ']'
			
			biol, tech = cmpd.get_rep_counts()

			LW.make_cmpd_graph(image_dir = pdf_path, 
				name = cmpd_name, 
				lc50 = lc50med, 
				lc50CI = lc50CI, 
				lcTrue = lcTrue,
				R2 = utils.format_lc_val(cmpd.curve_data.r2), 
				bio_reps = biol,
				tech_reps = tech)

		self.message = "Archiving data."
		self.progress = 96.0  #update progress for each compound
		self.save_archive()

		self.message = "Creating output csv."
		self.progress = 98.0  #update progress for each compound
		header = self.generate_csv_header()
		body = self.generate_csv_data_lines(header)
		if csv_outfile: self.save_csv(os.path.abspath(os.path.join(out_path, csv_outfile)), header, body)

		self.progress = 99.5  #update progress for each compound
		self.message = "Creating output pdf with dose-response curves."
		LW.write_file(out_path = os.path.abspath(os.path.join(out_path, pdf_outfile +".tex")), stat_lib = self.options )
		self.progress = 99.5
		self.message = "Statistical analysis complete."

