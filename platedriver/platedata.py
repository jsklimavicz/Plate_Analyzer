# platedriver/platedata.py
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
from datetime import datetime




### Contains treatment, concentration, date, etc for a plate. 

class PlateData:
	def __init__(self, dateval = None, plateID = 1, control_name = "VC", control_ID = "LM (T)"):
		self.date = dateval if dateval else datetime.today()
		self.wellID_data = {}
		self.initiate_data()
		self.plate_number = plateID
		self.control_name = control_name
		self.control_ID = control_ID

	def initiate_data(self):
		for row in range(ord("A"),ord("H")+1):
			for col in np.r_[1:11,12]: 
				well_dict = {}
				well_dict["Compound"] = ""
				well_dict["Ref ID"] = ""
				well_dict["ppm"] = ""
				well_dict["Rep"] = 1
				wellID = chr(row) + str(col)
				self.wellID_data[wellID] = well_dict

	def add_row_data(self, row, compound_ID = "", compound_name = "", max_conc = np.nan, rep_num = 1):
		for i in range(1,11):
			well_dict = {}
			well_dict["Compound"] = compound_name
			well_dict["Ref ID"] = compound_ID
			well_dict["ppm"] = max_conc / 2**(i-1)
			well_dict["Rep"] = rep_num
			wellID = row + str(i)
			self.wellID_data[wellID] = well_dict
		well_dict = {}
		well_dict["Compound"] = self.control_name 
		well_dict["Ref ID"] = self.control_ID
		well_dict["ppm"] = 0
		well_dict["Rep"] = 1
		wellID = row + "12"
		self.wellID_data[wellID] = well_dict

	def __str__(self):
		output = f"Date: {self.date}\n"
		output += f"PlateID: {self.plate_number}\n"
		output += f"Control Name: {self.control_name}\n"
		output += f"Control Code: {self.control_ID}\n"
		for i in self.wellID_data.keys():
			curr = self.wellID_data[i]
			output += f'Well: {i}; Ref ID: {curr["Ref ID"]}; Name: {curr["Compound"]}; Conc: {curr["ppm"]}; Rep: {curr["Rep"]}\n'
		return output