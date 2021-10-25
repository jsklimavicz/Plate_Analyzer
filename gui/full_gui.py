# gui/full_gui.py
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

from tkinter import *
from tkinter import ttk
import tkinter as tk
import os
import platform

import gui.utils as pdu
from gui.imageframe.analysisframe import AnalysisFrame as AnalysisFrame
from gui.statsframe.statsframe import StatsFrame as StatsFrame


#############################################################################
#                 Complete App and analysis driver
#############################################################################

class App(tk.Tk):
	'''
	Creates the overal GUI for both the Image AI and for the statistics. 
	These two interfaces are on different tabs.

	tk.Tk
	self.notebook
	|--self.image_frame tab #For image analysis
	|  |--AnalysisFrame
	|--self.stats_frame tab #UID selection for highlighting in plot
	|  |--StatsFrame

	'''
	def __init__(self):
		super().__init__()
		self.config = pdu.parse_config_file() #first look at configuration.
		self.scale = self.det_scaling()
		geom = [700* self.scale*1.5, 700* self.scale*1.0] 

		self.geometry(f'{geom[0]:n}x{geom[1]:n}')
		self.title('Merlin Bioassay Analyzer')

		# create a notebook
		style = ttk.Style()
		style.theme_settings("default", {"TNotebook.Tab": {"configure": {"padding": [50, 20]}}})
		self.notebook = ttk.Notebook(self)
		self.notebook.pack(pady=10, expand=True)

		# create frames
		self.image_frame = AnalysisFrame(self.notebook, 
							config =self.config, 
							scale =self.scale,
							width=geom[0], 
							height=geom[1])

		self.stats_frame = StatsFrame(self.notebook, 
							config =self.config, 
							scale =self.scale,
							width=geom[0], 
							height=geom[1])

		# add frames to notebook
		self.notebook.add(self.image_frame, text='Image Analysis')
		self.notebook.add(self.stats_frame, text='Statistics')

	def det_scaling(self):
		'''
		Apparently Ubuntu20.04 changed its display code such that scaling is wonky. 
		This makes the GUI on Ubuntu a more reasonable size, similar to what it would
		be on Windows.
		'''
		operating_sys = platform.system()
		scale = 1
		if operating_sys == 'Linux':
			version = os.uname()[3].split("~")[1].split("-")[0]
			if int(version.split(".")[0])>=18:
				scale = 2
		return scale


#############################################################################
#                                Run it!
#############################################################################

def main():
	app = App()
	app.mainloop()

if __name__ == "__main__": main()

