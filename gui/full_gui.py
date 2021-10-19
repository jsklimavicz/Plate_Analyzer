#!/usr/bin/env python3
#full_gui.py

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
	|--image_frame tab #For image analysis
	|  |--AnalysisFrame
	|--stats_frame tab #UID selection for highlighting in plot
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

