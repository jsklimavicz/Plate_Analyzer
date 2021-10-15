#!/usr/bin/env python3
#full_gui.py

from tkinter import *
from tkinter import ttk, filedialog, messagebox
from tkinter.ttk import Label, Style
import random
import os
from os import path, getcwd, mkdir
from os.path import exists
from datetime import datetime
import tkinter as tk
from gui.tooltip import Tooltip #as Tooltip
import gui.utils as pdu
import platform
import time
from stats.main import analyze_data
import multiprocessing
from gui.imageframe.analysisframe import AnalysisFrame as AnalysisFrame
from gui.statsframe.statsframe import StatsFrame as StatsFrame


#############################################################################
#                 Complete App and analysis driver
#############################################################################

class App(tk.Tk):
	def __init__(self):
		super().__init__()
		self.config = pdu.parse_config_file()
		self.scale = self.det_scaling()
		geom = [700* self.scale*1.5, 700* self.scale*1.0] 

		self.geometry(f'{geom[0]:n}x{geom[1]:n}')
		self.title('Notebook Demo')

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

