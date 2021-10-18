#!/usr/bin/env python3
#statsframe.py

from tkinter import *
from tkinter import ttk, filedialog, messagebox
from tkinter.ttk import Label, Style
from tkcalendar import Calendar 
import random
import os
from os import path, getcwd, mkdir
from os.path import exists
from datetime import datetime
import tkinter as tk
from gui.tooltip import Tooltip #as Tooltip
import threading
import platform
from threading import Thread
import time
from stats.main import analyze_data
from gui.statsframe.datatypeselectionframe import DataTypeSelectionFrame as DTSF
from gui.progmonitor import ProgMonitor



class StatsFrame(ttk.Frame):
	def __init__(self, container, config, scale = 1, **kwargs):
		super().__init__(container, **kwargs)
		self.container = container
		self.scale = scale
		self.config = config
		self.stats_obj = analyze_data(config_path = self.config["STATS_CONFIG_PATH"],
										MASK_RCNN_SAVE = 'True',
										MASK_RCNN_SAVE_NAME = 'larval_counts.csv')
		self.stats_obj.collect_data()

		self.__create_widgets()
		# self.auto_prog_bar_on = False
		# self.plate = None
		# self.__progress()
		
	def __create_widgets(self, scale = 1):

		self.font = ('Arial', 12*self.scale)
		self.selection_frame = DTSF(self, config = self.config, stats_obj = self.stats_obj, scale = self.scale)
		self.selection_frame.grid(column=0, row=0, padx=10, pady=20)

		# Button to run the analysis
		self.Statbutton = Button(self, text="Run Statistics", command = lambda: self.statistics_driver())
		self.Statbutton.grid(row=1, column=0, sticky=S)
		self.Statbutton.config(height = 2)
		Tooltip(self.Statbutton, text='Performs statistical analysis of larval count data, including dose-response curve fitting and LC calculations.')

		self.progmonitor = ProgMonitor(self, run_button = self.Statbutton)
		self.progmonitor.grid(row=2, column=0, sticky=S)

	def __stats_preprocessor(self):
		'''
		Simple preprocessor to update the disallowed UIDs for the different 
		compounds. 
		'''
		disallowed_uids = self.selection_frame.get_disallowed_uids()
		self.stats_obj.set_diallowed(disallowed_uids)

	def statistics_driver(self):
		'''
		Main driver for the Run Statistics button. Updates the disallowed list,
		and then performs the necessary statistics. 
		'''

		#Update stats object with diallowed compounds
		self.__stats_preprocessor()

		#Disable buttons to prevent weird things from happening. 
		self.Statbutton['state'] = tk.DISABLED 
		self.container.tab(0, state="disabled")
		self.selection_frame.disallowButton['state'] = tk.DISABLED 
		self.selection_frame.allowButton['state'] = tk.DISABLED 
		self.selection_frame.clearButton['state'] = tk.DISABLED 
		
		#create process thread for computing stats
		stats_thread=Thread(target = self.__stats_processor)

		#make a daemon to kill the plate thread if the GUI is closed
		stats_thread.daemon = True
		stats_thread.start() #start the thread
		time.sleep(0.5) #brief pause to allow monitoring thread to chill
		#continue monitoring the thread
		self.progmonitor.monitor(stats_thread, self.stats_obj, speed = 0.22)

	def __stats_processor(self):
		'''
		Performs statistical computations on computing threads
		'''	
		self.progmonitor.auto_prog_bar_on = True 
		self.stats_obj.full_process(new_datapath = self.config["OUT_DIR"])
		self.progmonitor.prog_bar.stop()
		self.selection_frame.disallowButton['state'] = tk.NORMAL 
		self.selection_frame.allowButton['state'] = tk.NORMAL 
		self.selection_frame.clearButton['state'] = tk.NORMAL 

