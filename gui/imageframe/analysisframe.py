#!/usr/bin/env python3
#analysisframe.py

from tkinter import *
import tkinter as tk
from tkinter import ttk
from tkinter.ttk import Label, Style

import os
import time
from threading import Thread

from gui.tooltip import Tooltip #as Tooltip
import gui.utils as pdu
from gui.imageframe.ioframe import IOFrame as IOFrame
from gui.imageframe.cmpdframe import CmpdFrame as CmpdFrame
from gui.progmonitor import ProgMonitor


class AnalysisFrame(ttk.Frame):
	def __init__(self, container, config, scale = 1, **kwargs):
		super().__init__(container, **kwargs)
		self.container = container
		self.scale = scale
		self.config = config
		self.__create_widgets()
		# self.auto_prog_bar_on = False
		# self.plate = None
		# self.__progress()

	def __create_widgets(self, scale = 1):
		# create the input frame
		self.font = ('Arial', 12*self.scale)
		self.input_frame = IOFrame(self, config = self.config, scale = self.scale)
		self.input_frame.grid(column=0, row=0, padx=10, pady=20)

		# create the button frame
		font = ('Arial', 12*self.scale)
		self.cmpd_frame = CmpdFrame(self, config = self.config, scale = self.scale)
		self.cmpd_frame.grid(column=0, row=1, padx=10, pady=20)
		style = ttk.Style()
		style.configure("correct_font", font = font)
		self.cmpd_frame.configure(borderwidth = 1, relief = 'raised', text = "Plate Configuration")

		# Button to run the analysis of images
		self.Convertbutton = Button(self, text="Run Image Analysis", command = lambda: self.threaded_plate_driver())
		self.Convertbutton.grid(row=2, column=0, sticky=S)
		self.Convertbutton.config(height = 2)
		Tooltip(self.Convertbutton, text='Analyzes the plate images and saves the data when complete.')

		self.progmonitor = ProgMonitor(self, run_button = self.Convertbutton)
		self.progmonitor.grid(row=3, column=0, sticky=S)

	def __plate_setup(self):
		# from platedriver.plate import Plate 
		from platedriver.platedata import PlateData
		from platedriver.plate import Plate 
		#update entry values
		# self.input_frame.update()
		self.cmpd_frame.update()
		# try : 
		#make a plate data object
		plate_data = PlateData(dateval = self.input_frame.config["DATE"], 
			plateID = self.input_frame.PlateIDButton.get(), 
			control_name = self.input_frame.config["CTRL_NAME"], 
			control_ID = self.input_frame.config["CTRL_CODE"])
		for i in range(len(self.cmpd_frame.param.row_names)):
			plate_data.add_row_data(row = self.cmpd_frame.param.row_names[i],
						compound_ID = self.cmpd_frame.param.cmpd_codes[i].get(), 
						compound_name = self.cmpd_frame.param.cmpd_names[i].get(), 
						max_conc = float(self.cmpd_frame.param.concentrations[i].get()), 
						rep_num = self.cmpd_frame.param.reps[i])

		
		self.plate = Plate(config = self.input_frame.config,
				save_dir = self.out_dir,
				plate_data=plate_data)

	def __plate_driver(self):
		
		if self.plate.error_message: return
		self.progmonitor.auto_prog_bar_on = True 
		self.plate.run_images(init_only=False)
		self.plate.detect_larvae()
		self.progmonitor.auto_prog_bar_on = False
		self.plate.save_csv(filename = self.input_frame.config["CSV_NAME"], \
						verbose_csv = self.input_frame.config["VERBOSE_CSV"], \
						open_csv=self.input_frame.config["OPEN_CSV"])
		self.progmonitor.prog_bar.stop()

	def threaded_plate_driver(self):
		#check to see if the image input folder exists to be read from:
		self.input_frame.update()

		self.input_frame.config['IMG_DIR_ID'] = pdu.refresh_output_dir(self.input_frame.config['IMG_DIR'])
		self.out_dir = os.path.join(self.input_frame.config['OUT_DIR'], self.input_frame.config['IMG_DIR_ID'])
		if not os.path.exists(self.out_dir): os.mkdir(self.out_dir)

		if not os.path.isdir(self.input_frame.config['IMG_DIR']):
			text = f'"{self.input_frame.config["IMG_DIR"]}" is not a directory. \nPlease check the image folder path specified above.'
			self.prog_label.config(text=text)
			return
		#check to see if the image output folder exists to be written to:
		if not os.path.isdir(self.input_frame.config['OUT_DIR']):
			text = f'{self.input_frame.config["OUT_DIR"]} is not a directory. \nPlease check the image output folder path specified above.'
			self.prog_label.config(text=text)
			return
		#check to see if the CSV file exists to be written to:
		file_path = os.path.join(self.out_dir, self.input_frame.config["CSV_NAME"])
		if not os.path.exists(file_path):
			try: 
				file = open(file_path, 'a')
				file.close()
				if os.stat(file_path).st_size == 0: os.remove(file_path)
			except:
				text = f'"{file_path}" cannot be written to. \nPlease check the output filename specified above.'
				self.prog_label.config(text=text)
				return
		elif os.stat(file_path).st_size == 0: os.remove(file_path)

		self.input_frame.reset_manual_indicators() #reset manual override halt on autoupdating date/csv
		#disable conert button
		self.Convertbutton['state'] = tk.DISABLED
		self.container.tab(1, state="disabled")
		#make progress bar
		
		self.__plate_setup()
		
		#make thread for the plate driver
		plate_thread=Thread(target=self.__plate_driver)
		#make a daemon to kill the plate thread if the GUI is closed
		plate_thread.daemon = True
		#start the thread
		plate_thread.start()
		time.sleep(0.5)
		#continue monitoring the thread
		# self.progmonitor.set_obj(obj = self.plate)
		self.progmonitor.monitor(plate_thread,obj = self.plate)
