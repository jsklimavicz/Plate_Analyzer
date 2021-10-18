#!/usr/bin/env python3
#ioframe.py

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
import gui.utils as pdu
import platform
import time
from stats.main import analyze_data

#############################################################################
#                 For entering background plate info
#############################################################################

class IOFrame(ttk.Frame):
	def __init__(self, container, config = None, scale = 1):
		self.scale = scale
		self.font = ('Arial', 12*self.scale)
		self.config = config
		self.plateID = self.config['PLATE_NUMBER']
		print(self.plateID)
		self.textbox_width=60
		super().__init__(container)
		# setup the grid layout manager
		self.columnconfigure(0, weight=1)
		self.columnconfigure(1, weight=2)
		self.columnconfigure(2, weight=2)
		self.columnconfigure(3, weight=2)
		self.columnconfigure(4, weight=2)
		self.columnconfigure(5, weight=2)
		self.row_num = 0
		self.csv_name_manual_mod = False
		self.date_manual_mod = False
		self.__create_widgets()
		

	def __create_widgets(self):
		s = ttk.Style()
		s.configure("my.TButton", font = self.font)
		full_width = 6
		#Row: Input Image Folder Selector
		self.imageFolderLabel = Label(self, text="Image folder path:")
		self.imageFolderLabel.grid(row=self.row_num, column=0, sticky=E)
		self.imagefolder_var = StringVar()
		self.imagefolder_var.set(self.config['MOST_RECENT_IMG_DIR'])
		self.imageFolderEntry = Entry(self, textvariable= self.imagefolder_var, width=self.textbox_width * self.scale, font = self.font)
		self.imageFolderEntry.grid(row=self.row_num, column=1, columnspan = full_width,  sticky=EW)
		self.imageBrowseButton = ttk.Button(self, text ='Browse', command = lambda: self.image_directory_chooser(), style = "my.TButton")
		self.imageBrowseButton.grid(row=self.row_num, column=full_width+1, sticky=(S,W))
		self.row_num += 1
		Tooltip(self.imageFolderLabel, text='Select the folder containing the Biotek platereader images.')
		Tooltip(self.imageBrowseButton, text='Select the folder containing the Biotek platereader images.')

		#Row: Image output folder selector
		self.outputFolderLabel = Label(self, text="Output parent directory:")
		self.outputFolderLabel.grid(row=self.row_num, column=0, sticky=E)
		self.outputfolder_var = StringVar()
		self.outputfolder_var.set(self.config['OUT_DIR'])
		self.outputFolderEntry = Entry(self, textvariable=self.outputfolder_var, width=self.textbox_width* self.scale, font = self.font)
		self.outputFolderEntry.grid(row=self.row_num, column=1, columnspan = full_width, sticky=EW)
		self.outputBrowseButton = ttk.Button(self, text ='Browse', command = lambda: self.image_output_chooser(), style = "my.TButton")
		self.outputBrowseButton.grid(row=self.row_num, column=full_width+1, sticky=(S,W))
		self.row_num += 1
		msg = 'Select a folder for saving all output images.'
		Tooltip(self.outputFolderLabel, text = msg)
		Tooltip(self.outputBrowseButton, text = msg)

		#Row: Image saving selector
		self.ImageChooserLabel = Label(self, text="Select output options:")
		self.ImageChooserLabel.grid(row=self.row_num, column=0, sticky=E)
		Tooltip(self.ImageChooserLabel, text='Controls which composite images are made.')

		self.INTMDLabel = Label(self, text="INTMD images:")
		self.INTMDLabel.grid(row=self.row_num, column=1, sticky=W)
		self.INTMD_var = IntVar()
		self.INTMD_var.set(self.config["SAVE_INTMD"])
		self.INTMDBox = Checkbutton(self, variable=self.INTMD_var)
		self.INTMDBox.grid(row=self.row_num, column=2,  sticky=W)
		msg = '''Check whis box if you would like to save the RGB images showing motion of larvae over the three frames.'''
		Tooltip(self.INTMDLabel, text=msg)

		self.MASKLabel = Label(self, text="Mask:")
		self.MASKLabel.grid(row=self.row_num, column=3, sticky=W)
		self.MASK_var = IntVar()
		self.MASK_var.set(self.config["SAVE_MASK"])
		self.MASKBox = Checkbutton(self, variable=self.MASK_var)
		self.MASKBox.grid(row=self.row_num, column=4,  sticky=W)
		Tooltip(self.MASKBox, text='Select this option to produce an image wherein identified and classified larvae are'+
			'filled in based on where Mask R-CNN identifies them. Masks are color-coded for live, dead, moribund, or egg.')

		self.BBOXImages = Label(self, text="Bounding Box:")
		self.BBOXImages.grid(row=self.row_num, column=5, sticky=W)
		self.BBOXImage_var = IntVar()
		self.BBOXImage_var.set(self.config["SAVE_BBOX"])
		self.outputVerboseEntry = Checkbutton(self, variable=self.BBOXImage_var)
		self.outputVerboseEntry.grid(row=self.row_num, column=6,  sticky=W)
		Tooltip(self.BBOXImages, text='Select this option to produce an image wherein identified and classified larvae are'+
			'have a box placed around them. Boxes are color-coded for alive, dead, moribund, or egg.')

		self.row_num += 1

		#Row: Date and Plate
		self.DateLabel = Label(self, text="Date:")
		self.DateLabel.grid(row=self.row_num, column=0, sticky=E)
		self.config["DATE"] = self.config["PREV_DATE"]
		date_string = self.config["DATE"].strftime(self.config["PRETTY_DATE_FORMAT"])
		self.DateEntry = Entry(self, textvariable=date_string, width=self.textbox_width* self.scale, font = self.font)
		self.DateEntry.grid(row=self.row_num, column=1, columnspan = full_width, sticky=EW)
		self.DateEntry.insert(END, date_string)
		self.DateButton = ttk.Button(self, text ='Pick date', command = lambda: self.date_picker(), style = "my.TButton")
		self.DateButton.grid(row=self.row_num, column=full_width+1, sticky=(S,W))
		self.row_num += 1
		msg = "Enter the date that the plate was started. Dates may be entered in " + \
			"MM/DD/YY format or chosen with the date selector. May be left blank."
		Tooltip(self.DateLabel, text=msg)
		Tooltip(self.DateButton, text=msg)

		self.PlateIDLabel = Label(self, text="Plate ID:")
		self.PlateIDLabel.grid(row=self.row_num, column=0, sticky=E)
		vlist = [1,2,3,4,5,6,7,8,9,10]
		self.PlateIDButton = ttk.Spinbox(self, values = vlist, wrap = True, format = '%2.0f', width=2, font = self.font)
		self.PlateIDButton.grid(row=self.row_num, column=1, sticky = (S,W))
		self.PlateIDButton.insert(0,self.plateID)
		self.row_num += 1
		Tooltip(self.PlateIDLabel, text="ID of the plate. Set to 1 unless more than one plate is run in a day.")

	def reset_manual_indicators(self):
		self.csv_name_manual_mod = False
		self.date_manual_mod = False

	def date_picker(self, scale = 1):
		w = 400 * self.scale
		h = 400 * self.scale
		win = Toplevel()
		win.wm_title("Select Assay Date")
		pdu.win_center(win, w, h)
		cal = Calendar(win, selectmode = 'day', 
				   font = 'Arial 14 bold',
				   foreground = 'Black',
				   selectforeground = 'Red',
				   firstweekday = "sunday",
				   showweeknumbers = False,
				   date_pattern = 'mm/dd/yyyy') 
		cal.pack(pady = 20 * self.scale) 
		def grad_date(): 
			self.config["DATE"] = cal.get_date()
			self.update_date()
			win.destroy()
		Button(win, text = "Select Date", command = grad_date).pack(pady = 20 * self.scale)

	def update_date(self):
		self.DateEntry.delete(0,END)
		self.DateEntry.insert(0,datetime.strptime(self.config["DATE"], '%m/%d/%Y').strftime(self.config["PRETTY_DATE_FORMAT"]))
		self.date_manual_mod = True

	def image_directory_chooser(self):
		new_name = filedialog.askdirectory(initialdir=self.config['IMG_DIR'])
		if new_name:
			self.config['IMG_DIR'] = new_name
			self.imageFolderEntry.delete(0,END)
			self.imageFolderEntry.insert(0,self.config['IMG_DIR'])
			#update predicted date and csv filename 
			self.config['MOST_RECENT_IMG_DIR'] = new_name 
			csv_name, currdate, dirdate, yesterday = pdu.output_filename_formater(self.config)
			# if not self.csv_name_manual_mod:
			#     self.config['CSV_NAME'] = csv_name
			#     self.outputCSVEntry.delete(0,END)
			#     self.outputCSVEntry.insert(0,self.config['CSV_NAME'])
			if not self.date_manual_mod:
				self.config["DIR_DATE"] = dirdate
				self.DateEntry.delete(0,END)
				self.DateEntry.insert(0,self.config["DIR_DATE"].strftime(self.config["PRETTY_DATE_FORMAT"]))


	def image_output_chooser(self):
		new_name = filedialog.askdirectory(initialdir=self.config['OUT_DIR'])
		if new_name:
			self.config['OUT_DIR'] = new_name
			self.outputFolderEntry.delete(0,END)
			self.outputFolderEntry.insert(0,self.config['OUT_DIR'])

	def __check_files__(self):
		if not path.exists(self.config['IMG_DIR']): 
			pdu.message_window(msg="Input file does not exist. Please select a valid input file.")
			return False
		try: 
			file_test = open(os.path.join(self.config['OUT_DIR'], self.config['CSV_NAME']), 'a')
			file_test.close()
		except:
			pdu.message_window(msg="Output file cannot be written to the specified path. Please try again.")
			return False
		return True

	def update(self):
		# self.config["VERBOSE_CSV"] = self.verbose_var.get()
		self.config["SAVE_BBOX"] = self.BBOXImage_var.get()
		self.config["SAVE_INTMD"] = self.INTMD_var.get()
		self.config["SAVE_MASK"] = self.MASK_var.get()
		self.config["IMG_DIR"] = self.imagefolder_var.get()
		# self.config["CSV_NAME"] = self.CSV_Var.get()
		self.config["OUT_DIR"] = self.outputfolder_var.get()