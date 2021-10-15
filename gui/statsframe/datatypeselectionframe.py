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
from platedriver.utils import Tooltip #as Tooltip
import threading
import platedriver.utils as pdu
import platform
from threading import Thread
import time
from stats.main import analyze_data

class DataTypeSelectionFrame(ttk.Frame):
	def __init__(self, container, config, selection_ID, scale = 1, **kwargs):
		super().__init__(container, **kwargs)
		self.scale = scale
		self.config = config
		self.textbox_width=60
		self.font = ('Arial', 12*self.scale)
		super().__init__(container)
		self.selection_ID = selection_ID
		# setup the grid layout manager
		self.columnconfigure(0, weight=1)
		self.columnconfigure(1, weight=2)
		self.columnconfigure(2, weight=2)
		self.__create_widgets()

	def __create_widgets(self):
		s = ttk.Style()
		s.configure("my.TButton", font = self.font)

		full_width = 2
		#Row: Input Image Folder Selector
		self.exclusionlabel = Label(self, text="Exclude runs by:")
		self.exclusionlabel.grid(row=0, column=0, sticky=E)
		self.exclusiontype = StringVar()
		self.exclusiontype.set(self.selection_ID[0])
		self.exclusionEntry = OptionMenu(self, 
							self.exclusiontype, 
							self.selection_ID[0], 
							*self.selection_ID)
		self.exclusionEntry.grid(row=0, column=1, columnspan = full_width,  sticky=EW)
