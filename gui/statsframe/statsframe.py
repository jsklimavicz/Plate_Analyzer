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
from gui.statsframe.dataselectionframe import DataSelectionFrame as DSF
from gui.statsframe.datatypeselectionframe import DataTypeSelectionFrame as DTSF

class StatsFrame(ttk.Frame):
	def __init__(self, container, config, scale = 1, **kwargs):
		super().__init__(container, **kwargs)
		self.scale = scale
		self.config = config
		self.selection_options = ["Compound", "Reference ID", "Date", "Plate ID", "Row ID"]
		self.__create_widgets()
		# self.auto_prog_bar_on = False
		# self.plate = None
		# self.__progress()
		
	def __create_widgets(self, scale = 1):

		self.font = ('Arial', 12*self.scale)
		self.selection_frame = DTSF(self, config = self.config, selection_ID = self.selection_options, scale = self.scale)
		self.selection_frame.grid(column=0, row=0, padx=10, pady=20)