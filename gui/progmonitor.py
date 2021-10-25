# gui/progmonitor.py
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
from gui.imageframe.ioframe import IOFrame as IOFrame
from gui.imageframe.cmpdframe import CmpdFrame as CmpdFrame
from threading import Thread

class ProgMonitor(ttk.Frame):
	def __init__(self, container, run_button, scale = 1, **kwargs):
		super().__init__(container, **kwargs)
		self.container = container
		self.scale = scale
		self.run_button = run_button
		self.__create_widgets()

	def __create_widgets(self):
		'''
		Creates the widgets for the progress monitoring class. 

		self.container
		self
		|--self.prog_bar (Progress Bar) #Shows the progress of the process
		|			with automatic movement included.
		|--self.prog_label (Label) #Provides messages from the class being
		|			monitored (obj.message).
		|--self.prog_percent (Label) #Provides a percent completion. This
					comes from the class being monitored (obj.progress)
		'''
		self.auto_prog_bar_on = False
		self.prog_bar = ttk.Progressbar(
			self, orient="horizontal",
			length=600, mode="determinate")
		self.prog_bar.grid(column=0, row=5, padx=20, pady=0)
		self.prog_label = Label(self, text="")
		self.prog_label.grid(row=4, column=0, padx=20, pady=0, sticky=EW)
		self.prog_percent = Label(self, text="")
		self.prog_percent.grid(row=6, column=0, padx=20, pady=0, sticky=EW)

	def monitor(self, thread, obj, **kwargs):
		if thread.is_alive():
			#update every 0.25 seconds
			self.after(250, lambda: self.monitor(thread, obj, **kwargs))
			self.__update_progress(obj, **kwargs)
			
		else:
			if obj.error_message:
				self.prog_label.config(text=obj.error_message)
				self.prog_bar['value'] == 0
			else:
				msg = "Done! You may close this window or process another "+\
							"plate."
				self.prog_label.config(text=msg)
				self.prog_bar['value'] = 100
			self.run_button['state'] = tk.NORMAL
			self.container.container.tab(0, state="normal")
			self.container.container.tab(1, state="normal")
			self.prog_percent.config(text="")

	def __update_progress(self, obj, speed = 0.15, **kwargs):
		'''
		Updates the progress by monitoring the object.
		obj: class object to be monitored.
		speed: rate at which to automatically move the progress bar. 
		'''
		self.prog_label.config(text=obj.message)
		if self.auto_prog_bar_on: 
			if self.prog_bar['value'] < obj.progress - 15:
				self.prog_bar['value'] = obj.progress - 10
			elif self.prog_bar['value'] < obj.progress:
				max_step = min(speed, obj.progress - self.prog_bar['value'])
				self.prog_bar.step(max_step)
			else:
				self.prog_bar['value'] = obj.progress
		self.prog_percent.config(text=f"{round(self.prog_bar['value']):n}%")