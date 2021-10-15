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

class DataTypeSelectionFrame(ttk.Frame):
	def __init__(self, container, config, stats_obj, scale = 1, **kwargs):
		super().__init__(container, **kwargs)
		self.scale = scale
		self.config = config
		self.textbox_width=60
		self.font = ('Arial', 12*self.scale)
		super().__init__(container)
		self.stats_obj = stats_obj

		self.uid_list = self.stats_obj.get_uid_list()

		self.make_var_lists()
		# setup the grid layout manager
		self.columnconfigure(0, weight=2)
		self.columnconfigure(1, weight=2)
		self.columnconfigure(2, weight=2)
		self.columnconfigure(3, weight=2)
		self.columnconfigure(4, weight=2)
		self.columnconfigure(5, weight=2)

		self.selection_options = ["Compound", "Reference ID", "Date", "Plate ID", "Row ID"]
		self.forbidden_list = []

		self.__create_widgets()

	def __create_widgets(self):
		s = ttk.Style()
		s.configure("my.TButton", font = self.font)

		full_width = 2
		#Row: Input Image Folder Selector
		self.exclusionlabel = Label(self, text="Exclude runs by:")
		self.exclusionlabel.grid(row=0, column=0, sticky=E)
		self.exclusiontype = StringVar()
		self.exclusiontype.set(self.selection_options[0])
		
		self.exclusionEntry = OptionMenu(self, 
							self.exclusiontype, 
							*self.selection_options,
							command=lambda _: self.list_update())
		self.exclusionEntry.grid(row=0, column=1, columnspan = full_width,  sticky=EW, padx=10, pady=20)

		self.allowed_list = Listbox(self, selectmode='multiple',exportselection=0, height=10)
		self.allowed_list.grid(row=1, column=0, columnspan = 2, rowspan = 4, sticky=EW, padx=10, pady=20)
		self.disallowed_list = Listbox(self, selectmode='multiple',exportselection=0, height=10)
		self.disallowed_list.grid(row=1, column=4, columnspan = 2, rowspan = 4, sticky=EW, padx=10, pady=20)
		
		self.disallowButton = Button(self, text="Remove from sample >", command = lambda: self.disallow())
		self.disallowButton.grid(row=2, column=2, sticky=S, padx=10, pady=20)
		self.disallowButton.config(height = 2)

		self.allowButton = Button(self, text="< Add to sample", command = lambda: self.allow())
		self.allowButton.grid(row=3, column=2, sticky=S, padx=10, pady=20)
		self.allowButton.config(height = 2)

		self.list_update()

	def disallow(self): 
		for i in self.allowed_list.curselection():
			self.forbidden_list.append(self.allowed_ids[i])
		self.list_update()

	def allow(self): 
		for i in self.disallowed_list.curselection():
			self.forbidden_list.remove(self.disallowed_ids[i])
		self.list_update()

	def list_update(self):
		var = self.exclusiontype.get()
		if var == "Compound": self.var_list = self.by_name
		elif var == "Reference ID": self.var_list = self.by_id
		elif var == "Date": self.var_list = self.by_date
		elif var == "Plate ID": self.var_list = self.by_plate
		elif var == "Row ID": self.var_list = self.by_row

		#remove all data
		self.allowed_ids = []
		self.allowed_list.delete(0,END)
		self.disallowed_ids = [] 
		self.disallowed_list.delete(0,END)

		#make allowed/disallowed lists
		for name in self.var_list:
			if name in self.forbidden_list:
				self.disallowed_ids.append(name)
			else:
				self.allowed_ids.append(name)
		#populate listBoxes
		for ind, item in enumerate(self.allowed_ids):
			self.allowed_list.insert(ind, item)
		for ind, item in enumerate(self.disallowed_ids):
			self.disallowed_list.insert(ind, item)


	def make_var_lists(self):
		UID_breakdown = [uid.split("_") for uid in self.uid_list]
		self.by_name = sorted(list(set([a[0]for a in UID_breakdown])))
		self.by_date = sorted(list(set([a[1]for a in UID_breakdown])))
		self.by_id = sorted(list(set([a[4]for a in UID_breakdown])))
		self.by_plate = sorted(list(set([f'{a[1]}_Plate_{a[2]}' for a in UID_breakdown])))
		self.by_row = sorted(list(set([f'{a[1]}_Plate_{a[2]}_Row_{a[3]}' for a in UID_breakdown])))

	def get_disallowed_uids(self):
		disallowed = []
		for uid in self.uid_list:
			for item in self.forbidden_list:
				if "Row" in item:
					item = item.split("_")
					item = f"{item[0]}_{item[2]}_{item[4]}"
				elif "Plate" in item:
					item = item.split("_")
					item = f"{item[0]}_{item[2]}"
				if item in uid: disallowed.append(uid)
		return disallowed

