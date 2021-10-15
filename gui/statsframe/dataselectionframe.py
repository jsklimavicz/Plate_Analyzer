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

class DataSelectionFrame(ttk.Frame):
	def __init__(self, container, config, scale = 1, **kwargs):
		super().__init__(container, **kwargs)
		self.scale = scale
		self.config = config
		self.textbox_width=40
		self.font = ('Arial', 12*self.scale)
		super().__init__(container)
		self.selection_ID = selection_ID
		# setup the grid layout manager
		self.columnconfigure(0, weight=3)
		self.columnconfigure(1, weight=1)
		self.columnconfigure(2, weight=3)
		self.__create_widgets()

	def __create_widgets(self):
		s = ttk.Style()
		s.configure("my.TButton", font = self.font)
		full_width = 2

		self.allowed_list = Listbox(self,selectmode='multiple',exportselection=0, height=10)
		# selected_text_list = [listbox.get(i) for i in listbox.curselection()]