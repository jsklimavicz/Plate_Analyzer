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

class DataSelectionFrame(ttk.Frame):
	def __init__(self, container, config, scale = 1, **kwargs):
		super().__init__(container, **kwargs)
		self.scale = scale
		self.config = config
		self.__create_widgets()
		self.auto_prog_bar_on = False
		self.plate = None
		self.__progress()