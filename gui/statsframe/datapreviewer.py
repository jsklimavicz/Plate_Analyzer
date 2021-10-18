#!/usr/bin/env python3
#datapreviewer.py

import tkinter as tk
from tkinter import *
from tkinter import ttk, filedialog, messagebox
from tkinter.ttk import Label, Style
import random
import os
from os import path, getcwd, mkdir
from os.path import exists
from datetime import datetime

from gui.tooltip import Tooltip #as Tooltip
from stats.main import analyze_data
import matplotlib

class DataPreviewer():
	def __init__(title="", msg=""):
    w = 600
    h = 600
    win = Toplevel()
    win.wm_title(title)
    win_center(win, w, h)
    l = ttk.Label(win, text=msg, wraplength=260)
    l.place(relx=.5, rely=.4, anchor="center")
    b = ttk.Button(win, text="Okay", command=win.destroy)
    b.place(relx=.5, rely=.65, anchor="center")