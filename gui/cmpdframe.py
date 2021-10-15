#!/usr/bin/env python3
#cmpdframe.py

from tkinter import *
from tkinter import ttk, filedialog, messagebox
from tkinter.ttk import Label, Style
import random
import os
from os import path, getcwd, mkdir
from os.path import exists
import platedriver.plate as pc
from datetime import datetime
import tkinter as tk
from platedriver.utils import Tooltip #as Tooltip
import platform
import time
from stats.main import analyze_data

#############################################################################
#                         For entering compound data
#############################################################################

class CmpdFrame(ttk.LabelFrame):
    class CmpdInfo:
        def __init__(self , config):
            self.row_names = []
            self.cmpd_names = []
            self.cmpd_codes = []
            self.concentrations = []
            self.reps = []

        def get_concentration(self, row):
            for i in range(0, len(self.row_names)):
                if self.row_names[i] == row: 
                    try: 
                        return self.concentrations[i].get()
                    except ValueError:
                        return float('NaN')
            return float('NaN')

        def get_name(self, row):
            for i in range(0, len(self.row_names)):
                if self.row_names[i] == row: return self.cmpd_names[i].get()
            return ""

        def get_code(self, row):
            for i in range(0, len(self.row_names)):
                if self.row_names[i] == row: return self.cmpd_codes[i].get()
            return ""

        def get_rep(self, row):
            for i in range(0, len(self.row_names)):
                if self.row_names[i] == row: return self.reps[i]
            return ""

    def __init__(self, container, scale = 1, config = None):
        self.scale = scale
        # container.tk.call('tk', 'scaling', '-displayof', '.', 2.0)
        super().__init__(container)
        self.config = config
        self.style = Style(container)
        self.font = ('Arial', 12*self.scale)
        self.style.configure("TLabel", font = self.font)
        # setup the grid layout manager
        self.columnconfigure(0, weight=1)
        self.grid()
        self.textBox = []
        self.cmpd_name =[]
        self.cmpd_code =[]
        self.concentration = []
        self.rep = []
        self.param = self.CmpdInfo(config)
        self.__create_widgets()

    def __create_widgets(self):
        #Column titles
        self.CodeLabel = ttk.Label(self, text="Compound Code")
        self.CodeLabel.grid(row=0, column=1, sticky = S)
        Tooltip(self.CodeLabel, text="Code name/number of compounds. May be left blank.")
        self.NameLabel = ttk.Label(self, text="Compound Name")
        self.NameLabel.grid(row=0, column=2, sticky = S)
        Tooltip(self.NameLabel, text="Name of treatment. May be left blank.")
        self.MaxConcLabel = ttk.Label(self, text="Max Conc.")
        self.MaxConcLabel.grid(row=0, column=3, sticky = S)
        Tooltip(self.MaxConcLabel, text="Highest well concentration of treatment, in ppm. May be left blank")
        self.RepLabel = ttk.Label(self, text="Rep.")
        self.RepLabel.grid(row=0, column=4, sticky = (S,W))

        for i in range(1,9): self.__create_cmpd_entry__(i)
       
        self.ControlLabel = ttk.Label(self, text=f"Control:")
        self.ControlLabel.grid(row=9, column=0, sticky = E)

        self.ControlCodeEntry = ttk.Entry(self, textvariable=self.config["CTRL_CODE"], width=15 * self.scale, font = self.font)
        self.ControlCodeEntry.grid(row=9, column=1, sticky = W)
        self.ControlCodeEntry.insert(END, self.config["CTRL_CODE"])
        Tooltip(self.ControlCodeEntry, text='Name of liquid media type (e.g. LM 2.0).')

        self.ControlNameEntry = ttk.Entry(self, textvariable=self.config["CTRL_NAME"], width=25 * self.scale, font = self.font)
        self.ControlNameEntry.grid(row=9, column=2, sticky = W)
        self.ControlNameEntry.insert(END, self.config["CTRL_NAME"])
        Tooltip(self.ControlNameEntry, text='Name of the control. This should almost always just be "control".')

        # self.button2 = ttk.Button(self, text="Convert", command = lambda: basic_info.convert()).grid(row=10, column=3, sticky=W)


    def __create_cmpd_entry__(self, row_number):
        vlist = [1,2,3,4,5,6,7,8,9,10]
        row_name = chr(row_number + 64)
        self.param.row_names.append(row_name)
        self.param.cmpd_names.append(StringVar())
        self.param.cmpd_codes.append(StringVar())
        self.param.reps.append(1)
        self.param.concentrations.append(DoubleVar())
        self.textBox.append(ttk.Label(self, text=f"Row {row_name:s}:").grid(row=row_number, column=0, sticky = E))
        self.cmpd_code.append(ttk.Entry(self, textvariable=self.param.cmpd_codes[row_number-1], width=15 * self.scale, font = self.font))
        self.cmpd_code[row_number-1].grid(row=row_number, column=1, sticky = W)
        self.cmpd_name.append(ttk.Entry(self, textvariable=self.param.cmpd_names[row_number-1], width=25 * self.scale, font = self.font))
        self.cmpd_name[row_number-1].grid(row=row_number, column=2, sticky = W)
        self.concentration.append(ttk.Entry(self, textvariable=self.param.concentrations[row_number-1], width=10 * self.scale, font = self.font))
        self.concentration[row_number-1].grid(row=row_number, column=3, sticky = W)
        self.rep.append(ttk.Spinbox(self, values = vlist, wrap = True, format = '%2.0f', width=2, font = self.font))
        self.rep[row_number-1].grid(row=row_number, column=4, sticky = W)
        self.rep[row_number-1].insert(0,1)

    def update(self):
        for i in range(0, len(self.param.reps)):
            self.param.reps[i] = int(self.rep[i].get())
        self.config["CTRL_NAME"] = self.ControlNameEntry.get()
        self.config["CTRL_CODE"] = self.ControlCodeEntry.get()
