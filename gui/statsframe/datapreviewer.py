#!/usr/bin/env python3
#datapreviewer.py


from tkinter import *
from tkinter import ttk, filedialog, messagebox
from tkinter.ttk import Label, Style
import random
import os
from os import path, getcwd, mkdir
from os.path import exists
from datetime import datetime
import numpy as np
import tkinter as tk

from gui.tooltip import Tooltip #as Tooltip
from stats.main import analyze_data
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from gui.utils import win_center
from itertools import compress

class DataPreviewer(tk.Frame):
	def __init__(self, parent, merlin_stats_obj, scale = 1):
		super(DataPreviewer, self).__init__()
		self.scale = scale
		w = 1000
		h = 1000
		self.top = Toplevel(parent)
		self.top.wm_title("Data Previewer")
		win_center(self.top, w, h)
		
		self.stats_obj = merlin_stats_obj

		self.columnconfigure(0, weight=2)
		self.columnconfigure(1, weight=2)
		self.columnconfigure(2, weight=2)
		self.__create_widgets()

	def __create_widgets(self, scale = 1):
		self.font = ('Arial', 12*self.scale)
		

		#cmpd List
		cmpd_frm = Frame(self.top)
		cmpd_frm.grid(row=0, column=0, rowspan = 5, columnspan = 2, sticky=N+S, padx=10, pady=10)
		exclusionlabel = Label(cmpd_frm, text="Select Compound.")
		exclusionlabel.grid(row=0, column=0, sticky=E)
		
		cmpd_list_frm = Frame(cmpd_frm)
		cmpd_list_frm.grid(row= 1, column=0, rowspan = 5, columnspan = 2, sticky=N+S, padx=10, pady=10)
		scrollbar = Scrollbar(cmpd_list_frm, orient="vertical")
		scrollbar.pack(side=RIGHT, fill=Y)
		self.cmpd_listbox = Listbox(cmpd_list_frm, selectmode='single',exportselection=0,yscrollcommand=scrollbar.set, width = 30)
		self.cmpd_listbox.pack(expand = True, fill = Y)
		self.cmpd_listbox.bind('<<ListboxSelect>>', self.plotselect)
		scrollbar.config(command=self.cmpd_listbox.yview)
		self.cmpd_list = sorted(list(self.stats_obj.cmpd_data.keys()))
		for ind, item in enumerate(self.cmpd_list):
			self.cmpd_listbox.insert(ind, item)
		Tooltip(self.cmpd_listbox, text="Select a compound to show the data that will be used for analysis.")
		
		#UID List
		uid_frm = Frame(self.top)
		uid_frm.grid(row=5, column=0, rowspan = 5, columnspan = 2, sticky=N+S, padx=10, pady=10)
		exclusionlabel = Label(uid_frm, text="Identifiers")
		exclusionlabel.grid(row=0, column=0, sticky=E)
		uid_list_frm = Frame(uid_frm)
		uid_list_frm.grid(row=1, column=0, rowspan = 5, columnspan = 2, sticky=N+S, padx=10, pady=10)
		scrollbar = Scrollbar(uid_list_frm, orient="vertical")
		scrollbar.pack(side=RIGHT, fill=Y)
		self.uid_listbox = Listbox(uid_list_frm, selectmode='single',exportselection=0,yscrollcommand=scrollbar.set, width = 30)
		self.uid_listbox.pack(expand = True, fill = Y)
		self.uid_listbox.bind('<<ListboxSelect>>', self.uidselect)
		scrollbar.config(command=self.uid_listbox.yview)
		Tooltip(self.uid_listbox, text="Identifier is of the format Cmpd_ID_Date_Plate_Row. Select an identifier to highlight its data in the plot.")

		#Figure frame
		#Making the figure and canvas here is important. When made in the plot function, a
		#memory leak results due to improper clearing of the plot. 
		self.fig_frame = Frame(self.top)
		self.fig_frame.grid(row=0, column=0, rowspan = 9, columnspan = 9, sticky=N+S)
		self.fig = Figure(figsize = (5,5), dpi = 150)
		self.canvas = FigureCanvasTkAgg(self.fig, master=self.top)
		self.canvas.get_tk_widget().grid(row=0, column=2, rowspan = 8, columnspan = 8, sticky=N+S)
		self.plot1 = self.fig.add_subplot(111)

		#initialize a plot with the first compound. 
		self.plotselect()

		# b = ttk.Button(self.top, text="Okay", command=self.top.destroy)
		# b.grid(column=0, row = 6)

	def plotselect(self, event = None):
		# Note here that Tkinter passes an event object to onselect()
		# See https://stackoverflow.com/a/12936031/8075803
		if event is not None:
			w = event.widget
			self.current_cmpd = self.cmpd_list[int(w.curselection()[0])]
		else:
			self.current_cmpd = self.cmpd_list[0]
		self.include_now = self.stats_obj.cmpd_data[self.current_cmpd].data['include_now']
		self.uid_list = sorted(list(self.stats_obj.cmpd_data[self.current_cmpd].data['unique_plate_ids']))
		self.uid_list = list(set(compress(self.uid_list, self.include_now)))
		
		#clear previous list
		self.uid_listbox.delete(0,END)
		for ind, item in enumerate(self.uid_list):
			self.uid_listbox.insert(ind, item)

		self.plot()

	def uidselect(self, event = None):
		self.current_uid = self.uid_list[int(event.widget.curselection()[0])]
		self.update_plot()


	def update_plot(self):
		self.plot1.clear()

		background_conc = []
		background_probs = []
		uid_conc = []
		uid_probs = []

		for ind, uid in enumerate(self.stats_obj.cmpd_data[self.current_cmpd].data['unique_plate_ids']):
			if uid == self.current_uid:
				uid_conc.append(self.conc[ind])
				uid_probs.append(self.probs[ind])
			else:
				background_conc.append(self.conc[ind])
				background_probs.append(self.probs[ind])

		self.plot1.plot(background_conc, 
				background_probs, 
				marker = '.', 
				mew = 0.0, 
				mfc = 'black', 
				ls = 'None')
		self.plot1.plot(uid_conc, 
				uid_probs, 
				marker = 'o', 
				mew = 0.0, 
				mfc = 'red', 
				ls = 'None')
		self.set_labels()
		self.canvas.draw()

	def plot(self):
		#Plots the data for the selected compound. No curve is drawn. 

		self.plot1.clear()
		self.curr_data = self.stats_obj.cmpd_data[self.current_cmpd]
		self.conc = self.curr_data.data["conc"].copy()
		self.probs = self.curr_data.data["live_count"]/(self.curr_data.data["live_count"] + self.curr_data.data["dead_count"])
		if self.stats_obj.options["JITTER"]: self.conc += np.random.uniform(-self.stats_obj.options["JITTER_FACTOR"], 
													self.stats_obj.options["JITTER_FACTOR"], 
													len(self.conc))
		self.plot1.plot(self.conc, 
				self.probs, 
				marker = '.', 
				mew = 0.0, 
				mfc = 'black', 
				ls = 'None')
		self.set_labels()
		self.canvas.draw()

	def set_labels(self):
		'''
		Set axes labels, ticks, title, etc. 
		'''

		def calc_x_ticks():
			'''
			Calculate the x-ticks for the graph based on the range of concentrations in the data. 
			'''
			lb, ub = round(min(self.conc)), round(max(self.conc))
			xticks = np.array(range(lb-1, ub+2, 1))
			xticklabels = [round(2**i) if i >=0 else 2.**i for i in xticks]
			for i, label in enumerate(xticklabels):
				if label < 1./64.: xticklabels[i] = "{:.1e}".format(label).replace("e-0", "e-")
				elif label < 1./4.: xticklabels[i] = "{:.3f}".format(label)
				else: xticklabels[i] = str(label) 
			return xticks, xticklabels

		self.plot1.set_xlabel('Concentration (ppm)')
		self.plot1.set_title(label = self.current_cmpd)
		self.plot1.set_ylabel('Percent Survival')
		xticks, xticklabels = calc_x_ticks()
		self.plot1.set_xticks(xticks, minor=False)
		self.plot1.set_xticklabels(xticklabels, rotation=90, fontsize=8)
		# plot1.xticks(ticks=xticks, rotation=90, fontsize=8, labels=xticklabels)
		# plot1.yticks(ticks=np.array(range(0, 5, 1))/4. , labels=np.array(range(0, 101, 25)) )
		self.plot1.set_yticks(np.array(range(0, 5, 1))/4., minor=False)
		self.plot1.set_yticklabels(np.array(range(0, 101, 25)), rotation=90, fontsize=8)
		self.plot1.grid(b=True, alpha = 0.5)
		self.plot1.set_xlim([min(xticks), max(xticks)])
		self.plot1.set_ylim([0,1])

		self.plot1.set_position([.125, .16, .8,.73])
		plt.subplots_adjust(bottom=0.75)

