#!/usr/bin/env python3
#datapreviewer.py

from tkinter import *
from tkinter.ttk import Label
import numpy as np
import tkinter as tk

from gui.tooltip import Tooltip
from stats.main import analyze_data
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from gui.utils import win_center
from itertools import compress
from stats.curvell import CI_finder

import math
from scipy.optimize import minimize
from scipy.special import gamma

class DataPreviewer(tk.Frame):
	def __init__(self, parent, merlin_stats_obj, scale = 1):
		super(DataPreviewer, self).__init__() #initialize root
		self.scale = scale
		w = 1000
		h = 1000
		self.top = Toplevel(parent)
		self.top.wm_title("Data Previewer")
		win_center(self.top, w, h)

		#pass on statistics object
		self.stats_obj = merlin_stats_obj

		self.columnconfigure(0, weight=2)
		self.columnconfigure(1, weight=2)
		self.columnconfigure(2, weight=2)
		self.__create_widgets()

	def __create_widgets(self, scale = 1):
		'''
		Makes the widget layout for the data previewer: 

		self.top 
		|--cmpd_frm #Compound selection for plotting
		|  |--cmpdlabel
		|  |--cmpd_list_frm w/ scrollbar #Inner frame 
		|  |  |--self.cmpd_listbox
		|--uid_frm #UID selection for highlighting in plot
		|  |--uidlabel
		|  |--uid_list_frm w/ scrollbar #Inner frame 
		|  |  |--self.uid_listbox
		|--self.fig_frame #The plot
		|  |--self.canvas
		|  |  |--self.fig

		'''

		#Compound Selection
		cmpd_frm = Frame(self.top)
		cmpd_frm.grid(row=0, column=0, rowspan = 5, columnspan = 2, 
						sticky=N+S, padx=10, pady=10)
		
		cmpdlabel = Label(cmpd_frm, text="Select Compound.") #Label
		cmpdlabel.grid(row=0, column=0, sticky=E)
		
		cmpd_list_frm = Frame(cmpd_frm) #inner frame
		cmpd_list_frm.grid(row= 1, column=0, rowspan = 5, columnspan = 2, 
						sticky=N+S, padx=10, pady=10)
		#Scrollbar for inner frame (will be bound to the listbox)
		scrollbar = Scrollbar(cmpd_list_frm, orient="vertical")
		scrollbar.pack(side=RIGHT, fill=Y)
		
		self.cmpd_listbox = Listbox(cmpd_list_frm, #Compound listbox
						selectmode='single',
						exportselection=0,
						yscrollcommand=scrollbar.set, 
						width = 30)
		self.cmpd_listbox.pack(expand = True, fill = Y) #pack/grid format gives tighter layout
		
		#bindings
		self.cmpd_listbox.bind('<<ListboxSelect>>', self.plotselect)
		self.cmpd_listbox.bind("<Down>", self.plotselect_scroll)
		self.cmpd_listbox.bind("<Up>", self.plotselect_scroll)
		scrollbar.config(command=self.cmpd_listbox.yview) #configure scrollbar to listbox
		
		#Populate the listbox with compound names.
		self.cmpd_list = sorted(list(self.stats_obj.cmpd_data.keys()))
		for ind, item in enumerate(self.cmpd_list):
			self.cmpd_listbox.insert(ind, item)
		Tooltip(self.cmpd_listbox, text="Select a compound to show the data that will be used for analysis.")
		
		#UID Selection
		uid_frm = Frame(self.top)
		uid_frm.grid(row=5, column=0, rowspan = 5, columnspan = 2, 
						sticky=N+S, padx=10, pady=10)
		
		uidlabel = Label(uid_frm, text="Identifiers") #Label
		uidlabel.grid(row=0, column=0, sticky=E)

		uid_list_frm = Frame(uid_frm) #Inner frame
		uid_list_frm.grid(row=1, column=0, rowspan = 5, columnspan = 2, 
						sticky=N+S, padx=10, pady=10)
		#Scrollbar for inner frame (will be bound to the listbox)
		scrollbar = Scrollbar(uid_list_frm, orient="vertical")
		scrollbar.pack(side=RIGHT, fill=Y)
		
		#UID listbox
		self.uid_listbox = Listbox(uid_list_frm, 
						selectmode='single',
						exportselection=0,
						yscrollcommand=scrollbar.set, 
						width = 30)
		self.uid_listbox.pack(expand = True, fill = Y)
		
		#bindings
		self.uid_listbox.bind('<<ListboxSelect>>', self.uidselect)
		self.uid_listbox.bind("<Down>", self.uidselect_scroll)
		self.uid_listbox.bind("<Up>", self.uidselect_scroll)
		scrollbar.config(command=self.uid_listbox.yview) #bind scrollbar
		Tooltip(self.uid_listbox, text="Identifier is of the format <date>_P<plate>_R<row>_ID. " +\
			"Select an identifier to highlight its data in the plot.")

		'''
		Figure frame
		Making the figure and canvas here is important. When made in the plot function, a
		memory leak results due to improper clearing of the plot. 
		'''
		self.fig_frame = Frame(self.top)
		self.fig_frame.grid(row=0, column=3, rowspan = 9, columnspan = 9, sticky=N+S)
		self.fig = Figure(figsize = (5,5), dpi = 150)
		self.canvas = FigureCanvasTkAgg(self.fig, master=self.fig_frame)
		self.canvas.get_tk_widget().grid(row=0, column=2, rowspan = 8, columnspan = 8, sticky=N+S)
		self.plot1 = self.fig.add_subplot(111)

		#initialize a plot with the first compound. 
		self.plotselect()

	def scroll(self, event):
		'''
		Based on https://stackoverflow.com/a/62711725/8075803
		Controls actions on up/down scrolling. 
		'''
		selection = event.widget.curselection()[0]
		if event.keysym == 'Up': selection += -1
		if event.keysym == 'Down': selection += 1
		if 0 <= selection < event.widget.size():
			event.widget.selection_clear(0, tk.END)
			event.widget.select_set(selection)

	def uidselect_scroll(self, event):
		#For arrow-based scrolling in the UID list.
		self.scroll(event)
		self.uidselect(event)

	def plotselect_scroll(self, event):
		#For arrow-based scrolling in the Compound list.
		self.scroll(event)
		self.plotselect(event)

	def plotselect(self, event = None):
		# Based on https://stackoverflow.com/a/12936031/8075803
		if event is not None:
			w = event.widget
			self.current_cmpd = self.cmpd_list[int(w.curselection()[0])]
		else:
			self.current_cmpd = self.cmpd_list[0]
		self.include_now = self.stats_obj.cmpd_data[self.current_cmpd].data['include_now'] #Get the allowed data list
		self.uid_list = list(self.stats_obj.cmpd_data[self.current_cmpd].data['unique_plate_ids']) #Get list of permitted IDs
		#Rescructure the UID to be moreuser-friendly
		self.uid_list = [s.split("_") for s in self.uid_list] 
		self.uid_list = [f"{s[4]}_Plate{s[1]}_Row{s[2]}_{s[3]}" for s in self.uid_list]
		self.uid_list = sorted(list(set(compress(self.uid_list, self.include_now))))
		
		#clear previous list
		self.uid_listbox.delete(0,END)
		for ind, item in enumerate(self.uid_list):
			self.uid_listbox.insert(ind, item)
		#update the plot
		self.plot()

	def uidselect(self, event = None):
		#Highlights the data from the currently-selected data. 
		self.current_uid = self.uid_list[int(event.widget.curselection()[0])]
		self.update_plot()

	def update_plot(self):
		'''
		Uses the current conc/probs from the selected compound to highlight the data from the 
		currently-selected UID. Note that this method actually creates a completely new plot;
		clearing the old plot is performed to prevent any memory leaks. 
		'''
		self.plot1.clear()

		#Split the concentrations and probabilities into two lists. 
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
		#Plot non-highlighted data
		self.plot1.plot(background_conc, 
				background_probs, 
				marker = '.', 
				mew = 0.0, 
				mfc = 'black', 
				ls = 'None')
		#plot highlighted data
		self.plot1.plot(uid_conc, 
				uid_probs, 
				marker = 'o', 
				mew = 0.0, 
				mfc = 'red', 
				ls = 'None')
		self.set_labels() #update axes, etc.
		self.canvas.draw() #draw the plot

	def plot(self):
		#Plots the data for the selected compound. Initially, no data is highlighted. 
		self.plot1.clear()
		self.curr_data = self.stats_obj.cmpd_data[self.current_cmpd]
		self.conc = self.curr_data.data["conc"].copy() #copied so that we don't change the original data. 
		#calculate survival probability. 
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

	# def make_curve(self):


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