#!/usr/bin/env python3
#datapreviewer.py

from tkinter import *
import tkinter as tk
from tkinter.ttk import Label

import math
from scipy.optimize import minimize, least_squares
from itertools import compress
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)

from stats.curvell import CI_finder
from gui.tooltip import Tooltip
from gui.utils import win_center

class DataPreviewer(tk.Frame):
	'''
	Provides an interactive interface to preview the data that will be 
	included in the statistical anlysis.
	'''
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
		|  |  |--self.curve_data_label

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
		#pack/grid format gives tighter layout
		self.cmpd_listbox.pack(expand = True, fill = Y) 
		
		#bindings
		self.cmpd_listbox.bind('<<ListboxSelect>>', self.plotselect)
		self.cmpd_listbox.bind("<Down>", self.plotselect_scroll)
		self.cmpd_listbox.bind("<Up>", self.plotselect_scroll)
		#configure scrollbar to listbox
		scrollbar.config(command=self.cmpd_listbox.yview) 
		
		#Populate the listbox with compound names.
		self.cmpd_list = sorted(list(self.stats_obj.cmpd_data.keys()))
		for ind, item in enumerate(self.cmpd_list):
			self.cmpd_listbox.insert(ind, item)
		msg = "Select a compound to preview data selected for analysis."
		Tooltip(self.cmpd_listbox, text=msg)
		
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
		msg = "Identifier is of the format <date>_P<plate>_R<row>_ID. " + \
				"Select an identifier to highlight its data in the plot."
		Tooltip(self.uid_listbox, text= msg)

		'''
		Figure frame
		Making the figure and canvas here is important. When made in the plot
		function, a memory leak results due to improper clearing of the plot. 
		'''
		self.fig_frame = Frame(self.top)
		self.fig_frame.grid(row=0, column=3, rowspan = 9, columnspan = 9, 
					sticky=N+S)
		self.fig = Figure(figsize = (5,5), dpi = 150)
		self.canvas = FigureCanvasTkAgg(self.fig, master=self.fig_frame)
		self.canvas.get_tk_widget().grid(row=0, column=3, rowspan = 8, 
					columnspan = 8, sticky=N+S)
		self.plot1 = self.fig.add_subplot(111)
		self.curve_data_label = Label(self.fig_frame, text="") #Label
		self.curve_data_label.grid(row=9, column=3, sticky=E)
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
			return True
		return False

	def uidselect_scroll(self, event):
		#For arrow-based scrolling in the UID list.
		reselect = self.scroll(event)
		if reselect: self.uidselect(event)

	def plotselect_scroll(self, event):
		#For arrow-based scrolling in the Compound list.
		reselect = self.scroll(event)
		if reselect: self.plotselect(event)

	def plotselect(self, event = None):
		# Based on https://stackoverflow.com/a/12936031/8075803
		if event is not None:
			w = event.widget
			self.current_cmpd = self.cmpd_list[int(w.curselection()[0])]
		else:
			self.current_cmpd = self.cmpd_list[0]
		#Get the allowed data list
		self.include_now = self.stats_obj.cmpd_data[
							self.current_cmpd].data['include_now']
		#Get list of permitted IDs 
		self.uid_list = list(self.stats_obj.cmpd_data[
							self.current_cmpd].data['unique_plate_ids']) 
		self.uid_list = sorted(list(set(compress(
							self.uid_list, self.include_now))))
		
		#clear previous list
		self.uid_listbox.delete(0,END)
		for ind, item in enumerate(self.uid_list):
			#Rescructure the UID to be moreuser-friendly
			s = item.split("_")
			item = f"{s[4]}_{s[1]}_Plate{s[2]}_Row{s[3]}"
			self.uid_listbox.insert(ind, item)
		#update the plot
		self.plot()

	def uidselect(self, event = None):
		#Highlights the data from the currently-selected data. 
		self.current_uid = self.uid_list[int(event.widget.curselection()[0])]
		self.update_plot()

	def update_plot(self):
		'''
		Uses the current conc/probs from the selected compound to highlight 
		the data from the currently-selected UID. Note that this method 
		actually creates a completely new plot; clearing the old plot is 
		performed to prevent any memory leaks. 
		'''
		self.plot1.clear()

		#Split the concentrations and probabilities into two lists. 
		background_conc = []
		background_probs = []
		uid_conc = []
		uid_probs = []
		for ind, uid in enumerate(self.stats_obj.cmpd_data[
						self.current_cmpd].data['unique_plate_ids']):
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
		self.plot_curve()
		self.canvas.draw() #draw the plot

	def plot(self):
		'''
		Plots the data for the selected compound. Initially, no data is 
		highlighted. 
		'''
		self.plot1.clear()
		self.curr_data = self.stats_obj.cmpd_data[self.current_cmpd]
		#copy concentrations so that we don't change the original data. 
		self.conc_orig = self.curr_data.data["conc"].copy() 
		self.conc = self.conc_orig.copy() 
		#calculate survival probability. 
		self.probs = np.array(self.curr_data.data["live_count"]/ \
					(self.curr_data.data["live_count"] + \
					self.curr_data.data["dead_count"]))
		if self.stats_obj.options["JITTER"]: 
			self.conc += np.random.uniform(
					-self.stats_obj.options["JITTER_FACTOR"],
					self.stats_obj.options["JITTER_FACTOR"], 
					len(self.conc))
		self.plot1.plot(self.conc, 
				self.probs, 
				marker = '.', 
				mew = 0.0, 
				mfc = 'black', 
				ls = 'None')
		self.set_labels()
		self.make_curve()
		self.canvas.draw()
		

	def make_curve(self):
		'''
		Generates the x and y values to produce a best-fit ll3 curve. The
		curve is fit using the maximum likelihood function in CI_Finder.ll3.
		The curve is then used to generate an approximate LC50 and an R2. 
		'''
		lb, ub = round(min(self.conc_orig)), round(max(self.conc_orig))
		num_points = 10 * (ub - lb + 2) + 1 #10 x-points (including both ends)
		self.x = np.linspace(lb-1, ub+1, num_points)
		
		#Find the best fit based on the number of parameters and curve-fitting
		#method specified in the options. 
		curve_type = self.stats_obj.options["CURVE_TYPE"].lower()
		params = 2 if "2" in curve_type else 3
		if "ls" in curve_type:
			meth = least_squares 
			func = CI_finder.least_squares_fit
		else:
			meth = minimize
			func = CI_finder.ll2 if params == 2 else CI_finder.ll3

		b = CI_finder.estimate_initial_b(self.conc_orig, self.probs, params=params)
		res = meth(func, b, args=(1-self.probs, self.conc_orig))

		fit = CI_finder.loglogit3 if len(res.x) == 3 else CI_finder.loglogit2
		self.y = fit(b = res.x, conc = self.x)
		self.plot_curve()
		r2 = CI_finder.find_r2(self.x, self.y, self.conc_orig, self.probs)
		lc50 = 2**(-res.x[0]/res.x[1])
		max_allowed = 2**(ub+2)
		min_allowed = 2**(lb-1)
		
		if lc50 < min_allowed: 
			lc50 = f"ND"
		elif lc50 > max_allowed: 
			lc50 = f">{max_allowed} ppm"
		else: 
			lc50 = f"{lc50.round(3)} ppm"
		self.curve_data_label.config(
					text=f"Approximate LC50: {lc50}  R2: {r2.round(3)}")

	def plot_curve(self):
		self.plot1.plot(self.x, self.y, ls = '-', c = 'blue')

	def set_labels(self):
		'''
		Set axes labels, ticks, title, etc. 
		'''

		def calc_x_ticks():
			'''
			Calculate the x-ticks for the graph based on the range of 
			concentrations in the data. 
			'''
			lb, ub = round(min(self.conc)), round(max(self.conc))
			xticks = np.array(range(lb-1, ub+2, 1))
			xticklabels = [round(2**i) if i >=0 else 2.**i for i in xticks]
			for i, label in enumerate(xticklabels):
				if label < 1./64.: xticklabels[i] = "{:.1e}".format(
							label).replace("e-0", "e-")
				elif label < 1./4.: xticklabels[i] = "{:.3f}".format(label)
				else: xticklabels[i] = str(label) 
			return xticks, xticklabels

		self.plot1.set_xlabel('Concentration (ppm)')
		self.plot1.set_title(label = self.current_cmpd)
		self.plot1.set_ylabel('Percent Survival')
		xticks, xticklabels = calc_x_ticks()
		self.plot1.set_xticks(xticks, minor=False)
		self.plot1.set_xticklabels(xticklabels, rotation=90, fontsize=8)
		self.plot1.set_yticks(np.array(range(0, 5, 1))/4., minor=False)
		self.plot1.set_yticklabels(np.array(range(0, 101, 25)), 
					rotation=90, fontsize=8)
		self.plot1.grid(b=True, alpha = 0.5)
		self.plot1.set_xlim([min(xticks), max(xticks)])
		self.plot1.set_ylim([0,1])

		self.plot1.set_position([.125, .16, .8,.73])
		plt.subplots_adjust(bottom=0.75)