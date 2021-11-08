# gui/statsframe/datapreviewer.py
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
import tkinter as tk
from tkinter.ttk import Label
from tkinter.messagebox import askyesno

import math
from scipy.optimize import minimize, least_squares
from itertools import compress
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from copy import copy, deepcopy

from stats.functionfit import FunctionFit
from stats.curvell import CI_finder
from gui.tooltip import Tooltip
from gui.utils import win_center



class DataPreviewer(tk.Frame):
	'''
	Provides an interactive interface to preview the data that will be 
	included in the statistical anlysis.
	'''
	def __init__(self, parent):
		super(DataPreviewer, self).__init__() #initialize root
		self.parent = parent
		self.scale = parent.scale
		w = 1200
		h = 1000
		self.top = Toplevel(parent)
		self.top.wm_title("Data Previewer")
		win_center(self.top, w, h)

		#get the statistics object
		self.stats_obj = deepcopy(self.parent.stats_obj)

		#make a deep copy to pass pack to parent
		self.UIDH = deepcopy(self.parent.UIDH)

		'''
		Monitors whether there was a net change to the UID list. A UID is
		included in this list if it switched from allowed to disallowed or
		vice versa, and removed from the lsit if it's switched back. If the
		length of this list is zero, then no net change as occurred. 
		'''
		self.UIDH_change_list = []

		self.columnconfigure(0, weight=2)
		self.columnconfigure(1, weight=2)
		self.columnconfigure(2, weight=2)
		self.__create_widgets()

	def wait(self, parent):
		'''
		Necessary function to ensure that the data selection window that 
		spawned this previewer waits until this window is closed, and then
		updates allowed/disallowed lists based on any changes to the UID
		handler. 
		'''
		parent.wait_window(self.top)
		return

	def __create_widgets(self, scale = 1):
		'''
		Makes the widget layout for the data previewer: 

		self.top 
		|--display_frm
		|  |--data_frm
		|  |  |--cmpd_frm #Compound selection for plotting
		|  |  |  |--cmpdlabel
		|  |  |  |--cmpd_list_frm w/ scrollbar #Inner frame 
		|  |  |  |  |--self.cmpd_listbox
		|  |  |--allowed_frm #allowed uid selection for highlighting in plot
		|  |  |  |--allowed_label
		|  |  |  |--allowed_list_frm w/ scrollbar #Inner frame 
		|  |  |  |  |--self.allowed_listbox
		|  |  |--excludeButton
		|  |  |--includeButton
		|  |  |--disallowed_frm #disallowed uid selection data highlighting
		|  |  |  |--disallowed_label
		|  |  |  |--disallowed_list_frm w/ scrollbar #Inner frame 
		|  |  |  |  |--self.disallowed_listbox
		|  |--self.fig_frame #The plot
		|  |  |--self.canvas
		|  |  |  |--self.fig
		|  |  |  |--self.curve_data_label
		|--bttn_frm
		|  |--accept_button #keep changes to UIDH
		|  |--cancel_button #cancel changes to UIDH

		'''
		display_frm = Frame(self.top)
		display_frm.grid(row=0, column=0, rowspan = 10, columnspan = 5, 
				sticky=N+S, padx=10, pady=10)
		data_frm = Frame(display_frm)
		data_frm.grid(row=0, column=0, rowspan = 9, columnspan = 2, 
						sticky=N+S, padx=10, pady=10)
		#Compound Selection
		cmpd_frm = Frame(data_frm)
		cmpd_frm.grid(row=0, column=0, rowspan = 2, columnspan = 2, 
						sticky=N+S, padx=10, pady=10)
		
		cmpdlabel = Label(cmpd_frm, text="Select Compound.") #Label
		cmpdlabel.grid(row=0, column=0, sticky=E)
		
		cmpd_list_frm = Frame(cmpd_frm) #inner frame
		cmpd_list_frm.grid(row= 1, column=0, rowspan = 2, columnspan = 2, 
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
		self.cmpd_listbox.bind("<Down>", 
			lambda event, func = self.plotselect: self.box_scroll(event, func))
		self.cmpd_listbox.bind("<Up>", 
			lambda event, func = self.plotselect: self.box_scroll(event, func))
		#configure scrollbar to listbox
		scrollbar.config(command=self.cmpd_listbox.yview) 
		#Populate the listbox with compound names.
		self.cmpd_list = sorted(list(self.stats_obj.cmpd_data.keys()))
		for ind, item in enumerate(self.cmpd_list):
			self.cmpd_listbox.insert(ind, item)
		msg = "Select a compound to preview data selected for analysis."
		Tooltip(self.cmpd_listbox, text=msg)
		

		#Allowed Selection
		allowed_frm = Frame(data_frm)
		allowed_frm.grid(row=2, column=0, rowspan = 2, columnspan = 2, 
						sticky=N+S, padx=10, pady=10)
		
		allowed_label = Label(allowed_frm, text="Included Replicates") #Label
		allowed_label.grid(row=0, column=0, sticky=E)

		allowed_list_frm = Frame(allowed_frm) #Inner frame
		allowed_list_frm.grid(row=1, column=0, rowspan = 2, columnspan = 2, 
						sticky=N+S, padx=10, pady=10)
		#Scrollbar for inner frame (will be bound to the listbox)
		scrollbar = Scrollbar(allowed_list_frm, orient="vertical")
		scrollbar.pack(side=RIGHT, fill=Y)
		msg = "These UIDs are included in the calculation of the dose-"+\
			"response curve shown. If left in this box, these UIDs will"+\
			" also be included in the final dose-response calculations. "+\
			"Data from these UIDs are shown in larger black circles in "+\
			"plot, and any currently-selected data is highlighted in red."
		Tooltip(allowed_label, text=msg)
		
		#UID listbox
		self.allowed_listbox = Listbox(allowed_list_frm, 
						selectmode='single',
						exportselection=0,
						yscrollcommand=scrollbar.set, 
						width = 30)
		self.allowed_listbox.pack(expand = True, fill = Y)
		
		#bindings
		self.allowed_listbox.bind('<<ListboxSelect>>', self.allow_select)
		self.allowed_listbox.bind("<Down>", 
			lambda event, func = self.allow_select: self.box_scroll(event, func))
		self.allowed_listbox.bind("<Up>", 
			lambda event, func = self.allow_select: self.box_scroll(event, func))
		scrollbar.config(command=self.allowed_listbox.yview) #bind scrollbar
		msg = "Identifier is of the format <date>_P<plate>_R<row>_ID. " + \
				"Select an identifier to highlight its data in the plot."
		Tooltip(self.allowed_listbox, text= msg)


		#Include button
		self.includeButton = Button(data_frm, 
					text="Include selection^",
					command = lambda: self.allow())
		self.includeButton.grid(row=4, column=0, 
					sticky=S, padx=10, pady=20)
		self.includeButton.config(height = 2)
		msg = "Move the UID from the disallowed list up to the allowed list."
		Tooltip(self.includeButton, text=msg)	

		#Exclude button
		self.excludeButton = Button(data_frm, 
					text="vExclude selection", 
					command = lambda: self.disallow())
		self.excludeButton.grid(row=4, column=1, 
					sticky=S, padx=10, pady=20)
		self.excludeButton.config(height = 2)
		msg = "Move the UID from the allowed list down to the disallowed list."
		Tooltip(self.excludeButton, text=msg)

		#Disallowed Selection
		disallowed_frm = Frame(data_frm)
		disallowed_frm.grid(row=5, column=0, rowspan = 2, columnspan = 2, 
						sticky=N+S, padx=10, pady=10)
		
		disallowed_label = Label(disallowed_frm, text="Excluded Replicates") #Label
		disallowed_label.grid(row=0, column=0, sticky=E)

		disallowed_list_frm = Frame(disallowed_frm) #Inner frame
		disallowed_list_frm.grid(row=1, column=0, rowspan = 2, columnspan = 2, 
						sticky=N+S, padx=10, pady=10)
		#Scrollbar for inner frame (will be bound to the listbox)
		scrollbar = Scrollbar(disallowed_list_frm, orient="vertical")
		scrollbar.pack(side=RIGHT, fill=Y)
		msg = "These UIDs are not included in the calculation of the dose-"+\
			"response curve shown. If left in this box, these UIDs will"+\
			" also not be included in the final dose-response calculations."+\
			"Data from these UIDs are shown in smalled gray circles in "+\
			"plot, and any currently-selected data is highlighted in red."
		Tooltip(disallowed_label, text=msg)

		#UID listbox
		self.disallowed_listbox = Listbox(disallowed_list_frm, 
						selectmode='single',
						exportselection=0,
						yscrollcommand=scrollbar.set, 
						width = 30)
		self.disallowed_listbox.pack(expand = True, fill = Y)
		
		#bindings
		self.disallowed_listbox.bind('<<ListboxSelect>>', self.disallow_select)
		self.disallowed_listbox.bind("<Down>", 
			lambda event, func = self.disallow_select: self.box_scroll(event, func))
		self.disallowed_listbox.bind("<Up>", 
			lambda event, func = self.disallow_select: self.box_scroll(event, func))
		scrollbar.config(command=self.disallowed_listbox.yview) #bind scrollbar
		msg = "Identifier is of the format <date>_P<plate>_R<row>_ID. " + \
				"Select an identifier to highlight its data in the plot."
		Tooltip(self.disallowed_listbox, text= msg)

		'''
		Figure frame
		Making the figure and canvas here is important. When made in the plot
		function, a memory leak results due to improper clearing of the plot. 
		'''
		self.fig_frame = Frame(display_frm)
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


		bttn_frm = Frame(self.top)
		bttn_frm.grid(row=10, column=0, rowspan = 2, columnspan = 5, 
				sticky=N+S, padx=10, pady=10)
		accept_button = Button(bttn_frm, 
					text="Accept changes", 
					command = lambda: self.accept_UIDH())
		accept_button.grid(row=0, column=0, rowspan = 2, columnspan = 2, 
					sticky=N+S, padx=10, pady=10)
		accept_button.config(height = 2)
		msg = "Keep the inclusions/exclusions made during this interactive "+\
					"session."
		Tooltip(accept_button, text= msg)
		cancel_button = Button(bttn_frm, 
					text="Cancel changes", 
					command = lambda: self.cancel_UIDH())
		cancel_button.grid(row=0, column=3, rowspan = 2, columnspan = 2, 
					sticky=N+S, padx=10, pady=10)
		cancel_button.config(height = 2)
		msg = "Cancel the inclusions/exclusions made during this "+\
					"interactive session."
		Tooltip(cancel_button, text= msg)

	def accept_UIDH(self):
		if len(self.UIDH_change_list) > 0: 
			msg ="Would you like to save the changes to the inclusion/exclusion"+\
					" lists that were made during this interactive session?"
			answer = askyesno(title = "Accept changes?", message = msg)
			if answer: 
				#Save modified objects back to parent
				self.parent.UIDH = self.UIDH
				self.parent.stats_obj = self.stats_obj
				self.top.destroy()
		else:
			self.top.destroy()

	def cancel_UIDH(self):
		if len(self.UIDH_change_list) > 0: 
			msg = "Would you like to cancel any changes to the inclusion/"+\
					"exclusion lists that were made during this interactive "+\
					"session?"
			answer = askyesno(title = "Cancel changes?", message = msg)
			if answer: 
				self.top.destroy()
		else:
			self.top.destroy()

	def update_change_list(self, uid):
		'''
		Updates the list of UIDs that have had a net change between the 
		include/exclude lists.  
		'''
		if uid in self.UIDH_change_list:
			self.UIDH_change_list.remove(uid)
		else:
			self.UIDH_change_list.append(uid)

	def select_first(self, list_type = "allowed"):
		'''
		Resets the selection to the first apprpriate ID after moving allowed
		or disallowed UIDs between the lists. 
		'''
		self.plotselect()

		if list_type == "allowed":
			curr_list = self.allowed_listbox
			uid_list = self.included_uid
		else:
			curr_list = self.disallowed_listbox
			uid_list = self.excluded_uid

		curr_list.selection_clear(0, "end")
		current_uid = uid_list[0]
		self.update_plot(current_uid)
		curr_list.selection_set(0)

	def allow(self):
		'''
		Moves the UID from the disallowed list to the allowed list, and updates the 
		UID handler and the actual data. 
		'''
		try:
			uid = self.excluded_uid[int(self.disallowed_listbox.curselection()[0])].uid
			self.UIDH.allow_if('uid', [uid])
			self.update_disallowed()
			try:
				self.select_first(list_type = "disallowed")
			except IndexError:
				self.select_first(list_type = "allowed")
			self.update_change_list(uid)
		except IndexError:
			pass

	def disallow(self):
		'''
		Moves the UID from the allowed list to the disallowed list, and updates the 
		UID handler and the actual data. 
		'''
		try:
			uid = self.included_uid[int(self.allowed_listbox.curselection()[0])].uid
			self.UIDH.disallow_if('uid', [uid])
			self.update_disallowed()
			try:
				self.select_first(list_type = "allowed")
			except IndexError:
				self.select_first(list_type = "disallowed")
			self.update_change_list(uid)
		except IndexError:
			pass

	def box_scroll(self, event, func):
		'''
		Adapted from https://stackoverflow.com/a/62711725/8075803
		Controls actions on up/down scrolling. 
		'''
		selection = event.widget.curselection()[0]
		if event.keysym == 'Up': selection += -1
		if event.keysym == 'Down': selection += 1
		if 0 <= selection < event.widget.size():
			event.widget.selection_clear(0, tk.END)
			event.widget.select_set(selection)
			func(event)

	def plotselect(self, event = None):
		'''
		Generates plots of data based on the data of the compound selected. 
		Adapted from https://stackoverflow.com/a/12936031/8075803
		'''
		if event is not None:
			ind = int(event.widget.curselection()[0])
		else:
			try:
				ind = int(self.cmpd_listbox.curselection()[0])
			except IndexError:
				ind = 0
		self.current_cmpd = self.cmpd_list[ind]
		#Get the allowed data list

		self.included_uid = self.UIDH.get_allowed_uids()
		self.included_uid = [uid for uid in self.included_uid if \
							uid.dict["Compound"]==self.current_cmpd]
		self.included_uid.sort(key=lambda x: str(x))
		self.excluded_uid = self.UIDH.get_disallowed_uids()
		self.excluded_uid = [uid for uid in self.excluded_uid if \
							uid.dict["Compound"]==self.current_cmpd]
		self.excluded_uid.sort(key=lambda x: str(x))
		#clear previous list
		self.allowed_listbox.delete(0,END)
		for ind, uid in enumerate(self.included_uid):
			self.allowed_listbox.insert(ind, str(uid))

		self.disallowed_listbox.delete(0,END)
		for ind, uid in enumerate(self.excluded_uid):
			self.disallowed_listbox.insert(ind, str(uid))
		#update the plot
		self.plot()
		self.update_plot()

	def update_disallowed(self):
		'''
		gets the list of all the uids from the diallowed list, and then sets the 
		include bool in the actual stats object for stats.
		'''
		disallowed = [uid.uid for uid in self.UIDH.get_disallowed_uids()]
		self.stats_obj.set_diallowed(disallowed)

	def change_select(self, event, listbox, uid_list):
		'''
		Highlights a UID in the desired listbox and removes any selection
		from the other listbox.
		'''
		listbox.selection_clear(0, 'end')
		try:
			if event is not None:
				w = event.widget
				current_uid = uid_list[int(w.curselection()[0])]
			else:
				current_uid = uid_list[0]
			self.update_plot(current_uid)
		except IndexError:
			pass

	def allow_select(self, event = None):
		'''
		Highlights the data from the currently-selected data in the allow list
		in the plot, and removes any highlight from the disallow box
		'''
		self.change_select(event, self.disallowed_listbox, self.included_uid)

	def disallow_select(self, event = None):
		'''
		Highlights the data from the currently-selected data in the disallow
		list in the plot, and removes any highlight from the allow box
		'''
		self.change_select(event, self.allowed_listbox, self.excluded_uid)

	def update_plot(self, current_uid = None):
		'''
		Uses the current conc/probs from the selected compound to highlight 
		the data from the currently-selected UID. Note that this method 
		actually creates a completely new plot; clearing the old plot is 
		performed to prevent any memory leaks. 
		'''
		self.plot1.clear()

		#plot active and inactive
		self.plot_data(self.conc_active, self.probs_active, style = "active")
		self.plot_data(self.conc_inactive, self.probs_inactive, 
														style = "inactive")

		if current_uid is not None:
			uid_conc = []
			uid_probs = []
			for ind, uid in enumerate(self.stats_obj.cmpd_data[
							self.current_cmpd].data['unique_plate_ids']):
				if uid == current_uid.uid:
					uid_conc.append(self.conc[ind])
					uid_probs.append(self.probs[ind])
			#plot highlighted data
			self.plot_data(uid_conc, uid_probs, style = "highlight")


		self.set_labels() #update axes, etc.
		if len(self.conc_active) > 2: 
			self.plot_curve()
		else:
			self.null_curve()
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
		self.probs = np.array(self.curr_data.data["dead_count"]/ \
					(self.curr_data.data["live_count"] + \
					self.curr_data.data["dead_count"]))

		if self.stats_obj.options["JITTER"]: 
			self.conc += np.random.uniform(
					-self.stats_obj.options["JITTER_FACTOR"],
					self.stats_obj.options["JITTER_FACTOR"], 
					len(self.conc))
		include_now = self.curr_data.data['include_now']
		self.calc_conc = np.array(list(compress(self.conc_orig, include_now)))
		self.conc_active = list(compress(self.conc, include_now))
		self.probs_active = np.array(list(compress(self.probs, include_now)))
		self.conc_inactive = list(compress(self.conc, 
								[not elem for elem in include_now]))
		self.probs_inactive = list(compress(self.probs, 
								[not elem for elem in include_now]))
		
		self.plot_data(self.conc_active, self.probs_active, style = "active")
		self.plot_data(self.conc_inactive, self.probs_inactive, 
														style = "inactive")

		self.set_labels()
		if len(self.conc_active) > 2:
			self.make_curve()
		else:
			self.null_curve()
		# else:

		self.canvas.draw()
		
	def null_curve(self):
		'''
		Provide a suitable comment for when no data is selected for plotting.
		'''
		self.curve_data_label.config(
					text="Insufficent data selected to fit a curve.\n"+\
					"Compound will not show up in data analysis.")

	def make_curve(self):
		'''
		Generates the x and y values to produce a best fit curve. The fit 
		curve is equivalent to the kind specified in the analysis_config.txt
		The curve is then used to generate an approximate LC50 and an R2. 
		'''
		lb, ub = round(min(self.conc_orig)), round(max(self.conc_orig))
		num_points = 10 * (ub - lb + 2) + 1 #10 x-points (including both ends)
		self.x = np.linspace(lb-1, ub+1, num_points)
		
		#Find the best fit based on the number of parameters and curve-fitting
		#method specified in the options. 
		ff = FunctionFit(**self.stats_obj.options)
		switch = self.stats_obj.options["CURVE_TYPE"].lower()
		b = ff.switch_fitter(switch, self.calc_conc, 
							self.probs_active, rev = False)
		self.y = ff.loglogit3(b = b, conc = self.x)

		self.plot_curve()
		r2 = CI_finder.find_r2(self.x, self.y, 
							self.calc_conc, 1-self.probs_active)

		if b[2] <= 1.0:
			lc50 = 2**(-b[0]/b[1])
		else:
			lc50 = 2**((math.log(2*b[2] - 1) - b[0])/b[1])

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
		#Does the actual curve plotting
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

	def plot_data(self, conc, probs, style = "background"):
		'''
		Plots the data points of a graph. This is included for consistency of 
		plots between plotting and replotting with highlighted data. 
		'''
		probs  = 1.0 - np.array(probs)
		if style == "active":
			self.plot1.plot(conc, 
					probs, 
					marker = '.', 
					mew = 0.0, 
					mfc = 'black', 
					ls = 'None')

		if style in ['background', 'inactive']:
			self.plot1.plot(conc, 
					probs,  
					marker = '.', 
					mew = 0.0, 
					ms = 4,
					mfc = 'gray', 
					ls = 'None')

		if style in ['highlight']:
			self.plot1.plot(conc, 
					probs, 
					marker = 'o', 
					mew = 1.0, 
					mec = 'red', 
					mfc = 'none',
					ms = 5.0,
					ls = 'None')
