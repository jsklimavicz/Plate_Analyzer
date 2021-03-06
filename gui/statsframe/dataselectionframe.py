# gui/statsframe/dataselectionframe.py
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
from tkinter import ttk
from tkinter.messagebox import askyesno, askokcancel, showinfo
from tkinter.ttk import Label, Style

import pickle
import hashlib
import hmac
import os

from gui.tooltip import Tooltip #as Tooltip
from stats.main import analyze_data
from gui.statsframe.uidhandler import UIDHandler
from gui.statsframe.datapreviewer import DataPreviewer
from stats.functionfit import FunctionFit as FF 

class DataSelectionFrame(ttk.Frame):
	cache_path = os.path.abspath('./stats/cache')
	sha_key = b'james is awesome'
	forb_hash = ".forbiddenhash"
	forb_pick = "forbidden." + "pickle"
	def __init__(self, container, config, stats_obj, scale = 1, **kwargs):
		super().__init__(container, **kwargs)
		self.container = container
		self.scale = scale
		self.config = config
		self.textbox_width=60
		self.font = ('Arial', 12*self.scale)
		self.stats_obj = stats_obj

		self.uid_list = list(set(self.stats_obj.get_uid_list()))

		self.make_var_lists()
		# setup the grid layout manager
		self.columnconfigure(0, weight=2)
		self.columnconfigure(1, weight=2)
		self.columnconfigure(2, weight=2)
		self.columnconfigure(3, weight=2)
		self.columnconfigure(4, weight=2)
		self.columnconfigure(5, weight=2)

		selection_dict = {	"Compound": self.by_name, 
							"Reference ID":self.by_id,
							"Date": self.by_date,
							"Plate ID": self.by_plate,
							"Row ID": self.by_row }

		self.selection_options = list(selection_dict.keys())

		self.UIDH = self.load_UIDH()
		if self.UIDH is None: 
			self.UIDH = UIDHandler(self.uid_list, selection_dict)
		else:
			self.UIDH = self.UIDH.update(self.uid_list, selection_dict)

		self.__create_widgets()

		ff = FF(**self.stats_obj.options)
		print("C Library running: ", ff.check_Clib_running())

	def __create_widgets(self):
		'''
		Creates the interactive features necessary for excluding data from 
		analysis. 

		self.container
		self (ttk.Frame) 
		|--self.exclusionlabel (Label) 
		|--self.exclusionEntry (OptionMenu) #Select category to remove data by.
		|				Saves entry to self.exclusiontype
		|--allow_disallow_frame (Frame)
		|  |--self.allowedlabel (label)
		|  |--self.disallowlabel (Label)
		|  |--afrm (Frame) #Frame to give allowed list a scroll bar
		|  |  |--self.allowed_list w/ scrollbar #allowed list
		|  |--dfrm (Frame) #Frame to give disallowed list a scroll bar
		|  |  |--self.disallowed_list w/ scrollbar #disallowed list 
		|  |--self.disallowButton #to move items from the allowed list to 
		|  |			the forbidden list
		|  |--self.allowButton #to move items from the forbidden list to 
		|  |			the allowed list
		|  |--self.PreviewData #to preview the data
		|  |--self.clearButton #to clear the forbidden list
		|  |--self.cacheButton #to clear the cache

		'''

		s = ttk.Style()
		s.configure("my.TButton", font = self.font)

		full_width = 2
		curr_row = 0
		#Data type selection for exclusion (date, plate, row, ID, and compound)
		self.exclusionlabel = Label(self, text="Exclude runs by:")
		self.exclusionlabel.grid(row=curr_row, column=0, sticky=E)
		self.exclusiontype = StringVar()
		self.exclusiontype.set(self.selection_options[0])
		self.exclusionEntry = OptionMenu(self, 
							self.exclusiontype, 
							*self.selection_options,
							command=lambda _: self.__list_update())
		self.exclusionEntry.grid(row=curr_row, column=1, 
				columnspan = full_width,  sticky=EW, padx=10, pady=20)
		self.exclusionEntry.config(height = 2)
		msg = "Select the data type for excluding data from the analysis."
		Tooltip(self.exclusionlabel, text=msg)
		Tooltip(self.exclusionEntry, text=msg)


		allow_disallow_frame= Frame(self)
		allow_disallow_frame.grid(row=1, column=0, rowspan = 4, 
					columnspan = 6, sticky=N+S)
		#List labels
		self.allowedlabel = Label(allow_disallow_frame, text="Permitted List")
		self.allowedlabel.grid(row=0, column=0, sticky=W)
		allow_msg = "Permitted categories. Data may be included only if it"+\
				" does not match any group in the Excluded List."
		Tooltip(self.allowedlabel, text=allow_msg)
		self.disallowlabel = Label(allow_disallow_frame, text="Excluded List")
		self.disallowlabel.grid(row=0, column=4, sticky=W)
		disallow_msg = "Excluded categories. Data matching any category in "+\
				"this list will be excluded. This list only shows the "+\
				"groups in the category selected by the dropdown menu above."
		Tooltip(self.disallowlabel, text=disallow_msg)

		#Lists
		# Scroll bar help from https://stackoverflow.com/a/24656407/8075803
		#allowed List
		afrm = Frame(allow_disallow_frame)
		afrm.grid(row=1, column=0, rowspan = 4, columnspan = 2, 
					sticky=N+S)
		scrollbar = Scrollbar(afrm, orient="vertical")
		scrollbar.pack(side=RIGHT, fill=Y)
		self.allowed_list = Listbox(afrm, selectmode='multiple',
					exportselection=0, yscrollcommand=scrollbar.set, 
					width = self.textbox_width)
		self.allowed_list.pack(expand = True, fill = Y)
		scrollbar.config(command=self.allowed_list.yview)
		Tooltip(self.allowed_list, text=allow_msg)
		
		#disallowed List
		dfrm = Frame(allow_disallow_frame)
		dfrm.grid(row=1, column=4, rowspan = 4, columnspan = 2, 
					sticky=N+S)
		scrollbar = Scrollbar(dfrm, orient="vertical")
		scrollbar.pack(side=RIGHT, fill=Y)
		self.disallowed_list = Listbox(dfrm, selectmode='multiple',
					exportselection=0, yscrollcommand=scrollbar.set, 
					width = self.textbox_width)
		self.disallowed_list.pack(expand = True, fill = Y)
		scrollbar.config(command=self.disallowed_list.yview)
		Tooltip(self.disallowed_list, text=disallow_msg)

		#Disallow button
		self.disallowButton = Button(allow_disallow_frame, 
					text="Remove from sample >",
					command = lambda: self.disallow())
		self.disallowButton.grid(row=2, column=2, 
					sticky=S, padx=10, pady=20)
		self.disallowButton.config(height = 2)
		msg = "Move selected group(s) from the Permitted List to the "+\
				"Excluded List"
		Tooltip(self.disallowButton, text=msg)	

		#Allow button
		self.allowButton = Button(allow_disallow_frame, 
					text="< Add to sample", 
					command = lambda: self.allow())
		self.allowButton.grid(row=3, column=2, 
					sticky=S, padx=10, pady=20)
		self.allowButton.config(height = 2)
		msg = "Move selected group(s) from the Excluded List to the "+\
				"Permitted List"
		Tooltip(self.allowButton, text=msg)

		#Preview Button
		self.PreviewData = Button(allow_disallow_frame, text="Preview Data", 
				command = lambda: self.__preview_data())
		self.PreviewData.grid(row=5, column=0, sticky=E+W+S, 
				columnspan = 2, padx=10, pady=20)
		self.PreviewData.config(height = 2)
		msg = 'Opens a window to allow a preview of the data by compound.'
		Tooltip(self.PreviewData, text=msg)

		#Button to clear excluded list. 
		self.clearButton = Button(allow_disallow_frame, 
				text="Clear Excluded List",
				command = lambda: self.__clear_list())
		self.clearButton.grid(row = 5, column=4, 
				sticky=E+W+S, padx=10, pady=20)
		self.clearButton.config(height = 2)
		Tooltip(self.clearButton, 
				text="Clear all groups from the Excluded List.")

		#Button to clear cache. 
		self.cacheButton = Button(allow_disallow_frame, 
				text="Clear Calculated Stats", 
				command = lambda: self.__clear_stats())
		self.cacheButton.grid(row = 5, column=5, 
				sticky=E+W+S, padx=10, pady=20)
		self.cacheButton.config(height = 2)
		msg ="Clear the calculated stats. This removes all past calculated "+\
				"dose-response curve data, but does not remove any larval "+\
				"count files."
		Tooltip(self.clearButton, text=msg)

		self.__list_update()

	def disallow(self): 
		#Remove data from the allowed list.
		remove_list = []
		var = self.exclusiontype.get()
		for i in self.allowed_list.curselection():
			remove_list.append(self.UIDH.allowed[var][i])
		self.UIDH.disallow_if(var, remove_list)
		self.__list_update()

	def allow(self): 
		#Move group from the disallowed list back to the allowed list.
		add_list = []
		var = self.exclusiontype.get()
		for i in self.disallowed_list.curselection():
			add_list.append(self.UIDH.disallowed[var][i])
		self.UIDH.allow_if(var, add_list)
		self.__list_update()

	def __list_update(self):
		var = self.exclusiontype.get()
		#remove all data
		self.allowed_list.delete(0,END)
		self.disallowed_list.delete(0,END)

		#populate listBoxes
		for ind, item in enumerate(self.UIDH.allowed[var]):
			self.allowed_list.insert(ind, item)
		for ind, item in enumerate(self.UIDH.disallowed[var]):
			self.disallowed_list.insert(ind, item)

	def make_var_lists(self):
		UID_breakdown = [uid.split("_") for uid in self.uid_list]
		self.by_name = sorted(list(set([a[0]for a in UID_breakdown])))
		self.by_date = sorted(list(set([a[1]for a in UID_breakdown])))
		self.by_id = sorted(list(set([a[4]for a in UID_breakdown])))
		self.by_plate = sorted(list(set([f'{a[1]}_Plate_{a[2]}' for a in UID_breakdown])))
		self.by_row = sorted(list(set([f'{a[1]}_Plate_{a[2]}_Row_{a[3]}' for a in UID_breakdown])))

	def get_disallowed_uids(self):
		disallowed = self.UIDH.get_disallowed_uids()
		self.save_UIDH()
		return [uid.uid for uid in disallowed]

	def save_UIDH(self):
		if not os.path.exists(self.cache_path): os.makedirs(self.cache_path)
		pickle_data = pickle.dumps(self.UIDH)
		digest =  hmac.new(self.sha_key, pickle_data, hashlib.sha1).hexdigest()
		header = '%s' % (digest)
		filepath = '.'
		with open(os.path.join(self.cache_path, self.forb_hash), 'w') as file:
			file.write(header)
		with open(os.path.join(self.cache_path, self.forb_pick), 'wb') as file:
			file.write(pickle_data)

	def load_UIDH(self):
		if os.path.exists(os.path.join(self.cache_path, self.forb_hash)):
			with open(os.path.join(self.cache_path, self.forb_hash), 'r') as file:
				pickle_hash = file.read().strip()
			with open(os.path.join(self.cache_path, self.forb_pick), 'rb') as file:
				pickled_data = file.read()
			digest = hmac.new(self.sha_key, pickled_data, 
						hashlib.sha1).hexdigest()
			if pickle_hash == digest:
				unpickled_UIDH = pickle.loads(pickled_data)
				return unpickled_UIDH
		else:
			return None

	def __clear_list(self): 
		'''
		Ask to remove all groups from the forbidden list if confirmed.
		'''
		if len(self.UIDH.get_disallowed_uids()) > 0:
			msg = "Are you sure you want to clear all values from the "+\
				"excluded list? Currently, this list contains replicates from "
			cmpd_list = self.UIDH.disallowed["Compound"]
			if len(cmpd_list) > 2:
				msg+=f"{', '.join(cmpd_list[0:-1])}" + \
					f", and {cmpd_list[-1]}"
			if len(cmpd_list) == 2:
				msg+=f"{cmpd_list[0]} and {cmpd_list[-1]} "
			elif len(cmpd_list) == 1: 
				msg += f"{cmpd_list[0]}"
			
			answer = askyesno(title = "Clear Excluded List?", message = msg)
			if answer: 
				self.UIDH.allow_all()
				self.__list_update()
		else: 
			showinfo(title = "Excluded list is empty", 
				message = "There are currently no replicates being excluded.")

	def __clear_stats(self): 
		'''
		Deletes the entire cache. Obviously this is a permanent thing, but we
		need to make sure the user knows this. 
		'''
		msg = "WARNING! Deleting the calculated stats will require " +\
			"recalculating all dose-response curve data. This action " +\
			"does NOT remove any larval count files. Are you sure you " +\
			"wish to procede?" 

		answer = askokcancel(title = "WARNING!", message = msg)
		if answer: 
			for cmpd_name, cmpd in self.stats_obj.cmpd_data.items():
				self.stats_obj.cmpd_data[cmpd_name] = \
							self.stats_obj.cmpd_data[cmpd_name].reset_curves()
			msg = "You have deleted the calculated statistics. Do you also "+\
					"wish to delete the cached data?\n\nWARNING: This "+\
					"action cannot be undone. This will remove cached "+\
					"statistical data, but will not remove any larval count "+\
					"files. This program will also be closed. Do you still "+\
					"wish to proceed?"
			rm_cache=askokcancel(title="Calculated statistics cleared",
				message = msg)
			if rm_cache:
				self.__delete_cache()
				self.container.quit()

	def __preview_data(self):
		'''
		Simple preprocessor to update the disallowed UIDs for the different 
		compounds. 
		'''
		#update forbidden list
		disallowed_uids = self.get_disallowed_uids()
		self.stats_obj.set_diallowed(disallowed_uids)
		DP = DataPreviewer(self)
		DP.wait(self)
		self.__list_update()


	def __delete_cache(self): 
		cache_path = self.stats_obj.cache_path
		
		remove_files = [
				#remove the picked larval count data and hash
				os.path.join(cache_path, self.stats_obj.archivefilename),
				os.path.join(cache_path, self.stats_obj.picklesha1hash),
				#remove the pickled forbidden list and hash
				os.path.join(self.cache_path, self.forb_hash),
				os.path.join(self.cache_path, self.forb_pick)]
		for file in remove_files:
			try: os.remove(file) 
			except FileNotFoundError: pass