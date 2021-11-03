# gui/statsframe/uidhandler.py
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>

class UIDHandler:
	def __init__(self, UID_in, exclude_types):
		'''
		UID_in is a list of all uids.
		exclude_types is a dictionary layed out as follows:
			{"Compound": [list of all compound names],
			 "Reference ID": [list of all reference ids],
			 "Date": [List of all dates],
			 "Plate ID": [List of all plate IDs],
			 "Row ID": [List of all row IDs]}
		'''
		self.UIDs = []
		for uid in UID_in:
			self.UIDs.append(UID(uid))

		self.all_groups = exclude_types
		self.disallowed = {}
		self.allowed = exclude_types.copy()
		for key in self.all_groups:
			if key not in self.allowed.keys(): self.allowed[key] = []
			if key not in self.disallowed.keys(): self.disallowed[key] = []

	def update(self, UID_update, exclude_update):
		disallowed_list = self.get_disallowed_uids()
		new_UIDH = UIDHandler(UID_update, exclude_update)
		for disallowed_uid in disallowed_list:
			for uid in new_UIDH.UIDs:
				uid.disallow_if("uid", disallowed_uid.uid)
		new_UIDH.update_lists()
		return new_UIDH

	def reset_DA_lists(self):
		for key in self.all_groups:
			self.allowed[key] = []
			self.disallowed[key] = []

	def disallow_if(self, var, vals):
		for val in vals:
			for uid in self.UIDs:
				uid.disallow_if(var, val)
		self.update_lists()

	def allow_if(self, var, vals):
		for val in vals:
			for uid in self.UIDs:
				uid.allow_if(var, val)
		self.update_lists()

	def update_lists(self):
		self.reset_DA_lists()
		for group, name_list in self.all_groups.items():
			for item in name_list:
				self.set_allowed_status(group, item)

	def get_allowed_uids(self):
		return list(set([uid for uid in self.UIDs if uid.allowed]))

	def get_disallowed_uids(self):
		return list(set([uid for uid in self.UIDs if not uid.allowed]))

	def set_allowed_status(self, var, val):
		allowed_list = self.get_allowed_uids()
		disallowed_list = self.get_disallowed_uids()
		for uid in allowed_list:
			if uid.dict[var] == val and uid.allowed:
				self.allowed[var].append(val)
				break
		for uid in disallowed_list:
			if uid.dict[var] == val and not uid.allowed:
				self.disallowed[var].append(val)
				break
		self.disallowed[var] = sorted(list(set(self.disallowed[var])))

class UID:
	def __init__(self, uid):
		self.uid = uid
		UID_breakdown = uid.split("_")
		self.dict = {}
		self.dict['Compound'] = UID_breakdown[0]
		self.dict['Reference ID'] = UID_breakdown[4]
		self.dict['Date'] = UID_breakdown[1]
		self.dict['plate'] = UID_breakdown[2]
		self.dict['Plate ID'] = f'{UID_breakdown[1]}_Plate_{UID_breakdown[2]}'
		self.dict['row'] = UID_breakdown[3]
		self.dict['Row ID'] = f'{UID_breakdown[1]}_Plate_{UID_breakdown[2]}_Row_{UID_breakdown[3]}'
		self.allowed = True

	def disallow(self): self.allowed = False

	def allow(self): self.allowed = True

	def disallow_if(self, var, val): self.change_allowed(var, val, action = self.disallow)

	def allow_if(self, var, val): self.change_allowed(var, val, action = self.allow)

	def change_allowed(self, var, val, action = None):
		if var.lower() == 'uid':
			if val == self.uid: action()
		if var.lower() in ['cmpd', 'compound']:
			if val == self.dict['Compound']: action()
		if var.lower() in ['refid', 'ref', 'id', 'reference id', 'ref id']:
			if val == self.dict['Reference ID']: action()
		if var.lower() in ['date', 'day']:
			if val == self.dict['Date']: action()
		if var.lower() in ['plate', 'plate id']:
			if val == self.dict['Plate ID']: action()
			else:
				val_split = val.split("_")
				if (val_split[0] == self.dict['Date']) and \
						(val_split[1] == self.dict['plate']): action()
		if var.lower() in ['row', 'row id']:
			if val == self.dict['Row ID']: action()
			else:
				val_split = val.split("_")
				if (val_split[0] == self.dict['Date']) and \
						(val_split[1] == self.dict['plate']) and\
						(val_split[2] == self.dict['row']): action()


	def __str__(self):
		return f"{str(self.dict)}, Allowed = {self.allowed}"