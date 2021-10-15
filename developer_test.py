import pickle
import hashlib
import hmac
import stats.utils as utils
import pandas as pd
from stats.curvell import CI_finder
from stats.compound import Compound
import os

##for testing purposes only. 
if __name__ == "__main__":
	filepath = "./stats/merlin_bioassay_archive_data.pickle"
	with open(filepath, 'rb') as file:
			pickled_data = file.read()
	cmpd_data = pickle.loads(pickled_data)

	by_name = sorted(list(cmpd_data.keys()))
	uid_list=[]
	for k, v in cmpd_data.items(): 
		uid_list += list(v.data['unique_plate_ids'])

	UID_breakdown = [uid.split("_") for uid in uid_list]
	by_date = sorted(list(set([a[1]for a in UID_breakdown])))
	by_id = sorted(list(set([a[4]for a in UID_breakdown])))
	by_plate = sorted([f'{a[1]}_Plate_{a[2]}' for a in UID_breakdown])
	by_row = sorted([f'{a[1]}_Plate_{a[2]}_Row_{a[3]}' for a in UID_breakdown])

	exclude = ["96P"]

	include_uid = [True for a in uid_list]
	for ind, a in enumerate(uid_list):
		for item in exclude:
			if item in a:
				include_uid[ind] = False
				break

	for a in zip(uid_list,include_uid): print(a)