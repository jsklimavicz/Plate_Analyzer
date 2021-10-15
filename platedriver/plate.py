#!/usr/bin/env python3
#plate.py
import math
import numpy as np
import glob
from os.path import basename, dirname, join, exists
from os import path, makedirs, listdir, mkdir
import platform
import pandas as pd
from skimage.io import imsave, imread
from platedriver.platedata import PlateData
from joblib import load, Parallel, delayed
from skimage.color import label2rgb, rgb2gray
from skimage.exposure import rescale_intensity
import skimage
import platedriver.well_prep as wpp
import wormAI.worm_detector as worm_detector
from wormAI.mrcnn import model as modellib, utils
from pathlib import Path
from datetime import datetime, timedelta
import gui.utils as pdu

#for parallel call to provide self pointer as arg
def unwrap_image_process(arg, **kwarg):
    return Plate.process_images_driver(*arg, **kwarg)

class Plate:
	suffix = 'Reflected_001.tif'

	def __init__(self, 
				config,
				save_dir,
				plate_data=PlateData(), 
				crop_dim = 800, 
				debug = False):
		self.config = config
		self.save_dir = save_dir
		self.debug = debug
		self.image_folder = config['IMG_DIR']
		self.message = "Starting..."
		self.progress = 0 # used for gui progressbar
		self.processed_well_count = 0 # used for gui progressbar
		self.error_message = None # used to pass errors to the gui
		self.well_ID_list = self.__get_well_IDs()
		self.crop_dim = crop_dim
		if len(self.well_ID_list) == 0: return # no images found; return to gui
		self.total_wells = len(self.well_ID_list)
		self.plate_data = plate_data #set the actual plate data. 
		self.plate_layout = [] #array for wells in plate bounding box
		self.__make_plate_layout() #makes plate layout and bounding box
		# print(self.plate_layout)
		self.well_list = {} #dictionary with wellID as keys, and well info as values
		self.verbose_plate_data = pd.DataFrame() 

	def run_images(self, init_only=False):
		self.init_only = init_only
		# self.__initiate_images() #make empty (black) images of each plate
		self.__process_images() #process all the images. This is where the fun happens
		
	def __check_for_subdir(self):
		'''This function checks to see if the provided directory contains image folders, or if it instead contains
		a subdirectory that contains the image folders. This is here because the Gen5 software creates a weird 
		directory system wherein images are two folders deep (likely to handle multiple plates in one save file.
		Adding to the issue, tkinter select directory methods don't show non-directory items in a folder, so it's
		difficult for the user to know whether the chosen directory actually contains images. 
		If the provided directory was in fact the superdirectory that doesn't contian images, but it does contain
		a subdirectory that does have files of the appropriate name format, then we will change the image directory path.'''
		sub_dirs = listdir(path=self.image_folder)
		if len(sub_dirs) == 0:
			self.error_message = "No images found in the folder selected. Please check the image folder path."
			return None
		head, tail = path.split(self.image_folder)
		if tail == "":
			head, tail = path.split(head)
			print(head, tail)
		for new_dir in sub_dirs:
			if new_dir[0:13] == tail[0:13]:
				#appropriate subdirectory potentially found. 
				well_ID_list = [] 
				for filename in glob.glob(path.join(self.image_folder, new_dir, ('*' + self.suffix))):
					image = basename(filename) # get just the filename portion
					well_ID_list.append(image.split("_")[0]) # seperate out the well ID
				if len(well_ID_list) >=3 :
					self.image_folder = path.join(self.image_folder, new_dir)
					return well_ID_list
		return None


	def __get_well_IDs(self):
		well_ID_list = [] #Loads wellIDs based on image names
		#Search for all eligible image files in the folder.
		for filename in glob.glob(path.join(self.image_folder, ('*' + self.suffix))):
			image = basename(filename) # get just the filename portion
			well_ID_list.append(image.split("_")[0]) # seperate out the well ID
		if len(well_ID_list) == 0:
			well_ID_list = self.__check_for_subdir()
			if not well_ID_list or len(well_ID_list) == 0:
				#Make we found something; else, log an error. 
				self.error_message = "No images found in the folder selected. Please check the image folder path."
		#get unique IDs
		id_list = list(set(well_ID_list))
		id_list.sort()
		return id_list

	#ADJUSTED	
	def __well_bounds(self):
		#Returns upper and lower bounds for the well array in the plate. 
		#Gaps between wells are allowed, but missing rows/columns at the edges of plates are not included. 
		#The goal of this is to keep image files as small as possible.
		min_char = 128
		max_char = 0
		min_num = 13
		max_num = 0
		#find bounding box of plate based on well IDs. Letters stored as ASCII for predictable well-ordering, 
		#and for easier row bounds setting in making the plate layout.
		for well in self.well_ID_list:
			row = ord(well[0:1]) #convert letter to number
			col = int(well[1:]) #convert str num to int
			min_char = min(row, min_char)
			max_char = max(row, max_char)
			min_num = min(col, min_num)
			max_num = max(col, max_num)
		return min_char, max_char, min_num, max_num

	#ADJUSTED
	def __make_plate_layout(self):
		# make a list of lists with well ID in each appropriate position.
		min_char, max_char, min_num, max_num = self.__well_bounds()
		ind = 0
		#Make 2D array of well IDs where they should be on the plate. 
		for ascii in range(min_char,max_char+1):
			self.plate_layout.append([])
			row = chr(ascii)
			for col in range(min_num,max_num+1): self.plate_layout[ind].append(row + str(col))
			ind += 1 
		self.plate_layout = np.array(self.plate_layout)

		#FOR NOW, ONLY MAKES IMAGES
	def __process_images(self):
		images = {}
		#Start the image processing step after updating parameters
		self.message = "Starting multiprocessing pool..."
		self.progress = 2
		#Use joblib's parallel to perform num_cpu jobs. Note that because we're in a class, and joblib's parallel is wonky, 
		# we must call and outside function (unwrap_image_process()) to get Parallel to work as expected. Here, we wrap the 
		# "self" point for the arg and the well IDs for the **kwargs Sharing memory is on to allow persistance of class 
		# variables, e.g. number of wells counted and the progress made for the progres bar.
		dict_list = Parallel(n_jobs=-1, require='sharedmem')(delayed(unwrap_image_process)(well) for well in zip([self]*len(self.well_ID_list), self.well_ID_list))

		self.message = "Calculating processed images..."
		for line in dict_list:
			for wellID in line:
				self.well_list[wellID] = line[wellID]

	def process_images_driver(self, wellID):
		try:
			# Open the before and after images
			image_paths = glob.glob(join(self.image_folder, (wellID + '_*' + self.suffix)))
			if not image_paths: return None
			# print(image_paths)
			image_paths.sort()
			imageA = self.__well_image_opener(image_paths[0])
			imageB = self.__well_image_opener(image_paths[1])
			imageC = self.__well_image_opener(image_paths[2])

			# increment the number of wells counted. 
			self.processed_well_count +=1 
			# change progress bar message to show what row we're on in the plate.
			self.message = f"Performing preprocessing on well {(wellID+'.'):3} This will take a few minutes; go grab some tea!"
			# increment progress bar value before and after the time-consuming read well process.
			# A larger chunk of completion time is allocated *before* the process to allow the 
			# progress bar to take samll, slow steps for most of the process to allow the user to know
			# that the process is still running.
			self.progress += 30 * (0.3 / self.total_wells)
			print(f"Processing {wellID}.", end = "\r")

			comp = self.__read_well(wellID, imageA,imageB,imageC)
			# self.batch_image_save(wellID, comp)

			self.progress += 30 * (0.7 / self.total_wells)
			# make a tuple of the values that will be returned for use in making the composite class images
			image_dict = {}
			image_dict[wellID] = {'wellID': wellID, 'composite': comp}
			#return the stuff!
			return image_dict
		# this is for when a well is included in the well array but images were not taken, e.g. column 11 of the plate. 
		# Also useful if we want to do some analyzing of the plate but only on some wells.
		except IOError: pass

	def __well_image_opener(self, image_path): return imread(image_path)

	def __read_well(self, wellID, imageA, imageB, imageC):
		#perform preprocessing
		IPP = wpp.ImagePreprocessor(wellID, blur_sigma = 1, gamma = 1.4, crop_dim = self.crop_dim, debug=self.debug)
		IPP.well_preprocess(imageA,imageB,imageC) if not self.init_only else IPP.quick_gray_rescale(imageA)
		return IPP.get_composite()

	def normalize8(self, I): return (((I - I.min())/(I.max()-I.min())) * 255).astype(np.uint8)

	def detect_larvae(self):
		
		self.message = f"Starting larval neural network..."
		config = worm_detector.WormConfig(self.config)
		model = modellib.MaskRCNN(mode="inference", config=config, model_dir="./wormAI/mrcnn/log/")
		config.display()
		WEIGHTS_PATH = worm_detector.get_weights_path()
		model.load_weights(WEIGHTS_PATH, by_name=True)
		well_count = 0

		if not exists(self.save_dir): mkdir(self.save_dir)
		if self.config["SAVE_MASK"]: 
			self.mask_dir = join(self.save_dir, "mask")
			if not exists(self.mask_dir): mkdir(self.mask_dir)
		if self.config["SAVE_BBOX"]: 
			self.bbox_dir = join(self.save_dir, "bbox")
			if not exists(self.bbox_dir): mkdir(self.bbox_dir)
		if self.config["SAVE_INTMD"]: 
			self.comp_dir = join(self.save_dir, "comp")
			if not exists(self.comp_dir): mkdir(self.comp_dir)

		start_time = datetime.utcnow()
		for wellID in self.well_ID_list:

			im = self.well_list[wellID]['composite']
			r = model.detect([im], verbose=0)[0]

			self.message = f"Counting larvae in well {wellID}."
			self.progress += 65 * (0.7 / self.total_wells)
			if well_count < 5:
				self.message += " Estimating time remaining..."
			else:
				now = datetime.utcnow()
				time_per_well = (now - start_time).total_seconds()/(well_count+1)
				timeleft = round((self.total_wells - well_count) * time_per_well)
				if timeleft > 10:
					hours, remainder = divmod(timeleft, 3600)
					mins, secs = divmod(remainder, 60)
					secs = f"{secs}s " if (mins >0 or hours>0) else f"{secs} seconds"
					hours = f"{hours}h " if hours>0 else ""
					mins = f"{mins}m " if mins >0 else ""

					time_format = f"{hours}{mins}{secs}"

					self.message += f" Estimated time remaining: {time_format}"
				else:
					self.message += " Finishing up..."

			well_count += 1

			self.well_list[wellID]['rois'] = r["rois"]
			self.well_list[wellID]['class_ids'] = (r["class_ids"])
			self.well_list[wellID]['scores'] = (r["scores"])
			if self.config["SAVE_BBOX"] or self.config["SAVE_MASK"] or self.config["SAVE_INTMD"]:
				gray = (rgb2gray(im[:,:,0:3])*255).astype(np.uint8)
				splash_img, bbox_img = worm_detector.color_splash_and_bbox(gray,r)
				if self.config["SAVE_MASK"]: 
					filename = f'{wellID}_splash.png'
					imsave(join(self.mask_dir, filename), self.normalize8(splash_img))
				if self.config["SAVE_BBOX"]: 
					filename = f'{wellID}_bbox.png'
					imsave(join(self.bbox_dir, filename), self.normalize8(bbox_img)) 
				if self.config["SAVE_INTMD"]: 
					rgb_im = (np.dstack((gray, im[:,:,4], im[:,:,5]))*255).astype(np.uint8)
					filename = f'{wellID}_comp.png'
					imsave(join(self.comp_dir, filename), self.normalize8(rgb_im)) 
			self.progress += 65 * (0.3 / self.total_wells)

	def save_csv(self, filename = None, verbose_csv = False, open_csv=False):
		if not exists(self.save_dir): mkdir(self.save_dir)
		self.message = "Saving the csv file(s)..."
		# print(verbose_csv)
		data, verb_data, ctrl_warning = self.__make_summary(verbose_csv = verbose_csv)
		if self.config["VERBOSE_CSV"]: 
			verb_filename = f"{self.config['CSV_NAME'][:-4]}_objectdata.csv"
			verb_fp = join(self.save_dir, verb_filename)
			with open(verb_fp,'w') as file: file.write("\n".join(verb_data))
		with open(join(self.save_dir, self.config['CSV_NAME']),'w') as file: file.write("\n".join(data))
		if self.config["OPEN_CSV"] and "Windows" in platform.system(): 
			from os import startfile
			startfile(join(self.save_dir, filename))
		if ctrl_warning: pdu.message_window(title=ctrl_warning['title'], msg=ctrl_warning['msg'])

	def __make_summary(self, verbose_csv = False):
		# print(verbose_csv)
		inv_worm_classes = {"alive": 1,
						"dead": 2,
						"moribund": 3,
						"egg": 4, 
						"L2": 5, 
						"aged_dead": 6, 
						"artifact": 7 }
		worm_classes = {v: k for k, v in inv_worm_classes.items()}
		self.message = "Summarizing data results ..."
		if verbose_csv: 
			verb_columns = ["Ref ID", "Compound", "Date", "ppm", "Rep", "Plate", 'Row', 'Column', "Media", "Object ID", "Bounding Box", "Class", "Class Name", "Score"]
			verb_header = ",".join(verb_columns)
			verb_data = [verb_header]
		else:
			verb_data = None
		columns = ["Ref ID", "Date", "ppm", "Rep", "Plate", 'Row', 'Column', "Total", "Dead", "Alive", "Media", "Notes", "% Mortality"]
		header = ",".join(columns)
		data = [header]

		ctrl_dead = 0
		ctrl_total = 0
		for row in range(0,len(self.plate_layout)):
			for col in range(0,len(self.plate_layout[row])):
				wellID = self.plate_layout[row][col]	
				try:
					curr_row = {**self.plate_data.wellID_data[wellID]}
					curr_row['Row'] = row + 1
					curr_row['Column'] = col + 1
					curr_row['Date'] = self.plate_data.date.strftime(self.config["PRETTY_DATE_FORMAT"])
					curr_row['Plate'] = self.plate_data.plate_number
					curr_row['Media'] = self.plate_data.control_name
					classes = self.well_list[wellID]['class_ids'].tolist()
					curr_row['Alive'] = classes.count(inv_worm_classes["alive"])
					curr_row['Dead'] = classes.count(inv_worm_classes["dead"]) + classes.count(inv_worm_classes["moribund"])
					curr_row['Total'] = curr_row['Alive'] + curr_row['Dead']
					if col == 11: #control
						ctrl_dead += curr_row['Dead']
						ctrl_total += curr_row['Total']

					mort = round(100. * float(curr_row['Dead']) / (1e-7 + float(curr_row['Total'])))
					curr_row['% Mortality'] = f"{mort}%"

					#Autopopulate the notes column with anything potentially of note
					curr_row['Notes'] = ""
					if (curr_row['Column'] == 12 and mort >= 25): curr_row['Notes'] += "High control mortality. "
					if (curr_row['Total'] > 25): curr_row['Notes'] += "High worm count. "
					if (curr_row['Total'] == 0): 
						curr_row['Notes'] += "NO LARVAE DETECTED. "
					elif (curr_row['Total'] < 5): 
						curr_row['Notes'] += "Low worm count. "
					egg_count = classes.count(inv_worm_classes["egg"])
					L2_count = classes.count(inv_worm_classes["L2"])
					aged_dead_count = classes.count(inv_worm_classes["aged_dead"])
					artifact_count = classes.count(inv_worm_classes["artifact"])
					uncounted_count = egg_count + L2_count + aged_dead_count
					if (uncounted_count >= 5) or (artifact_count >= 5) or (uncounted_count + artifact_count > 8): 
						art_str = []
						if egg_count > 0: art_str.append(f"{egg_count} egg{'s' if egg_count > 1 else ''}")
						if L2_count > 0: art_str.append(f"{L2_count} L2 larva{'e' if L2_count > 1 else ''}")
						if aged_dead_count > 0: art_str.append(f"{aged_dead_count} aged dead larva{'e' if aged_dead_count > 1 else ''}")
						if artifact_count > 0: art_str.append(f"{artifact_count} artifact{'s' if artifact_count > 1 else ''}")
						curr_row['Notes'] += "; ".join(art_str)


					curr_line = "" 
					for item in columns:
						curr_line += f"{curr_row[item]},"
					data.append(curr_line[0:-1])
					if verbose_csv:
						objID=0
						scores = self.well_list[wellID]['scores']
						bboxes = self.well_list[wellID]['rois'] 
						for obj in classes: 
							verb_row = curr_row.copy()
							verb_row["Bounding Box"] = f'"{bboxes[objID].tolist()}"'
							verb_row["Score"] = f"{float(scores[objID]):.4f}"
							verb_row["Class"] = f"{obj}"
							verb_row["Class Name"] = f"{worm_classes[obj]}"
							objID+=1
							verb_row["Object ID"] = f"{objID}"
							verb_line = ""
							for item in verb_columns:
								verb_line += f"{verb_row[item]},"
							verb_data.append(verb_line[0:-1])
					# print(verb_data)
				except ValueError:
					# exit()
					self.error_message = "Cannot properly form output file. Alert James to the issue."
				except KeyError:
					# print(f"Skipping well {wellID}.")
					# exit()
					continue
		ctrl_warning = None
		ctrl_mort = float(ctrl_dead)/float(1.e-7 + ctrl_total)
		if  ctrl_mort> 0.15:
			ctrl_warning = {'title': "High Control Mortality", 'msg': f"Average control mortality was {round(ctrl_mort*100)}%. Review data to assess its utility."}
		return data, verb_data, ctrl_warning
