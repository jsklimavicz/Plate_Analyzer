#!/usr/bin/env python3
#plate.py
import math
import numpy as np
import glob
from os.path import basename, dirname, join
from os import path, makedirs, listdir
import platform
import pandas as pd
from skimage.io import imsave, imread
from platedriver.platedata import PlateData
from joblib import load, Parallel, delayed
from skimage.color import label2rgb
from skimage.exposure import rescale_intensity
import skimage
import platedriver.well_prep as wpp
import wormAI.worm_detector as worm_detector
from wormAI.mrcnn import model as modellib, utils

#for parallel call to provide self pointer as arg
def unwrap_image_process(arg, **kwarg):
    return Plate.process_images_driver(*arg, **kwarg)

#for parallel call to provide self pointer as arg
def unwrap_image_save(arg, **kwarg):
    return Plate.batch_image_save(*arg, **kwarg)

class Plate:
	suffix = 'Reflected_001.tif'

	def __init__(self, image_folder, save_dir, plate_data=PlateData(), crop_dim = 800, debug = False):
		self.debug = debug
		self.image_folder = image_folder
		self.save_dir = save_dir
		self.plate_message = "Starting..."
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

	def run_images(self, init_only=False, intmd_save=True):
		self.init_only = init_only
		# self.__prepare_classifier() #load classifiers
		# self.__initiate_images() #make empty (black) images of each plate
		self.intmd_save = intmd_save
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
		self.rgb_images = {}
		#Start the image processing step after updating parameters
		self.plate_message = "Starting multiprocessing pool..."
		self.progress = 2
		print(self.plate_message)
		#Use joblib's parallel to perform num_cpu jobs. Note that because we're in a class, and joblib's parallel is wonky, we must call and outside 
		#function (unwrap_image_process()) to get Parallel to work as expected. Here, we wrap the "self" point for the arg and the well IDs for the **kwargs
		#Sharing memory is on to allow persistance of class variables, e.g. number of wells counted and the progress made for the progres bar.
		self.plate_message = "Processing raw images..."
		dict_list = Parallel(n_jobs=-1, require='sharedmem')(delayed(unwrap_image_process)(well) for well in zip([self]*len(self.well_ID_list), self.well_ID_list))
		for well, img, uid in dict_list: self.rgb_images[well]={'well': well, 'img': img, 'uid': uid}

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
			self.plate_message = f"Performing preprocessing on well images (currently on row {wellID[0:1]}; go grab some tea or coffee)..."
			# increment progress bar value before and after the time-consuming read well process.
			# A larger chunk of completion time is allocated *before* the process to allow the 
			# progress bar to take samll, slow steps for most of the process to allow the user to know
			# that the process is still running.
			self.progress += 30 * (0.7 / self.total_wells)
			print(f"Processing {wellID}.", end = "\r")

			image = self.__read_well(wellID, imageA,imageB,imageC)

			path = basename(self.image_folder)
			path = path.split("_")
			path = f'{path[0]}_{path[1]}'
			uid = f'{path}_{wellID}'

			# filename = f'{uid}.png'
			# print()
			# image = imread(join(self.save_dir, filename))

			if self.intmd_save: 
				filename = f'{uid}.png'
				imsave(join(self.save_dir, filename), self.normalize8(image)) 
			self.progress += 30 * (0.3 / self.total_wells)
			# make a tuple of the values that will be returned for use in making the composite class images
			image_tuple = (wellID, image, uid)
			return image_tuple

		except IOError: 
			# this is for when a well is included in the well array but images were not taken, e.g.
			# column 11 of the plate. Also useful if we want to do some analyzing of the plate but
			# only on some wells.
			pass

	def __well_image_opener(self, image_path):
		image = imread(image_path)
		return image

	def __read_well(self, wellID, imageA, imageB, imageC):
		#perform preprocessing
		IPP = wpp.ImagePreprocessor(wellID, blur_sigma = 1, gamma = 1.4, crop_dim = self.crop_dim, debug=self.debug)
		IPP.well_preprocess(imageA,imageB,imageC) if not self.init_only else IPP.quick_gray_rescale(imageA)
		return IPP.get_image()

	def normalize8(self, I):
		mn = I.min()
		mx = I.max()
		I = ((I - mn)/(mx-mn)) * 255
		return I.astype(np.uint8)

	def detect_larvae(self, splash = False, bbox = False):
		self.plate_message = f"Starting larval neural network..."
		config = worm_detector.WormConfig()
		model = modellib.MaskRCNN(mode="inference", config=config, model_dir="./wormAI/mrcnn/log/")
		config.display()
		WEIGHTS_PATH = worm_detector.get_weights_path()
		model.load_weights(WEIGHTS_PATH, by_name=True)
		well_count = 0
		for wellID in self.well_ID_list:
			self.plate_message = f"Counting larvae in well {wellID}."
			self.progress += 65 * (0.7 / self.total_wells)
			if well_count < self.total_wells/2: 
				self.plate_message += " This will take a few minutes. Go grab some coffee or tea."
			elif well_count < 0.85*self.total_wells:
				self.plate_message += " You still have a few minutes to grab some tea!"
			else:
				self.plate_message += " Almost done..."
			well_count += 1
			im = skimage.img_as_ubyte(self.rgb_images[wellID]['img'])
			r = model.detect([im], verbose=0)[0]
			self.rgb_images[wellID]['rois'] = r["rois"]
			self.rgb_images[wellID]['class_ids'] = (r["class_ids"])
			self.rgb_images[wellID]['scores'] = (r["scores"])
			if bbox or splash:
				path = basename(self.image_folder)
				path = path.split("_")
				path = f'{path[0]}_{path[1]}'
				splash_img, bbox_img = worm_detector.color_splash_and_bbox(im,r)
				if splash: 
					filename = f'{path}_{wellID}_splash.png'
					imsave(join(self.save_dir, filename), self.normalize8(splash_img))
				if bbox: 
					filename = f'{path}_{wellID}_bbox.png'
					imsave(join(self.save_dir, filename), self.normalize8(bbox_img)) 
			self.progress += 65 * (0.3 / self.total_wells)


	def save_csv(self, filename = None, verbose_csv = False, open_csv=False):
		self.plate_message = "Saving the csv file(s)..."
		# print(verbose_csv)
		data, verb_data = self.__make_summary(verbose_csv = verbose_csv)
		if verbose_csv: 
			verb_filename = f"{filename[:-4]}_objectdata.csv"
			verb_fp = join(self.save_dir, verb_filename)
			with open(verb_fp,'w') as file: file.write("\n".join(verb_data))
		with open(join(self.save_dir, filename),'w') as file: file.write("\n".join(data))
		if open_csv and "Windows" in platform.system(): 
			from os import startfile
			startfile(join(self.save_dir, filename))

	def __make_summary(self, verbose_csv = False):
		# print(verbose_csv)
		worm_classes = {1: "alive", \
						2: "dead", \
						3: "moribund", \
						4: "egg"}
		inv_worm_classes = {v: k for k, v in worm_classes.items()}
		self.plate_message = "Summarizing data results ..."
		if verbose_csv: 
			verb_columns = ["Ref ID", "Date", "ppm", "Rep", "Plate", 'Row', 'Column', "Media", "Object ID", "Bounding Box", "Class", "Class Name", "Score"]
			verb_header = ",".join(verb_columns)
			verb_data = [verb_header]
		else:
			verb_data = None
		columns = ["Ref ID", "Date", "ppm", "Rep", "Plate", 'Row', 'Column', "Total", "Dead", "Alive", "Media"]
		header = ",".join(columns)
		data = [header]
		for row in range(0,len(self.plate_layout)):
			for col in range(0,len(self.plate_layout[row])):
				wellID = self.plate_layout[row][col]
				try:
					curr_row = {**self.plate_data.wellID_data[wellID]}
					curr_row['Row'] = row + 1
					curr_row['Column'] = col + 1
					curr_row['Date'] = self.plate_data.date
					curr_row['Plate'] = self.plate_data.plate_number
					curr_row['Media'] = self.plate_data.control_ID 
					classes = self.rgb_images[wellID]['class_ids'].tolist()
					curr_row['Alive'] = classes.count(inv_worm_classes["alive"])
					curr_row['Dead'] = classes.count(inv_worm_classes["dead"]) + classes.count(inv_worm_classes["moribund"])
					curr_row['Total'] = curr_row['Alive'] + curr_row['Dead']
					curr_line = ""
					for item in columns:
						curr_line += f"{curr_row[item]},"
					data.append(curr_line[0:-1])
					if verbose_csv:
						objID=0
						scores = self.rgb_images[wellID]['scores']
						bboxes = self.rgb_images[wellID]['rois'] 
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
		return data, verb_data

	# #ADJUSTED
	# def __initiate_composite_image(self):
	# 	#makes empty images the size of the plate
	# 	self.plate_progress = 1 #Move progress bar a wee bit
	# 	self.plate_message = "Initiating images..." #Update messaging for gui
	# 	#calculate dimensions of the plate based on wells and image size
	# 	dim1 = self.crop_dim*len(self.plate_layout) 
	# 	dim2 = self.crop_dim*len(self.plate_layout[0])
	# 	#Make image arrays of zeros of the appropriate size
	# 	self.full_plate = np.zeros((dim1, dim2, 3), dtype=np.float32) if not self.init_only else np.zeros((dim1, dim2), dtype=np.float32)