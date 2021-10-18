#!/usr/bin/env python3
#gui/utils.py

import os.path 
from os import walk
from datetime import datetime, timedelta
import re

from tkinter import *
from tkinter import ttk
import tkinter as tk

def win_center(window, w, h):
    window.withdraw()
    window.update_idletasks()  # Update "requested size" from geometry manager
    window.geometry('%dx%d' % (w, h))
    x = (window.winfo_screenwidth() - w) / 2
    y = (window.winfo_screenheight() - h) / 2
    window.geometry('+%d+%d' % (x, y))
    # This seems to draw the window frame immediately, so only call deiconify() after setting correct window position
    window.deiconify()

def message_window(title="", msg=""):
    w = 300
    h = 200
    win = Toplevel()
    win.wm_title(title)
    win_center(win, w, h)
    l = ttk.Label(win, text=msg, wraplength=260)
    l.place(relx=.5, rely=.4, anchor="center")
    b = ttk.Button(win, text="Okay", command=win.destroy)
    b.place(relx=.5, rely=.65, anchor="center")

def refresh_output_dir(img_path):
    dir_basename = os.path.basename(img_path).split("_")
    return f"{dir_basename[0]}_{dir_basename[1]}"

def output_filename_formater(config_dict):
    csv_name = config_dict['ORIG_CSV_NAME'] if "ORIG_CSV_NAME" in config_dict else'larval_counts_$f.csv'
    currdate = datetime.today()
    currdate_str = currdate.strftime(config_dict["FILE_DATE_FORMAT"])
    path = os.path.basename(config_dict["MOST_RECENT_IMG_DIR"])
    dirdate = datetime.strptime(path.split("_")[0], '%y%m%d')
    dirdate_str = dirdate.strftime(config_dict["FILE_DATE_FORMAT"])
    if config_dict["MOST_RECENT_IMG_DIR"]:
        yesterday = dirdate - timedelta(days=1)
    else:
        yesterday = datetime.today() - timedelta(days=1)
    yesterday_str = yesterday.strftime(config_dict["FILE_DATE_FORMAT"])

    if '$d' in csv_name: csv_name = csv_name.replace('$d', currdate_str)
    if '$f' in csv_name: csv_name = csv_name.replace('$f', dirdate_str)
    if '$y' in csv_name: csv_name = csv_name.replace('$y', yesterday_str)
    if csv_name[-4:] != ".csv" : csv_name = csv_name + ".csv"
    return csv_name, currdate, dirdate, yesterday

def parse_config_file(verbose = 0):
    root = os.path.abspath('.')
    config_file_path = os.path .join(root, 'config/config.txt')
    #defaults
    config_dict = {'IMG_DIR': '.', 
                'STATS_CONFIG_PATH': './analysis_config.txt',
                'OUT_DIR': './out', 
                'SAVE_INTMD': False, 
                'SAVE_MASK': True, 
                'SAVE_BBOX': True, 
                'VERBOSE_CSV': True, 
                'CTRL_NAME': 'control', 
                'CTRL_CODE': 'VC',
                'MEDIA_NAME': 'LM (T)', 
                'OPEN_CSV': False,
                'PLATE_NUMBER': 1, 
                'DIR_DATE_FORMAT': '%y%m%d_%H%M%S', 
                'FILE_DATE_FORMAT': '%y%m%d', 
                'PRETTY_DATE_FORMAT': '%m/%d/%Y' ,
                'CSV_NAME': 'larval_counts_$f.csv',
                'DETECTION_NMS_THRESHOLD': 0.6, 
                'DETECTION_MIN_CONFIDENCE': 0.8}
    try:
        with open(config_file_path, 'r') as file:
            line_count = 0
            for line in file:
                line_count += 1
                if line[0]=="#" or not line.strip(): continue
                (key, val) = line.split('=')
                key = key.strip()
                val = val.strip().strip('"').strip("'")
                if key in ['SAVE_INTMD', 'SAVE_MASK', 'SAVE_BBOX', 'VERBOSE_CSV']: val = True if "true" in val.lower() else False
                config_dict[key] = val
    except IOError:
        msg = f'''Configuration file was not found. The config.txt file should be present in the \
same directory as the main executable for this program. Continuing with default values.'''
        message_window(title="Configuration File Problem", msg=msg)
    except ValueError:
        msg = '''Configuration file was found but is not formatted correctly. Please be sure that each line is empty, \
or contains exactly one equal sign, or starts with # (comments). Check line {line_count} of the config file.\
Continuing with default values.'''
        message_window(title="Configuration File Problem", msg=msg)



    config_dict['ORIG_CSV_NAME'] = config_dict['CSV_NAME']
    config_dict["MOST_RECENT_IMG_DIR"] = find_most_recent_img_dir(config_dict)
    dir_basename = os.path.basename(config_dict["MOST_RECENT_IMG_DIR"]).split("_")
    config_dict['IMG_DIR_ID'] = f"{dir_basename[0]}_{dir_basename[1]}"
    config_dict["CSV_NAME"], config_dict["CURR_DATE"], config_dict["DIR_DATE"], config_dict["PREV_DATE"] = output_filename_formater(config_dict)
    
    '''
    Extract a plate number from the file path name if specified. Note that the 
    user-named folder is *not* the deepest directory, and is instead the 
    penultimate directory. We thus need to strip off the final directory. 
    '''
    penultimate_directory = os.path.dirname(config_dict["MOST_RECENT_IMG_DIR"])
    penultimate_directory = os.path.basename(penultimate_directory).lower()
    if "plate" in penultimate_directory:
        print(penultimate_directory)
        plate_split = penultimate_directory.split("plate")
        if len(plate_split) >=1: 
            plate_num = re.match(r'\d', plate_split[1]) #find first digit after the plate split.
            if plate_num : config_dict['PLATE_NUMBER'] = plate_num.group(0)

    if verbose > 2: print(config_dict)

    return config_dict

def find_most_recent_img_dir(config_dict):
    if "IMG_DIR" in config_dict: 
        #find most recent folder that matches Cytation file structure
        dirs = [x[0] for x in walk(config_dict["IMG_DIR"])]
        recent_time = 0
        image_folder = dirs[0]
        for folder in dirs[1:]:
            if os.path.getmtime(folder) > recent_time :
                # print(folder, os.path.getmtime(folder))
                recent_time = os.path.getmtime(folder)
                image_folder = folder
        return image_folder
    return None

