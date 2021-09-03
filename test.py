#!/usr/bin/env python3
#analyze_output.py

from platedriver.plate import Plate
import os
from os import path
from time import time as time
import platedriver.utils as pdu


def parse_config_file():
    root = path.abspath('.')
    config_file_path = path.join(root, 'config.txt')
    config_dict = {}
    try:
        with open(config_file_path, 'r') as file:
            for line in file:
                if line[0]=="#" or not line.strip(): continue
                (key, val) = line.split('=')
                config_dict[key.strip()] = val.strip().strip('"').strip("'")
                # print(key, ",", val)
    except IOError:
        msg = '''Configuration file was not found. The config.txt file should be present in the \
same directory as the main executable for this program.'''
        message_window(title="Configuration File Problem", msg=msg)
    except ValueError:
        msg = '''Configuration file was found but is not formatted correctly. Please be sure that each line is empty, \
or contains exactly one equal sign, or starts with # (comments). '''
        message_window(title="Configuration File Problem", msg=msg)
    return config_dict

if __name__ == '__main__':
    config_dict = pdu.parse_config_file(verbose = 0)
    # if "IMG_DIR" in config_dict: 
    #     #find most recent folder that matches Cytation file structure
    #     dirs = [x[0] for x in os.walk(config_dict["IMG_DIR"])]
    #     recent_time = 0
    #     image_folder = dirs[0]
    #     for folder in dirs[1:]:
    #         if os.path.getmtime(folder) > recent_time :
    #             recent_time = os.path.getmtime(folder)
    #             image_folder = folder
    #     print(image_folder)
        

    CSVOuputName = config_dict['CSV_NAME']
    save_folder = config_dict["OUT_DIR"] if "OUT_DIR" in config_dict else ""
    log_folder = config_dict["LOG_DIR"] if "LOG_DIR" in config_dict else ""


    # plate = Plate("/Users/James/Documents/platereader/")
    plate = Plate(config_dict['MOST_RECENT_IMG_DIR'])

    plate.plate_data.add_row_data("A", compound_ID = "JSK-9999Z", compound_name = "JSK Cmpd", max_conc = 128, rep_num = 1)
    # plate.plate_data.add_row_data("B", compound_ID = "test2", compound_name = "test2", max_conc = 128, rep_num = 1)
    # plate.plate_data.add_row_data("C", compound_ID = "test3", compound_name = "test3", max_conc = 256, rep_num = 1)
    # plate.plate_data.add_row_data("D", compound_ID = "test4", compound_name = "test4", max_conc = 128, rep_num = 1)
    # plate.plate_data.addnd_ID = "test6", compound_name = "test6", max_conc = 8, rep_num = 1)
    # plate.plate_data.add_row_data("G", compound_ID = "test7", compound_name = "test7", max_conc = 128, rep_num = 1)
    # plate.plate_data.add_row_data("H", compound_ID = "test8", compound_name = "test8", max_conc = 128, rep_num = 1)

    plate.run_images(init_only=False, intmd_save=True, save_dir=log_folder)
    plate.detect_larvae()
    plate.save_csv(out_dir = save_folder, filename = CSVOuputName, verbose_csv = config_dict['VERBOSE_CSV'], open_csv=True)