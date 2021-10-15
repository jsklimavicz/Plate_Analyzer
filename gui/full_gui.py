#!/usr/bin/env python3
#full_gui.py

from tkinter import *
from tkinter import ttk, filedialog, messagebox
from tkinter.ttk import Label, Style
from tkcalendar import Calendar 
import random
import os
from os import path, getcwd, mkdir
from os.path import exists
import platedriver.plate as pc
from datetime import datetime
import tkinter as tk
from platedriver.utils import Tooltip #as Tooltip
import threading
from platedriver.plate import Plate 
from platedriver.platedata import PlateData
import platedriver.utils as pdu
import platform
from threading import *
import time
from stats.main import analyze_data
import multiprocessing
from gui.ioframe import IOFrame as IOFrame
from gui.cmpdframe import CmpdFrame as CmpdFrame


#############################################################################
#                 Complete App and analysis driver
#############################################################################

class App(tk.Tk):
    def __init__(self):
        self.config = pdu.parse_config_file()
        self.scale = self.det_scaling()
        super().__init__()
        self.title('Merlin Bioassay Image Analyzer')
        geom = [700, 700] 
        self.geometry(f'{geom[0]* self.scale*1.5:n}x{geom[1]* self.scale*1:n}')
        # self.resizable(0, 0)
        # windows only (remove the minimize/maximize button)
        # self.attributes('-toolwindow', True)

        # layout on the root window
        self.columnconfigure(0, weight=1)
        self.__create_widgets()
        self.auto_prog_bar_on = False
        self.plate = None
        self.__progress()

    def det_scaling(self):
        operating_sys = platform.system()
        scale = 1
        if operating_sys == 'Linux':
            version = os.uname()[3].split("~")[1].split("-")[0]
            if int(version.split(".")[0])>=18:
                scale = 2
        return scale


    def __create_widgets(self, scale = 1):
        # create the input frame
        self.font = ('Arial', 12*self.scale)
        self.input_frame = IOFrame(self, config = self.config, scale = self.scale)
        self.input_frame.grid(column=0, row=0, padx=10, pady=20)

        # create the button frame
        font = ('Arial', 12*self.scale)
        self.cmpd_frame = CmpdFrame(self, config = self.config, scale = self.scale)
        self.cmpd_frame.grid(column=0, row=1, padx=10, pady=20)
        style = ttk.Style()
        style.configure("correct_font", font = font)
        self.cmpd_frame.configure(borderwidth = 1, relief = 'raised', text = "Plate Configuration")

        # Button to run the analysis of images
        self.Convertbutton = Button(self, text="Run Image Analysis", command = lambda: self.threaded_plate_driver())
        self.Convertbutton.grid(row=2, column=0, sticky=S)
        self.Convertbutton.config(height = 2)
        Tooltip(self.Convertbutton, text='Analyzes the plate images and saves the data when complete.')

        # Button to run the analysis
        self.Statbutton = Button(self, text="Run Statistics", command = lambda: self.statistics_driver())
        self.Statbutton.grid(row=3, column=0, sticky=S)
        self.Statbutton.config(height = 2)
        Tooltip(self.Statbutton, text='Performs statistical analysis of larval count data, including dose-response curve fitting and LC calculations.')

    def __plate_setup(self):
        #update entry values
        # self.input_frame.update()
        self.cmpd_frame.update()
        # try : 
        #make a plate data object
        plate_data = PlateData(dateval = self.input_frame.config["DATE"], 
            plateID = self.input_frame.PlateIDButton.get(), 
            control_name = self.input_frame.config["CTRL_NAME"], 
            control_ID = self.input_frame.config["CTRL_CODE"])
        for i in range(len(self.cmpd_frame.param.row_names)):
            plate_data.add_row_data(row = self.cmpd_frame.param.row_names[i],
                        compound_ID = self.cmpd_frame.param.cmpd_codes[i].get(), 
                        compound_name = self.cmpd_frame.param.cmpd_names[i].get(), 
                        max_conc = float(self.cmpd_frame.param.concentrations[i].get()), 
                        rep_num = self.cmpd_frame.param.reps[i])

        
        self.plate = Plate(config = self.input_frame.config,
                        save_dir = self.out_dir,
                        plate_data=plate_data)

    def __plate_driver(self):
        if self.plate.error_message: return
        self.auto_prog_bar_on = True 
        self.plate.run_images(init_only=False)
        self.plate.detect_larvae()
        self.auto_prog_bar_on = False
        self.plate.save_csv(filename = self.input_frame.config["CSV_NAME"], \
                        verbose_csv = self.input_frame.config["VERBOSE_CSV"], \
                        open_csv=self.input_frame.config["OPEN_CSV"])
        self.prog_bar.stop()

    def threaded_plate_driver(self):
        #check to see if the image input folder exists to be read from:
        self.input_frame.update()

        self.input_frame.config['IMG_DIR_ID'] = pdu.refresh_output_dir(self.input_frame.config['IMG_DIR'])
        self.out_dir = os.path.join(self.input_frame.config['OUT_DIR'], self.input_frame.config['IMG_DIR_ID'])
        if not exists(self.out_dir): mkdir(self.out_dir)

        if not path.isdir(self.input_frame.config['IMG_DIR']):
            text = f'"{self.input_frame.config["IMG_DIR"]}" is not a directory. \nPlease check the image folder path specified above.'
            self.prog_label.config(text=text)
            return
        #check to see if the image output folder exists to be written to:
        if not path.isdir(self.input_frame.config['OUT_DIR']):
            text = f'{self.input_frame.config["OUT_DIR"]} is not a directory. \nPlease check the image output folder path specified above.'
            self.prog_label.config(text=text)
            return
        #check to see if the CSV file exists to be written to:
        file_path = os.path.join(self.out_dir, self.input_frame.config["CSV_NAME"])
        if not path.exists(file_path):
            try: 
                file = open(file_path, 'a')
                file.close()
                if os.stat(file_path).st_size == 0: os.remove(file_path)
            except:
                text = f'"{file_path}" cannot be written to. \nPlease check the output filename specified above.'
                self.prog_label.config(text=text)
                return
        elif os.stat(file_path).st_size == 0: os.remove(file_path)

        self.input_frame.reset_manual_indicators() #reset manual override halt on autoupdating date/csv
        #disable conert button
        self.Convertbutton['state'] = tk.DISABLED
        self.Statbutton['state'] = tk.DISABLED
        #make progress bar
        
        self.__plate_setup()
        
        #make thread for the plate driver
        plate_thread=Thread(target=self.__plate_driver)
        #make a daemon to kill the plate thread if the GUI is closed
        plate_thread.daemon = True
        #start the thread
        plate_thread.start()
        time.sleep(0.5)
        #continue monitoring the thread
        self.monitor(plate_thread, obj = self.plate)

    def statistics_driver(self):
        self.input_frame.update()
        self.Convertbutton['state'] = tk.DISABLED
        self.Statbutton['state'] = tk.DISABLED 
        self.stats_obj = analyze_data(config_path = self.input_frame.config["STATS_CONFIG_PATH"],
                                        MASK_RCNN_SAVE = 'True',
                                        MASK_RCNN_SAVE_NAME = 'larval_counts.csv')
        # stats_thread = multiprocessing.Process(target=self.__stats_processor )
        stats_thread=Thread(target = self.__stats_processor)
        #make a daemon to kill the plate thread if the GUI is closed
        stats_thread.daemon = True
        #start the thread
        stats_thread.start()
        time.sleep(0.5)
        #continue monitoring the thread
        self.monitor(stats_thread, self.stats_obj, speed = 0.25)

    def __stats_processor(self):
        self.auto_prog_bar_on = True 
        self.stats_obj.full_process(new_datapath = self.input_frame.config["OUT_DIR"])
        self.prog_bar.stop()

    def monitor(self, thread, obj, **kwargs):
        if thread.is_alive():
            #update every 0.25 seconds
            self.after(250, lambda: self.monitor(thread, obj, **kwargs))
            self.__update_progress(obj, **kwargs)
            # print(obj.progress)
            # print(obj.message)
            
        else:
            if isinstance(obj, Plate) and obj.error_message:
                self.prog_label.config(text=obj.error_message)
                self.prog_bar['value'] == 0
            else:
                self.prog_label.config(text="Done! You may close this window or process another plate.")
                self.prog_bar['value'] = 100
            self.Convertbutton['state'] = tk.NORMAL
            self.Statbutton['state'] = tk.NORMAL
            self.prog_percent.config(text="")

    def __progress(self):
        self.prog_bar = ttk.Progressbar(
            self.master, orient="horizontal",
            length=600, mode="determinate")
        self.prog_bar.grid(column=0, row=5, padx=20, pady=0)
        self.prog_label = Label(self, text="")
        self.prog_label.grid(row=4, column=0, padx=20, pady=0, sticky=EW)
        self.prog_percent = Label(self, text="")
        self.prog_percent.grid(row=6, column=0, padx=20, pady=0, sticky=EW)

    def __update_progress(self, obj, speed = 0.15, **kwargs):
        self.prog_label.config(text=obj.message)
        # print(self.plate.progress)
        if self.auto_prog_bar_on: 
            if self.prog_bar['value'] < obj.progress - 15:
                self.prog_bar['value'] = obj.progress - 10
            elif self.prog_bar['value'] < obj.progress:
                # print(speed, obj.progress, self.prog_bar['value'])
                max_step = min(speed, obj.progress - self.prog_bar['value'])
                self.prog_bar.step(max_step)
            else:
                self.prog_bar['value'] = obj.progress
        # elif self.prog_bar['value'] < self.nextval:
        #     max_step = min(speed, self.nextval - self.prog_bar['value'])
        #     self.prog_bar.step(max_step)
        self.prog_percent.config(text=f"{round(self.prog_bar['value']):n}%")


#############################################################################
#                               Run that shizz
#############################################################################

def main():
    app = App()
    app.mainloop()

if __name__ == "__main__": main()

