#!/usr/bin/env python3
#platereader_gui.py

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
import threading
from threading import *
import time
from pathlib import Path


#############################################################################
#                 For entering background plate info
#############################################################################

class IOFrame(ttk.Frame):
    def __init__(self, container, config = None, scale = 1):
        self.scale = scale
        self.font = ('Arial', 12*self.scale)
        self.config = config
        self.plateID = 1
        self.textbox_width=60
        super().__init__(container)
        # setup the grid layout manager
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=2)
        self.columnconfigure(2, weight=2)
        self.columnconfigure(3, weight=2)
        self.columnconfigure(4, weight=2)
        self.columnconfigure(5, weight=2)
        self.row_num = 0
        self.csv_name_manual_mod = False
        self.date_manual_mod = False
        self.__create_widgets()
        

    def __create_widgets(self):
        s = ttk.Style()
        s.configure("my.TButton", font = self.font)
        full_width = 6
        #Row: Input Image Folder Selector
        self.imageFolderLabel = Label(self, text="Image folder path:")
        self.imageFolderLabel.grid(row=self.row_num, column=0, sticky=E)
        self.imagefolder_var = StringVar()
        self.imagefolder_var.set(self.config['MOST_RECENT_IMG_DIR'])
        self.imageFolderEntry = Entry(self, textvariable= self.imagefolder_var, width=self.textbox_width * self.scale, font = self.font)
        self.imageFolderEntry.grid(row=self.row_num, column=1, columnspan = full_width,  sticky=EW)
        self.imageBrowseButton = ttk.Button(self, text ='Browse', command = lambda: self.image_directory_chooser(), style = "my.TButton")
        self.imageBrowseButton.grid(row=self.row_num, column=full_width+1, sticky=(S,W))
        self.row_num += 1
        Tooltip(self.imageFolderLabel, text='Select the folder containing the Biotek platereader images.')
        Tooltip(self.imageBrowseButton, text='Select the folder containing the Biotek platereader images.')

        #Row: Image output folder selector
        self.outputFolderLabel = Label(self, text="Output parent directory:")
        self.outputFolderLabel.grid(row=self.row_num, column=0, sticky=E)
        self.outputfolder_var = StringVar()
        self.outputfolder_var.set(self.config['OUT_DIR'])
        self.outputFolderEntry = Entry(self, textvariable=self.outputfolder_var, width=self.textbox_width* self.scale, font = self.font)
        self.outputFolderEntry.grid(row=self.row_num, column=1, columnspan = full_width, sticky=EW)
        self.outputBrowseButton = ttk.Button(self, text ='Browse', command = lambda: self.image_output_chooser(), style = "my.TButton")
        self.outputBrowseButton.grid(row=self.row_num, column=full_width+1, sticky=(S,W))
        self.row_num += 1
        msg = 'Select a folder for saving all output images.'
        Tooltip(self.outputFolderLabel, text = msg)
        Tooltip(self.outputBrowseButton, text = msg)

        #Row: Image saving selector
        self.ImageChooserLabel = Label(self, text="Select output options:")
        self.ImageChooserLabel.grid(row=self.row_num, column=0, sticky=E)
        Tooltip(self.ImageChooserLabel, text='Controls which composite images are made.')

        self.INTMDLabel = Label(self, text="INTMD images:")
        self.INTMDLabel.grid(row=self.row_num, column=1, sticky=W)
        self.INTMD_var = IntVar()
        self.INTMD_var.set(self.config["SAVE_INTMD"])
        self.INTMDBox = Checkbutton(self, variable=self.INTMD_var)
        self.INTMDBox.grid(row=self.row_num, column=2,  sticky=W)
        msg = '''Check whis box if you would like to save the RGB images showing motion of larvae over the three frames.'''
        Tooltip(self.INTMDLabel, text=msg)

        self.MASKLabel = Label(self, text="Mask:")
        self.MASKLabel.grid(row=self.row_num, column=3, sticky=W)
        self.MASK_var = IntVar()
        self.MASK_var.set(self.config["SAVE_MASK"])
        self.MASKBox = Checkbutton(self, variable=self.MASK_var)
        self.MASKBox.grid(row=self.row_num, column=4,  sticky=W)
        Tooltip(self.MASKBox, text='Select this option to produce an image wherein identified and classified larvae are'+
            'filled in based on where Mask R-CNN identifies them. Masks are color-coded for live, dead, moribund, or egg.')

        self.BBOXImages = Label(self, text="Bounding Box:")
        self.BBOXImages.grid(row=self.row_num, column=5, sticky=W)
        self.BBOXImage_var = IntVar()
        self.BBOXImage_var.set(self.config["SAVE_BBOX"])
        self.outputVerboseEntry = Checkbutton(self, variable=self.BBOXImage_var)
        self.outputVerboseEntry.grid(row=self.row_num, column=6,  sticky=W)
        Tooltip(self.BBOXImages, text='Select this option to produce an image wherein identified and classified larvae are'+
            'have a box placed around them. Boxes are color-coded for alive, dead, moribund, or egg.')

        self.row_num += 1


        # #Row: Output filename
        # self.outputCSVLabel = Label(self, text="Output filename:")
        # self.outputCSVLabel.grid(row=self.row_num, column=0, sticky=E)
        # self.CSV_Var = StringVar()
        # self.CSV_Var.set(self.config["CSV_NAME"])
        # self.csv_name_manual_mod = False
        # self.outputCSVEntry = Entry(self, textvariable=self.CSV_Var, width=40, font = self.font)
        # self.outputCSVEntry.grid(row=self.row_num, column=1, columnspan = full_width, sticky=EW)
        # # self.outputCSVEntry.insert(END, self.CSV_Var.get())
        # # self.outputCSVBrowseButton = ttk.Button(self, text ='Browse', command = lambda: self.csv_output(), style = "my.TButton")
        # # self.outputCSVBrowseButton.grid(row=self.row_num, column=full_width+1, sticky=(S,W))
        # self.row_num += 1
        # Tooltip(self.outputCSVLabel, text="Select a filename and location for the csv file " +
        #     "of the live/dead larval counts.")
        # # Tooltip(self.outputCSVBrowseButton, text="Select a filename and location for the csv file " +
        # #     "of the live/dead larval counts.")

        # #Row: Verbose output 
        # self.outputVerboseLabel = Label(self, text="Save a verbose csv:")
        # self.outputVerboseLabel.grid(row=self.row_num, column=1, columnspan=3, sticky=W)
        # self.verbose_var = IntVar()
        # self.verbose_var.set(self.config["VERBOSE_CSV"])
        # self.outputVerboseEntry = Checkbutton(self, variable=self.verbose_var)
        # self.outputVerboseEntry.grid(row=self.row_num, column=4,  sticky=W)
        # self.row_num += 1
        # Tooltip(self.outputVerboseLabel, text="Prints a .csv output file containing live/dead estimates for "+
        #     "each object in every well. This option is used mostly for debugging purposes. Saves to a "+
        #     "csv file in the same folder as the normal output but with ""_verbose"" added to the basename.")

        #Row: Date and Plate
        self.DateLabel = Label(self, text="Date:")
        self.DateLabel.grid(row=self.row_num, column=0, sticky=E)
        self.config["DATE"] = self.config["PREV_DATE"]
        date_string = self.config["DATE"].strftime(self.config["PRETTY_DATE_FORMAT"])
        self.DateEntry = Entry(self, textvariable=date_string, width=self.textbox_width* self.scale, font = self.font)
        self.DateEntry.grid(row=self.row_num, column=1, columnspan = full_width, sticky=EW)
        self.DateEntry.insert(END, date_string)
        self.DateButton = ttk.Button(self, text ='Pick date', command = lambda: self.date_picker(), style = "my.TButton")
        self.DateButton.grid(row=self.row_num, column=full_width+1, sticky=(S,W))
        self.row_num += 1
        msg = "Enter the date that the plate was started. Dates may be entered in " + \
            "MM/DD/YY format or chosen with the date selector. May be left blank."
        Tooltip(self.DateLabel, text=msg)
        Tooltip(self.DateButton, text=msg)

        self.PlateIDLabel = Label(self, text="Plate ID:")
        self.PlateIDLabel.grid(row=self.row_num, column=0, sticky=E)
        vlist = [1,2,3,4,5,6,7,8,9,10]
        self.PlateIDButton = ttk.Spinbox(self, values = vlist, wrap = True, format = '%2.0f', width=2, font = self.font)
        self.PlateIDButton.grid(row=self.row_num, column=1, sticky = (S,W))
        self.PlateIDButton.insert(0,1)
        self.row_num += 1
        Tooltip(self.PlateIDLabel, text="ID of the plate. Set to 1 unless more than one plate is run in a day.")

    def reset_manual_indicators(self):
        self.csv_name_manual_mod = False
        self.date_manual_mod = False

    def date_picker(self, scale = 1):
        w = 400 * self.scale
        h = 400 * self.scale
        win = Toplevel()
        win.wm_title("Select Assay Date")
        pdu.win_center(win, w, h)
        cal = Calendar(win, selectmode = 'day', 
                   font = 'Arial 14 bold',
                   foreground = 'Black',
                   selectforeground = 'Red',
                   firstweekday = "sunday",
                   showweeknumbers = False,
                   date_pattern = 'mm/dd/yyyy') 
        cal.pack(pady = 20 * self.scale) 
        def grad_date(): 
            self.config["DATE"] = cal.get_date()
            self.update_date()
            win.destroy()
        Button(win, text = "Select Date", command = grad_date).pack(pady = 20 * self.scale)

    def update_date(self):
        self.DateEntry.delete(0,END)
        self.DateEntry.insert(0,datetime.strptime(self.config["DATE"], '%m/%d/%Y').strftime(self.config["PRETTY_DATE_FORMAT"]))
        self.date_manual_mod = True

    def image_directory_chooser(self):
        new_name = filedialog.askdirectory(initialdir=self.config['IMG_DIR'])
        if new_name:
            self.config['IMG_DIR'] = new_name
            self.imageFolderEntry.delete(0,END)
            self.imageFolderEntry.insert(0,self.config['IMG_DIR'])
            #update predicted date and csv filename 
            self.config['MOST_RECENT_IMG_DIR'] = new_name 
            csv_name, currdate, dirdate, yesterday = pdu.output_filename_formater(self.config)
            # if not self.csv_name_manual_mod:
            #     self.config['CSV_NAME'] = csv_name
            #     self.outputCSVEntry.delete(0,END)
            #     self.outputCSVEntry.insert(0,self.config['CSV_NAME'])
            if not self.date_manual_mod:
                self.config["DIR_DATE"] = dirdate
                self.DateEntry.delete(0,END)
                self.DateEntry.insert(0,self.config["DIR_DATE"].strftime(self.config["PRETTY_DATE_FORMAT"]))


    # def csv_output(self):
    #     new_name = self.param.CSVOuputName = filedialog.asksaveasfilename(initialdir=self.config['OUT_DIR'], 
    #         initialfile = self.config['CSV_NAME'],
    #         filetypes =[('CSV Files', '*.csv')])
    #     if new_name:
    #         self.config['CSV_NAME'] = os.path.basename(new_name)
    #         self.outputCSVEntry.delete(0,END)
    #         self.outputCSVEntry.insert(0,self.config['CSV_NAME'])

    def image_output_chooser(self):
        new_name = filedialog.askdirectory(initialdir=self.config['OUT_DIR'])
        if new_name:
            self.config['OUT_DIR'] = new_name
            self.outputFolderEntry.delete(0,END)
            self.outputFolderEntry.insert(0,self.config['OUT_DIR'])

    def __check_files__(self):
        if not path.exists(self.config['IMG_DIR']): 
            pdu.message_window(msg="Input file does not exist. Please select a valid input file.")
            return False
        try: 
            file_test = open(os.path.join(self.config['OUT_DIR'], self.config['CSV_NAME']), 'a')
            file_test.close()
        except:
            pdu.message_window(msg="Output file cannot be written to the specified path. Please try again.")
            return False
        return True

    def update(self):
        # self.config["VERBOSE_CSV"] = self.verbose_var.get()
        self.config["SAVE_BBOX"] = self.BBOXImage_var.get()
        self.config["SAVE_INTMD"] = self.INTMD_var.get()
        self.config["SAVE_MASK"] = self.MASK_var.get()
        self.config["IMG_DIR"] = self.imagefolder_var.get()
        # self.config["CSV_NAME"] = self.CSV_Var.get()
        self.config["OUT_DIR"] = self.outputfolder_var.get()


#############################################################################
#                         For entering compound data
#############################################################################

class CmpdFrame(ttk.LabelFrame):
    class CmpdInfo:
        def __init__(self , config):
            self.row_names = []
            self.cmpd_names = []
            self.cmpd_codes = []
            self.concentrations = []
            self.reps = []

        def get_concentration(self, row):
            for i in range(0, len(self.row_names)):
                if self.row_names[i] == row: 
                    try: 
                        return self.concentrations[i].get()
                    except ValueError:
                        return float('NaN')
            return float('NaN')

        def get_name(self, row):
            for i in range(0, len(self.row_names)):
                if self.row_names[i] == row: return self.cmpd_names[i].get()
            return ""

        def get_code(self, row):
            for i in range(0, len(self.row_names)):
                if self.row_names[i] == row: return self.cmpd_codes[i].get()
            return ""

        def get_rep(self, row):
            for i in range(0, len(self.row_names)):
                if self.row_names[i] == row: return self.reps[i]
            return ""

    def __init__(self, container, scale = 1, config = None):
        self.scale = scale
        # container.tk.call('tk', 'scaling', '-displayof', '.', 2.0)
        super().__init__(container)
        self.config = config
        self.style = Style(container)
        self.font = ('Arial', 12*self.scale)
        self.style.configure("TLabel", font = self.font)
        # setup the grid layout manager
        self.columnconfigure(0, weight=1)
        self.grid()
        self.textBox = []
        self.cmpd_name =[]
        self.cmpd_code =[]
        self.concentration = []
        self.rep = []
        self.param = self.CmpdInfo(config)
        self.__create_widgets()

    def __create_widgets(self):
        #Column titles
        self.CodeLabel = ttk.Label(self, text="Compound Code")
        self.CodeLabel.grid(row=0, column=1, sticky = S)
        Tooltip(self.CodeLabel, text="Code name/number of compounds. May be left blank.")
        self.NameLabel = ttk.Label(self, text="Compound Name")
        self.NameLabel.grid(row=0, column=2, sticky = S)
        Tooltip(self.NameLabel, text="Name of treatment. May be left blank.")
        self.MaxConcLabel = ttk.Label(self, text="Max Conc.")
        self.MaxConcLabel.grid(row=0, column=3, sticky = S)
        Tooltip(self.MaxConcLabel, text="Highest well concentration of treatment, in ppm. May be left blank")
        self.RepLabel = ttk.Label(self, text="Rep.")
        self.RepLabel.grid(row=0, column=4, sticky = (S,W))

        for i in range(1,9): self.__create_cmpd_entry__(i)
       
        self.ControlLabel = ttk.Label(self, text=f"Control:")
        self.ControlLabel.grid(row=9, column=0, sticky = E)

        self.ControlCodeEntry = ttk.Entry(self, textvariable=self.config["CTRL_CODE"], width=15 * self.scale, font = self.font)
        self.ControlCodeEntry.grid(row=9, column=1, sticky = W)
        self.ControlCodeEntry.insert(END, self.config["CTRL_CODE"])
        Tooltip(self.ControlCodeEntry, text='Name of liquid media type (e.g. LM 2.0).')

        self.ControlNameEntry = ttk.Entry(self, textvariable=self.config["CTRL_NAME"], width=25 * self.scale, font = self.font)
        self.ControlNameEntry.grid(row=9, column=2, sticky = W)
        self.ControlNameEntry.insert(END, self.config["CTRL_NAME"])
        Tooltip(self.ControlNameEntry, text='Name of the control. This should almost always just be "control".')

        # self.button2 = ttk.Button(self, text="Convert", command = lambda: basic_info.convert()).grid(row=10, column=3, sticky=W)


    def __create_cmpd_entry__(self, row_number):
        vlist = [1,2,3,4,5,6,7,8,9,10]
        row_name = chr(row_number + 64)
        self.param.row_names.append(row_name)
        self.param.cmpd_names.append(StringVar())
        self.param.cmpd_codes.append(StringVar())
        self.param.reps.append(1)
        self.param.concentrations.append(DoubleVar())
        self.textBox.append(ttk.Label(self, text=f"Row {row_name:s}:").grid(row=row_number, column=0, sticky = E))
        self.cmpd_code.append(ttk.Entry(self, textvariable=self.param.cmpd_codes[row_number-1], width=15 * self.scale, font = self.font))
        self.cmpd_code[row_number-1].grid(row=row_number, column=1, sticky = W)
        self.cmpd_name.append(ttk.Entry(self, textvariable=self.param.cmpd_names[row_number-1], width=25 * self.scale, font = self.font))
        self.cmpd_name[row_number-1].grid(row=row_number, column=2, sticky = W)
        self.concentration.append(ttk.Entry(self, textvariable=self.param.concentrations[row_number-1], width=10 * self.scale, font = self.font))
        self.concentration[row_number-1].grid(row=row_number, column=3, sticky = W)
        self.rep.append(ttk.Spinbox(self, values = vlist, wrap = True, format = '%2.0f', width=2, font = self.font))
        self.rep[row_number-1].grid(row=row_number, column=4, sticky = W)
        self.rep[row_number-1].insert(0,1)

    def update(self):
        for i in range(0, len(self.param.reps)):
            self.param.reps[i] = int(self.rep[i].get())
        self.config["CTRL_NAME"] = self.ControlNameEntry.get()
        self.config["CTRL_CODE"] = self.ControlCodeEntry.get()


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
        self.nextval = 0
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

        # Button to run the analysis,
        self.Convertbutton = Button(self, text="Run Analysis", command = lambda: self.threaded_plate_driver())
        self.Convertbutton.grid(row=2, column=0, sticky=S)
        self.Convertbutton.config(height = 2)
        Tooltip(self.Convertbutton, text='Analyzes the plate and saves the data when complete.')

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
        # else: 
        #     msg = f'{file_path} already exists. Please select a different name.'
        #     pdu.message_window(title = "CSV file already exists", msg = msg)
        #     return

        self.input_frame.reset_manual_indicators() #reset manual override halt on autoupdating date/csv
        #disable conert button
        self.Convertbutton['state'] = tk.DISABLED
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
        self.monitor(plate_thread)

    def monitor(self, thread):
        if thread.is_alive():
            #update every 0.25 seconds
            self.after(250, lambda: self.monitor(thread))
            self.__update_progress()
            
        else:
            if self.plate.error_message:
                self.prog_label.config(text=self.plate.error_message)
                self.prog_bar['value'] == 0
            else:
                self.prog_label.config(text="Done! You may close this window or process another plate.")
                self.prog_bar['value'] = 100
            self.Convertbutton['state'] = tk.NORMAL
            self.prog_percent.config(text="")

    def __progress(self):
        self.prog_bar = ttk.Progressbar(
            self.master, orient="horizontal",
            length=600, mode="determinate")
        self.prog_bar.grid(column=0, row=4, padx=20, pady=0)
        self.prog_label = Label(self, text="")
        self.prog_label.grid(row=3, column=0, padx=20, pady=0, sticky=EW)
        self.prog_percent = Label(self, text="")
        self.prog_percent.grid(row=5, column=0, padx=20, pady=0, sticky=EW)

    def __update_progress(self):
        self.prog_label.config(text=self.plate.plate_message)
        # print(self.plate.progress)
        if self.auto_prog_bar_on: 
            # if self.plate.progress < 8 and \
            #         self.prog_bar['value'] < 10 and \
            #         self.prog_bar['value'] < self.plate.progress:
            #     self.prog_bar.step(0.15)
            if self.prog_bar['value'] < self.plate.progress - 15:
                self.prog_bar['value'] = self.plate.progress - 10
            elif self.prog_bar['value'] < self.plate.progress:
                max_step = min(0.15, self.plate.progress - self.prog_bar['value'])
                self.prog_bar.step(max_step)
            else:
                self.prog_bar['value'] = self.plate.progress
        elif self.prog_bar['value'] < self.nextval:
            max_step = min(0.15, self.nextval - self.prog_bar['value'])
            self.prog_bar.step(max_step)
        self.prog_percent.config(text=f"{round(self.prog_bar['value']):n}%")


#############################################################################
#                               Run that shizz
#############################################################################

def main():
    app = App()
    app.mainloop()

if __name__ == "__main__": main()

