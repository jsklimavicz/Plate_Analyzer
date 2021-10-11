from tkinter import *
from tkinter import ttk
import os.path 
from os import walk
from datetime import datetime, timedelta
import tkinter as tk


def parse_config_file(verbose = 0):
    root = os.path.abspath('.')
    config_file_path = os.path .join(root, 'config.txt')
    #defaults
    config_dict = {'IMG_DIR': '.', \
                'STATS_CONFIG_PATH': './analysis_config.txt',
                'OUT_DIR': './out', \
                'SAVE_INTMD': False, \
                'SAVE_MASK': True, \
                'SAVE_BBOX': True, \
                'VERBOSE_CSV': True, \
                'CTRL_NAME': 'control', \
                'CTRL_CODE': 'VC',\
                'MEDIA_NAME': 'LM (T)', \
                'OPEN_CSV': False,\
                'DIR_DATE_FORMAT': '%y%m%d_%H%M%S', \
                'FILE_DATE_FORMAT': '%y%m%d', \
                'PRETTY_DATE_FORMAT': '%m/%d/%Y' ,\
                'CSV_NAME': 'larval_counts_$f.csv',
                'DETECTION_NMS_THRESHOLD': 0.6, \
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
                recent_time = os.path.getmtime(folder)
                image_folder = folder
        return image_folder
    return None

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

def normalize8(self, I):
        mn = I.min()
        mx = I.max()
        I = ((I - mn)/(mx-mn)) * 255
        return I.astype(np.uint8)

class Tooltip:
    '''
    It creates a tooltip for a given widget as the mouse goes on it.

    see:

    http://stackoverflow.com/questions/3221956/
           what-is-the-simplest-way-to-make-tooltips-
           in-tkinter/36221216#36221216

    http://www.daniweb.com/programming/software-development/
           code/484591/a-tooltip-class-for-tkinter

    - Originally written by vegaseat on 2014.09.09.

    - Modified to include a delay time by Victor Zaccardo on 2016.03.25.

    - Modified
        - to correct extreme right and extreme bottom behavior,
        - to stay inside the screen whenever the tooltip might go out on
          the top but still the screen is higher than the tooltip,
        - to use the more flexible mouse positioning,
        - to add customizable background color, padding, waittime and
          wraplength on creation
      by Alberto Vassena on 2016.11.05.

      Tested on Ubuntu 16.04/16.10, running Python 3.5.2
      Tested on Ubuntu 20.04/16.10, running Python 3.8.8
      Tested on Windows 10, running Python 3.8.8

    TODO: themes styles support
    '''

    def __init__(self, widget,
                 *,
                 bg='#FFFFEA',
                 pad=(5, 3, 5, 3),
                 text='widget info',
                 waittime=400,
                 wraplength=250):

        self.waittime = waittime  # in miliseconds, originally 500
        self.wraplength = wraplength  # in pixels, originally 180
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.onEnter)
        self.widget.bind("<Leave>", self.onLeave)
        self.widget.bind("<ButtonPress>", self.onLeave)
        self.bg = bg
        self.pad = pad
        self.id = None
        self.tw = None

    def onEnter(self, event=None):
        self.schedule()

    def onLeave(self, event=None):
        self.unschedule()
        self.hide()

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(self.waittime, self.show)

    def unschedule(self):
        id_ = self.id
        self.id = None
        if id_:
            self.widget.after_cancel(id_)

    def show(self):
        def tip_pos_calculator(widget, label,
                               *,
                               tip_delta=(10, 5), pad=(5, 3, 5, 3)):
            w = widget
            s_width, s_height = w.winfo_screenwidth(), w.winfo_screenheight()
            width, height = (pad[0] + label.winfo_reqwidth() + pad[2],
                             pad[1] + label.winfo_reqheight() + pad[3])

            mouse_x, mouse_y = w.winfo_pointerxy()
            x1, y1 = mouse_x + tip_delta[0], mouse_y + tip_delta[1]
            x2, y2 = x1 + width, y1 + height
            x_delta = x2 - s_width
            if x_delta < 0:
                x_delta = 0
            y_delta = y2 - s_height
            if y_delta < 0:
                y_delta = 0
            offscreen = (x_delta, y_delta) != (0, 0)
            if offscreen:
                if x_delta:
                    x1 = mouse_x - tip_delta[0] - width
                if y_delta:
                    y1 = mouse_y - tip_delta[1] - height
            offscreen_again = y1 < 0  # out on the top
            if offscreen_again:
                # No further checks will be done.
                # TIP: A further mod might automagically augment the wraplength when the tooltip is too high to be
                # kept inside the screen.
                y1 = 0
            return x1, y1

        bg = self.bg
        pad = self.pad
        widget = self.widget

        # creates a toplevel window
        self.tw = tk.Toplevel(widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)

        win = tk.Frame(self.tw,
                       background=bg,
                       borderwidth=0)
        label = tk.Label(win,
                          text=self.text,
                          justify=tk.LEFT,
                          background=bg,
                          relief=tk.SOLID,
                          borderwidth=0,
                          wraplength=self.wraplength)
        label.grid(padx=(pad[0], pad[2]),
                   pady=(pad[1], pad[3]),
                   sticky=tk.NSEW)
        win.grid()
        x, y = tip_pos_calculator(widget, label)
        self.tw.wm_geometry("+%d+%d" % (x, y))

    def hide(self):
        tw = self.tw
        if tw:
            tw.destroy()
        self.tw = None
