# Merlin Bioassay Data Analyzer

The purpose of this program is to analyze the images taken from the Biotek Cytation5 platereader. The program also has the ability to perform dose-response curve analysis and calculate LC values.

## Use
For Windows 10, the recommended method is to use the included batch file to invoke the program. The PYTHON_HOME variable should be set to the computer's path to the python executable. 

Alternatively, one may invoke the program directly with `main.py`, which takes no command line arguments. 

All variables used by the program are set in either `./config.txt` or in `./stats/analysis_config.txt`; the former contains arguments for the GUI and general program options and filepaths, while the latter contains options related to data analysis, statistical options, and filepaths related to the statistical analysis output. The `./docs/stats_config` file describes the permitted variables and their actions for the statistical analysis and plotting options, while `./docs/general_config` describes the permitted options for the general program. 

## Requirements

This program requires the use of Python 3. The testing was performed with Python 3.8. The Anaconda distribution is recommended for ease; otherwise, scikit-image, numpy, scipy, tensorflow, tkinter, and other packages may be needed. 

The GUI does require the addition of tkcalendar, which is easily installed using pip. 

Testing has been performed on Ubuntu 20.04 and Windows 10 using Python 3.8.8.

## Graphical User Interface


## Data Collection and Analysis

### Data Collection

The pipeline from image collection to data output is as follows:
1. The platereader takes three images of the well. Images are spaced approximately 0.2 seconds apart, and are saved locally. 
2. 

