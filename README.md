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

Tensorflow is used for analyzing the images with the Mask R-CNN algorithm. The process is highly parallelized, and runs fastest with a CUDA-enabled NVIDIA graphics card; however, images can also be processed on multiple CPU cores. 

## Graphical User Interface


## Data Collection and Analysis Pipeline

The pipeline from image collection to data output is as follows:

1. The platereader takes three images of the well. Images are spaced approximately 0.2 seconds apart, and are saved locally. 
2. For each well, these three images are used to make a composite image with six channels as follows:
	1. Red channel of the first image
	2. Green channel of the first image
	3. Blue channel of the first image
	4. Prewitt edge filter of guassian blurred grayscale version of the first image
	5. Grayscale of the second image
	6. Grayscale of the third image
3. The grayscale of the first image is also subjected to a Circle Hough Transform to determine the edge of the well. The area outside the well is set to mean value of the inside of the well for each image channel. The image is also downsized to 800x800 pixels and single precision. 
4. The trained Mask R-CNN algorithm processes each image to classify larvae as dead, alive, moribund, L2 larvae, egg, long-dead, or artifact. 
5. The number of each class is counted, and, if the user has selected the option, composite images are made. In images with bounding boxes or masks shown, the following colors are used:
	- alive: green
	- dead: red
	- moribund: yellow
	- L2 larvae: purple
	- egg: blue
	- long-dead: brown
	- artifact: aqua
6. A .csv file containing the count data for each well is produced. Only the live, dead, and moribund larvae are counted; the moribund larvae are included with the dead larvae. The other four groups are not included in the count data; however, if there are more than five total objects not classified as live, dead, or moribund, a note is included in the .csv file stating the number of eggs, L2 larvae, long-dead larvae, and artifacts in the well. 

## Training Set

The Mask R-CNN algorith was trained using a set of 800 hand-classified wells. 