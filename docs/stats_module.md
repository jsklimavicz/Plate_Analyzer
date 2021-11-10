# Statistics Module

## Overview 

The statistics module fits dose-response live/dead count data to the curve *s*(*x*) = *b*<sub>2</sub>/(1 + exp(*b*<sub>0</sub> + *b*<sub>1</sub>*x*). 

## Setup and Installation

### C Library

The statistic module fits bootstrapped data to log-logistic models. Functions to perform curve-fitting are written in python for this module with no further configuration necessary. However, because python is relatively slow, the functions to perform the curve-fitting have also been rewritten in C with use of the [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/). While compiling this library is not required, it provides a substantial speed-up in the statistical analysis of data. 

#### Linux
GSL 2.7 was downloaded directly from the main GNU ftp site, configured with the `./configure` file, and installed with `make`/`make install`. The C functions for this module were compiled into a shared object library `cloglik.so` with the command `gcc -fPIC -pedantic -Wall -Werror -shared -o cloglik.so cloglik.c ll_mod.c multimin.c ls.c -lgsl -lgslcblas -lm`. 

#### Mac OS
GSL was installed using [Homebrew](https://brew.sh/) using the command `brew install gsl`. The C functions for this module were compiled into a shared object library `cloglik.so` with the command `gcc -fPIC -pedantic -Wall -Werror -shared -o cloglik.so cloglik.c ll_mod.c multimin.c ls.c -lgsl -lgslcblas -lm`. 

#### Windows 10
Briefly, Mingw-w64 was used to port gcc, and MSYS2 was used to directly install GSL via `pacman -S mingw-w64-x86_64-gsl`. The dynamic link library file `cloglik.dll` was compiled with the command `gcc cloglik.c ll_mod.c multimin.c ls.c -O -pedantic -s -shared -I"C:\path\to\msys64\ming64\include" -L"C:\path\to\msys64\ming64\bin" -lgsl -lgslcblas -lm -o cloglik.dll`. The GSL .dll library must be in the path; alternatively, it can be added manually in the file `./stats/functionfit.py` with `os.add_dll_directory` (search for 'libgsl.dll').

### Output graphics

The plots of all the graphs are saved as both .png graphics and .pdf vector graphics in the output path. The .pdf graphics are used to generate a LaTeX file that can be compiled into a .pdf document is `pdflatex` is in the path; `pdflatex` may be acquired through an suitable [TeX distribution](https://www.latex-project.org/get/).

## Configuration

The statistical analysis and curve-fitting portion of this program can be controlled through the use of numerous keywords described in the [statistical analysis configuration options](./stats_config.md).

## Previewing Data

The statistics tab of the program contains the button `Preview Data`, found under the `Permitted List`. This button opens up an interactive tab that allows the user to see the survival data for each compound for each replicate for both included and excluded data. The user maye also move specific replicates between the include and exclude lists, and a preview dose-response curve of the variety specified in the the `analysis_config.txt` file is generated on the fly. Included data is shown as larger black dots, while excluded data is shown as gray smaller dots; the actively elected replicate is outlined in red as a slightly larger dot. The user is reminded that while this feature is very useful for identifying replicates that have a clear problem (e.g., outlier replicates), the feature should not be used to exclude data without cause to make dose-response curve fit better or look tidier. To keep changes to included/excluded replicates, the user must click `Accept changes`, and a confirmation dialog will confirm that the changes should be made.