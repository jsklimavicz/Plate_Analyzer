# stats/latex_writer.py
# This file is a component of Plate_Analyzer, which can count Drosophila
# L1 larvae, classify them as alive or dead, and determine dose-repsonse
# information based on live/dead count data. 

# Copyright (C) 2021 James Klimavicz

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import sys
import argparse
import os, errno
import glob
import re
import time
import numpy as np
import math
import pandas as pd
import time



class LatexWriter:
	'''
	Class for making a LaTeX file and corresponding .pdf from all of the graphed data. 
	'''

	#control number of plots on a page.
	gridx = 2
	gridy = 3

	def __init__(self, img_folder, title = "Merlin Bioassay Results"):
		self.image_folder = img_folder
		self.make_header()
		self.make_title(title)
		self.cmpd_graph_list = []
		self.graph_count = 0
		self.make_end()

	def make_header(self):
		'''
		Header section for the LaTeX file.
		'''
		self.header = ["\\documentclass{article}"]
		self.header.append("\\usepackage[letterpaper,left=0.75in, right=0.75in, top=0.75in, bottom=0.75in]{geometry}")
		self.header.append("\\usepackage{graphicx}")
		self.header.append("\\usepackage[justification=centering]{subcaption}")
		self.header.append("\\captionsetup[subfigure]{format=hang,justification=raggedright,singlelinecheck=false}")
		self.header.append("\\usepackage{amsmath,amsthm,amsfonts,amssymb,mathtools, bm}")
		self.header.append("\\usepackage[mmddyyyy,HHmmss]{datetime}")
		self.header.append(f"\\graphicspath{{{{{self.image_folder:s}}}}}")
		self.header.append("\\setcounter{topnumber}{8}")
		self.header.append("\\setcounter{bottomnumber}{8}")
		self.header.append("\\setcounter{totalnumber}{8}")
		self.header.append("\\usepackage{fancyhdr}")
		self.header.append("\\pagestyle{fancy}")
		self.header.append("\\fancyhf{}")
		self.header.append("\\fancyhead[L]{Merlin Bioassay Results}")
		self.header.append("\\fancyfoot[C]{\\thepage}")
		self.header.append("\\fancyhead[R]{Compiled on \\today\\ at \\currenttime}")
		self.header.append("\\renewcommand{\\headrulewidth}{0pt}")
		self.header.append("\n\n")
		self.header.append("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")


	def make_title(self, title):
		'''
		Title section of document 
		'''
		self.title = [f"\\title{{{title:s}}}"]
		self.title.append("\\author{James Klimavicz}")
		self.title.append("\\date{\\today}")
		self.title.append("\\begin{document}")
		self.title.append("\n\n")


	def make_end(self):
		'''
		Portion to close out the LaTeX document. 
		'''
		self.end = ["\n\n"]
		self.end.append("\\end{document}")

	def make_cmpd_caption(self, name, lc50, lc50CI, R2, bio_reps, tech_reps, lcTrue = True):
		'''
		Creates the caption for each image. Retruns the caption as a string. 
		**<Compound Name>** LC50: <LC50> ppm <CI>
		<bio_reps> biological reps; R^2: <R2> 
		'''
		caption = f'\\textbf{{{name}}}' 
		if lcTrue: caption += ' LC$_{50}$: ' + f'{lc50} ppm {lc50CI} \\\\ \n'
		else: caption += ' LC$_{50}$ cannot be estimated with the given data. \\\\ \n'

		if bio_reps == 1: caption += f'{bio_reps} biol. rep; '
		else: caption += f'{bio_reps} biol. reps; '
		if tech_reps == 1: caption += f'{tech_reps} tech. rep; '
		else: caption += f'{tech_reps} tech. reps; '
		caption += f'R$^2$: {R2}'
		#TODO: add reference compound rel. potency
		return caption

	def make_cmpd_graph(self, image_dir, name, lc50, lc50CI, R2, bio_reps, tech_reps, lcTrue = True):
		'''
		Portion that the graphic and caption for a given compound. Input values are 
		from the compound info from curve fitting. 
		'''
		graph_lines = []
		image_width = 1/self.gridx
		img_line = f"\\includegraphics[width = {{0.95\\textwidth}}]{{{image_dir:s}}}"
		name, alpha_name = self.latexify_name(name)
		# print(alpha_name, name)
		caption = f'\\caption*{{{self.make_cmpd_caption(name, lc50, lc50CI, R2, bio_reps, tech_reps, lcTrue = lcTrue)}}}'
		graph_lines.append(f"   \\begin{{subfigure}}{{{image_width:.3f}\\textwidth}}")
		graph_lines.append("      \\centering")
		graph_lines.append("      "+img_line)
		graph_lines.append("      \\vspace{-0.05cm}")
		graph_lines.append("      "+caption)
		graph_lines.append("      \\vspace{0.1cm}")
		graph_lines.append("   \\end{subfigure}%")
		self.graph_count += 1
		self.cmpd_graph_list.append((alpha_name, graph_lines))


	def make_body(self):
		'''
		Driver for including all figures/captions for the output file. 
		'''
		self.body = [] #array of lines. 
		self.graph_count = 0
		row_end = False

		#sort by alphabetizing name, determined by latexify_name
		self.cmpd_graph_list.sort(key=lambda x:x[0]) 

		#include a figure of subfigures with subcaptions. 
		#Layout determined by self.gridx and self.gridy
		for index, (alpha_name,graph_lines) in enumerate(self.cmpd_graph_list):
			if self.graph_count % (self.gridx * self.gridy) == 0:  
				self.body.append("\\begin{figure}[thp!]")
			for i in graph_lines: self.body.append(i)
			self.graph_count += 1
			if self.graph_count % (self.gridx * self.gridy) == 0: 
				self.body.append("\\end{figure}")
				self.graph_count = 0
				self.body.append("\\clearpage")
				self.body.append("\\pagebreak")
			if self.graph_count % self.gridx == 0: 
				self.body.append("\\vspace{-0.1cm}")
		if self.graph_count % (self.gridx * self.gridy) != 0: 
			self.body.append("\\end{figure}")
		self.body.append("\\pagebreak")
		self.body.append("\n\n")

	def write_file(self, out_path, stat_lib):
		'''
		Creates the actual LaTeX file at out_path
		'''
		self.make_body()
		file = open(out_path, 'w')
		file.write("\n".join(self.header))
		file.write("\n".join(self.title))
		file.write("\n".join(self.body))
		self.make_writeup(stat_lib)
		file.write("\n".join(self.writeup))
		file.write("\n".join(self.end))
		file.close()
		self.make(out_path)

	def make(self, out_path):
		'''
		Convert the LaTeX file out_path to a pdf in the same directory. 
		'''
		pwd = os.getcwd()
		os.chdir(os.path.dirname(out_path))
		os.system(f'pdflatex -interaction=nonstopmode {out_path:s} '+\
				'--shell-escape  --enable-write18  > pdf_tex.log 2>&1')
		os.chdir(pwd)

	def latexify_name(self, name):
		'''
		Stylizes the name of a compound to include LaTeX-style letters, e.g. 
		greek letters, italicized stereochem, etc. Currently, R- and S- 
		stereochem designators and lowercase greek letters are implemented. 
		Only looks at the start of the name to ensure that middles of names 
		are not replaced, though this behavior may not be desirable later. 
		Returns the LaTeX-style name and a name without these designators for
		the purpose of alphabetizing.
		'''
		greek_letters = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
							'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', 
							'nu','xi', 'omicron','pi', 'rho', 'sigma','tau', 
							'upsilon','phi', 'chi','psi','omega']
		mod = False
		alpha_name = name

		#stereochem designators
		if name[:2].lower() == "s-": 
			alpha_name = name[2:]
			name = "\\textit{S}-" + alpha_name
			mod = True
		elif name[:2].lower() == "r-": 
			alpha_name = name[2:]
			name = "\\textit{R}-" + alpha_name
			mod = True
		else:
			#Greek letters
			for letter in greek_letters:
				if name[:len(letter)+1].lower() == letter+'-': 
					alpha_name = name[len(letter)+1:]
					name = "$\\" + letter + '$-' + alpha_name.capitalize()
					mod = True
					break
		#Return name and alphabetizing name 
		if len(name)< 5:
			return name.upper(), alpha_name.lower()
		elif mod:
			return name, alpha_name.lower()
		else:
			return name[0].capitalize() + name[1:], alpha_name.lower()


	def make_writeup(self, stat_lib):
		from stats.functionfit import FunctionFit

		def get_method():
			ff = FunctionFit(**stat_lib)
			c_used = ff.use_C_lib
			method = "ls" if "ls" in stat_lib["CURVE_TYPE"].lower() else "ll"
			GSL = "using a C interface to the GNU Scientific Library. "
			py = "using the scipy module. "
			if method == "ls": 
				if c_used:
					solve_id = stat_lib['LS_METHOD']
					if solve_id == 0: m = "Levenberg-Marquardt algorithm"
					elif solve_id == 1: m = "Levenberg-Marquardt algorithm with geodesic acceleration"
					elif solve_id == 2: m = "dogleg algorithm"
					elif solve_id == 3: m = "double dogleg algorithm"
					else: m = "2D subspace algorithm"
				else: m = "trust region reflective algorithm"
				meth_str = f" Least-squares fitting was performed using the {m} algorithm "
			else:
				if c_used:
					solve_id = stat_lib['OPTIM_METHOD']
					if solve_id == 0: m = "Fletcher-Reeves conjugate gradient algorithm"
					elif solve_id == 1: m = "Polak-Ribi\\`{e}re conjugate gradient algorithm"
					elif solve_id == 2: m = "vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm"
					elif solve_id == 3: m = "steepest descent algorithm"
					elif solve_id == 4: m = "Nelder-Mead simplex algorithm"
					elif solve_id == 5: m = "optimized vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm"
					elif solve_id == 6: m = "optimized Nelder-Mead simplex algorithm"
					elif solve_id == 7: m = "optimized Nelder-Mead simplex algorithm with random initialization"
					else: m = "optimized vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm"
				else: m = "vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm"
				meth_str = f" Optimization of the log-likelihood function was performed using the {m} "
			meth_str += GSL if c_used else py
			return meth_str

		def func_str(b2_val):
			return "$$\\frac{" + f"{b2_val}" + "}{\\left(1 + \\exp(b_0+b_1x) \\right)}$$"

		def beta_dist(beta_param):
			if abs(beta_param)<1e-6:
				return "Heldane's prior, the improper prior $\\text{Beta}(0,0)$"
			elif abs(beta_param - 0.5)<1e-6:
				return "Jeffreys' prior, the distribution $\\text{Beta}(0.5,0.5)$"	
			elif abs(beta_param - 1)<1e-6: 	
				return "the Bayes-Laplace prior, the distribution $\\text{Beta}(1,1)$"
			else: return "the distribution $\\text{Beta}"+f"({beta_param},{beta_param})$"


		self.writeup = []
		summ = "Data analysis was performed using the statistics module for" +\
				" the Merlin Data Analysis program. " 
		self.writeup.append(summ)

		boot_summ = "Live/dead counts from the bioassay were used to generate"+\
			" new survival probabilities using a Beta prior. The user-"+\
			f"specified prior is {beta_dist(stat_lib['BETA_PRIOR']):s}, and " +\
			f"{stat_lib['BOOTSTRAP_ITERS']} bootstrap iterations were used. "
		boot_summ += "When either the live count or dead count was equal to 0,"+\
			f"the prior {beta_dist(stat_lib['BETA_PRIOR_0']):s} was used to"+\
			" avoid the sunrise problem. "
		boot_summ += "Correlation between wells in a replicate was modelled"+\
			" by generating multivariate normal random variables with "+\
			"correlation $\\rho = " + f"{stat_lib['RHO']}" + "$, which "+\
			"were then converted to quantiles, and then back-converted to "+\
			"probabilities in the appropriate beta distribution. "
		self.writeup.append(boot_summ)

		curve_summ = "\nEach iteration of bootstrapped dose-response data was fit to "
		prior_summ = ""
		switch = stat_lib["CURVE_TYPE"].lower()
		if switch == 'auto':
			curve_summ += f"the curve {func_str('b_2')}"+\
			" by maximizing the log-likelihood function, where $b_2$ was fixed at 1.0 "+\
			"when background mortailty is low. You may allow $b_2$ to vary freely by "+\
			"setting \\texttt{CURVE_TYPE = ls3} in the \\texttt{analysis\\_config.txt} file, or "+\
			"keep $b_2$ fixed at 1.0 by setting \\texttt{CURVE_TYPE = ls2}. "
			prior_summ += "Priors on parameters were $(b_0, b_1) \\sim \\mathcal{N}(\\bm 0,\\sigma \\mc I_2)$" +\
			"and $b_2 \\sim \\text{Beta}(\\alpha, \\beta)$, where "+\
			"$\\sigma = " + f"{stat_lib['LL_SIGMA']}$, " +\
			"$\\alpha = "+f"{stat_lib['LL_BETA1']}$, " + \
			"and $\\beta = "+f"{stat_lib['LL_BETA2']}$, " +\
			"as defined by \\texttt{LL\\_SIGMA}, \\texttt{LL\\_BETA1}, "+\
			"and \\texttt{LL\\_BETA2} "+\
			"in the \\texttt{analysis\\_config.txt}, respectively." 
		if switch in ["2", "ll2", 2]:
			curve_summ += f"the curve {func_str(1)}"+\
			" by maximizing the log-likelihood function. A scale parameter may be" +\
			" included by setting \\texttt{CURVE_TYPE = ls3} in the \\texttt{analysis\\_config.txt} file. "
			prior_summ += "Priors on parameters were $(b_0, b_1) \\sim "+\
			"\\mathcal{N}(\\bm 0,\\sigma \\mc I_2)$," +\
			" where $\\sigma = " + f"{stat_lib['LL_SIGMA']}$, as defined "+\
			"by \\texttt{LL_SIGMA} in the \\texttt{analysis\\_config.txt} file." 
		elif switch in ["3", "ll3", 3]:
			curve_summ += f"the curve {func_str('b_2')}"+\
			" by maximizing the log-likelihood function. "
			prior_summ += "Priors on parameters were $(b_0, b_1) \\sim \\mathcal{N}(\\bm 0,\\sigma \\mc I_2)$" +\
			"and $b_2 \\sim \\text{Beta}(\\alpha, \\beta)$, where "+\
			"$\\sigma = " + f"{stat_lib['LL_SIGMA']}$, " +\
			"$\\alpha = "+f"{stat_lib['LL_BETA1']}$, " + \
			"and $\\beta = "+f"{stat_lib['LL_BETA2']}$, " +\
			"as defined by \\texttt{LL\\_SIGMA}, \\texttt{LL\\_BETA1}, and \\texttt{LL\\_BETA2} "+\
			"in the \\texttt{analysis\\_config.txt} file, respectively. " 
		elif switch in ["ls3"]:
			curve_summ += f"the curve {func_str('b_2')}"+\
			" using the least squares approach. "
		elif switch in ["ls2"]:
			curve_summ += f"the curve {func_str(1)}"+\
			" using the least squares approach. "
		elif switch in ["best", "aic"]:
			curve_summ += f"the curve {func_str('b_2')}"+\
			" by maximizing the log-likelihood function, where $b_2$ was first fixed at 1.0 "+\
			"(two free parameters) and then allowed to vary freely (three parameters)."+\
			" The Akaike information criterion of the two fits were "+\
			"compared, and the model with the lower AIC was chosen for the iteration. "

		self.writeup.append(curve_summ + prior_summ + get_method() + "\n")

		graph_summ = ""
		if stat_lib['PLOT_ERROR_BARS']:
			graph_summ += "Credible intervals for the data points are shown "+\
					f"at the {round(stat_lib['ERROR_BAR_CI']*100)}\\% level when "+\
					f"fewer than {stat_lib['ERROR_BAR_CUTOFF']} replicates "+\
					"are used. "
		if stat_lib['PLOT_LINE']:
			graph_summ += "The best-fit line is calculated as the median "+\
					"value of all fitted curves at a given concentration. "
		if stat_lib['PLOT_CURVE_ERROR']:
			graph_summ += "The error region for the curve respresents a "+\
					f"{round(stat_lib['CURVE_CI']*100)}\\% "
			graph_summ += "confidence " if "ls" not in switch else "credible "
			graph_summ += "region, as determined by quantiles of predicted"+\
					" survivals at each concentration."
		self.writeup.append(graph_summ)

