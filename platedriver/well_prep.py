# platedriver/well_prep.py
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


import numpy as np
from skimage.filters import sobel, gaussian, threshold_triangle,threshold_otsu, roberts, prewitt, meijering, sato, frangi, hessian
from skimage.exposure import adjust_gamma
from skimage.transform import hough_circle, hough_circle_peaks, resize
from skimage.draw import disk
from skimage.color import rgb2gray
from skimage.util import invert
from skimage.morphology import binary_erosion, square
from skimage.feature import canny
from scipy import ndimage as ndi
import math


import matplotlib.pyplot as plt
import skimage.io


class ImagePreprocessor:


	def __init__(self, wellID, blur_sigma = 3, gamma = 1.2, crop_dim = 800, debug = False, background_level = 0.6):
		self.wellID = wellID
		self.blur_sigma = blur_sigma
		self.gamma = gamma
		self.crop_dim = crop_dim
		self.pad_val = 50
		self.unfocused = False
		self.debug = debug
		self.default_background = background_level

	def true_normalize(self, im):
		ave_val = im[1,1,:]
		for (i,ave) in zip(range(3), ave_val): im[:,:,i] = im[:,:,i] - ave
		im = im-np.min(im)
		im = im/np.max(im)
		ave = im[0,0,0]
		scale = self.default_background/ave
		a = (1 - self.default_background)/(1-ave)
		b = 1-a 
		prev_min = im.min()
		prev_max = im.max()
		# im = scale * im if ave > self.default_background else a * im + b
		gamma = math.log(self.default_background)/math.log(ave)
		im = np.power(im, gamma)
		if self.debug: print(f"Before: {ave} After: {im[0,0,0]}, prev_max: {prev_max}, max: {im.max()}, prev_min: {prev_min}, min: {im.min()}")
		return(im)

	def blur_gamma(self, im):
		im = gaussian(im, sigma = self.blur_sigma)
		im = adjust_gamma(im, gamma=self.gamma, gain=self.gamma)
		return im

	def find_well_edge(self):
		test = -0.4*self.imageA[:,:,1]+1.1*self.imageA[:,:,2]+0.5*self.imageA[:,:,0]
		k = resize(prewitt(gaussian(gaussian(test,4), sigma = 2)),(1024,1024))

		thresh = threshold_triangle(k) #set a threshold
		edges = k > thresh

		hough_radii = np.array([824, 826, 828, 830])/2#set the radii for Hough transform
		hough_res = hough_circle(edges, hough_radii)
		accum, cx, cy, radii  = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=5)

		# print(self.wellID, accum, cx*2, cy*2, radii*2)

		return cy[0]*2, cx[0]*2, radii[0]*2,accum[0]

	def circle_crop(self, im, mask, center_y, center_x, adj):
		miny = int(center_y-adj)
		maxy = int(center_y+adj)
		minx = int(center_x-adj)
		maxx = int(center_x+adj)
		ave_background_red = np.sum(np.where(mask == 1, im[:,:,0],0))/np.sum(mask)
		ave_background_green = np.sum(np.where(mask == 1, im[:,:,1],0))/np.sum(mask)
		ave_background_blue = np.sum(np.where(mask == 1, im[:,:,2],0))/np.sum(mask)
		im[:,:,0] = np.where(mask == 1, im[:,:,0],ave_background_red)
		im[:,:,1] = np.where(mask == 1, im[:,:,1],ave_background_green)
		im[:,:,2] = np.where(mask == 1, im[:,:,2],ave_background_blue)
		# im = np.pad(im, pad_width=((self.pad_val,self.pad_val),(self.pad_val,self.pad_val),(0,0)), mode = 'edge' )
		#if self.debug: print("Pre-crop:", self.wellID, self.imageA.shape, miny, maxy, minx, maxx, self.crop_dim, center_y, center_x, rad, accum)
		im = im[miny:maxy, minx:maxx, :]
		return im

	def crop_and_remove_background(self):
		self.init_dim = self.imageA.shape
		self.tol = self.init_dim[0] + self.pad_val
		center_y, center_x, rad, accum = self.find_well_edge()
		if (center_y+rad > self.tol or center_y-rad < -self.pad_val or center_x+rad > self.tol or center_x-rad < -self.pad_val):
			center_y = 1100
			center_x = 1100
			rad = 850
			self.unfocused = True
		mask = np.zeros(self.imageA.shape[0:2])
		circy, circx = disk((center_y, center_x), rad, shape=mask.shape)
		mask[circy, circx] = 1
		adj = 850

		self.imageA = self.circle_crop(self.imageA, mask, center_y, center_x, adj)
		self.imageB = self.circle_crop(self.imageB, mask, center_y, center_x, adj)
		self.imageC = self.circle_crop(self.imageC, mask, center_y, center_x, adj)

		if self.debug: print("After normalization:", self.wellID, self.imageA.shape)

	def well_preprocess(self, imageA, imageB, imageC):


		self.imageA = imageA.astype(float)
		self.imageB = imageB.astype(float)
		self.imageC = imageC.astype(float)
		self.crop_and_remove_background()
		dims = (self.crop_dim,self.crop_dim,3)

		self.imageA = resize(self.imageA, dims)
		self.imageB = resize(self.imageB, dims)
		self.imageC = resize(self.imageC, dims)


		self.imageA = self.true_normalize(self.imageA)
		self.imageB = self.true_normalize(self.imageB)
		self.imageC = self.true_normalize(self.imageC)

		a = rgb2gray(self.imageA)
		b = rgb2gray(self.imageB)
		c = rgb2gray(self.imageC)

		self.composite = np.dstack((b+c, a+c, a+b))/2
		self.composite = self.true_normalize(self.composite)

	def quick_gray_rescale(self, input_image):
		self.image = rgb2gray(input_image)
		self.crop_and_remove_background()
		self.rescale_factor = self.rescale_val(self.image.shape[0])
		self.image = rescale(self.image, (self.rescale_factor,self.rescale_factor))

	def get_images(self):
		return (self.imageA*255).astype(np.uint8), (self.composite*255).astype(np.uint8)

	def get_composite(self):
		a = self.composite[:,:,1] + self.composite[:,:,2] - self.composite[:,:,0] 
		b = self.composite[:,:,0] + self.composite[:,:,2] - self.composite[:,:,1] 
		c = self.composite[:,:,1] + self.composite[:,:,0] - self.composite[:,:,2] 
		k = prewitt(gaussian(a,1.5))
		k = np.power((k -k.min())/ (k.max() -k.min() ),0.5)
		self.multilayer = np.dstack((self.imageA, k, b, c))
		return (self.multilayer*255).astype(np.uint8)

# #well_rev.py
# import numpy as np
# from skimage.filters import sobel, gaussian, threshold_triangle,threshold_otsu, roberts, prewitt, meijering, sato, frangi, hessian
# from skimage.exposure import adjust_gamma
# from skimage.transform import hough_circle, hough_circle_peaks, rescale
# from skimage.draw import disk
# from skimage.color import rgb2gray
# from skimage.util import invert


# import matplotlib.pyplot as plt
# import skimage.io

# class ImagePreprocessor:

# 	def __init__(self, wellID, blur_sigma = 3, gamma = 1.4, crop_dim = 800, debug = False):
# 		self.wellID = wellID
# 		self.blur_sigma = blur_sigma
# 		self.gamma = gamma
# 		self.crop_dim = crop_dim
# 		self.pad_val = 50
# 		self.unfocused = False
# 		self.debug = debug

# 	def true_normalize(self):
# 		# for i in range(3):
# 		# 	self.image[:,:,i] = self.image[:,:,i]-self.image[:,:,i].min()
# 		# 	self.image[:,:,i] = self.image[:,:,i]/self.image[:,:,i].max()
# 		self.image = self.image-self.image.min()
# 		self.image = self.image/self.image.max()

# 	def blur_gamma(self, im):
# 		im = gaussian(im, sigma = self.blur_sigma)
# 		im = adjust_gamma(im, gamma=self.gamma, gain=self.gamma)
# 		return im

# 	def find_well_edge(self, image_channel):
# 		#apply blur and gamma correction
# 		image_channel = self.blur_gamma(image_channel)
# 		# edges=sobel(image_channel) #use edge detection method
# 		# edges=hessian(image_channel, black_ridges=False) #use edge detection method
# 		edges=prewitt(image_channel) #use edge detection method
# 		thresh = threshold_triangle(edges) #set a threshold
# 		edges = edges > thresh
# 		# fig, ax = plt.subplots(figsize=(12, 9))
# 		# ax.imshow(edges)
# 		# # ax.imshow(invert(edges))
# 		# plt.show()
# 		# exit()
# 		hough_radii = np.array([820, 823])#set the radii for Hough transform
# 		hough_res = hough_circle(edges, hough_radii)
# 		accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=1)
# 		# print(accums[0], cx[0], cy[0], radii[0])
# 		edges = edges/edges.max() #normalize the image
# 		center_y = cy[0] 
# 		center_x = cx[0]
# 		radius = radii[0]
# 		accum = accums[0]
# 		return center_y, center_x, radius, accum

# 	def crop_and_remove_background(self):
# 		if self.image.ndim == 3 :
# 			image_channel = self.image[:,:,0]
# 		else:
# 			image_channel = self.image

# 		self.init_dim = image_channel.shape
# 		self.tol = self.init_dim[0] + self.pad_val
# 		center_y, center_x, rad, accum = self.find_well_edge(image_channel)
# 		if (center_y+rad > self.tol or center_y-rad < -self.pad_val or center_x+rad > self.tol or center_x-rad < -self.pad_val):
# 			center_y = 1100
# 			center_x = 1100
# 			rad = 850
# 			self.unfocused = True
# 		mask = np.zeros(image_channel.shape)
# 		circy, circx = disk((center_y, center_x), rad, shape=mask.shape)
# 		mask[circy, circx] = 1

# 		adj = 900
# 		miny = int(center_y-adj)
# 		maxy = int(center_y+adj)
# 		minx = int(center_x-adj)
# 		maxx = int(center_x+adj)

# 		if self.image.ndim == 3 :
# 			ave_background_red = np.sum(np.where(mask == 1, self.image[:,:,0],0))/np.sum(mask)
# 			ave_background_green = np.sum(np.where(mask == 1, self.image[:,:,1],0))/np.sum(mask)
# 			ave_background_blue = np.sum(np.where(mask == 1, self.image[:,:,2],0))/np.sum(mask)
# 			self.image[:,:,0] = np.where(mask == 1, self.image[:,:,0],ave_background_red)
# 			self.image[:,:,1] = np.where(mask == 1, self.image[:,:,1],ave_background_green)
# 			self.image[:,:,2] = np.where(mask == 1, self.image[:,:,2],ave_background_blue)
# 			#pad image
# 			self.image = np.pad(self.image, pad_width=((self.pad_val,self.pad_val),(self.pad_val,self.pad_val),(0,0)), mode = 'edge' )
# 			if self.debug: print("Pre-crop:", self.wellID, self.image.shape, miny, maxy, minx, maxx, self.crop_dim, center_y, center_x, rad, accum)
# 			self.image = self.image[miny:maxy, minx:maxx, :]
# 		else:
# 			ave_background = np.sum(np.where(mask == 1, self.image[:,:],0))/np.sum(mask)
# 			self.image = np.where(mask == 1, self.image[:,:],ave_background)
# 			self.image = np.pad(self.image, pad_width=((self.pad_val,self.pad_val),(self.pad_val,self.pad_val)), mode = 'edge' )
# 			if self.debug: print("Pre-crop:", self.wellID, self.image.shape, miny, maxy, minx, maxx, self.crop_dim, center_y, center_x, rad, accum)
# 			self.image = self.image[miny:maxy, minx:maxx]

# 		self.true_normalize()
# 		if self.debug: print("After normalization:", self.wellID, self.image.shape)

# 	def well_preprocess(self, imageA, imageB, imageC):
# 		self.image = np.zeros((imageA.shape[0],imageA.shape[1] ,3))
# 		#convert to grayscale if needed
# 		# print('converting to grayscale')
# 		if imageA.ndim==3: imageA = rgb2gray(imageA)
# 		if imageB.ndim==3: imageB = rgb2gray(imageB)
# 		if imageC.ndim==3: imageC = rgb2gray(imageC)

# 		self.image[:,:,0] = imageB/2 + imageC/2
# 		self.image[:,:,1] = imageA/2 + imageC/2
# 		self.image[:,:,2] = imageA/2 + imageB/2

# 		#crop to around the well and apply uniform background
# 		# print('ccroping and removing bkgrd')
# 		self.crop_and_remove_background()
# 		#resize image
# 		# print('Rescaling')
# 		self.rescale_factor = self.rescale_val(self.image.shape[0])
# 		self.image = rescale(self.image, (self.rescale_factor,self.rescale_factor,1))
# 		if self.debug: print("After rescale:", self.wellID, self.image.shape)
# 		# print('done processing')

# 	def rescale_val(self, image_dim):
# 		return self.crop_dim/image_dim

# 	def quick_gray_rescale(self, input_image):
# 		self.image = rgb2gray(input_image)
# 		self.crop_and_remove_background()
# 		self.rescale_factor = self.rescale_val(self.image.shape[0])
# 		self.image = rescale(self.image, (self.rescale_factor,self.rescale_factor))

# 	def get_image(self):
# 		return self.image