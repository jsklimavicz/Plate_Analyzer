# misc/image_augment.py
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

import os
import sys
import math
import random
import numpy as np
import cv2
import json
import skimage.io
import datetime
import matplotlib.pyplot as plt
import random

from joblib import load, Parallel, delayed

import matplotlib.pyplot as plt
import skimage.io
from skimage.util import img_as_ubyte
from skimage.draw import polygon_perimeter
from skimage.exposure import adjust_gamma 
from skimage.transform import rotate, hough_circle, hough_circle_peaks
from skimage.filters import gaussian, prewitt, threshold_triangle
from skimage.color import rgb2gray
from skimage.draw import disk
import copy
from datetime import datetime

import statsmodels.api as sm 
from scipy import stats

#for parallel call to provide self pointer as arg
def unwrap_image_process(arg, **kwarg):
    return DatasetAugmenter.process_images_driver(*arg, **kwarg)

class DatasetAugmenter():
    def __init__(self, args):
        self.args=args
        self.extract_json_info()
        self.args.json_name = os.path.basename(os.path.abspath(self.args.json_path))
        self.processed_count = 0

    def extract_json_info(self):
        if self.args.verbose > 1: print(f"Extracting json annotations from {os.path.abspath(self.args.json_path)}")
        self.annotations = json.load(open(os.path.abspath(self.args.json_path)))
        self.annotations = dict(self.annotations)
        self.img_metadata = self.annotations["_via_img_metadata"]#["_via_img_metadata"]
        self.img_metadata  = list(self.img_metadata .values())  # don't need the dict keys
        self.img_metadata  = [a for a in self.img_metadata if a['regions']] #only keep image dictionaries with polygons

    def make_new_json(self):
        if self.args.verbose > 1: print("Creating new json output")
        dict_keys = [f"{a['filename']}{a['size']}" for a in self.train_metadata]
        self.annotations["_via_img_metadata"] = dict(zip(dict_keys,self.train_metadata))
        self.annotations["_via_image_id_list"] = dict_keys
        train_json = json.dumps(self.annotations, separators=(',', ':'))

        dict_keys = [f"{a['filename']}{a['size']}" for a in self.val_metadata]
        self.annotations["_via_img_metadata"] = dict(zip(dict_keys,self.val_metadata))
        self.annotations["_via_image_id_list"] = dict_keys
        val_json = json.dumps(self.annotations, separators=(',', ':'))
        return train_json, val_json

    def augment_dataset(self):
        self.train_metadata = []
        self.val_metadata = []
        if self.args.verbose > 0: 
            self.print_label_counts()
            print(f"Augmenting each image {self.args.train_number + self.args.val_number + self.args.omit_originals - 1} times...")
        self.total_augmentations = (self.args.train_number + self.args.val_number + self.args.omit_originals - 1) * len(self.img_metadata)
        self.curr_aug = 0
        if self.args.mode=="count": exit()
        self.starttime = datetime.utcnow()
        if self.args.verbose > 0: print(f"Starting augmentation of {len(self.img_metadata)} images.", end = "\r")
        par_list = Parallel(n_jobs=-1, require='sharedmem')(delayed(unwrap_image_process)(a) for a in zip([self]*len(self.img_metadata), self.img_metadata))
        if self.args.verbose > 0:  print(f"Augmentation complete.{'    '*40}")
            
        for item in par_list:
            for annotations in item[0]: self.train_metadata.append(annotations)
            for annotations in item[1]: self.val_metadata.append(annotations)
        if self.args.verbose > 0: print(f"Image processing complete. Creating json file....")
        train_json, val_json = self.make_new_json()
        with open(os.path.join(self.args.train_dir,self.args.json_name), 'w') as file: file.write(train_json)
        with open(os.path.join(self.args.val_dir,self.args.json_name), 'w') as file: file.write(val_json)

        if self.args.verbose > 0: print(f"Done!")

    def print_label_counts(self): 
        def param_est(counts):
            ave = np.mean(counts)
            sd = np.std(counts)
            return ave, sd 

        region_dict = {}
        worm_counts = []
        min_val = 100
        max_val = 0
        total_points = 0
        for anno in self.img_metadata:
            curr_len = len(anno["regions"])
            if curr_len < min_val: min_val = curr_len
            if curr_len > max_val: max_val = curr_len
            worms = 0
            for region in anno["regions"]:
                status = region['region_attributes']['type']
                total_points += len(region['shape_attributes']['all_points_x'])
                if status in ["alive", "dead", "moribund", "egg", "L2"]: worms +=1
                if status in region_dict:
                    region_dict[status]["count"] = region_dict[status]["count"] + 1
                    region_dict[status]["well_status"] = True 
                else:
                    region_dict[status] = {}
                    region_dict[status]["count"] = 1
                    region_dict[status]["well_status"] = True 
                    region_dict[status]["well_count"] = 0
            for status in region_dict.keys():
                if region_dict[status]["well_status"]:
                    region_dict[status]["well_count"] += 1
                    region_dict[status]["well_status"] = False 

            worm_counts.append(worms)
        total_regions = sum([region_dict[a]["count"] for a in region_dict.keys()])
        total_worms = sum([region_dict[a]["count"] for a in ['alive','dead','moribund']])

        ave, sd = param_est(worm_counts)
        model = sm.NegativeBinomial(endog=worm_counts, exog = np.ones(len(worm_counts)), loglike_method="nb2")
        result = model.fit(start_params=[np.log(ave+1), .5], retall=True)
        x = np.arange(min_val-1, max_val+1)
        mu = np.exp(result.params[0]) 
        alpha = result.params[1]
        var = mu + alpha*mu**2 
        p = mu/var 
        n = mu**2/(var - mu)
        print(f"Negative binomial parameters: n={round(n,4)}, p={round(p,4)}, mu={round(mu,4)}")

        fig, ax = plt.subplots(figsize=(12, 9))
        ax.hist(worm_counts, bins = (max_val - min_val -1))
        ax.plot(x, len(self.img_metadata)*stats.nbinom.pmf(x, n, p), 'ro')
        plt.show()
        print(f"{len(self.img_metadata)} images found; {total_regions} total labels found; {total_worms} total worms found. {total_points} total vertices.") 
        print(f"{' '*5}{'Type'.ljust(10)}{'Count'.rjust(5)}{'Wells'.rjust(6)}")
        for key in region_dict.keys(): print(f"{' '*5}{(key + ':').ljust(10)}{str(region_dict[key]['count']).rjust(5)}",
            f"{str(region_dict[key]['well_count']).rjust(5)}")


    def process_images_driver(self, attr):
        new_train_metadata = []
        new_val_metadata = []
        if self.args.verbose > 2: print(attr)
        if not self.args.omit_originals: 
            attr['filename'] = attr['filename'][:-4] + ".npy"
            new_val_metadata.append(attr) #copy original metadata
        base_filename = attr['filename'][:-4]
        image_path = os.path.join(self.args.input_dir, f"{base_filename}.npy")
        image = np.load(image_path)
        self.processed_count += 1
        if not self.args.omit_originals: np.save(os.path.join(self.args.val_dir,f"{base_filename}.npy"), \
                    (image*255).astype(np.uint8)) #copy image
        IM = ImageLabelAugmenter(image, attr, self.args)
        for i in range(self.args.train_number):
            new_im, new_json, tag = IM.augment_image()
            output_img_path = os.path.join(self.args.train_dir, new_json['filename'])
            np.save(output_img_path, (new_im*255).astype(np.uint8))
            new_json["size"] = os.path.getsize(output_img_path)
            new_train_metadata.append(new_json)
            if self.args.verbose > 0: self.get_status()

            # fig, ax = plt.subplots(figsize=(12, 9))
            # ax.imshow(draw_polygons((new_im[:,:,:3]*255).astype(np.uint8),new_json))
            # plt.show()

            # fig, ax = plt.subplots(figsize=(12, 9))
            # ax.imshow(draw_polygons((new_im[:,:,3:]*255).astype(np.uint8),new_json))
            # plt.show()

        for i in range(self.args.val_number + self.args.omit_originals - 1):
            new_im, new_json, tag = IM.augment_image()
            output_img_path = os.path.join(self.args.val_dir, new_json['filename'])
            np.save(output_img_path, (new_im*255).astype(np.uint8))
            new_json["size"] = os.path.getsize(output_img_path)
            new_val_metadata.append(new_json)
            if self.args.verbose > 0: self.get_status()
        # return "new_train_metadata", "new_val_metadata"
        return new_train_metadata, new_val_metadata

    def get_status(self):
        current_image_msg = f"On image {self.processed_count} of {len(self.img_metadata)}. {' '*5}"
        self.curr_aug += 1
        if self.curr_aug % 10 != 1: return
        curr_percent = round(100*self.curr_aug/self.total_augmentations,1)
        n_space = 30
        n_star = math.floor(curr_percent/100*n_space)
        n_space -= n_star
        elapsed_time = (datetime.utcnow() - self.starttime).total_seconds()
        meantime = elapsed_time/self.curr_aug
        timeleft = round((self.total_augmentations - self.curr_aug) * meantime)

        hours, remainder = divmod(timeleft, 3600)
        mins, secs = divmod(remainder, 60)
        secs = f"{secs}s " if (mins >0 or hours>0) else f"{secs} seconds"
        hours = f"{hours}h " if hours>0 else ""
        mins = f"{mins}m " if mins >0 else ""

        time_format = f"{hours}{mins}{secs}"
        print(f"[{'='*n_star}>{'.'*(n_space-1)}] {curr_percent}%  {current_image_msg} Estimated time remaining: {time_format}{' '*20}", end = "\r")



class ImageLabelAugmenter():
    '''
    Helps to augment machine learning datasets by generating variations of labelled images. 
    Images are subjected to a gamma adjustment, a rotation, a rescaling, potential flipping, and potential 
    gaussian bluring. The labelling polygons are subjected to the same roation and rescaling. 
    '''
    def __init__(self, image, json_annotation, args):
        '''
        image: a (n, m, 3) or (n, m) numpy array of an image.
        json_annotation: The json annotation for the file that contains all region/shape information. 
        '''
        self.image = image
        self.json = json_annotation
        self.args = args
        self.__initialize()

    def __initialize(self):
        k = prewitt(gaussian(rgb2gray(self.image[:,:,:3]), sigma = 2))
        thresh = threshold_triangle(k) #set a threshold
        k = k > thresh
        
        hough_radii = np.array([824, 826, 828, 830])*800/1700 #set the radii for Hough transform
        hough_res = hough_circle(k, hough_radii)
        accum, cx, cy, radii  = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=1)

        self.mask = np.zeros(self.image.shape)
        circy, circx = disk((cy[0], cx[0]), radii[0], shape=self.mask.shape)
        self.mask[circy, circx, :] = 1


    def __augment_parameters(self):
        '''Internal method that generates the random variables for the image augmentation.'''
        theta = (random.random() - .5) * self.args.theta #gives a value between in [-180, 180)
        scale = np.random.normal(1, self.args.scale)#gives a value between 0.95 and 1.05 for resizing
        flip = True if random.random() < self.args.flip else False
        gamma = self.args.gamma[0] + random.random() * (self.args.gamma[1] - self.args.gamma[0])  #gives a value between 0.7 and 1.6 for gamma/gain adjustment
        blur = True if random.random() < self.args.blur else False # blurs an image 7% of the time.
        GB_flip = True if random.random() < self.args.GB_flip else False #flips the Green and Blue channels 20% of the time.
        add_noise = True if random.random() < self.args.add_noise else False #makes a noisy image 2% of the time
        if self.args.verbose > 1: print(f"theta: {theta} scale: {scale} flip: {flip} gamma: {gamma} blur: {blur}, GB_flip: {GB_flip}")
        return theta, scale, flip, gamma, blur, GB_flip, add_noise

    def __get_gaussian_value(self):
        '''Internal method for determining a sigma value to use when performing a gaussian blur.'''
        gauss_val = random.random()
        if gauss_val < 0.1:
            sigma = 0.5
        elif gauss_val <0.3:
            sigma = 1
        elif gauss_val<0.55:
            sigma = 1.25
        elif gauss_val <0.72:
            sigma = 1.5
        elif gauss_val <0.85:
            sigma = 1.75
        elif gauss_val <0.92:
            sigma = 2
        elif gauss_val <0.97:
            sigma = 2.25
        else:
            sigma = 2.5
        if self.args.verbose > 1: print(f"sigma: {sigma}")
        return sigma

    def resize_and_crop(self, image, scale):
        '''Resizes an image using a scale parameter, and then either center-crops or pads the image to be the original size.'''
        old_size = image.shape #save old size
        resize_scale = (scale, scale, 1) 
        image = skimage.transform.rescale(image, scale)
        new_size = image.shape
        size_diffx = abs(new_size[0] - old_size[0])
        crop1Dx = (math.floor(size_diffx/2), math.ceil(size_diffx/2))
        size_diffy = abs(new_size[1] - old_size[1])
        crop1Dy = (math.floor(size_diffy/2), math.ceil(size_diffy/2))
        crop_width = (crop1Dx, crop1Dy, (0,0)) if image.ndim == 3 else (crop1Dx, crop1Dy)
        if self.args.verbose > 2: print(f"Resized image dimensions: {new_size}")
        if scale > 1:
            image = skimage.util.crop(image, crop_width, copy=True, order='K') 
            if self.args.verbose > 2: print(f"Cropping with {crop_width} for crop margins")
        else: 
            if self.args.verbose > 2: print(f"Padding with {crop_width} for pad margins")
            image = np.pad(image, crop_width, mode= 'edge') 
        return image

    def gen_affine_mat(self, theta, scale, flip):
        #create affine transformation matrix
        theta_rad = theta /180*math.pi
        if flip:
            affine = scale * np.array([[-math.cos(theta_rad), math.sin(theta_rad)],
                [math.sin(theta_rad), math.cos(theta_rad)]])
        else:
            affine = scale * np.array([[math.cos(theta_rad), -math.sin(theta_rad)],
                [math.sin(theta_rad), math.cos(theta_rad)]])
        if self.args.verbose > 2: print(f"Affine matrix for polygon transformations: {affine}")
        return affine

    def affine_polygon(self, polygon, affine_mat):
        #transform polygon coordinates to a matrix
        coords = np.row_stack((polygon['all_points_y'], polygon['all_points_x']))
        if self.args.verbose > 2: print(f"Old polygon coordinates: {coords}")
        # print(coords)
        #subtract the center from each set of coordinate (rotation happens about the center)
        height, width = self.image.shape[:2]
        center = [height/2, width/2]
        center_coords = np.repeat(np.array([[center[0]],[center[1]]]), coords.shape[1] ,axis = 1) 
        edge_coords = np.repeat(np.array([[height-1],[width-1]]), coords.shape[1] ,axis = 1) 
        nil_coords = np.repeat(np.array([[0],[0]]), coords.shape[1] ,axis = 1) 
        coords = coords - center_coords 
        new_coords = np.matmul(affine_mat, coords) 
        new_coords = new_coords + center_coords 
        #check to make sure that new coordinates do not go out of bounds of image. 
        new_coords = np.where(new_coords >= edge_coords, edge_coords-1, new_coords)
        new_coords = np.where(new_coords < nil_coords, nil_coords, new_coords)

        if self.args.verbose > 2: print(f"New polygon coordinates: {new_coords}")
        return new_coords

    def generate_output_json(self, theta, scale, flip):
        new_json = copy.deepcopy(self.json)
        mat = self.gen_affine_mat(theta, scale, flip)

        # polygons = [r['shape_attributes'] for r in new_json['regions']]
        for r in new_json['regions']:
            p = r['shape_attributes'] 
            new_coords = np.round(self.affine_polygon(p , mat))
            p['all_points_y'] = [int(x) for x in new_coords[0,:].tolist()]
            p['all_points_x'] = [int(x) for x in new_coords[1,:].tolist()]
            r['shape_attributes'] = p
        return new_json

    def img_transform_label(self, theta, scale, flip, gamma, sigma1=None, sigma2=None, GB_flip=False):
        '''Makes a string describing the transformed image.
        String format is: 
        R<clockwise rotation in degrees>
        S<scale factor * 1000>
        G<gamma adjustment * 1000>
        [B<sigma for gaussian blur>] (if performed)
        [F] (if image is flipped)
        '''
        img_transform = f"R{round(theta%360):03d}S{round(scale*1000):04d}G{round(gamma*1000):04d}"
        if sigma1: img_transform = f"{img_transform}B{sigma1*100:n}"
        if sigma2: img_transform = f"{img_transform}N{sigma2*100:n}"
        if flip: img_transform = f"{img_transform}F"
        if GB_flip: img_transform = f"{img_transform}RBG"
        if self.args.verbose > 1: print(f"Image transformation summary: {img_transform}")
        return img_transform


    def augment_image(self):
        '''Performs dataset augmentation.
        Returns:
            new_image: The transformed image.
            new_json: a json description that has updated polygon information and a new filename.
                Updated file name contains the old filename with the transform info appended.
            img_transform: a string description of the transformation performed. 
                See self.img_transform_label for description.
        '''

        theta, scale, flip, gamma, blur, GB_flip, add_noise = self.__augment_parameters() 

        new_image = np.where(self.mask, self.image**gamma, self.image) #skimage.exposure.adjust_gamma(self.image, gamma=gamma)
        new_image = skimage.transform.rotate(new_image, theta, mode= 'edge') 

        # fig, ax = plt.subplots(figsize=(12, 9))
        # ax.imshow(self.image)
        # plt.show()
        if flip: new_image = new_image[::-1, ...]
        sigma1 = None
        sigma2 = None
        if GB_flip: new_image = new_image[:,:,[0,1,2,3,5,4]]
        if blur: 
            sigma1 = self.__get_gaussian_value()
            new_image = np.where(self.mask, gaussian(new_image, sigma = sigma1, multichannel=True), new_image)
        if add_noise:
            sigma2 = self.__get_gaussian_value()
            new_image = np.where(self.mask, np.random.normal(loc = new_image, scale = sigma2/255, size = new_image.shape), new_image)
            new_image = np.where(new_image<0,0, np.where(new_image>1, 1, new_image))
        new_image = self.resize_and_crop(new_image, scale) 
        new_json = self.generate_output_json(theta, scale, flip)
        img_transform = self.img_transform_label(theta, scale, flip, gamma, sigma1, sigma2,GB_flip)
        new_json["filename"] = f"{self.json['filename'][:-4]}_{img_transform}.npy"#update filename to include transformation

        return new_image, new_json, img_transform


def draw_polygons(image, json):
    for region in json["regions"]:
        poly = region['shape_attributes']
        rr, cc = polygon_perimeter(poly['all_points_y'], poly['all_points_x'], shape = image.shape, clip = True)
        image[rr, cc, :] = [1, 0, 0]
    return image


if __name__ == '__main__':
    import argparse

     # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Perform affine transformations on larval images and polygons.')
    parser.add_argument('-d', '--input_dir', required=True,
                        default = os.path.abspath("."),
                        metavar = "/path/to/image/dir",
                        help='Directory of the worm dataset')
    parser.add_argument('-t', '--train_dir', required=True, 
                        default = None,
                        metavar = "/path/to/train/dir",
                        help='Full path to the training directory.')
    parser.add_argument('-v', '--val_dir', required=True, 
                        default = None,
                        metavar = "/path/to/val/dir",
                        help='Full path to the validation directory.')
    parser.add_argument('--mode', required=False, 
                        default = 'count',
                        help='Counting performed for classes in json file for mode "count". Actual augmenting with "augment".')
    parser.add_argument('-j', '--json_path', required=False, 
                        default = './annotations/annotations/worm_annotations.json',
                        metavar = "json_filename.json",
                        help='Full path to the json annotation file.')
    parser.add_argument('-n', '--train_number', required=False,  type=int,
                        default=50,
                        help='Number of augmented images to add to training set')
    parser.add_argument('-m', '--val_number', required=False,  type=int,
                        default=20,
                        help='Number of training images to add to the validation set. If original images are included, they are placed in the validation set.')
    parser.add_argument('--random_seed', required=False,  type=int,
                        default=random.random(),
                        help='Random seed for consistent results.')
    parser.add_argument('--verbose', required=False,  type=int,
                        default=1,
                        help='Determines how much information is printed to stdout.')
    parser.add_argument('--omit_originals', default = False, action='store_true',
                        help='Include this flag if you do not wish for the original images to be included in the validation set.')
    parser.add_argument('--theta', required=False,  type=float,
                        default=180,
                        help='Determines range of rotations. Permissible rotations uniform on are +/- theta degrees (0<theta<=180).')
    parser.add_argument('--scale', required=False,  type=float,
                        default=0.02,
                        help='Determines potential scaling. Scaling is then distributed normally ~N(1,scale).')
    parser.add_argument('--flip', required=False,  type=float,
                        default=0.5,
                        help='Probability of flipping about the y-axis.')
    parser.add_argument('--gamma', required=False,  type=float, nargs=2,
                        default=[0.8,1.2],
                        help='Determines upper and lower limits on gamma correction.')
    parser.add_argument('--blur', required=False,  type=float,
                        default=0.07,
                        help='Determines probability of blurring the image. Blurring is done based on a gaussian blur with sigma defined in the code.')
    parser.add_argument('--GB_flip', required=False,  type=float,
                        default=0,
                        help='Determines probability of swapping the G and B channels of the image.')
    parser.add_argument('--add_noise', required=False,  type=float,
                        default=0.3,
                        help='Determines probability adding gaussian noise to the image. Blurring is done based on a gaussian blur with sigma defined in the code.')
    args = parser.parse_args()


    random.seed(a=args.random_seed)

    dataset_dir = os.path.abspath(args.input_dir)
    train_dir = os.path.abspath(args.train_dir)
    val_dir = os.path.abspath(args.val_dir)

    print(dataset_dir)
    # print(os.path.join(dataset_dir, args.json_name))
    print(train_dir)


    DA =  DatasetAugmenter(args=args)
    DA.augment_dataset()
