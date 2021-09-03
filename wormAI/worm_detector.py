"""
Mask R-CNN detector for Merlin bioassay.

"""
import os
import sys
import math
import random
import numpy as np
import cv2
import json
import skimage
import datetime
import matplotlib.pyplot as plt
from wormAI.mrcnn.config import Config
from wormAI.mrcnn import model as modellib, utils

def get_weights_path():
    # Path to trained weights
    return os.path.abspath("./wormAI/library/worm_mask-r-cnn.h5")


class WormConfig(Config):
    """Configuration for training on the toy  dataset.
    Derives from the base Config class and overrides some values.
    """
    # Give the configuration a recognizable name
    NAME = "worms"

    IMAGE_CHANNEL_COUNT = 6

    # Image mean (RGB): all channels have approximately equal color except channel 3
    MEAN_PIXEL = np.array([153.0, 152.8, 152.9, 30.0, 152.8, 152.9])

    # We use a GPU with 12GB memory, which can fit two images.
    # Adjust down if you use a smaller GPU.
    IMAGES_PER_GPU = 1

    # Number of classes (including background)
    NUM_CLASSES = 1 + 7  # Background + (live + dead + moribund + egg + L2 + long dead + artifact)

    # Number of training steps per epoch
    STEPS_PER_EPOCH =1000


    # LOSS_WEIGHTS = {'rpn_class_loss': 50.0, 'rpn_bbox_loss': 5.0, 'mrcnn_class_loss': 5.0, 'mrcnn_bbox_loss': 5.0, 'mrcnn_mask_loss': 1.0}
    # LOSS_WEIGHTS = {'rpn_class_loss': 10.0, 'rpn_bbox_loss': 5.0, 'mrcnn_class_loss': 2.0, 'mrcnn_bbox_loss': 2.0, 'mrcnn_mask_loss': 1.0}
    LOSS_WEIGHTS = {'rpn_class_loss': 1.0, 'rpn_bbox_loss': 10.0, 'mrcnn_class_loss': 1.0, 'mrcnn_bbox_loss': 1.0, 'mrcnn_mask_loss': 1.0}

    LEARNING_MOMENTUM = 0.9

    VALIDATION_STEPS=50

    BACKBONE = "resnet101"

    def __init__(self, args):
        super().__init__()
        # Non-maximum suppression threshold for detection
        self.DETECTION_NMS_THRESHOLD = float(args["DETECTION_NMS_THRESHOLD"])
        # Skip detections with < 80% confidence
        self.DETECTION_MIN_CONFIDENCE = float(args["DETECTION_MIN_CONFIDENCE"])


def color_splash_and_bbox(image, detect):
    """Apply color splash effect.
    image: RGB image [height, width, 3]
    mask: instance segmentation mask [height, width, instance count]
    Returns result image.
    """
    # Make a grayscale copy of the image. The grayscale copy still
    # has 3 RGB channels, though.
    gray_splash = image.copy()
    gray_bbox = skimage.color.gray2rgb(image.copy())*255
    name_dict = {"alive": 1,"dead": 2,"moribund": 3,"egg": 4, "L2": 5, "aged_dead": 6, "artifact": 7 }
    color_list = ['green', 'red', 'yellow', 'blue', 'purple', 'brown', 'cyan']
    max_val = 2**round(math.log2(gray_bbox.max()))-1
    color_vals = [  np.array([0,1,0]), 
                    np.array([1,0,0]), 
                    np.array([1,1,0]), 
                    np.array([0,0,1]), 
                    np.array([0.5,0,1]), 
                    np.array([0.65, 0.27, 0.07]),
                    np.array([0,1,1])
                  ]
    color_vals = [a*255 for a in color_vals]
    color_dict = dict(zip(color_list, color_vals))
    #first make splash image.
    labels = np.zeros(gray_splash.shape)
    mask = detect['masks']
    obj_labels = detect["class_ids"]
    for i in range(mask.shape[-1]): labels = np.where(mask[...,i], obj_labels[i],labels)
    labels = labels - 1
    for i in range(len(color_vals)): labels[0,i] = i
    splash = skimage.color.label2rgb(labels, skimage.color.gray2rgb(gray_splash), alpha = 0.2, bg_label=-1,\
        colors = color_list )
    #now make bbox image
    rois = detect['rois']
    bboxes=0
    for roi in rois:
        obj_label = obj_labels[bboxes]
        bboxes +=1
        if obj_label not in [name_dict['artifact']]:
            color = color_list[obj_label-1]
            start = (roi[0],roi[1])
            end = (roi[2],roi[3])
            rr,cc = skimage.draw.rectangle_perimeter(start, end, shape=gray_bbox.shape)
            skimage.draw.set_color(gray_bbox, (rr,cc), color_dict[color])
    return splash, gray_bbox


############################################################
#  Startup
############################################################

if __name__ == '__main__':
    config = WormConfig()
    model = modellib.MaskRCNN(mode="inference", config=config,
                                  model_dir=args.logs)
    model.load_weights(WEIGHTS_PATH, by_name=True)
    detect_and_color_splash(model, image_path=args.image)