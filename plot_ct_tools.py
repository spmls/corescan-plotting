"""
Created on Friday March 15 16:22 2019

tools to work with tiffs from the GeoTek RXCT scanner

@author: SeanPaul La Selle
"""

import os
import sys
import glob
import re
import tkinter
from tkinter import filedialog
import numpy as np
import xml.etree.ElementTree
from skimage.external import tifffile
from skimage.transform import downscale_local_mean
import matplotlib as matplotlib
import matplotlib.pyplot as plt
import exifread
import cv2


###############################################################################
def ct_xml(filename=''):
    """
    read in GeoTek RXCT xml data to a dictionary
    """
    ## Get directory if not specified in function call
    if not filename:
        filename = filedialog.askopenfilename()
        if not filename:
            sys.exit()
    fname = filename
    ## Parse the xml file
    tree = xml.etree.ElementTree.parse(fname)
    root = tree.getroot()
    ## Add element tags and attributes to a dictionary
    dic = {}
    for elem in root.iter():
        if isinteger(elem.text):
            dic[elem.tag] = int(elem.text)
        elif isfloat(elem.text):
            dic[elem.tag] = float(elem.text)
        else:
            dic[elem.tag] = elem.text
    return dic

###############################################################################
def ct_in(filename='',xml_fname=''):
    """
    read in RXCT data from reconstructed orthogonal views (tiffs)
    """
    ## Get filename if not specified in function call
    if not filename:
        filename = filedialog.askopenfilename()
        if not filename:
            sys.exit()
    im = tifffile.imread(filename)
    # Determine the directory of the file
    directory = os.path.dirname(filename)
    ## Read xml file
    if not xml_fname:
        xml_fname = glob.glob(os.path.splitext(filename)[0]+'*.xml')[0]
    xml_dic = ct_xml(xml_fname)
    return im, xml_dic
###############################################################################
def ct_crop_rotate(ct_data,ct_xml,thresh_val, plot=False):
    """
    crop out "air" in a ct image, and rotate to 90 degrees vertical
    As of 4/22/2019, reading in image as 8 bit to perform thresholding,
    apply rotation and cropping to 32 bit image (cv2 doesn't handle 16 bit well)
    Need to apply a threshold value (check plot to make sure whole core is
    being selected and only background is masked out)
    """
    image_16bit = ct_data.astype('uint16')
    image_8bit = (image_16bit/256).astype('uint8')
    image_32bit = ct_data.astype('uint32')
    # blur to reduce noise
    blur = cv2.GaussianBlur(image_8bit,(3,3),0)
    # Find edges of core using thresholding
    ret, thresh = cv2.threshold(blur, thresh_val, 256, 1)
    thresh = 255-thresh
    # Find contours of threshold image
    contours, hierarchy = cv2.findContours(thresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    # Find the index of the largest contour
    areas = [cv2.contourArea(c) for c in contours]
    max_index = np.argmax(areas)
    cnt=contours[max_index]
    # Find the rotated rectangle that encloses this contour
    rect = cv2.minAreaRect(cnt)
    box = cv2.boxPoints(rect)
    box=np.int0(box)
    #crop image inside bounding box
    scale = 1  # cropping margin, 1 == no margin
    W = rect[1][0]
    H = rect[1][1]

    Xs = [i[0] for i in box]
    Ys = [i[1] for i in box]
    x1 = min(Xs)
    x2 = max(Xs)
    y1 = min(Ys)
    y2 = max(Ys)

    angle = rect[2]
    rotated = False
    if angle < -45:
        angle += 90
        rotated = True

    center = (int((x1+x2)/2), int((y1+y2)/2))
    size = (int(scale*(x2-x1)), int(scale*(y2-y1)))

    M = cv2.getRotationMatrix2D((size[0]/2, size[1]/2), angle, 1.0)
    cropped = cv2.getRectSubPix(np.float32(image_32bit), size, center) # converted to 32 bit
    cropped = cv2.warpAffine(cropped, M, size)

    croppedW = W if not rotated else H
    croppedH = H if not rotated else W

    image = cv2.getRectSubPix(
        cropped, (int(croppedW*scale), int(croppedH*scale)),
                        (size[0]/2, size[1]/2))

    # Update xml so that new images scale correctly, list of dictionaries
    xml = ct_xml.copy()
    xml['physical-width'] = image.shape[1]/xml['pixels-per-CM']
    xml['physical-height'] = image.shape[0]/xml['pixels-per-CM']
    xml['pixel-width'] = image.shape[1]
    xml['scan-lines'] = image.shape[0]

    if plot is True:
        fig = plt.figure(figsize=(11,17))
        ax1=plt.subplot(121)
        ax1.imshow(thresh,cmap = matplotlib.cm.gray)
        ax1.add_patch(matplotlib.patches.Polygon(xy=box,edgecolor='red',
                                            linewidth=1.0, fill=False))
        ax1.set_title('Thresholded image with bounding rectangle')
        ax2=plt.subplot(122)
        ax2.imshow(image,cmap = matplotlib.cm.gray)
        ax2.set_title('Cropped and rotated')

    return image, xml

###############################################################################
def ct_crop_rotate_multiple(ct_data,ct_xml,thresh_val,top_depths,plot=False):
    """
    crops out "air" in a ct image, and rotates to 90 degrees vertical for
    a scan with multiple core segments (2-3 Russian cores, for example)

    As of 4/22/2019, reading in image as 8 bit to perform thresholding,
    apply rotation and cropping to 32 bit image (cv2 doesn't handle 16 bit well)

    ct_data: ct image object read in by "ct_in"
    ct_xml: xml data associated with ct_data
    thresh_val: threshold value used to classify pixel intensities in grayscale
    top_depths: depth of the upper parts of each 'segment' in the scan
    plot: output plots of image processing steps, default is False
    """
    n = np.size(top_depths)
    image_16bit = ct_data.astype('uint16')
    image_8bit = (image_16bit/256).astype('uint8')
    image_32bit = ct_data.astype('uint32')
    # blur to reduce noise
    blur = cv2.GaussianBlur(image_8bit,(5,5),0)

    # Find edges of core using thresholding
    ret, thresh = cv2.threshold(blur, thresh_val, 256, 1)
    thresh = 255-thresh
    # Find contours of threshold image
    contours, hierarchy = cv2.findContours(thresh,cv2.RETR_TREE,
                                            cv2.CHAIN_APPROX_SIMPLE)

    # Sort contours from top to bottom
    boundingBoxes = [cv2.boundingRect(c) for c in contours]
    (contours, boundingBoxes) = zip(*sorted(zip(contours, boundingBoxes),
                                key=lambda b:b[1][1],reverse=False))

    # Find the indices of the largest n contours
    areas = np.array([cv2.contourArea(c) for c in contours])
    max_indices = np.argpartition(areas,-n)[-n:]

    cropped_rotated_images = []
    cropped_rotated_xml = []
    boxes = []
    for i,m in enumerate(max_indices):
        cnt=contours[m]
        # Find the rotated rectangle that encloses this contour
        rect = cv2.minAreaRect(cnt)
        box = cv2.boxPoints(rect)
        box = np.int0(box)
        boxes.append(box)
        #crop image inside bounding box
        scale = 1  # cropping margin, 1 == no margin
        W = rect[1][0]
        H = rect[1][1]

        Xs = [i[0] for i in box]
        Ys = [i[1] for i in box]
        x1 = min(Xs)
        x2 = max(Xs)
        y1 = min(Ys)
        y2 = max(Ys)

        angle = rect[2]
        rotated = False
        if angle < -45:
            angle += 90
            rotated = True

        center = (int((x1+x2)/2), int((y1+y2)/2))
        size = (int(scale*(x2-x1)), int(scale*(y2-y1)))

        M = cv2.getRotationMatrix2D((size[0]/2, size[1]/2), angle, 1.0)
        cropped = cv2.getRectSubPix(np.float32(image_32bit), size, center) # converted to 32 bit
        cropped = cv2.warpAffine(cropped, M, size)

        croppedW = W if not rotated else H
        croppedH = H if not rotated else W

        image = cv2.getRectSubPix(
            cropped, (int(croppedW*scale), int(croppedH*scale)),
                        (size[0]/2, size[1]/2))

        # Add to list of cropped and rotated images
        cropped_rotated_images.append(image)

        # Update xml so that new images scale correctly, list of dictionaries
        xml = ct_xml.copy()
        xml['physical-width'] = image.shape[1]/xml['pixels-per-CM']
        xml['physical-height'] = image.shape[0]/xml['pixels-per-CM']
        xml['pixel-width'] = image.shape[1]
        xml['scan-lines'] = image.shape[0]
        xml['physical-top'] = top_depths[i]*100.
        xml['coreID'] = ct_xml['coreID']+str(" %d-%d cm"
                                %(xml['physical-top']/100,
                                xml['physical-top']/100+xml['physical-height']))
        cropped_rotated_xml.append(xml)

    # Plot steps if plot=True
    if plot is True:
        ## Get screen size for figure
        # root = tkinter.Tk()
        # pix2in = root.winfo_fpixels('1i')
        # screen_width = root.winfo_screenwidth()/pix2in
        # screen_height = root.winfo_screenheight()/pix2in
        # image_h2w = round(ct_xml['physical-height']/ct_xml['physical-width'])
        # fig = plt.figure(figsize=(screen_height/image_h2w, screen_height))
        fig = plt.figure(figsize=(11,17))
        ax1=plt.subplot(131)
        ax1.imshow(ct_data,cmap = matplotlib.cm.gray)
        ax1.set_title('original image')

        ax2=plt.subplot(132)
        ax2.imshow(blur,cmap = matplotlib.cm.gray)
        ax2.set_title('blur')

        ax3=plt.subplot(133)
        ax3.imshow(thresh,cmap=matplotlib.cm.gray)
        for b in boxes:
            ax3.add_patch(matplotlib.patches.Polygon(xy=b,edgecolor='red',
                                                linewidth=1.0,
                                                fill=False))
        ax3.set_title('thresholded and contoured')
    return cropped_rotated_images, cropped_rotated_xml

###############################################################################
def ct_plot(ct_data, ct_xml,vmin=0, vmax=50000):
    """
    plot a CT scan image, using xml data to generate scale

    -dpi sets the saved image resolution at
    """
    ## Load the ct image file and xml data
    im = ct_data
    ct_xml = ct_xml
    ## Get screen size for figure
    root = tkinter.Tk()
    pix2in = root.winfo_fpixels('1i')
    screen_width = root.winfo_screenwidth()/pix2in
    screen_height = root.winfo_screenheight()/pix2in
    image_h2w = round(ct_xml['physical-height']/ct_xml['physical-width'])
    fig = plt.figure(figsize=(screen_height/image_h2w, screen_height))
    ## Plot images
    aspect = 'equal'
    ax = plt.subplot(1, 1, 1)
    plt.imshow(im, aspect=aspect, extent=(0,ct_xml['physical-width'],\
                        ct_xml['physical-top']/100+ct_xml['physical-height'],\
                        ct_xml['physical-top']/100),
                        vmin=vmin,vmax=vmax)
    ax.set_title(ct_xml['coreID'])

    return fig

###############################################################################
def dpi_ct_plot(fig,ct_xml):
    """
    calculate a dpi to preserve the resolution of the original image
    """
    ## Calculate the new dpi from the original resolution and new image size
    orig_h_pixels = ct_xml['scan-lines']
    orig_w_pixels = ct_xml['pixel-width']
    bbox = fig.axes[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    new_w_in, new_h_in = bbox.width, bbox.height
    ## Average the new dpis
    new_dpi = np.round(np.average([orig_w_pixels/new_w_in,
                                    orig_h_pixels/new_h_in]))
    ## Calculate
    return new_dpi

###############################################################################
def ct_histogram(ct_data):
    """
    plot a histogram of the CT scan intensities
    """
    ## Downsample data to about 1 million pixels (typical 2x2 image at 100
    ## micron resolution has about 15 million pixels)
    scale_factor = int(round(np.size(ct_data)/10e5))
    ds_ct_data = downscale_local_mean(ct_data,(scale_factor,scale_factor))
    ## Calculate histogram, uses 500 bins
    hist, bins = np.histogram(ds_ct_data,bins=100)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:])/2
    ## Plot histogram
    fig = plt.figure(figsize=(10,10))
    plt.bar(center, hist, align='center',width=width)
    plt.xlabel('Intensity')
    plt.ylabel('Count')
###############################################################################
def get_reconstructed_slope_offset(recon_slice_folder=''):
    """
    get the slope and offset used to calculate the intensity in a reconstructed
    slice.  this is useful for correcting reconstructed XZ views back to
    unscaled attenuation coefficients. input is the path to the folder
    containing reconstructed slices.
    """
    slice_folder = recon_slice_folder
    slice_image = next(os.path.join(slice_folder,f) for \
                    f in os.listdir(slice_folder)\
                    if os.path.isfile(os.path.join(slice_folder,f)))
    img = open(slice_image,'rb')
    tags = exifread.process_file(img)
    img.close()
    ## Parse the exif tags to get the slope and offset
    string = str(tags['Image ImageDescription'])
    match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
    final_list = [float(x) for x in re.findall(match_number,string)]
    slope, offset = final_list[0],final_list[1]
    return slope, offset
###############################################################################
def remove_reconstructed_intensity(ct_data,slope,offset):
    """
    use the slope and offset values from a reconstructed slice to revert
    CT data back to an unscaled range of attenuation coefficients
    """
    ct_data = ct_data*slope
    ct_data = ct_data+offset
    return ct_data
###############################################################################
def isfloat(x):
    """
    determine if a string can be converted to float
    """
    try:
        a = float(x)
    except ValueError:
        return False
    except TypeError:
        return False
    else:
        return True

###############################################################################
def isinteger(x):
    """
    determine if a string can be converted to an integer
    """
    try:
        a = int(x)
    except ValueError:
        return False
    except TypeError:
        return False
    else:
        return True

###############################################################################
def cm2inch(*tupl):
    """
    convert centimeters to inches
    """
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
###############################################################################
def get_ax_size(ax):
    """
    get the size of a matplotlib figure axis in inches
    """
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    return width,height
