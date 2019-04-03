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


###############################################################################
def ct_xml(directory=''):
    """
    read in GeoTek RXCT xml data to a dictionary
    """
    ## Get directory if not specified in function call
    if not directory:
        directory = filedialog.askdirectory()
        if not directory:
            sys.exit()
    dname = os.path.split(directory)[1]
    fname = glob.glob(directory+"/*.xml")[0]
    ## Parse the xml file
    tree = xml.etree.ElementTree.parse(fname)
    root = tree.getroot()
    ## Add element tags and attributes to a dictionary
    dic = {}
    for elem in root.iter():
        if isinteger(elem.text) is True:
            dic[elem.tag] = int(elem.text)
        elif isfloat(elem.text) is True:
            dic[elem.tag] = float(elem.text)
        else:
            dic[elem.tag] = elem.text
    return dic

###############################################################################
def ct_in(filename=''):
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
    xml_dic = ct_xml(directory=directory)
    return im, xml_dic

###############################################################################
def ct_plot(ct_file='',vmin=0, vmax=50000):
    """
    plot a CT scan image, using xml data to generate scale

    -dpi sets the saved image resolution at
    """
    ## Load the ct image file and xml data
    im, ct_xml = ct_in(ct_file)
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
