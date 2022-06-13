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
import tifffile
from skimage.transform import downscale_local_mean
import matplotlib as matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, MultiCursor
from matplotlib.ticker import MultipleLocator
import scipy
import exifread
import cv2


###############################################################################
def ct_xml(filename=''):
    """
    read in GeoTek RXCT xml data to a dictionary
    """
    # Get directory if not specified in function call
    if not filename:
        tk_root = tkinter.Tk()
        tk_root.wm_withdraw()
        filename = filedialog.askopenfilename()
        tk_root.destroy()
        if not filename:
            sys.exit()
    fname = filename
    # Parse the xml file
    tree = xml.etree.ElementTree.parse(fname)
    root = tree.getroot()
    # Add element tags and attributes to a dictionary
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


def ct_xml_write(new_filename='', orig_filename=''):
    """
    write to GeoTek RXCT xml
    """
    # Get directory if not specified in function call
    if not orig_filename:
        tk_root = tkinter.Tk()
        tk_root.wm_withdraw()
        filename = filedialog.askopenfilename()
        tk_root.destroy()
        if not filename:
            sys.exit()
    orig_fname = orig_filename
    if not new_filename:
        new_fname = orig_filename
    else:
        new_fname = new_filename
    # Parse the xml file
    tree = xml.etree.ElementTree.parse(orig_fname)
    root = tree.getroot()
    # Add element tags and attributes to a dictionary
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


def ct_in(filename='', xml_fname=''):
    """
    read in RXCT data from reconstructed orthogonal views (tiffs)
    """
    # Get filename if not specified in function call
    if not filename:
        tk_root = tkinter.Tk()
        tk_root.wm_withdraw()
        filename = filedialog.askopenfilename()
        tk_root.destroy()
        if not filename:
            sys.exit()
    im = tifffile.imread(filename)
    # Determine the directory of the file
    directory = os.path.dirname(filename)
    # Read xml file
    if not xml_fname:
        xml_fname = glob.glob(os.path.splitext(filename)[0]+'*.xml')[0]
    xml_dic = ct_xml(xml_fname)
    return im, xml_dic
###############################################################################


def auto_crop_rotate(ct_data, ct_xml, thresh_val, plot=False):
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
    blur = cv2.GaussianBlur(image_8bit, (3, 3), 0)
    # Find edges of core using thresholding
    ret, thresh = cv2.threshold(blur, thresh_val, 256, 1)
    thresh = 255-thresh
    # Find contours of threshold image
    image, contours, hierarchy = cv2.findContours(thresh, mode=cv2.RETR_TREE,
                                                  method=cv2.CHAIN_APPROX_SIMPLE)
    # Find the index of the largest contour
    areas = [cv2.contourArea(c) for c in contours]
    max_index = np.argmax(areas)
    cnt = contours[max_index]
    # Find the rotated rectangle that encloses this contour
    rect = cv2.minAreaRect(cnt)
    box = cv2.boxPoints(rect)
    box = np.int0(box)
    # crop image inside bounding box
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
    cropped = cv2.getRectSubPix(np.float32(image_32bit), size, center)  # converted to 32 bit
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
        fig = plt.figure(figsize=(11, 17))
        ax1 = plt.subplot(121)
        ax1.imshow(thresh, cmap=matplotlib.cm.gray)
        ax1.add_patch(matplotlib.patches.Polygon(xy=box, edgecolor='red',
                                                 linewidth=1.0, fill=False))
        ax1.set_title('Thresholded image with bounding rectangle')
        ax2 = plt.subplot(122)
        ax2.imshow(image, cmap=matplotlib.cm.gray)
        ax2.set_title('Cropped and rotated')

    return image, xml

###############################################################################


def ct_crop_rotate_multiple(ct_data, ct_xml, thresh_val, top_depths, plot=False):
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
    blur = cv2.GaussianBlur(image_8bit, (5, 5), 0)

    # Find edges of core using thresholding
    ret, thresh = cv2.threshold(blur, thresh_val, 256, 1)
    thresh = 255-thresh
    # Find contours of threshold image
    image, contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE,
                                                  cv2.CHAIN_APPROX_SIMPLE)

    # Sort contours from top to bottom
    boundingBoxes = [cv2.boundingRect(c) for c in contours]
    (contours, boundingBoxes) = zip(*sorted(zip(contours, boundingBoxes),
                                            key=lambda b: b[1][1], reverse=False))

    # Find the indices of the largest n contours
    areas = np.array([cv2.contourArea(c) for c in contours])
    max_indices = np.argpartition(areas, -n)[-n:]

    cropped_rotated_images = []
    cropped_rotated_xml = []
    boxes = []
    for i, m in enumerate(max_indices):
        cnt = contours[m]
        # Find the rotated rectangle that encloses this contour
        rect = cv2.minAreaRect(cnt)
        box = cv2.boxPoints(rect)
        box = np.int0(box)
        boxes.append(box)
        # crop image inside bounding box
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
        cropped = cv2.getRectSubPix(np.float32(image_32bit), size, center)  # converted to 32 bit
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
                                             % (xml['physical-top']/100,
                                                xml['physical-top']/100+xml['physical-height']))
        cropped_rotated_xml.append(xml)

    # Plot steps if plot=True
    if plot is True:
        # Get screen size for figure
        # root = tkinter.Tk()
        # pix2in = root.winfo_fpixels('1i')
        # screen_width = root.winfo_screenwidth()/pix2in
        # screen_height = root.winfo_screenheight()/pix2in
        # image_h2w = round(ct_xml['physical-height']/ct_xml['physical-width'])
        # fig = plt.figure(figsize=(screen_height/image_h2w, screen_height))
        fig = plt.figure(figsize=(11, 17))
        ax1 = plt.subplot(121)
        ax1.imshow(ct_data, cmap=matplotlib.cm.gray)
        ax1.set_title('original image')
        ax1.set_ylabel('pixels')

        ax2 = plt.subplot(122)
        ax2.imshow(thresh, cmap=matplotlib.cm.gray)
        for b in boxes:
            ax2.add_patch(matplotlib.patches.Polygon(xy=b, edgecolor='red',
                                                     linewidth=1.0,
                                                     fill=False))
        ax2.set_title('thresholded and contoured')
    return cropped_rotated_images, cropped_rotated_xml

################################################################################


def crop_custom(ct_data, ct_xml, units='cm', bbox=None, plot=False):
    """
    extract a subset of the input image using a bounding box defined by:
    [x0,x1,y0,y1] where x0,y0 are top left. x1,y1 are bottom right.
    by default, coordinates are in centimeters, but can be defined in pixels if
    units='pixels'.
    """

    if bbox == None:
        print('need to define bbox, see help')
        return
    else:
        x0, x1, y0, y1 = bbox[0], bbox[1], bbox[2], bbox[3]

    # convert from cm to pixels and visa versa if necessary
    cm2pix = ct_xml['pixels-per-CM']
    top = ct_xml['physical-top']/100.
    if units == 'cm':
        xp0, xp1 = int(x0*cm2pix), int(x1*cm2pix)
        yp0 = int(np.abs(top-y0)*cm2pix)
        yp1 = int(np.abs(top-y1)*cm2pix)
    elif units == 'pixels':
        xp0, xp1, yp0, yp1 = x0, x1, y0, y1
        x0, x1 = xp0/cm2pix, xp1/cm2pix
        y0, y1 = top + yp0/cm2pix, top + yp1/cm2pix

    # Extract
    ct_crop = ct_data[yp0-1:yp1-1, xp0-1:xp1-1]

    # Plot original
    if plot == True:
        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.imshow(ct_data, aspect='equal', extent=(0, ct_xml['physical-width'],
                                                    ct_xml['physical-top']/100 +
                                                    ct_xml['physical-height'],
                                                    ct_xml['physical-top']/100))
        ax1.plot([x0, x0, x1, x1, x0], [y0, y1, y1, y0, y0], 'ro-')
        ax1.set_title('original')

        ax2.imshow(ct_crop, aspect='equal', extent=(x0, x1, y1, y0))
        ax2.set_ylim(y1, y0)
        ax2.set_title('cropped')

    # Update xml so that new images scale correctly, list of dictionaries
    xml = ct_xml.copy()
    xml['physical-width'] = x1-x0
    xml['physical-height'] = y1-y0
    xml['pixel-width'] = ct_crop.shape[1]
    xml['scan-lines'] = ct_crop.shape[0]
    xml['physical-top'] = y0*100.
    xml['coreID'] = str("%s %d-%d cm"
                        % (xml['coreID'],
                           xml['physical-top']/100,
                           xml['physical-top']/100+xml['physical-height']))
    return ct_crop, xml

###############################################################################


def ct_plot(ct_data, ct_xml, vmin=0, vmax=50000, cmap=matplotlib.cm.inferno):
    """
    plot a CT scan image, using xml data to generate scale

    -dpi sets the saved image resolution at
    """
    # Load the ct image file and xml data
    im = ct_data
    ct_xml = ct_xml
    # Get screen size for figure
    root = tkinter.Tk()
    pix2in = root.winfo_fpixels('1i')
    screen_width = root.winfo_screenwidth()/pix2in
    screen_height = root.winfo_screenheight()/pix2in
    image_h2w = round(ct_xml['physical-height']/ct_xml['physical-width'])
    fig = plt.figure(figsize=(screen_height/image_h2w, screen_height))
    # Plot images
    aspect = 'equal'
    ax = plt.subplot(1, 1, 1)
    plt.imshow(im, aspect=aspect, extent=(0, ct_xml['physical-width'],
                                          ct_xml['physical-top']/100+ct_xml['physical-height'],
                                          ct_xml['physical-top']/100),
               vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_title(ct_xml['coreID'])
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(1))

    return fig

###############################################################################


def dpi_ct_plot(fig, ct_xml):
    """
    calculate a dpi to preserve the resolution of the original image
    """
    # Calculate the new dpi from the original resolution and new image size
    orig_h_pixels = ct_xml['scan-lines']
    orig_w_pixels = ct_xml['pixel-width']
    bbox = fig.axes[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    new_w_in, new_h_in = bbox.width, bbox.height
    # Average the new dpis
    new_dpi = np.round(np.average([orig_w_pixels/new_w_in,
                                   orig_h_pixels/new_h_in]))
    # Calculate
    return new_dpi

###############################################################################


def ct_histogram(ct_data):
    """
    plot a histogram of the CT scan intensities
    """
    # Downsample data to about 1 million pixels (typical 2x2 image at 100
    # micron resolution has about 15 million pixels)
    scale_factor = int(round(np.size(ct_data)/10e5))
    ds_ct_data = downscale_local_mean(ct_data, (scale_factor, scale_factor))
    # Calculate histogram, uses 500 bins
    hist, bins = np.histogram(ds_ct_data, bins=100)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:])/2
    # Plot histogram
    fig = plt.figure(figsize=(10, 10))
    plt.bar(center, hist, align='center', width=width)
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
    slice_image = next(os.path.join(slice_folder, f) for
                       f in os.listdir(slice_folder)
                       if os.path.isfile(os.path.join(slice_folder, f)))
    img = open(slice_image, 'rb')
    tags = exifread.process_file(img)
    img.close()
    # Parse the exif tags to get the slope and offset
    string = str(tags['Image ImageDescription'])
    match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
    final_list = [float(x) for x in re.findall(match_number, string)]
    slope, offset = final_list[0], final_list[1]
    return slope, offset
###############################################################################


def remove_reconstructed_intensity(ct_data, slope, offset):
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
    return width, height

################################################################################


def extract_profile(ct_data, ct_xml, transect=None):
    """
    extract pixel values from ct image.
    if transect is not specified, take a line down the center
    if transect is specified, should be of form: [x0,x1,y0,y1] (in cm)
    """
    cm2pix = ct_xml['pixels-per-CM']
    top = ct_xml['physical-top']/100
    if transect == None:
        x0, y0 = ct_xml['physical-width']/2, top
        x1, y1 = ct_xml['physical-width']/2, top+ct_xml['physical-height']
        length = ct_xml['scan-lines']
    else:
        x0, y0 = transect[0], transect[2]
        x1, y1 = transect[1], transect[3]
        length = ct_xml['scan-lines']

    xp0, xp1 = x0*cm2pix, x1*cm2pix
    yp0, yp1 = np.abs(top-y0)*cm2pix, np.abs(top-y1)*cm2pix
    x = np.linspace(xp0, xp1, length)
    y = np.linspace(yp0, yp1, length)
    zi = ct_data[y.astype(np.int)-1, x.astype(np.int)-1]
    y_cm = top+y/cm2pix

    return zi, y_cm, [x0, x1, y0, y1]


################################################################################
def lamina_fft_filter(signal, high_freq_thresh=0.03, plot=False):
    """
    calculates the fft of a given signal (e.g. profile of intensities down core)
    and filters out high frequencies above "high_freq_thresh". returns the
    filtered signal.
    5/11/2020 NEED TO FIGURE OUT IMAGINARY NUMBERS ISSUE
    """
    # Calculate FFT of signal and the power
    sig_fft = np.fft.fft(signal)
    power = np.abs(sig_fft)
    sample_freq = np.fft.fftfreq(signal.size, d=1)

    # Find peak Frequency
    pos_mask = np.where(sample_freq > 0)
    freqs = sample_freq[pos_mask]
    peak_freq = freqs[power[pos_mask].argmax()]

    # Remove high frequencies
    high_freq_fft = sig_fft.copy()
    high_freq_fft[np.abs(sample_freq) > high_freq_thresh] = 0
    filtered_sig = np.fft.ifft(high_freq_fft).real

    if plot == True:
        fig = plt.figure(figsize=(11, 8.5))
        ax1 = fig.add_subplot(2, 1, 1)
        ax1.plot(sample_freq, power)
        ax1.set_xlabel('Frequency [Hz]')
        ax1.set_ylabel('Power')

        ax2 = fig.add_subplot(2, 1, 2)
        ax2.plot(signal)
        ax2.plot(filtered_sig, 'r-')
        ax2.set_title('Peak Frequency')

    return filtered_sig

################################################################################


def lamina_peaks(signal, prominence=500, width=5):
    """
    a simple peak finder using scipy.signal.find_peaks.
    define the threshold prominence (in intensity units) and minimum peak width
    (in pixels) for identifying peaks.
    """
    peaks, peaks_prop = scipy.signal.find_peaks(signal,
                                                prominence=(prominence, None),
                                                width=width)

    valleys, valleys_prop = scipy.signal.find_peaks(-1.*signal,
                                                    prominence=(prominence, None),
                                                    width=width)
    return peaks, valleys

################################################################################


def autocount_laminae_fft(ct_data, ct_xml, layout='horizontal', transect=None):
    """
    automate counting of laminae
    if transect is not specified, it will extract a line from the center of the
    full length of the core. If specified, transect should be in the form:
    [x0,x1,y0,y1] (x0 and x1 will be equal for a straight vertical line)
    """
    # Load the ct image file and xml data
    im = ct_data
    ct_xml = ct_xml

    # Get screen size for figure
    root = tkinter.Tk()
    pix2in = root.winfo_fpixels('1i')
    screen_width = root.winfo_screenwidth()/pix2in
    screen_height = root.winfo_screenheight()/pix2in
    image_h2w = round(ct_xml['physical-height']/ct_xml['physical-width'])
    pix2cm = ct_xml['pixels-per-CM']
    if layout == 'horizontal':
        fig, (ax1, ax) = plt.subplots(nrows=2,
                                      figsize=(screen_width*0.8, screen_height*0.8),
                                      sharex=True)
        plt.subplots_adjust(bottom=0.25)
    else:
        fig, (ax1, ax) = plt.subplots(ncols=2,
                                      figsize=(screen_width*0.8, screen_height/image_h2w*10),
                                      sharey=True)

    # Extract transect
    zi, y_cm, profile_points = extract_profile(ct_data, ct_xml, transect=transect)
    x0_cm, x1_cm = profile_points[0], profile_points[1]
    y0_cm, y1_cm = profile_points[2], profile_points[3]

    # Calculate vmin and vmax based on the transect values
    vmin0 = np.min(zi)
    vmax0 = np.max(zi)
    vmin = np.min(ct_data)
    vmax = np.max(ct_data)

    # filter out high frequencies with fourier transform
    freq_thresh = 0.03
    zi_filter = lamina_fft_filter(zi, high_freq_thresh=freq_thresh, plot=False)

    # get peaks
    prominence = 500
    width = 5
    peaks, valleys = lamina_peaks(zi_filter, prominence=prominence,
                                  width=width)

    # Plot images
    aspect = 'equal'
    if layout == 'horizontal':
        ct_plot = ax.imshow(np.rot90(im, k=1, axes=(1, 0)), aspect=aspect,
                            extent=(ct_xml['physical-top']/100 +
                                    ct_xml['physical-height'],
                                    ct_xml['physical-top']/100,
                                    0,
                                    ct_xml['physical-width']),
                            vmin=vmin0, vmax=vmax0)
        transect = ax.plot([y0_cm, y1_cm], [x0_cm, x1_cm], 'ro-')
        ax.set_xlim(y0_cm, y1_cm)
        ax.set_ylabel('Width [cm]')
        ax.set_xlabel('Depth [cm]')
    else:
        ct_plot = ax.imshow(im, aspect=aspect,
                            extent=(0, ct_xml['physical-width'],
                                    ct_xml['physical-top']/100+ct_xml['physical-height'],
                                    ct_xml['physical-top']/100),
                            vmin=vmin0, vmax=vmax0)
        transect = ax.plot([x0_cm, x1_cm], [y0_cm, y1_cm], 'ro-')
        ax.set_ylim(y1_cm, y0_cm)
        ax.set_xlabel('Width [cm]')

    ax.set_title(ct_xml['coreID'])

    # Plot transect values
    if layout == 'horizontal':
        ax1.plot(y_cm, zi, 'k-', linewidth=1.0)
        ax1.set_xlabel('Depth [cm]')
        ax1.set_ylabel('Intensity')
        ax1.set_xlim(y0_cm, y1_cm)
        ax1.tick_params(labelbottom=True)
        peaks_plot, = ax1.plot(y_cm[peaks], zi_filter[peaks], 'rx')
        valleys_plot, = ax1.plot(y_cm[valleys], zi_filter[valleys], 'bx')
        fft_plot, = ax1.plot(y_cm, zi_filter, 'r--', linewidth=1.0)
    else:
        ax1.plot(zi, y_cm, 'k-', linewidth=1.0)
        ax1.set_ylabel('Depth [cm]')
        ax1.set_xlabel('Intensity')
        ax1.set_ylim(y1_cm, y0_cm)
        peaks_plot, = ax1.plot(zi_filter[peaks], y_cm[peaks], 'rx')
        valleys_plot, = ax1.plot(zi_filter[valleys], y_cm[valleys], 'bx')
        fft_plot, = ax1.plot(zi_filter, y_cm, 'r--', linewidth=1.0)

    # Vmin and vmax slider
    ax_max = plt.axes([0.25, 0.21, 0.65, 0.01])
    ax_min = plt.axes([0.25, 0.16, 0.65, 0.01])

    smin = Slider(ax_min, 'min intensity', vmin, vmax, valinit=vmin0,
                  dragging=True, valfmt='%i')
    smax = Slider(ax_max, 'max intensity', vmin, vmax, valinit=vmax0,
                  dragging=True, valfmt='%i')

    # Fourier slider
    ax_fft = plt.axes([0.25, 0.11, 0.65, 0.01])
    sfft = Slider(ax_fft, 'fft filter threshold',
                  0, 0.1, dragging=True, valfmt='%.2f',
                  valinit=freq_thresh)

    # Peak finding sliders
    ax_prom = plt.axes([0.25, 0.06, 0.65, 0.01])
    sprom = Slider(ax_prom, 'peak prominence',
                   0, 1000, dragging=True, valfmt='%i',
                   valinit=prominence)
    ax_wid = plt.axes([0.25, 0.01, 0.65, 0.01])
    swid = Slider(ax_wid, 'min peak width (pixels)',
                  0, 100, dragging=True, valfmt='%i',
                  valinit=width)

    def update(val):
        # update contrast sliders
        vmin, vmax = smin.val, smax.val
        ct_plot.set_clim(vmin=vmin, vmax=vmax)
        ax1.set_ylim(vmin, vmax)

        # update fft threshold slider
        freq_thresh = sfft.val
        zi_filter = lamina_fft_filter(zi, high_freq_thresh=freq_thresh, plot=False)
        prominence = sprom.val
        width = swid.val
        peaks, valleys = lamina_peaks(zi_filter, prominence=prominence,
                                      width=width)
        fft_plot.set_data(y_cm, zi_filter)
        peaks_plot.set_data(y_cm[peaks], zi_filter[peaks])
        valleys_plot.set_data(y_cm[valleys], zi_filter[valleys])

        fig.canvas.draw_idle()

    smin.on_changed(update)
    smax.on_changed(update)
    sfft.on_changed(update)
    sprom.on_changed(update)
    swid.on_changed(update)

    sliders = [smin, smax, sfft, sprom, swid]
    # # ## Horizontal line across both plots (multicursor)
    multi = MultiCursor(fig.canvas, (ax, ax1), color='r', lw=0.5,
                        horizOn=False, vertOn=True)

    return fig, multi, sliders

################################################################################


def pick_laminae(ct_data, ct_xml, layout='horizontal'):
    """
    interactive plot. If you click on a plot, it will write the depth at that
    pixel to a list.
    """
    # Load the ct image file and xml data
    im = ct_data
    ct_xml = ct_xml

    # Get screen size for figure
    root = tkinter.Tk()
    pix2in = root.winfo_fpixels('1i')
    screen_width = root.winfo_screenwidth()/pix2in
    screen_height = root.winfo_screenheight()/pix2in
    image_h2w = round(ct_xml['physical-height']/ct_xml['physical-width'])
    pix2cm = ct_xml['pixels-per-CM']
    if layout == 'horizontal':
        fig, (ax1, ax) = plt.subplots(nrows=2,
                                      figsize=(screen_width*0.8, screen_height*0.8),
                                      sharex=True)
        plt.subplots_adjust(bottom=0.25)
    else:
        fig, (ax1, ax) = plt.subplots(ncols=2,
                                      figsize=(screen_width*0.8, screen_height/image_h2w*10),
                                      sharey=True)

    # Extract transect
    zi, y_cm, profile_points = extract_profile(ct_data, ct_xml)
    x0_cm, x1_cm = profile_points[0], profile_points[1]
    y0_cm, y1_cm = profile_points[2], profile_points[3]

    # Calculate vmin and vmax based on the transect values
    vmin0 = np.min(zi)
    vmax0 = np.max(zi)
    vmin = np.min(ct_data)
    vmax = np.max(ct_data)

    # Plot images
    aspect = 'equal'
    if layout == 'horizontal':
        ct_plot = ax.imshow(np.rot90(im, k=1, axes=(1, 0)), aspect=aspect,
                            extent=(ct_xml['physical-top']/100 +
                                    ct_xml['physical-height'],
                                    ct_xml['physical-top']/100,
                                    0,
                                    ct_xml['physical-width']),
                            vmin=vmin0, vmax=vmax0)
        transect = ax.plot([y0_cm, y1_cm], [x0_cm, x1_cm], 'ro-')
        ax.set_xlim(y0_cm, y1_cm)
        ax.set_ylabel('Width [cm]')
        ax.set_xlabel('Depth [cm]')
    else:
        ct_plot = ax.imshow(im, aspect=aspect,
                            extent=(0, ct_xml['physical-width'],
                                    ct_xml['physical-top']/100+ct_xml['physical-height'],
                                    ct_xml['physical-top']/100),
                            vmin=vmin0, vmax=vmax0)
        transect = ax.plot([x0_cm, x1_cm], [y0_cm, y1_cm], 'ro-')
        ax.set_ylim(y1_cm, y0_cm)
        ax.set_xlabel('Width [cm]')

    ax.set_title(ct_xml['coreID'])

    contact_depths = []

    def onclick(event):
        x, y = event.xdata, event.ydata
        if layout == 'horizontal':
            if event.inaxes == ax:
                ax.plot(x, y, 'rx')
                contact_depths.append(x)
                ax.vlines(x, ymin=0, ymax=ct_xml['physical-width'], ls='--', lw=0.5)
                ax1.vlines(x, ymin=vmin0, ymax=vmax0, ls='--', lw=0.5)
            if event.inaxes == ax1:
                ax1.plot(x, y, 'rx')
                contact_depths.append(x)
                ax.vlines(x, ymin=0, ymax=ct_xml['physical-width'], ls='--', lw=0.5)
                ax1.vlines(x, ymin=vmin0, ymax=vmax0, ls='--', lw=0.5)
        fig.canvas.draw()
    pick = fig.canvas.mpl_connect('button_press_event', onclick)

    # Plot transect values
    if layout == 'horizontal':
        ax1.plot(y_cm, zi, 'k-', linewidth=1.0)
        ax1.set_xlabel('Depth [cm]')
        ax1.set_ylabel('Intensity')
        ax1.set_xlim(y0_cm, y1_cm)
        ax1.tick_params(labelbottom=True)
    else:
        ax1.plot(zi, y_cm, 'k-', linewidth=1.0)
        ax1.set_ylabel('Depth [cm]')
        ax1.set_xlabel('Intensity')
        ax1.set_ylim(y1_cm, y0_cm)

    # Vmin and vmax slider
    ax_max = plt.axes([0.25, 0.21, 0.65, 0.01])
    ax_min = plt.axes([0.25, 0.16, 0.65, 0.01])

    smin = Slider(ax_min, 'min intensity', vmin, vmax, valinit=vmin0,
                  dragging=True, valfmt='%i')
    smax = Slider(ax_max, 'max intensity', vmin, vmax, valinit=vmax0,
                  dragging=True, valfmt='%i')

    def update(val):
        # update contrast sliders
        vmin0, vmax0 = smin.val, smax.val
        ct_plot.set_clim(vmin=vmin0, vmax=vmax0)
        ax1.set_ylim(vmin0, vmax0)
        fig.canvas.draw_idle()

    smin.on_changed(update)
    smax.on_changed(update)

    sliders = [smin, smax]

    # # ## Horizontal line across both plots (multicursor)
    multi = MultiCursor(fig.canvas, (ax, ax1), color='r', lw=0.5,
                        horizOn=False, vertOn=True)

    return fig, multi, sliders, contact_depths

################################################################################


def picks_to_excel(contact_depths, filename=''):
    """
    save list of picked depths to a csv file
    """
    depths = sorted(contact_depths)  # sort the depths
    depths = np.array(depths)

    if not filename:
        root = tkinter.Tk()
        root.wm_withdraw()
        filename = filedialog.asksaveasfilename()
        root.destroy()
        if not filename:
            sys.exit()
    fname = filename

    np.savetxt(fname, depths.T, delimiter=',', fmt='%f')
    print('Saved picks as %s' % fname)

################################################################################


def extract_average_ortho(ct, ct_xml, window_width=1):
    '''extract downcore ct intensity from the center of the core of orthogonal view.
        will average across the specified distance (centered on core midpoint).
        best to perform this on a cropped core image to avoid blank space.
    '''
    cm2pix = ct_xml['pixels-per-CM']
    width = ct_xml['physical-width']
    midpoint = width/2
    left = midpoint - window_width/2
    right = midpoint + window_width/2
    left_pix, right_pix = round(cm2pix*left), round(cm2pix*right)
    avg = [np.average(ct[i, left_pix:right_pix]) for i in range(np.shape(ct)[0])]
    return avg
