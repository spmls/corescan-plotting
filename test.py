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
import scipy
import matplotlib as matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, MultiCursor
import exifread
import cv2
%matplotlib qt5

#%%
################################################################################
def count_laminae(ct_data,ct_xml,layout='horizontal',transect=None):
    """
    automate counting of laminae
    if transect is not specified, it will extract a line from the center of the
    full length of the core. If specified, transect should be in the form:
    [x0,x1,y0,y1] (x0 and x1 will be equal for a straight vertical line)
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
    pix2cm = ct_xml['pixels-per-CM']
    if layout == 'horizontal':
        fig,(ax1,ax) = plt.subplots(nrows=2,
            figsize=(screen_width*0.8,screen_height*0.8),
            sharex=True)
        plt.subplots_adjust(bottom=0.2)
    else:
        fig,(ax1,ax) = plt.subplots(ncols=2,
            figsize=(screen_width*0.8,screen_height/image_h2w*10),
            sharey=True)
    ## Calculate vmin and vmax automatically
    vmin0 = np.min(ct_data)
    vmax0 = np.max(ct_data)

    ## Extract transect
    if transect == None:
        x0,y0 = ct_xml['pixel-width']/2, 0
        x1,y1 = ct_xml['pixel-width']/2,ct_xml['scan-lines']
        length = ct_xml['scan-lines']
        x = np.linspace(x0,x1,length)
        y = np.linspace(y0,y1-1,length)
        x0_cm = x0/pix2cm
        x1_cm = x0_cm
        y0_cm, y1_cm = y0/pix2cm, y1/pix2cm
    else:
        x0_cm,y0_cm = transect[0],transect[2]
        x1_cm,y1_cm = transect[1],transect[3]
        x0,x1 = x0_cm*pix2cm, x1_cm*pix2cm
        y0,y1 = y0_cm*pix2cm, y1_cm*pix2cm
        length = ct_xml['scan-lines']
        x = np.linspace(x0,x1,length)
        y = np.linspace(y0,y1-1,length)
    zi = ct_data[y.astype(np.int),x.astype(np.int)]
    y_cm = y/ct_xml['pixels-per-CM']

    # filter out high frequencies with fourier transform
    zi_filter=lamina_fft_filter(zi,high_freq_thresh=0.03,plot=False)

    ## Plot images
    aspect = 'equal'
    if layout == 'horizontal':
        ct_plot = ax.imshow(np.rot90(im,k=1,axes=(1,0)), aspect=aspect,
                             extent=(ct_xml['physical-top']/100+\
                                     ct_xml['physical-height'],
                                     ct_xml['physical-top']/100,
                                     0,
                                     ct_xml['physical-width']),
                            vmin=vmin0,vmax=vmax0)
        transect = ax.plot([y0_cm,y1_cm],[x0_cm,x1_cm],'ro-')
        ax.set_xlim(y0_cm,y1_cm)
        ax.set_ylabel('Width [cm]')
        ax.set_xlabel('Depth [cm]')
    else:
        ct_plot = ax.imshow(im, aspect=aspect,
                             extent=(0,ct_xml['physical-width'],\
                             ct_xml['physical-top']/100+ct_xml['physical-height'],\
                             ct_xml['physical-top']/100),
                             vmin=vmin0,vmax=vmax0)
        transect = ax.plot([x0_cm,x1_cm],[y0_cm,y1_cm],'ro-')
        ax.set_ylim(y1_cm,y0_cm)
        ax.set_xlabel('Width [cm]')

    ax.set_title(ct_xml['coreID'])

    ## Plot transect values
    if layout == 'horizontal':
        ax1.plot(y_cm,zi,'k-',linewidth=0.5)
        ax1.set_xlabel('Depth [cm]')
        ax1.set_ylabel('Intensity')
    else:
        ax1.plot(zi,y_cm,'k-',linewidth=0.5)
        ax1.set_ylabel('Depth [cm]')
        ax1.set_xlabel('Intensity')

    ## Peak finding
    peaks, peak_prop = scipy.signal.find_peaks(zi_filter,
                                               prominence=(1000,None))
    valleys, valley_pro = scipy.signal.find_peaks(-1.*zi_filter,
                                                  prominence=(1000,None))
    ax1.plot(y_cm[peaks],zi_filter[peaks],'rx')
    ax1.plot(y_cm[valleys],zi_filter[valleys],'bx')
    ax1.plot(y_cm,zi_filter,'r--',linewidth=0.3)


    ## Vmin and vmax slider
    ax_max = plt.axes([0.25, 0.05, 0.65, 0.01])
    ax_min = plt.axes([0.25, 0.1, 0.65, 0.01])

    smin = Slider(ax_min,'min intensity',vmin0,vmax0,valinit=vmin0,
                  dragging=True,valfmt='%i')
    smax = Slider(ax_max,'max intensity',vmin0,vmax0,valinit=vmax0,
                  dragging=True,valfmt='%i')

    def update(val):
        vmin,vmax = smin.val, smax.val
        ct_plot.set_clim(vmin=vmin,vmax=vmax)
        ax1.set_ylim(vmin,vmax)
        fig.canvas.draw_idle()

    smin.on_changed(update)
    smax.on_changed(update)

    # # ## Horizontal line across both plots (multicursor)
    multi = MultiCursor(fig.canvas, (ax,ax1), color='r', lw=0.5,
                        horizOn=False, vertOn=True)

    return fig,multi

################################################################################
def lamina_fft_filter(signal,high_freq_thresh=0.03,plot=False):

    # Calculate FFT of signal and the power
    sig_fft = np.fft.fft(signal)
    power = np.abs(sig_fft)
    sample_freq = np.fft.fftfreq(signal.size,d=1)

    # Find peak Frequency
    pos_mask = np.where(sample_freq>0)
    freqs = sample_freq[pos_mask]
    peak_freq = freqs[power[pos_mask].argmax()]

    # Remove high frequencies
    high_freq_fft = sig_fft.copy()
    high_freq_fft[np.abs(sample_freq) > high_freq_thresh] = 0
    filtered_sig = np.fft.ifft(high_freq_fft)

    if plot == True:
        fig = plt.figure(figsize=(11,8.5))
        ax1 = fig.add_subplot(2,1,1)
        ax1.plot(sample_freq, power)
        ax1.set_xlabel('Frequency [Hz]')
        ax1.set_ylabel('Power')

        ax2 = fig.add_subplot(2,1,2)
        ax2.plot(signal)
        ax2.plot(filtered_sig,'r-')
        ax2.set_title('Peak Frequency')

    return filtered_sig

# %% Load and test
import plot_ct_tools as pct
ct, xml = pct.ct_in('I:\Floras\Floras_RXCT\selected_orthogonal_views\VC09'\
                    '\FLO17-VC09-000-150cm\SliceViews_XZView2.TIF')
ct, xml = pct.ct_crop_rotate(ct,xml,70,plot=False)

# %%
if fig:
    plt.close()
fig,multi = count_laminae(ct,xml,transect=[3.1,3.1,60,100])




# %%
fig = plt.figure(figsize=(11,8.5))

# Calculate FFT of signal and the power and plot
sig_fft = np.fft.fft(signal)
power = np.abs(sig_fft)
sample_freq = np.fft.fftfreq(signal.size,d=1)
ax1 = fig.add_subplot(2,1,1)
ax1.plot(sample_freq, power)
ax1.set_xlabel('Frequency [Hz]')
ax1.set_ylabel('Power')

# Find peak Frequency
pos_mask = np.where(sample_freq>0)
freqs = sample_freq[pos_mask]
peak_freq = freqs[power[pos_mask].argmax()]

# Remove high frequencies
high_freq_fft = sig_fft.copy()
high_freq_fft[np.abs(sample_freq) > 0.03] = 0
filtered_sig = np.fft.ifft(high_freq_fft)
ax2 = fig.add_subplot(2,1,2)
ax2.plot(signal)
ax2.plot(filtered_sig,'r-')
ax2.set_title('Peak Frequency')



# %%

fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
ax1.imshow(np.rot90(ct,k=1,axes=(1,0)), aspect='equal',
           extent=(xml['physical-top']/100+\
           xml['physical-height'],
           xml['physical-top']/100,0,xml['physical-width']))
ax1.set_xlim(60,100)
ax2.plot(signal)
multi = MultiCursor(fig.canvas, (ax1, ax2), color='r', lw=1,
                    horizOn=False, vertOn=True)
plt.show()
