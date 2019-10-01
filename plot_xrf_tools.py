"""
Created on Wednesday Septebmer 25 17:07 2019

tools to work with XRF data from the Geotek MSCL (Olympus head)

@author: SeanPaul La Selle
"""

import os
import sys
import glob
import tkinter
from tkinter import filedialog
import numpy as np
import csv
import pandas
import matplotlib as matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import warnings
from corescan_plotting import plot_ct_tools, plot_linescan_tools

#%% TESTING
# Load XRF data
###############################################################################
def xrf_in(filename=''):
    """
    read in Geotek XRF data from from .out file
    """
    ## Get filename if not specified in function call
    if not filename:
        filename = filedialog.askopenfilename()
        if not filename:
            sys.exit()
    data = csv_xrf_parser(filename)

    # Determine the directory of the file
    directory = os.path.dirname(filename)
    ## Read other files
    # if not xml_fname:
    #     xml_fname = glob.glob(os.path.splitext(filename)[0]+'*.xml')[0]
    # xml_dic = linescan_xml(xml_fname)
    return data

def csv_xrf_parser(filename):
    """
    parses a Geotek XRF .out file, returns the elements and an array with
    depths, counts, ppm and errors
    """
    with open(filename) as csvfile:
        readcsv = csv.reader(csvfile,delimiter='\t')
        header=[]
        data = []
        for i,row in enumerate(readcsv): # Assume header is 9 rows
            header.append(row)
            if(i>=9):
                break
        elements = header[7][5::2] # Assume elements start on the 7th row
        for row in readcsv: # From here, csv should be data
            data.append([float(i) for i in row])
    for i,r in enumerate(data): # Need to pad rows with empty data
        if len(r) != len(max(data,key=len)):
            r = np.append(r,np.zeros((len(max(data,key=len))-len(r))))
            data[i] = r
    data = np.reshape(data,(np.shape(data)[0],len(max(data,key=len))))

    return data
# %%
filename="/Volumes/tsudisk/Cascadia/Floras Lake/Floras_XRF/VC22-667-817cm_archive/VC22-667-817cm_archive.out"
data = xrf_in(filename)
