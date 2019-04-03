from ct import plot_ct_tools
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
%matplotlib inline


# %% testing the ct_in function to read in RXCT tiffs and xml
ct_data, ct_xml = plot_ct_tools.ct_in('ct/test_images/CT/FLO18-VC29-420-570cm/FLO18-VC29-420-570cm_XZView.TIF')
ct_xml

#%% testing the ct_plot function
filename = 'ct/test_images/CT/FLO18-VC29-420-570cm/FLO18-VC29-420-570cm_XZView.TIF'
fig = plot_ct_tools.ct_plot(filename,vmin=18000,vmax=36000)
dpi = plot_ct_tools.dpi_ct_plot(fig,ct_xml)
fig.savefig('test1466.png',type='png',dpi=dpi,bbox_inches='tight',transparent=True)
#%% testing the ct_histogram function
plot_ct_tools.ct_histogram(ct_data)
