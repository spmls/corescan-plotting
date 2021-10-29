corescan_plotting
=======

`corescan_plotting` contains tools for analyzing and plotting outputs from various core scanning instruments located at the USGS Pacific Coastal and Marine Science Center.

Setup
-----
For USGS python users, I recommend using [Miniconda](https://docs.conda.io/en/latest/miniconda.html). 

Install the appropriate Python 3.x Miniconda installer. If on Windows, choose *Just Me* and for location, something like `C:\Users\USERNAME\Miniconda3`, but this is up to you. If this is your first tie working with conda environments, IOOS has a [useful guide](https://ioos.github.io/ioos_code_lab/content/ioos_installation_conda.html) on setting up conda on both macos and windows and is perhaps a good place to start before creating a new conda environment for the corescan_plotting package.

Create a conda environment (see: [conda.io](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands)) with python=3.7

Once created, activate the environment and install corescan_plotting with:
`conda install -c spmls corescan_plotting`

opencv is a tricky dependency to install, you may have to reinstall it with: `pip install opencv-python==4.2.0.32` if 
corescan_plotting fails to import in python.

Usage
-----
See examples.ipynb

```
from corescan_plotting import ct
import matplotlib.pyplot as plt

im, xml = ct.ct_in() # this will open up a window to load in a Geotek orthogonal view image and associated xml file
ct.ct_plot(im, xml) # this will plot the loaded ct file
plt.show()
```