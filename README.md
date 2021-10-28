corescan_plotting
=======

`corescan_plotting` contains tools for analyzing and plotting outputs from various core scanning instruments located at the USGS Pacific Coastal and Marine Science Center.

Setup
-----
For USGS python users, I recommend using [Miniconda](https://docs.conda.io/en/latest/miniconda.html). 

Install the appropriate Python 3.x Miniconda installer. If on Windows, choose *Just Me* and for location, something like `C:\Users\USERNAME\Miniconda3`, but this is up to you. 

Create a conda environment (see: [conda.io](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands)) with python=3.7

Once created, activate the environment and install corescan_plotting with:
`conda install -c spmls corescan_plotting`

opencv is a tricky dependency to install, you may have to reinstall it with: `pip install opencv-python==4.2.0.32` if 
corescan_plotting fails to import in python.