
| CI          | [![GitHub Workflow Status][github-ci-badge]][github-ci-link] [![anaconda latest release data][latest-release-badge][latest-release-link]                  		               |
| :---------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
| **Package** |                                                          [![Conda][conda-badge]][conda-link] [![PyPI][pypi-badge]][pypi-link]                                                          |
| **License** |                                                                         [![License][license-badge]][repo-link]									       |

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
This installs the latest stable version.

Usage
-----
See [examples.ipynb](https://code.usgs.gov/slaselle/corescan_plotting/-/blob/master/examples.ipynb) for more detailed examples.

```
from corescan_plotting import ct
import matplotlib.pyplot as plt

im, xml = ct.ct_in() # this will open up a window to load in a Geotek orthogonal view image and associated xml file
ct.ct_plot(im, xml) # this will plot the loaded ct file
plt.show()
```

[github-ci-badge]: https://img.shields.io/github/actions/workflow/status/xarray-contrib/datatree/main.yaml?branch=main&label=CI&logo=github
[github-ci-link]: https://github.com/spmls/corescan_plotting/actions?query=workflow%3ACI
[latest-release-badge]: https://anaconda.org/conda-forge/qutip/badges/latest_release_date.svg
[latest-rlease-link]: https://img.shields.io/github/release-date/spmls/corescan_plotting
[pypi-badge]: https://img.shields.io/pypi/v/xarray-datatree?logo=pypi
[pypi-link]: https://pypi.org/project/corescan-plotting/
[conda-badge]: https://img.shields.io/conda/vn/conda-forge/xarray-datatree?logo=anaconda
[conda-link]: https://anaconda.org/spmls/corescan_plotting
[license-badge]: https://anaconda.org/conda-forge/qutip/badges/license.svg
[repo-link]: https://github.com/spmls/corescan_plotting
