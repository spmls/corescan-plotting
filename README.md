
| CI          | [![GitHub Workflow Status][github-ci-badge]][github-ci-link] [![Code Coverage Status][codecov-badge]][codecov-link] [![pre-commit.ci status][pre-commit.ci-badge]][pre-commit.ci-link] |
| :---------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
| **Docs**    |                                                                     [![Documentation Status][rtd-badge]][rtd-link]                                                                     |
| **Package** |                                                          [![Conda][conda-badge]][conda-link] [![PyPI][pypi-badge]][pypi-link]                                                          |
| **License** |                                                                         [![License][license-badge]][repo-link]     

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
[github-ci-link]: https://github.com/xarray-contrib/datatree/actions?query=workflow%3ACI
[codecov-badge]: https://img.shields.io/codecov/c/github/xarray-contrib/datatree.svg?logo=codecov
[codecov-link]: https://codecov.io/gh/xarray-contrib/datatree
[rtd-badge]: https://img.shields.io/readthedocs/xarray-datatree/latest.svg
[rtd-link]: https://xarray-datatree.readthedocs.io/en/latest/?badge=latest
[pypi-badge]: https://img.shields.io/pypi/v/xarray-datatree?logo=pypi
[pypi-link]: https://pypi.org/project/xarray-datatree
[conda-badge]: https://img.shields.io/conda/vn/conda-forge/xarray-datatree?logo=anaconda
[conda-link]: https://anaconda.org/conda-forge/xarray-datatree
[license-badge]: https://img.shields.io/github/license/xarray-contrib/datatree
[repo-link]: https://github.com/xarray-contrib/datatree
[pre-commit.ci-badge]: https://results.pre-commit.ci/badge/github/xarray-contrib/datatree/main.svg
[pre-commit.ci-link]: https://results.pre-commit.ci/latest/github/xarray-contrib/datatree/main
