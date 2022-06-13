from setuptools import setup, find_packages, Extension
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = '0.9'
DESCRIPTION = 'Python tools for plotting scanned core data from the USGS Pacific Coastal and Marine Science Center'
LONG_DESCRIPTION = ''

# Setting up
setup(
    name="corescan_plotting",
    version=VERSION,
    author="SeanPaul La Selle",
    author_email="<slaselle@usgs.gov>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(),
    install_requires=[],
    keywords=['python', 'geology', 'core', 'ct', 'xrf', 'linescan'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)
