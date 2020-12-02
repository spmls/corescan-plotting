import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="corescan_plotting",
    version="0.0.2",
    author="SeanPaul La Selle",
    author_email="slaselle@usgs.gov",
    description="Tools for analyzing and plotting outputs from various core scanning instruments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://code.usgs.gov/slaselle/corescan_plotting",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
