OzFluxQC
========

Welcome to the repository for the OzFluxQC code.

OzFluxQC is a suite of Python scripts integrated into a single GUI that is designed to simplify and standardise the quality control, post-processing, gap filling and partitioning of data from flux towers.  OzFluxQC has been developed by the OzFlux (http://ozflux.org.au) community in Australia and is used by the community as the operational tool for processing data from the OzFlux network of flux towers.  OzFluxQC is not limited to Australia and can be used for flux tower data collected anywhere in the world.  Using OzFluxQC does not require any knowledge of Python (though we would always recommend people learn Python anyway!), all aspects of the processing can be controlled via the simple GUI and by editing text files.  OzFluxQC can read data from Excel workbooks and CSV files and uses netCDF files (http://www.unidata.ucar.edu/software/netcdf/) for storing intermediate and final output data.

The following documentation gives basic information on how to install and use OzFluxQC.

#Installation
There are 3 steps to installing OzFluxQC:
* 1. Install Python.
* 2. Install the "git" version control software.
* 3. Clone this repository using the "git" version control software.

##Installing Python
OzFluxQC is written for Python V2.7 and uses a number of standard and 3rd party Python modules.

OzFlux uses and recommends the Anaconda (https://www.continuum.io/) Python V2.7 distribution.  This Python distribution comes with all of the modules used by OzFluxQC and all except 1 are installed by default.  Adding the 1 required module that is not installed by default is very easy, thanks to the conda package manager, and is explained below.  Using the Anaconda distribution is not essential, just very convenient, and it is possible to use any Python V2.7 environment provided the required modules are installed.  There is a list of the required modules in the /docs folder of this repository.

To install the Anaconda Python V2.7 distribution, follow these steps:
* 1. Download the Anaconda Python V2.7 installer for your operating system from https://www.continuum.io/downloads.
* 2. Follow the instructions on the Anaconda web page to install the Anaconda Python V2.7 distribution.
* 3. Accept all the defaults during the installation, including having Anaconda append the path to this Python installation to your system PATH environment variable.
* 4. The default installation provides everything needed to run OzFluxQC with 1 exception, the module required to read and write netCDF files.  This can be installed as follows:
  1. Open a command line window or terminal session.
  2. At the command prompt, type "conda install netcdf4" and follow the instructions.  Accept all of the defaults and the netCDF module will be installed.
