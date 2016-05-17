## IceBridge Arctic sea ice topography scripts

Data processing and plotting used to assess Arctic sea ice topography using high-resolution IceBridge laser (ATM) data, from the following publication:

Petty, A. A., M. C. Tsamados, N. T. Kurtz, S. L. Farrell, T. Newman, J. Harbeck, D. L. Feltham, J. A. Richter-Menge (2016), Characterizing Arctic sea ice topography using high-resolution IceBridge data, The Cryosphere, doi:10.5194/tcd-9-6495-2015.

Note that individual descriptions should be included at the top of each script.

The 'calc' scripts are all used to process IceBridge ATM data into 'topography' datasets used by the 'plot' plotting scripts. 

calc_multi_atm.py is the primary processing script for generating surface feature data from the IceBridge ATM data (and posAV data).

The processed data can be downloaded from ADD ZENEDO LINK WITH DERIVED DATA.
Place this in the Data_output folder. Note that this data is in a binary Python format so data readers (included in IB_functions.py) are needed.

The following raw IceBridge datasets are needed before this script will run:
The IceBridge ATM data: https://nsidc.org/data/docs/daac/icebridge/ilatm1b/. 
The IceBridge DMS imagery: http://nsidc.org/data/iodms1b}. 
The IceBridge IDCSI4 and quick-look sea ice data: http://nsidcorg/data/docs/daac/icebridge/evaluation_products/sea ice-freeboard-snowdepth-thickness-quicklook-index.html and http://nsidc.org/data/idcsi4.html. 

For some of the extra processing/plotting, the following datasets are also required:
The daily OSI-SAF ice type data: http://saf.met.no/p/ice/
The nearest coastline proximity data: http://oceancolor.gsfc.nasa.gov/DOCS/DistFromCoast/.

Note also that Python 2.7 was used for all processing. I have not tested these scripts in Python 3.



