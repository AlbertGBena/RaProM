# ImproveProcessRawMRR - RaProM.py

RaProM is a MRR processing methodology, with enhanced spectra processing and Doppler dealiasing, that produces as output data a number of fields which include equivalent reflectivity (Ze), Doppler fall speed and derived parameters such as spectral width, skewness, and kurtosis, plus a simplified precipitation hydrometeor type classification (drizzle, rain, mixed, snow, and hail), and additional variables depending on the precipitation hydrometeor type. MRR stands for Micro Rain Radar, a Doppler vertically pointing radar manufactured by Metek GmbH. A description of the RaProM processing and examples is available at: <br/>
Garcia-Benadi A, Bech J, Gonzalez S, Udina M, Codina B, Georgis J-F. Precipitation Type Classification of Micro Rain Radar Data Using an Improved Doppler Spectral Processing Methodology. Remote Sens. 2020, 12, 4113. https://doi.org/10.3390/rs12244113<br/><br/>
**Note1: the scripts are designed to work with MRR-2.**<br/><br/>
**Note2: a different version similar to RaProM, called RaProM-Pro (https://github.com/AlbertGBena/RaProM-Pro), is available for MRR-Pro files. More information on RaProM-Pro is available at: <br/>
Garcia-Benadí A, Bech J, Gonzalez S, Udina M, Codina B. A New Methodology to Characterise the Radar Bright Band Using Doppler Spectral Moments from Vertically Pointing Radar Observations. Remote Sens. 2021, 13, 4323. https://doi.org/10.3390/rs13214323**<br/><br/>

An additional Python script CorrecRawFile is also included. If you know that your raw file has writing  errors (temporal inconsistencies), please execute this file before RaProM. You can execute this file from the command window or idle.

More information at: Garcia-Benadí et al (2020)
https://doi.org/10.3390/rs12244113

## Versions and dependences

The main script is called RaProM.py and it is available in python 2.7., 3.8. ,and 3.11. The following libraries are necessary:

For 2.7 and 3.8. python version

	numpy , version 1.14.5 or later until 1.19.

	miepython, version 1.3.0 or later (matplotlib is necessary for this library works)

	netCDF4, version 1.2.7 or later(cftime is necessary for this library works)

For 3.11 pyhton version

	numpy , version 1.21.6

	miepython, version 2.2.1 (matplotlib is necessary for this library works)

	netCDF4, version 1.7.2 or later(cftime is necessary for this library works)


The script works with the MRR raw archives.

The libraries can be installed with pip, using these sentences:

	pip install numpy
	pip install miepython
	pip install netCDF4
	pip install matplotlib
	pip install cftime

If you have already installed one of the libraries but need to change the version, you can use this syntaxis:

	pip install numpy~=1.21.6

## How to cite

If you use this script for your publication, please cite as:

Garcia-Benadi A, Bech J, Gonzalez S, Udina M, Codina B, Georgis J-F. Precipitation Type Classification of Micro Rain Radar Data Using an Improved Doppler Spectral Processing Methodology. Remote Sens. 2020, 12, 4113.DOI: 10.3390/rs12244113  

## Outputs
The script produces the following outputs from MRR raw data:<br />
**W:** Fall speed with aliasing correction (m s<sup>-1</sup>)<br />
**spectral width:** Spectral width of the dealiased velocity distribution (m s<sup>-1</sup>)<br />
**skewness:** Skewness of the dealiased velocity distribution<br />
**kurtosis:** Kurtosis of the dealiased velocity distribution<br />
**PIA:** Path Integrated Attenuation calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**PIA_all:** Path Integrated Attenuation calculated assuming all hydrometeors are in liquid phase regardless of hydrometeor type classification<br />
**Type:** Predominant hydrometeor type numerical value where possible values are: -20 (hail), -15 (mixed), -10 (snow), 0 (mixed), 5 (drizzle), 10 (rain) and 20 (unknown precipitation)<br />
**LWC:** Liquid Water Content (g m<sup>-3</sup>) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**RR:** Rain Rate (mm h<sup>-1</sup>) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**SR:** Snow Rate (mm h<sup>-1</sup>)<br />
**Z:** Radar reflectivity (dBZ) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**Za:** Attenuated radar reflectivity (dBZ) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**Ze:** Equivalent radar reflectivity (dBZ)<br />
**N(D):** Drop Size Distribution (log10(m<sup>-3</sup> mm<sup>-1</sup>)) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**N(D) in_function_of_time_and_height** Drop Size Distribution (log10(m<sup>-3</sup> mm<sup>-1</sup>)) in function of time and height calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**SNR:** Signal to noise ratio from signal without dealiasing (dB)<br />
**Noise:** Noise from spectra reflectivity  (m<sup>-1</sup>)<br />
**N<sub>w</sub>:** Intercept of the gamma distribution normalized to the Liquid Water Content (log10(m<sup>-3</sup> mm<sup>-1</sup>)) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**D<sub>m</sub>:** Mean mass-weighted raindrop diameter (mm) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**BB<sub>bottom</sub>:** Bright Band bottom height  (m) (above ground level)<br />
**BB<sub>top</sub>:** Bright Band top height (m) (above ground level)<br />
**TyPrecipi:** Precipitation regime numerical value where possible values are: 5 (convective), 0 (transition) and -5 (stratiform) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
<br />
Notice that PIA and PIA_all have 1 height bin more, because the first element is at 0 m height a.g.l. imposed by the manufacturer

## How to execute the script
The script can be executed from a command line at the system prompt (see MS-Windows example):<br />
<br />
![commandWindow](https://user-images.githubusercontent.com/35369817/67784656-64703d00-fa6c-11e9-94fa-0e616d703168.JPG)
<br />
at the directory where RaProM_XX.py has been copied, where XX is 27, 38, or 3-11 in function of your python version:
```
python RaProM_XX.py

```

The script has some additional command line execution options. Please note that their use may imply a substantial increase of the netcdf output file (see below). Command line options are (more than one is possible, in any order):<br /> 
<br /> 
**<i>-hxxx</i>**: this option forces the MRR antenna height to be at xxx meters above sea level which is important if the height was not correctly configured in the original raw data file. xxx can be a float or an integer value.<br />
<br /> 
**<i>-Myyy</i>**: this option modifies the MRR radar constant so it will affect Z, RR and other variables. M is the multiplicative bias calculated by comparing the MRR rainfall (RR_MRR) with a reference rainfall value such as a rain gauge (RR_REF), M=RR_MRR/RR_REF. M can be a float or an integer, typically close to 1.<br />
<br /> 

The syntax of these options are:
```
python RaProM_XX.py -h100.8

```
This example forces the antenna height to be at 100.8 m above sea level.<br />
```
python RaProM_XX.py -M0.78

```
This example assumes a multiplicative bias of 0.87 between MRR2 and reference rainfall.<br />

The script asks the directory where the raw files to be processed are located (it will process all the MRR raw files of the folder selected), for example:
```
C:\mrrdata\test\
```
**NOTE 1: the path must end with \\ in Windows or a / in Linux**<br />
**NOTE 2: Be careful to avoid using spaces and special characters in your file path.**<br />


The script asks for the integration time (in seconds, usually 60)

The script indicates the number of raw files in the folder and starts the process.

The result is stored in a netcdf file with the same name but finished "-processed"

## Do you have any problem with your data?
If so, your RAW files may be corrupted. There is a new script for this called CorrecRawFiles-py_XX.py .
This script analyses every line in the original RAW file and fixes it. If errors are found, a new file with 
the same name but finished as -corrected will be created.
To execute the script follow the same steps described above.

## Contact
If you have any question, please contact with Albert at albert.garcia@meteo.ub.edu  or   albert.garcia-benadi@upc.edu
