# ImproveProcessRawMRR - RaProM.py

RaProM is a novel MRR processing methodology, with enhanced spectra processing and Doppler dealiasing, that produces as output data a number of fields which include equivalent reflectivity (Ze), Doppler fall speed and derived parameters such as spectral width, skewness, and kurtosis, plus a simplified precipitation type classification (drizzle, rain, mixed, snow, and hail), and additional variables depending on the precipitation type. **Note: the scripts works for MRR-2.**
An additional file CorrecRawFile is added. If you know that your raw file has written errors, please execute this file before RaProM. You can execute this file with command window or idle.

More information at: Garcia-Benadí et al (2020)
https://doi.org/10.3390/rs12244113

## Versions and dependences

The main script is called RaProM.py and it is avalaible in python 2.7. and 3.8. The following libraries are necessary::

	numpy , version 1.14.5 or later

	miepython, version 1.3.0 or later

	netCDF4, version 1.2.7 or later

The script works with the MRR raw archives.


## How to cite

If you use this script for your publication, please cite as:

Garcia-Benadi, A.; Bech, J.; Gonzalez, S.; Udina, M.; Codina, B.; Georgis, J.-F. Precipitation Type Classification of Micro Rain Radar Data Using an Improved Doppler Spectral Processing Methodology. Remote Sens. 2020, 12, 4113.DOI: 10.3390/rs12244113  

## Outputs
The script produces the following outputs from MRR raw data:<br />
**W:** fall speed with aliasing correction (m s<sup>-1</sup>)<br />
**spectral width:** spectral width of the dealiased velocity distribution (m s<sup>-1</sup>)<br />
**skewness:** skewness of the dealiased velocity distribution<br />
**kurtosis:** kurtosis of the dealiased velocity distribution<br />
**PIA:** Path Integrated Attenuation<br />
**Type:** Hydrometeor type (unknown[20], rain [10], drizzle [5], snow [-10], mixed [-15] and hail [-20])<br />
**LWC:** Liquid water content (g m<sup>-3</sup>)<br />
**RR:** Rain rate (mm h<sup>-1</sup>)<br />
**SR:** Snow rate (mm h<sup>-1</sup>)<br />
**Z:** Reflectivity considering only liquid drops (dBZ)<br />
**Ze:** Equivalent Reflectivity (dBZ)<br />
**N(D):** Drop Size Distribution (log10(m<sup>-3</sup> mm<sup>-1</sup>))<br />
**SNR:** Signal noise relation from signal without deliasing (dB)<br />
**Noise:** Noise from spectra reflectivity (m<sup>-1</sup>)<br />
**N<sub>w</sub>:** Intercept of the gamma distribution normalized to the liquid water content (log10(m<sup>-3</sup> mm<sup>-1</sup>))<br />
**D<sub>m</sub>:** Mean mass-wighted raindrop diameter (mm)<br />
**BB<sub>bottom</sub>:** Bright Band bottom height  (m) (above MRR level)<br />
**BB<sub>top</sub>:** Bright Band top height (m) (above MRR level)<br />
**TyPrecipi:** Rainfall type where the value 5 is convective, 0 is transition and -5 is stratiform"<br />


## How to execute the script
The script can be executed from a command line at the system prompt (see MS-Windows example):<br />
<br />
![commandWindow](https://user-images.githubusercontent.com/35369817/67784656-64703d00-fa6c-11e9-94fa-0e616d703168.JPG)
<br />
at the directory where RaProM_XX.py has been copied, where XX is 27 or 38 in function of your python version:
```
python RaProM_XX.py

```

The script asks the directory where the raw files to be processed are located (it will process all the MRR raw files of the folder selected), for example:
```
c:\mrrdata\test\
```
**NOTE 1: the path must end with \\ in Windows or a / in Linux**<br />
**NOTE 2: Be careful to not have spaces and special characters in your file path**<br />

The script asks for the integration time (in seconds, usually 60)

The script indicates the number of raw files in the folder and starts the process.

The result is stored in a netcdf file with the same name but finished "-processed"

## Do you have any problem with your data?
If so, your RAW files may be corrupted. There is a new script for this called CorrecRawFiles-py_XX.py .
This script analyze every line in the origin RAW file and fixes it. If errors are found, a new file with 
the same name but finished as -corrected will be created.
To execute the script follow the same steps described above.

## Contact
If you have any question, please contact with Albert at albert.garcia@meteo.ub.edu  or   albert.garcia-benadi@upc.edu
