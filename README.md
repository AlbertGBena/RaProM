# ImproveProcessRawMrr - RaProM.py

RaProM is a novel MRR processing methodology, with enhanced spectra processing and Doppler dealiasing, that produces as output data a number of fields which include equivalent reflectivity (Ze), Doppler fall speed and derived parameters such as spectral width, skewness, and kurtosis, plus a simplified precipitation type classification (drizzle, rain, mixed, snow, and hail), plus additional variables depending on the precipitation type. 

More information at: [aquí va la ruta a l'article]

## Versions and dependences

The main script is RaProM.py avalaible in python 2.7. and 3.8. It is necessary the next libraries:

	numpy , version up or equal 1.14.5

	miepython, version up or equal 1.3.0

	netCDF4, version up or equal 1.2.7

The script works with the MRR raw archives.


## How to cite

If you use this script for your publication, please cite as:
[aquí va la cita oficial de l'article]


## Outputs
The script produces the following outputs from MRR raw data:<br />
**W:** fall speed with aliasing correction (m s<sup>-1</sup>)<br />
**spectral width:** spectral width of the spectral reflectivity  with dealiasing (m s<sup>-1</sup>)<br />
**skewness:** skewness of the spectral reflectivity with dealiasing<br />
**kurtosis:** kurtosis of the spectral reflectivity with dealiasing<br />
**PIA:** Path Integrated Attenuation<br />
**Type:** Type from hydrometeor (unknown[20], rain [10], drizzle [5], snow [-10], mixed [-15] and hail [-20])<br />
**LWC:** Liquid water content (g m<sup>-3</sup>)<br />
**RR:** Rain rate (mm hr<sup>-1</sup>)<br />
**SR:** Snow rate (mm hr<sup>-1</sup>)<br />
**Z:** Reflectivity considering only liquid drops (dBZ)<br />
**Ze:** Equivalent Reflectivity (dBZ)<br />
**N(D):** Drop Size Distribution (log10(m<sup>-3</sup> mm<sup>-1</sup>))<br />
**SNR:** Signal noise relation from signal without deliasing (dB)<br />
**Noise:** Noise from spectra reflectivity (m<sup>-1</sup>)<br />
**N<sub>w</sub>:** Intercept of the gamma distribution normalized to the liquid water content (log10(m<sup>-3</sup> mm<sup>-1</sup>))<br />
**D<sub>m</sub>:** Mean mass-wighted raindrop diameter (mm)<br />
**Fall speed variability:** Estimate of the fall speed variability (m s<sup>-1</sup>)<br />
**BB<sub>bottom</sub>:** height from the Bright Band bottom (m)<br />
**BB<sub>top</sub>:** height from the Bright Band top (m)<br />


## How to execute the script
To execute the script you must to open with window command (window prompt),<br />
<br />
![commandWindow](https://user-images.githubusercontent.com/35369817/67784656-64703d00-fa6c-11e9-94fa-0e616d703168.JPG)
<br />
and in the folder where is the RaProM_XX.py write, where XX is 27 or 38 in function of your python version:
```
python RaProM_XX.py
```

The script ask you where are the raw archieves to process (it will process all the MRR archives of the folder selected). You must write the correct root, for example:
```
c:\mrrdata\test\
```
**NOTE: the path must end with \\ in Windows or a / in Linux**<br />

The script ask you for the number for integration time (usually 60)

The script indicates the number of raw files in the folder and start the process.

The result is a netcdf file with the same name but finished "-processed"


## Contact
If you have any question, please contact with Albert at albert.garcia@meteo.ub.edu  or   albert.garcia-benadi@upc.edu
