# ImproveProcessRawMrr - RaProM.py

RaProM is a novel MRR processing methodology, with enhanced spectral processing and Doppler dealiasing, that produces as output data a number of fields which include equivalent reflectivity (Ze), Doppler fall speed and derived parameters such as spectral width, skewness, and kurtosis, plus a simplified precipitation type classification (drizzle, rain, snow, and hail), plus additional variables depending on the precipitation type. 

More information at: [aquí va la ruta a l'article]

## Versions and dependences

The main script is RaProM.py avalaible in python 2.7. and 3.8. It is necessary the next libraries:

	*numpy , version up or equal 1.14.5

	*miepython, version up or equal 1.3.0

	*netCDF4, version up or equal 1.2.7

The script works with the MRR raw archives.


## How to cite

If you use this script for your publication, please cite as:
[aquí va la cita oficial de l'article]


## Outputs
The script produces the following outputs from MRR raw data:"\n"
**W:** vertical speed with aliasing correction (m s-1)"\n"
**spectral width:** spectral width with aliasing (m s-1)
**skewness:** skewness of the spectral reflectivity with dealiasing
**kurtosis:** kurtosis of the spectral reflectivity with dealiasing
**PIA:** Path Integrated Attenuation
**Type:** Type from hydrometeor (unknown[20], rain [10], drizzle ]5], snow [-10], mixed [-15] and hail [-20])
**LWC:** Liquid water content (g m-3)
**RR:** Rain rate (mm hr-1)
**SR:** Snow rate (mm hr-1)
**Z:** Reflectivity considering only liquid drops (dBZ)
**Ze:** Equivalent Reflectivity (dBZ)
**Vmov:** Verical movement (+1 downward -1 upward)
**N(D):** Drop Size Distribution (log10(m-3 mm-1))
**SNR:** Signal noise relation from signal without deliasing (dB)
**Noise:** Noise (m-1)
**Nw:** Intercept of the gamma distribution normalized to the liquid water content (log10(mm-1 m-3))
**Dm:** Mean mass-wighted raindrop diameter (mm)
**Fall speed variability:** Estimate of the fall speed variability (m s-1)
**BB_bottom:** height from the Bright Band bottom (m)
**BB_top:** height from the Bright Band top (m)


## How to execute the script
To execute the script you must to open with window command (window prompt),
![commandWindow](https://user-images.githubusercontent.com/35369817/67784656-64703d00-fa6c-11e9-94fa-0e616d703168.JPG)

and in the folder where is the RaProM.py write:
```
python RaProM.py
```

The script ask you where are the raw archieves to process (it will process all the MRR archives of the folder selected). You must write the correct root, for example:
```
c:\mrrdata\test\
```
**NOTE: the path must end with a "\" in Windows or a "/" in Linux**

The script ask you for the number for integration time (usually 60)

The script indicates the number of raw files in the folder and start the process.

The result is a netcdf file with the same name but finished "-processed"


## Contact
If you have any question, please contact with Albert Benadí at albert.garcia@meteo.ub.edu  or   albert.garcia-benadi@upc.edu
