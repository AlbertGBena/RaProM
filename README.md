# ImproveProcessRawMrr
Script to process data from raw Mrr to netcdf

The main script is RaProM_v0.py and it is necessary the next libraries:

numpy , version 1.14.5

miepython, version 1.3.0

To execute the script you must to open with window command (window prompt), 

![commandWindow](https://user-images.githubusercontent.com/35369817/67784656-64703d00-fa6c-11e9-94fa-0e616d703168.JPG)

and in the folder where is the py write:

python RaProM_v0.py


The script ask you where is the raw data to process, you must write the correct root: c:\mrrdata\tes\

The script indicates the number of raw files in the folder and start the process.

The result is a netcdf file with the same name but finished "-processed"
