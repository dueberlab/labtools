# !/usr/bin/python

import os

### LCMS data from the single quad can be batch exported, but the CDF files are in subdirectories
### and all have the same name. This script will move them to the parent directory and rename the
### CDF files after the name of the directory in which they are contained.


target_directory = "/path/to/directory/"
cdf_filename = "/MSSIG01.CDF"

for dir in os.listdir(target_directory):
    old_name = target_directory + dir + cdf_filename
    new_name = target_directory + dir.replace(".aia", ".CDF")
    os.renames(old_name, new_name)