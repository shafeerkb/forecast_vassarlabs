#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 10:05:40 2020

@author: vassar
"""

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import pandas as pd
import xesmf as xe
import glob
import time
import os

import sys

os.chdir(Path().absolute())
print(Path().absolute())

sys.path.append("./src")
from src.WRF_3km_INDIA import WRF3km_fn,heatmap

start_time = time.time()

data_dir=Path("./data").absolute()
brk_rh= pd.read_csv(data_dir/"brk_rh.csv")["x"]

in_file_path="/home/vassar/Documents/forcastdata/WRF/"
out_file_path="/home/vassar/2020022600"
slice_hr=6
ref_date='2020-02-26 00:00:00'



WRF3km_fn(in_file_path,out_file_path,slice_hr,ref_date,data_dir)   

    
print("Total time=", time.time() - start_time,"Sec")


