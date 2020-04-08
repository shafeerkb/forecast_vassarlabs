#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 22:19:41 2020

@author: vassar
"""

import pandas as pd
import time
import os

start_time = time.time()


Ini_date='2020-04-06 00:00:00'



date_time=pd.date_range(start=Ini_date, periods=75,freq='H')
os.mkdir(date_time[0].strftime("%Y%m%d%H/"))
os.chdir(date_time[0].strftime("%Y%m%d%H/"))

for tm in date_time[1:]:
    nc_url=date_time[0].strftime("%Y%m%d%H/")+"WRF3km_INDIA-"+tm.strftime("%Y-%m-%d_%H.nc")
    YRL="axel -a -n 4 ftp://nwp:nwp@125.21.185.50/nwp-data/WRF_NETCDF/"+nc_url
#    print(YRL)
    os.system(YRL)
    
print(time.time() - start_time)