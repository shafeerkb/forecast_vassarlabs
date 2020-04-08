#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 17:56:20 2020

@author: vassar
"""


import xarray as xr
import xesmf as xe
import pandas as pd
import numpy as np


#mask = xr.open_rasterio('/home/vassar/Documents/Rahul/SuBBasin_Intersection /New_subbasin shapefile/subbasin.tif')
mask = xr.open_dataset('/home/vassar/Documents/Rahul/SuBBasin_Intersection /New_subbasin shapefile/subbasin.nc')

#mask = xr.open_dataset('/home/vassar/Documents/Rahul/Block&&SubBasin/Block_UUID (1)/Block.nc')



data=xr.open_dataset("/home/vassar/Documents/forcastdata/WRF/out/Total_precipitationto.nc")
#
#/home/vassar/Documents/forcastdata/GFS_RIMES/2020022600/GFS_Heat_Index_t2.nc

# ds_out = xr.Dataset({'lon': (['lon'], mask.lon),
#                      'lat': (['lat'], mask.lat),
#                     }
#                    )
# ds_outq




regridder = xe.Regridder(data, mask, 'bilinear',reuse_weights=True);
regridder  # print basic regridder information.

Heat_Index = regridder(data.Heat_Index)

mask_table=pd.read_csv("/home/vassar/Documents/forecast_vassarlabs/data/SUBBASIN_att_table.csv")

#mask_table=pd.read_csv("/home/vassar/Documents/Rahul/Block&&SubBasin/Block_UUID (1)/Block_uuid.csv")
import time
start_time = time.time()
mydf=mask_table[['OBJECTID', 'SUB_BASIN', 'UUID']]#   pd.DataFrame(columns=['OBJECTID_1',"dsad"])
mydf = mydf.reindex(mydf.columns.tolist() +  list(pd.to_datetime(data.Heat_Index.time.values).strftime("%Y-%b_%d_%H")), axis=1)  # version > 0.20.0

row=0
for ID in mask_table["OBJECTID"]:
    mydf.iloc[row,3:]=Heat_Index.where(mask.Band1 == ID).mean(dim=['lat','lon']).values
    row=row+1
    print(time.time() - start_time)
    print(ID)
    
mydf.to_csv("sample_table.csv")
