#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 15:42:41 2020

@author: vassar

"""

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
from matplotlib import colors
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from multiprocessing import Process
from pathlib import Path

#def WRF3km_INDIA(in_file_path,out_file_path,slice_hr=6):



def heatmap(input_nc_variable, col, brk):
    start_time = time.time()

    interp_methord='bicubic'#bilinear #
    proj = ccrs.epsg(3857)  #set projection
    dpi=200                 #image quality
    bounds=list(brk)
    cmap =colors.ListedColormap(list(col))    
    norm = colors.BoundaryNorm(bounds, cmap.N)

    img_extent = (min(input_nc_variable.lon.values), max(input_nc_variable.lon.values),
                  min(input_nc_variable.lat.values), max(input_nc_variable.lat.values))  
    os.mkdir(input_nc_variable.name)
    t=0
    for time_ in np.datetime_as_string(input_nc_variable.time.values,unit='h'):
        print(time_,"    ",input_nc_variable.name)  
        fig = plt.figure(frameon=False,dpi=dpi)
        ax = plt.axes(projection=proj)
        ax.set_extent(img_extent, ccrs.PlateCarree())
        ax.outline_patch.set_visible(False)
        ax.background_patch.set_visible(False)
        ax.imshow(input_nc_variable[t,:,:].values,interpolation=interp_methord,alpha=.8, cmap=cmap,origin = 'upper', norm=norm, extent=img_extent, transform=ccrs.PlateCarree())
        #ax.contour(grid1_lon,grid1_lat,grid1,alpha=1, levels=bounds, colors='k', transform=ccrs.PlateCarree(),linewidths=0.05)    
        #ax.contour(cntr,grid1_lon,grid1_lat,grid1, colors='black')
        #ax.clabel(cntr, inline=True, fontsize=8)    
        
        fig.savefig(input_nc_variable.name+"/"+time_+'.png',bbox_inches='tight',transparent=True,pad_inches = 0.0)
        plt.close()
        convert="convert "+input_nc_variable.name+"/"+time_+'.png '+ input_nc_variable.name+"/"+time_+'_4_bit.png'
        os.system(convert)
        os.remove(input_nc_variable.name+"/"+time_+'.png')
        t=t+1
    print(time.time() - start_time)

def nc2table(nc_vname,mask,mask_table,data_dir,weight_regrd):
    regridder = xe.Regridder(nc_vname, mask, 'bilinear',reuse_weights=True,filename=data_dir/weight_regrd)
    regridded_nc = regridder(nc_vname)
    A1=regridded_nc.values
    M1=mask.values
    A=np.append( [M1],A1, axis=0)
    arr = A.transpose(1,2,0)
    df = pd.concat([pd.DataFrame(x) for x in arr], keys=np.arange(max(arr.shape)))
    df.dropna(inplace=True)
    clm = list(pd.to_datetime(nc_vname.time.values).strftime("%Y-%b_%d_%H"))  # version > 0.20.0
    clm.insert(0, "Mask_id")
    df.columns=clm
    df=df.groupby("Mask_id").mean()
    DF=pd.merge(mask_table,df,on="Mask_id",how="inner")  
    DF.to_csv(nc_vname.name+".csv")
    print(nc_vname.name)
     

def nc2table_max(nc_vname,mask,mask_table,data_dir,weight_regrd):
    regridder = xe.Regridder(nc_vname, mask, 'bilinear',reuse_weights=True,filename=data_dir/weight_regrd)
    regridded_nc = regridder(nc_vname)
    A1=regridded_nc.values
    M1=mask.values
    A=np.append( [M1],A1, axis=0)
    arr = A.transpose(1,2,0)
    df = pd.concat([pd.DataFrame(x) for x in arr], keys=np.arange(max(arr.shape)))
    df.dropna(inplace=True)
    clm = list(pd.to_datetime(nc_vname.time.values).strftime("%Y-%b_%d_%H"))  # version > 0.20.0
    clm.insert(0, "Mask_id")
    df.columns=clm
    df=df.groupby("Mask_id").max()
    DF=pd.merge(mask_table,df,on="Mask_id",how="inner")  
    DF.to_csv(nc_vname.name+".csv")
    print(nc_vname.name)
     


def GFS_fn(in_file_path,out_file_path,slice_hr):
    os.mkdir(out_file_path)
    
    start_time = time.time()
    
    GFS_globe=xr.open_mfdataset(in_file_path+"*.grib",engine='cfgrib',compat='override',combine='by_coords')
    GFS_grib=GFS_globe.sel(longitude=slice(65,100), latitude=slice(40,5))
    
    ref_date=str(GFS_grib.time.values)#'2020-02-26 00:00:00'
    GFS_grib.to_netcdf(in_file_path+'GFS_INDIA_'+ref_date[:13]+'.nc')
    
    GFS = xr.open_dataset(in_file_path+'GFS_INDIA_'+ref_date[:13]+'.nc')
    
    t_bins = GFS.step[::2]
    
    
    GFS.tp.groupby_bins('step', bins=t_bins,include_lowest=True,right=False).groups
    tp=GFS.tp.groupby_bins('step', bins=t_bins,include_lowest=True,right=False).sum(dim=('step'))  #GFS.tp[:,1440,0].values
    tmax=GFS.tmax.groupby_bins('step', bins=t_bins,include_lowest=True,right=False).max(dim=('step')) #max
    r2=GFS.r2.groupby_bins('step', bins=t_bins,include_lowest=True,right=False).max(dim=('step'))     #max
    
    W_SPD=np.sqrt(GFS.u10**2+GFS.v10**2)
    W_speed =W_SPD.groupby_bins('step', bins=t_bins,include_lowest=True,right=False).max(dim=('step'))
    W_speed_index =W_SPD.groupby_bins('step', bins=t_bins,include_lowest=True,right=False).map(lambda x: np.argmax(x,axis=0))
    
    HI_tmax=mpcalc.heat_index(GFS.tmax.values *units.degK,GFS.r2.values/100 * units.dimensionless,mask_undefined=False).to(units.degC)
    HI_t2=mpcalc.heat_index(GFS.t2m.values *units.degK,GFS.r2.values/100 * units.dimensionless,mask_undefined=False).to(units.degC)
    Heatindex_tmax_hr=GFS.tmax.copy()*0+HI_tmax.magnitude
    Heatindex_t2_hr=GFS.t2m.copy()*0+HI_t2.magnitude
    Heat_Index_tmax=Heatindex_tmax_hr.groupby_bins('step', bins=t_bins,include_lowest=True,right=False).max(dim=('step'))   #max of Heat Index
    Heat_Index_t2  =Heatindex_t2_hr.groupby_bins('step', bins=t_bins,include_lowest=True,right=False).max(dim=('step'))   #max of Heat Index
    
    print(time.time() - start_time)
    U10=W_speed.values.copy()*0
    V10=U10.copy()*0
    for tm in range( W_speed.shape[0]):
        print(tm)
        a=W_speed_index[tm,:,:].values+(tm*2)
        m,n = a.shape
        I,J = np.ogrid[:m,:n]
        U10[tm,:,:] = GFS.u10.values[a, I, J]
        V10[tm,:,:] = GFS.v10.values[a, I, J]
      
    
    
    ds = xr.Dataset({'Total_precipitation': (['time','lat','lon' ],  tp.values),
                      'T_max': (['time','lat','lon'], tmax.values-273.15),
                      'Relative_Humidity_2m': (['time','lat','lon'], r2.values),
                      'U10': (['time','lat','lon'], U10),
                      'V10': (['time','lat','lon'], V10),
                      'wind_speed': (['time','lat','lon'], W_speed.values),
                      'Heat_Index_tmax': (['time','lat','lon'], Heat_Index_tmax.values),
                      'Heat_Index_t2': (['time','lat','lon'], Heat_Index_t2.values)},
                      
                     coords={'lon': (['lon'], GFS.longitude.values),
                             'lat': (['lat'], GFS.latitude.values),
                             'time': pd.date_range(ref_date, periods= W_speed.shape[0]+1,freq='6H')[1:],
                             'reference_time': pd.Timestamp(ref_date)})
    
    os.chdir(out_file_path)          
    ds.Total_precipitation.to_netcdf("GFS_Total_precipitationto.nc")
    ds.T_max.to_netcdf("GFS_T_max.nc")
    ds.Relative_Humidity_2m.to_netcdf("GFS_Relative_Humidity_2m.nc")
    ds.U10.to_netcdf("GFS_U10.nc")
    ds.V10.to_netcdf("GFS_V10.nc")
    ds.wind_speed.to_netcdf("GFS_wind_speed.nc")
    ds.Heat_Index_tmax.to_netcdf("GFS_Heat_Index_tmax.nc")
    ds.Heat_Index_t2.to_netcdf("GFS_Heat_Index_t2.nc")

    lat_new=np.linspace(ds['lat'].values.min(), 
                        ds['lat'].values.max(), num=int((ds['lat'].values.max()-ds['lat'].values.min())/0.25)+1, endpoint=True)
    lon_new=np.linspace(ds['lon'].values.min(), 
                            ds['lon'].values.max(), num=int((ds['lon'].values.max()-ds['lon'].values.min())/0.25)+1, endpoint=True)
    
    
    ds_out = xr.Dataset({'lat': (['lat'], lat_new),
                         'lon': (['lon'], lon_new),
                        }
                       )

    regridder = xe.Regridder(ds, ds_out, 'bilinear',reuse_weights=True);
    regridder  # print basic regridder information.
    
    U10p25 = regridder(ds.U10)
    V10p25 = regridder(ds.V10)
    
    U10p25.to_netcdf("GFS0p25_U10m.nc")
    V10p25.to_netcdf("GFS0p25_V10m.nc")

        
    print(time.time() - start_time)
   
    
    data_dir_raw=Path("/home/vassar/Documents/forecast_vassarlabs/data")


    brk_rh= pd.read_csv(data_dir_raw/"brk_rh.csv")["x"]
    col_rh= pd.read_csv(data_dir_raw/"col_rh.csv")["x"]
    
    
    brk_tem= pd.read_csv(data_dir_raw/"brk_tem.csv")["x"]
    col_tem= pd.read_csv(data_dir_raw/"col_tem.csv")["x"]
    
    brk_ws= pd.read_csv(data_dir_raw/"brk_ws.csv")["x"]
    col_ws= pd.read_csv(data_dir_raw/"col_ws.csv")["x"]
    
    brk_rf=[0,.1,2.5,15.6,64.5,115.6,204.5,8000]
    col_rf=["#FFFFFF","#C3FDCA","#01FF04","#048500","#FDC0CB","#FC0300","#610301"]
    
    mask_table=pd.read_csv("/home/vassar/Documents/Rahul/Block&&SubBasin/Block_UUID (1)/Block_uuid.csv")
    mask_table=mask_table[['OBJECTID_1',"State_Name",'district', 'block',  'UUID']]
    mask_table.columns=['Mask_id',"State_Name",'district', 'block',  'UUID']
    
    mask = xr.open_dataset('/home/vassar/Documents/Rahul/Block&&SubBasin/Block_UUID (1)/Block.nc')

    nc2table_max(ds.T_max,mask.Band1,mask_table,data_dir_raw,"bilinear_281x281_3032x2923.nc") 
    nc2table_max(ds.Relative_Humidity_2m,mask.Band1,mask_table,data_dir_raw,"bilinear_281x281_3032x2923.nc") 
    nc2table(ds.Total_precipitation,mask.Band1,mask_table,data_dir_raw,"bilinear_281x281_3032x2923.nc") 
    nc2table_max(ds.wind_speed,mask.Band1,mask_table,data_dir_raw,"bilinear_281x281_3032x2923.nc") 
    nc2table_max(ds.Heat_Index_t2,mask.Band1,mask_table,data_dir_raw,"bilinear_281x281_3032x2923.nc") 
   

    Process(target=heatmap, args=(ds.T_max,col_tem,brk_tem)).start()
    Process(target=heatmap, args=(ds.Heat_Index_t2,col_tem,brk_tem)).start()
    Process(target=heatmap, args=(ds.Relative_Humidity_2m,col_rh,brk_rh) ).start()
    Process(target=heatmap, args=(ds.wind_speed,col_ws,brk_ws) ).start()   
    Process(target=heatmap, args=(ds.Total_precipitation,col_rf,brk_rf) ).start()  
    
        
        
    
