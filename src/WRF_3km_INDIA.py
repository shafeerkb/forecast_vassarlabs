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


def download(Ini_date):
    start_time = time.time()
    Ini_date='2020-04-06 00:00:00'
    date_time=pd.date_range(start=Ini_date, periods=75,freq='H')
    os.mkdir(date_time[0].strftime("%Y%m%d%H/"))
    os.chdir(date_time[0].strftime("%Y%m%d%H/"))
    
    for tm in date_time[1:]:
        nc_url=date_time[0].strftime("%Y%m%d%H/")+"WRF3km_INDIA-"+tm.strftime("%Y-%m-%d_%H.nc")
        URL="axel -a -n 4 ftp://nwp:nwp@125.21.185.50/nwp-data/WRF_NETCDF/"+nc_url
    #    print(YRL)
        os.system(URL)        
    print("Downloded in :",round(time.time() - start_time,2),"Sec")


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
#        print(time_,"    ",input_nc_variable.name)  
        fig = plt.figure(frameon=False,dpi=dpi)
        ax = plt.axes(projection=proj)
        ax.set_extent(img_extent, ccrs.PlateCarree())
        ax.outline_patch.set_visible(False)
        ax.background_patch.set_visible(False)
        ax.imshow(input_nc_variable[t,:,:].values,interpolation=interp_methord, origin='lower',alpha=.8, cmap=cmap, norm=norm, extent=img_extent, transform=ccrs.PlateCarree())
        #ax.contour(grid1_lon,grid1_lat,grid1,alpha=1, levels=bounds, colors='k', transform=ccrs.PlateCarree(),linewidths=0.05)    
        #ax.contour(cntr,grid1_lon,grid1_lat,grid1, colors='black')
        #ax.clabel(cntr, inline=True, fontsize=8)    
        
        fig.savefig(input_nc_variable.name+"/"+time_+'.png',bbox_inches='tight',transparent=True,pad_inches = 0.0)
        plt.close()
        convert="convert "+input_nc_variable.name+"/"+time_+'.png '+ input_nc_variable.name+"/"+time_+'_4_bit.png'
        os.system(convert)
        os.remove(input_nc_variable.name+"/"+time_+'.png')
        t=t+1
    print(input_nc_variable.name,"Heatmap  :",round(time.time() - start_time,2),"Sec")

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
     
def FFWI():
    pass

def WRF3km_fn(in_file_path,out_file_path,slice_hr,ref_date,data_dir):
    os.mkdir(out_file_path)
    
    start_time = time.time()
    
    WRF3km=xr.open_mfdataset(in_file_path+"WRF3km_INDIA*.nc",concat_dim='time',combine='nested')
#    sun_=WRF3km.SWNETB*-1
    
    urls =glob.glob(in_file_path+"WRF3km_INDIA*.nc")
    urls.sort()
    lead_time_hr =list(map (lambda x: x[-16:-3], urls))
    
    hr_since=(pd.to_datetime(lead_time_hr, format="%Y-%m-%d_%H")-pd.to_datetime(ref_date, format="%Y-%m-%d %H:%M:%S"))// np.timedelta64(1, 'h')
    WRF3km['time']=hr_since-1
    
    HI=mpcalc.heat_index(WRF3km.T2m.values *units.degC,WRF3km.RH2m.values/100 * units.dimensionless,mask_undefined=False).to(units.degC)
    
    Heatindex_hr=WRF3km.T2m.copy()*0+HI.magnitude
    
    #t_bins = list(range(0,WRF3km.time.size,6))
    t_bins = list(range(0,hr_since[-1],6))
    
    W_SPD=np.sqrt(WRF3km.U10**2+WRF3km.V10**2)
    
    group=WRF3km.APCP.groupby_bins('time', bins=t_bins,include_lowest=True,right=False).groups
    
    APCP=WRF3km.APCP.groupby_bins('time', bins=t_bins,include_lowest=True,right=False).sum(dim=('time'))  #sum of rainfall
    T2m =WRF3km.T2m.groupby_bins('time', bins=t_bins,include_lowest=True,right=False).max(dim=('time'))   #max of temperature
    RH2m =WRF3km.RH2m.groupby_bins('time', bins=t_bins,include_lowest=True,right=False).max(dim=('time')) #max of Relative Humidity
    #SWNETB =WRF3km.SWNETB.groupby_bins('time', bins=t_bins,include_lowest=True,right=False).mean(dim=('time'))
    #LWNETB =WRF3km.LWNETB.groupby_bins('time', bins=t_bins,include_lowest=True,right=False).mean(dim=('time'))
    Heat_Index=Heatindex_hr.groupby_bins('time', bins=t_bins,include_lowest=True,right=False).max(dim=('time'))   #max of Heat Index
    
    #U10 =WRF3km.U10.groupby_bins('time', bins=t_bins,include_lowest=True,right=False).max(dim=('time'))
    #V10 =WRF3km.V10.groupby_bins('time', bins=t_bins,include_lowest=True,right=False).max(dim=('time'))
    
    W_speed =W_SPD.groupby_bins('time', bins=t_bins,include_lowest=True,right=False).max(dim=('time'))
    W_speed_index =W_SPD.groupby_bins('time', bins=t_bins,include_lowest=True,right=False).map(lambda x: np.argmax(x,axis=0))
    
    print(time.time() - start_time)
    U10=W_speed.values.copy()*0
    V10=W_speed.values.copy()*0
    for tm in range(U10.shape[0]):
        print(tm,group[list(group)[tm]][0],tm*6)
    #    a=W_speed_index[tm,:,:].values+(tm*6)
        a=W_speed_index[tm,:,:].values+ group[list(group)[tm]][0]
        m,n = a.shape
        I,J = np.ogrid[:m,:n]
        U10[tm,:,:] = WRF3km.U10.values[a, I, J]
        V10[tm,:,:] = WRF3km.V10.values[a, I, J]
        
        
    
    ds = xr.Dataset({'Total_precipitation': (['time','lat','lon' ],  APCP.values),
                      'Temperature_2m': (['time','lat','lon'], T2m.values),
                      'Relative_Humidity_2m': (['time','lat','lon'], RH2m.values),
                      'U10': (['time','lat','lon'], U10),
                      'V10': (['time','lat','lon'], V10),
                      'wind_speed': (['time','lat','lon'], W_speed.values),
                      'Heat_Index': (['time','lat','lon'], Heat_Index.values)},
                      
                     coords={'lon': (['lon'], WRF3km.lon.values),
                             'lat': (['lat'], WRF3km.lat.values),
                             'time': pd.date_range(ref_date, periods=len(t_bins),freq='6H')[1:],
                             'reference_time': pd.Timestamp(ref_date)})
    
    os.chdir(out_file_path)
    
    ds.Total_precipitation.to_netcdf("Total_precipitationto.nc")
    ds.Temperature_2m.to_netcdf("Temperature_2m.nc")
    ds.Relative_Humidity_2m.to_netcdf("Relative_Humidity_2m.nc")
    ds.U10.to_netcdf("U10.nc")
    ds.V10.to_netcdf("V10.nc")
    ds.wind_speed.to_netcdf("wind_speed.nc")
    ds.Heat_Index.to_netcdf("Heat_Index.nc")
    
    
    mask_table=pd.read_csv("/home/vassar/Documents/Rahul/Block&&SubBasin/Block_UUID (1)/Block_uuid.csv")
    mask_table=mask_table[['OBJECTID_1',"State_Name",'district', 'block',  'UUID']]
    mask_table.columns=['Mask_id',"State_Name",'district', 'block',  'UUID']
    
    mask = xr.open_dataset('/home/vassar/Documents/Rahul/Block&&SubBasin/Block_UUID (1)/Block.nc')

    nc2table_max(ds.Temperature_2m,mask.Band1,mask_table,data_dir,"bilinear_1284x1284_3032x2923.nc") 
    nc2table_max(ds.Relative_Humidity_2m,mask.Band1,mask_table,data_dir,"bilinear_1284x1284_3032x2923.nc") 
    nc2table(ds.Total_precipitation,mask.Band1,mask_table,data_dir,"bilinear_1284x1284_3032x2923.nc") 
    nc2table_max(ds.wind_speed,mask.Band1,mask_table,data_dir,"bilinear_1284x1284_3032x2923.nc") 
    nc2table_max(ds.Heat_Index,mask.Band1,mask_table,data_dir,"bilinear_1284x1284_3032x2923.nc") 

    print(time.time() - start_time)
    
    lat_new=np.linspace(WRF3km['lat'].values.min(), 
                            WRF3km['lat'].values.max(), num=int((WRF3km['lat'].values.max()-WRF3km['lat'].values.min())/0.25)+1, endpoint=True)
    lon_new=np.linspace(WRF3km['lon'].values.min(), 
                            WRF3km['lon'].values.max(), num=int((WRF3km['lon'].values.max()-WRF3km['lon'].values.min())/0.25)+1, endpoint=True)
    
    
    ds_out = xr.Dataset({'lat': (['lat'], lat_new),
                         'lon': (['lon'], lon_new),
                        }
                       )
    ds_out
    
    
    regridder = xe.Regridder(WRF3km, ds_out, 'bilinear',reuse_weights=True,filename=data_dir/"bilinear_1284x1284_141x141.nc");
    regridder  # print basic regridder information.
    
    U10p25 = regridder(ds.U10)
    V10p25 = regridder(ds.V10)
    
    U10p25.to_netcdf("WRF0p25_U10m.nc")
    V10p25.to_netcdf("WRF0p25_V10m.nc")
    
    wind=xr.merge([U10p25,V10p25])
    
    #wind.isel(time=range(2,72,6)).to_netcdf("WRF0p25WIND10m.nc")h
    wind.to_netcdf("WRF0p25WIND10m.nc")
    
    print(time.time() - start_time)
   
    
#    data_dir_raw=Path("/home/vassar/Documents/forecast_vassarlabs/data")


    brk_rh= pd.read_csv(data_dir/"brk_rh.csv")["x"]
    col_rh= pd.read_csv(data_dir/"col_rh.csv")["x"]
    
    
    brk_tem= pd.read_csv(data_dir/"brk_tem.csv")["x"]
    col_tem= pd.read_csv(data_dir/"col_tem.csv")["x"]
    
    brk_ws= pd.read_csv(data_dir/"brk_ws.csv")["x"]
    col_ws= pd.read_csv(data_dir/"col_ws.csv")["x"]
    
    brk_rf=[0,.1,2.5,15.6,64.5,115.6,204.5,8000]
    col_rf=["#FFFFFF","#C3FDCA","#01FF04","#048500","#FDC0CB","#FC0300","#610301"]
    

   

    Process(target=heatmap, args=(ds.Temperature_2m,col_tem,brk_tem)).start()
    Process(target=heatmap, args=(ds.Heat_Index,col_tem,brk_tem)).start()
    Process(target=heatmap, args=(ds.Relative_Humidity_2m,col_rh,brk_rh) ).start()
    Process(target=heatmap, args=(ds.wind_speed,col_ws,brk_ws) ).start()   
    Process(target=heatmap, args=(ds.Total_precipitation,col_rf,brk_rf) ).start()  
    
        
        
    
