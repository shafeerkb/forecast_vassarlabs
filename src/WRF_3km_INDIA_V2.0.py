#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 15:45:04 2020

@author: vassar
"""


import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import pandas as pd
import xesmf as xe
import glob
import sys
import time


import os
from matplotlib import colors
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from multiprocessing import Process


def NC_slice(nc_with_missing,slice_hr):
#    slice_hr=6
    start_time = time.time()
    t_bins = nc_with_missing.time[0::slice_hr]
    
    group=nc_with_missing.rf.groupby_bins('time', bins=t_bins,include_lowest=True,right=True).groups
    
    rf   =nc_with_missing.rf.groupby_bins('time', bins=t_bins,include_lowest=True,right=True).sum(dim=('time'))         #sum of rainfall
    temp =nc_with_missing.temp.groupby_bins('time', bins=t_bins,include_lowest=True,right=True).max(dim=('time'))   #max of temperature
    hum  =nc_with_missing.hum.groupby_bins('time', bins=t_bins,include_lowest=True,right=True).max(dim=('time'))      #max of Relative Humidity
    fl   =nc_with_missing.fl.groupby_bins('time', bins=t_bins,include_lowest=True,right=True).max(dim=('time'))   #max of Heat Index

    W_speed =nc_with_missing.wind_mag.groupby_bins('time', bins=t_bins,include_lowest=True,right=True).max(dim=('time'))
    W_speed_=nc_with_missing.wind_mag.copy()
    W_speed_=W_speed_.fillna(-1)    
#    W_speed_index =nc_with_missing.wind_mag.groupby_bins('time', bins=t_bins,include_lowest=True,right=True).map(lambda x: np.argmax(x,axis=0))
    W_speed_index =W_speed_.groupby_bins('time', bins=t_bins,include_lowest=True,right=True).map(lambda x: np.argmax(x,axis=0))
    
    U10=W_speed.values.copy()*0
    V10=W_speed.values.copy()*0
    for tm in range(U10.shape[0]):
        print(tm,group[list(group)[tm]][0],tm*slice_hr)
    #    a=W_speed_index[tm,:,:].values+(tm*6)
        a=W_speed_index[tm,:,:].values+ group[list(group)[tm]][0]
        m,n = a.shape
        I,J = np.ogrid[:m,:n]
        U10[tm,:,:] = nc_with_missing.wind_u.values[a, I, J]
        V10[tm,:,:] = nc_with_missing.wind_v.values[a, I, J]
    
          
    ds = xr.Dataset({'rf': (['time','lat','lon' ],  rf.values),
                      'temp': (['time','lat','lon'], temp.values),
                      'hum': (['time','lat','lon'], hum.values),
                      'wind_u': (['time','lat','lon'], U10),
                      'wind_v': (['time','lat','lon'], V10),
                      'wind_mag': (['time','lat','lon'], W_speed.values),
                      'fl': (['time','lat','lon'], fl.values)},
                      
                     coords={'lon': (['lon'], nc_with_missing.lon.values),
                             'lat': (['lat'], nc_with_missing.lat.values),
                             'time': pd.date_range(ref_date, periods=len(t_bins),freq='6H')[1:],
                             'reference_time': pd.Timestamp(ref_date)})  
               
    print("Slice  :",round(time.time() - start_time,2),"Sec")
    return ds

def NC_file_com(ref_file,input_folder):
        
#    ref_file=data_dir+"WRF3km_INDIA-2020-02-26_09.nc"
    
    start_time = time.time()
    ref_nc=xr.open_dataset(ref_file)
    urls =glob.glob(input_folder + "WRF3km_INDIA*.nc")
    urls.sort()
    
    # for file in urls:
    #     fl=xr.open_dataset(file)
    #     print(pd.to_datetime(fl.RH2m.attrs['valid_time'], format="%Y-%m-%d_%H"))   
    
    lead_time_hr =list(map (lambda x: x[-16:-3], urls))
    hr_since=(pd.to_datetime(lead_time_hr, format="%Y-%m-%d_%H")-pd.to_datetime(ref_date, format="%Y-%m-%d %H:%M:%S"))// np.timedelta64(1, 'h')

    if(hr_since.max()<72):
        step_tot=72
    else:
        step_tot=hr_since.max()

    dummy_data=np.full([step_tot+1,len(ref_nc.lat), len(ref_nc.lon)], np.nan,dtype="float32")
#    dummy_data=np.full([step_tot+1,len(ref_nc.lat), len(ref_nc.lon)], -1E38,dtype="float32")    
    nc_with_missing = xr.Dataset({'rf': (['time','lat','lon' ],  dummy_data.copy()),
                  'temp': (['time','lat','lon'], dummy_data.copy()),
                  'hum': (['time','lat','lon'], dummy_data.copy()),
                  'wind_u': (['time','lat','lon'], dummy_data.copy()),
                  'wind_v': (['time','lat','lon'], dummy_data.copy()),
                  'wind_mag': (['time','lat','lon'],dummy_data.copy()),
                  'fl': (['time','lat','lon'], dummy_data.copy())},
                  
                 coords={'lon': (['lon'], ref_nc.lon.values),
                         'lat': (['lat'], ref_nc.lat.values),
                         'time': pd.date_range(ref_date, periods=step_tot+1,freq='1H'),
                         'reference_time': pd.Timestamp(ref_date)})

    WRF3km=xr.open_mfdataset(input_folder + "WRF3km_INDIA*.nc",concat_dim='time',combine='nested')
    HI=mpcalc.heat_index(WRF3km.T2m.values *units.degC,WRF3km.RH2m.values/100 * units.dimensionless,mask_undefined=False).to(units.degC)
    W_SPD=np.sqrt(WRF3km.U10**2+WRF3km.V10**2)   
    
    nc_with_missing.rf[np.array(hr_since),:,:]=WRF3km.APCP.values
    nc_with_missing.temp[np.array(hr_since),:,:]=WRF3km.T2m.values
    nc_with_missing.hum[np.array(hr_since),:,:]=WRF3km.RH2m.values
    nc_with_missing.wind_u[np.array(hr_since),:,:]=WRF3km.U10.values
    nc_with_missing.wind_v[np.array(hr_since),:,:]=WRF3km.V10.values
    nc_with_missing.wind_mag[np.array(hr_since),:,:]=W_SPD.values
    nc_with_missing.fl[np.array(hr_since),:,:]=HI.magnitude

    print("NC_read  :",round(time.time() - start_time,2),"Sec")
    
    #wind=xr.merge([nc_with_missing.wind_u,nc_with_missing.wind_u]) 
    #wind.isel(time=range(3,72,6)).to_netcdf("WRF0p25WIND10m.nc")
    
    ds=NC_slice(nc_with_missing,slice_hr=6)
    return ds

def heatmap(input_nc_variable, col, brk,dpi):
    # input_nc_variable=ds.temp
    # col=col_tem
    # brk=brk_tem
    # dpi=200
    
    start_time = time.time()
    interp_methord='bicubic'#bilinear #
    proj = ccrs.epsg(3857)  #set projection
                     #image quality 
    bounds=list(brk)
    cmap =colors.ListedColormap(list(col))    
    norm = colors.BoundaryNorm(bounds, cmap.N)

    img_extent = (min(input_nc_variable.lon.values), max(input_nc_variable.lon.values),
                  min(input_nc_variable.lat.values), max(input_nc_variable.lat.values))  
    os.mkdir(input_nc_variable.name)
    t=0   
    for time_ in list(input_nc_variable.time.values.astype('uint64')):
#    for time_ in np.datetime_as_string(input_nc_variable.time.values,unit='h'):
        time_=str(int(time_/1000000))
#        print(time_,"    ",input_nc_variable.name,t,input_nc_variable.time.values[t])  
        fig = plt.figure(frameon=False,dpi=dpi)
        ax = plt.axes(projection=proj)
        ax.set_extent(img_extent, ccrs.PlateCarree())
        ax.outline_patch.set_visible(False)
        ax.background_patch.set_visible(False)
        ax.imshow(input_nc_variable[t,:,:].values,interpolation=interp_methord, origin='lower',alpha=.8, cmap=cmap, norm=norm, extent=img_extent, transform=ccrs.PlateCarree())
        #ax.contour(grid1_lon,grid1_lat,grid1,alpha=1, levels=bounds, colors='k', transform=ccrs.PlateCarree(),linewidths=0.05)    
        #ax.contour(cntr,grid1_lon,grid1_lat,grid1, colors='black')
        #ax.clabel(cntr, inline=True, fontsize=8)    
        
        fig.savefig(input_nc_variable.name+"/"+time_+'temp.png',bbox_inches='tight',transparent=True,pad_inches = 0.0)
        plt.close()
        convert="convert "+input_nc_variable.name+"/"+time_+'temp.png '+ input_nc_variable.name+"/"+time_+'.png'
        os.system(convert)
        os.remove(input_nc_variable.name+"/"+time_+'temp.png')
        t=t+1
    print(input_nc_variable.name,"Heatmap  :",round(time.time() - start_time,2),"Sec")


def nc2table(nc_vname,mask,mask_table,weight_regrd):

    start_time = time.time()
    regridder = xe.Regridder(nc_vname, mask, 'bilinear',reuse_weights=True,filename=data_dir+weight_regrd)
    regridded_nc = regridder(nc_vname)
    A1=regridded_nc.values
    M1=mask.values
    A=np.append( [M1],A1, axis=0)
    arr = A.transpose(1,2,0)
    df = pd.concat([pd.DataFrame(x) for x in arr], keys=np.arange(max(arr.shape)))
    time_=list(nc_vname.time.values.astype('uint64')/1000000)
    clm = [str(int(i)) for i in time_]  
#    clm = list(pd.to_datetime(nc_vname.time.values).strftime("%Y-%b_%d_%H"))  # version > 0.20.0
    clm.insert(0, "Mask_id")
    df.columns=clm
    df.dropna(subset=['Mask_id'],inplace=True)
    df=df.groupby("Mask_id").mean()
    
    DF=pd.merge(mask_table,df,on="Mask_id",how="inner")  
    DF.to_csv(output_folder+nc_vname.name+".csv")
    print(nc_vname.name+" NC 2 tABLE (MEAN):",round(time.time() - start_time,2),"Sec")
     

def nc2table_max(nc_vname,mask,mask_table,weight_regrd):

    start_time = time.time()
    regridder = xe.Regridder(nc_vname, mask, 'bilinear',reuse_weights=True,filename=data_dir+weight_regrd)
    regridded_nc = regridder(nc_vname)
    A1=regridded_nc.values
    M1=mask.values
    A=np.append( [M1],A1, axis=0)
    arr = A.transpose(1,2,0)
    df = pd.concat([pd.DataFrame(x) for x in arr], keys=np.arange(max(arr.shape)))
    time_=list(nc_vname.time.values.astype('uint64')/1000000)
    clm = [str(int(i)) for i in time_]  
#    clm = list(pd.to_datetime(nc_vname.time.values).strftime("%Y-%b_%d_%H"))  # version > 0.20.0
    clm.insert(0, "Mask_id")
    df.columns=clm
    df.dropna(subset=['Mask_id'],inplace=True)
    df=df.groupby("Mask_id").max()
    
    DF=pd.merge(mask_table,df,on="Mask_id",how="inner")  
    DF.to_csv(output_folder+nc_vname.name+".csv")
    print(nc_vname.name+" NC 2 tABLE (MAX) :",round(time.time() - start_time,2),"Sec")


def Save_nc(ds):

    start_time = time.time()
    ds.rf.to_netcdf(output_folder + "rf.nc")
    ds.temp.to_netcdf(output_folder + "temp.nc")
    ds.hum.to_netcdf(output_folder + "hum.nc")
    ds.wind_u.to_netcdf(output_folder + "wind_u.nc")
    ds.wind_v.to_netcdf(output_folder + "wind_v.nc")
    ds.wind_mag.to_netcdf(output_folder + "wind_mag.nc")
    ds.fl.to_netcdf(output_folder + "fl.nc")
    print("NC file saved:",round(time.time() - start_time,2),"Sec")

def Regridp25(ds):

    lat_new=np.linspace(ds['lat'].values.min(), 
                        ds['lat'].values.max(), num=int((ds['lat'].values.max()-ds['lat'].values.min())/0.25)+1, endpoint=True)
    lon_new=np.linspace(ds['lon'].values.min(), 
                            ds['lon'].values.max(), num=int((ds['lon'].values.max()-ds['lon'].values.min())/0.25)+1, endpoint=True)
    ds_out = xr.Dataset({'lat': (['lat'], lat_new),
                         'lon': (['lon'], lon_new),
                        })
#    ds_out
    regridder = xe.Regridder(ds, ds_out, 'bilinear',reuse_weights=True);
    regridder  # print basic regridder information.
    
    U10p25 = regridder(ds.wind_u)
    V10p25 = regridder(ds.wind_v)
    
    U10p25.to_netcdf(wind_output_folder + "wind_u.nc")
    V10p25.to_netcdf(wind_output_folder + "wind_v.nc")
    
    wind=xr.merge([U10p25,V10p25])   

    wind.to_netcdf(wind_output_folder + "wind.nc")
    
def Save_table(ds):

    mask_table=pd.read_csv(data_dir+"Block_uuid.csv")
    mask_table=mask_table[['OBJECTID_1',"State_Name",'district', 'block',  'UUID']]
    mask_table.columns=['Mask_id',"State_Name",'district', 'block',  'UUID']
    
    mask = xr.open_dataset(data_dir+'Block.nc')
    
    nc2table_max(ds.temp,mask.Band1,mask_table,"bilinear_1284x1284_3032x2923.nc") 
    nc2table_max(ds.hum,mask.Band1,mask_table,"bilinear_1284x1284_3032x2923.nc") 
    nc2table(ds.rf,mask.Band1,mask_table,"bilinear_1284x1284_3032x2923.nc") 
    nc2table_max(ds.wind_mag,mask.Band1,mask_table,"bilinear_1284x1284_3032x2923.nc") 
    nc2table_max(ds.fl,mask.Band1,mask_table,"bilinear_1284x1284_3032x2923.nc") 
        
def Save_heatmaps(ds):
#    global data_dir
#    global output_folder
    brk_rh= pd.read_csv(data_dir+"brk_rh.csv")["x"]
    col_rh= pd.read_csv(data_dir+"col_rh.csv")["x"]   
    
    brk_tem= pd.read_csv(data_dir+"brk_tem.csv")["x"]
    col_tem= pd.read_csv(data_dir+"col_tem.csv")["x"]
    
    brk_ws= pd.read_csv(data_dir+"brk_ws.csv")["x"]
    col_ws= pd.read_csv(data_dir+"col_ws.csv")["x"]
    
    brk_rf=[0,.1,2.5,15.6,64.5,115.6,204.5,8000]
    col_rf=["#FFFFFF","#C3FDCA","#01FF04","#048500","#FDC0CB","#FC0300","#610301"]
    
    os.chdir(output_folder)
    #heatmap(ds.temp,col_tem,brk_tem,dpi=200)
    
    Process(target=heatmap, args=(ds.temp,col_tem,brk_tem,200)).start()
    Process(target=heatmap, args=(ds.fl,col_tem,brk_tem,200)).start()
    Process(target=heatmap, args=(ds.hum,col_rh,brk_rh,200) ).start()
    Process(target=heatmap, args=(ds.wind_mag,col_ws,brk_ws,150) ).start()   
    Process(target=heatmap, args=(ds.rf,col_rf,brk_rf,200) ).start()
def Download(Ini_date):
    global input_folder
    start_time = time.time()

    date_time=pd.date_range(start=Ini_date, periods=75,freq='H')
    input_folder=input_folder+date_time[0].strftime("%Y%m%d%H/")
    os.mkdir(input_folder)
    os.chdir(input_folder)
    
    for tm in date_time[1:]:
        nc_url=date_time[0].strftime("%Y%m%d%H/")+"WRF3km_INDIA-"+tm.strftime("%Y-%m-%d_%H.nc")
        URL="axel -a -n 4 ftp://nwp:nwp@125.21.185.50/nwp-data/WRF_NETCDF/"+nc_url
    #    print(YRL)
        os.system(URL)        
    print("Downloded in :",round(time.time() - start_time,2),"Sec")
    

if __name__ == '__main__':
    
    '''
    python combine_and_generate_6_hour_nc_with_HI.py '2020-02-26 00:00:00' "/home/vassar/Documents/forcastdata/WRF/" "/home/vassar/temp/" "/home/vassar/temp/wind/" "/home/vassar/Documents/forecast_vassarlabs/data/"
    
    argv1=reference date and time
    argv2=input folder
    argv3=output folder
    argv4=meta data folder
    
    '''
    
    start_time = time.time()
    ref_date=sys.argv[1]  #'2020-02-26 00:00:00'
    input_folder= sys.argv[2]
    output_folder = sys.argv[3]
    wind_output_folder = sys.argv[4]
    data_dir = sys.argv[5]
    
    
    # ref_date='2020-02-26 00:00:00' 
    # input_folder="/home/vassar/Documents/forcastdata/WRF/" 
    # output_folder="temp/" 
    # wind_output_folder="temp/"
    # data_dir= "/home/vassar/Documents/forecast_vassarlabs/data/"
    
    #Download(ref_date)   #downloading

    ds=NC_file_com(data_dir+"WRF3km_INDIA-2020-02-26_09.nc",input_folder)
    Save_nc(ds) 
    Regridp25(ds) 
    Save_table(ds) 
    Save_heatmaps(ds)  
    print(time.time() - start_time)








