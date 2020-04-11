#!/usr/bin/env python
# coding: utf-8

# In[1]:


from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import tqdm
import os
import re
import copy
import datetime
import numba
from numba import jit,prange


# ### Setup

# In[2]:


infolder = "D:/data/rgAgCFSR/regrid/"

calname = "Sacks_ZARC_fill_fill_120d"

outfolder = "rgagcfsr_computed_partial/"+calname+"/"

masterpref = "AgCFSR_"
tempsuf = "_tavg.nc4"
tmaxsuf = "_tmax.nc4"
tminsuf = "_tmin.nc4"

calfolder = "bcagcfsr_calendars/"+calname+"/"
calsuf = ".crop.calendar.fill.nc"
planvar = "plant.start"
harvvar = "harvest.end"

mskfname = infolder+"seamask.nc"

crops = ["Soybeans","Maize","Cotton"]
# crops = ["Maize","Cotton"]
# crops = ["Soybeans"]

alims = {"Soybeans":10.0, "Maize":10.0, "Cotton":15.0}
blims = {"Soybeans":30.0, "Maize":29.0, "Cotton":32.0}
clims = {"Soybeans":31.4, "Maize":29.9, "Cotton":34.3}
flims = {"Soybeans":38.0, "Maize":38.0, "Cotton":38.0}

# The basis of calculation will be harvest years
# hyears = [2002]
# hyears = [2003]
# hyears = [2009]
# hyears = [2009,2010,2011,2012,2013,2014]
# hyears = [2010,2011,2012,2013,2014]
hyears = list(range(2002,2008+1))
# hyears = list(range(2002,2003+1))
print(hyears)

deltats = [0,1,2,3,4,5]


# ### Defining functions

# In[3]:


# Opens and concatenates a harvest year with the equivalent planting year
# def concat_clim(infolder,climpref,hyear):
def concat_clim(infolder,masterpref,climsuf,hyear):
    pyear = hyear - 1

#     climharr = xr.open_dataarray(infolder+climpref+str(hyear)+".nc")
#     climparr = xr.open_dataarray(infolder+climpref+str(pyear)+".nc")
#     climharr = xr.open_dataarray(infolder+masterpref+str(hyear)+climsuf+".nc")
#     climparr = xr.open_dataarray(infolder+masterpref+str(pyear)+climsuf+".nc")
    climharr = xr.open_dataarray(infolder+masterpref+str(hyear)+climsuf, decode_times = False)
    climparr = xr.open_dataarray(infolder+masterpref+str(pyear)+climsuf, decode_times = False)

    climarr = xr.concat([climparr,climharr], dim = "time")
    return climarr


# In[4]:


# Calculates AGDD for a single day (ascending numbers, constrained in both sides)
@jit(nopython=True)
def calc_agdd_day(Tmin,Tmax,Tlo,Thi):
    res = 0.005 #Resolution (dt, in days) on which to evaluate the T sine curve

    t = np.arange(0,1,res)
    nt = t.shape[0]

    Tamp = (Tmax-Tmin)/2.0
    Tmed = (Tmax+Tmin)/2.0

    T = Tmed + Tamp*np.sin(t*(2.0*np.pi/1))

    gdd = np.sum(np.invert(np.isnan(
        np.where((T>=Tlo)&(T<=Thi),T,np.nan)
    ))*(T-Tlo))/nt

    return gdd

@jit(nopython=True)
def calc_agdd_point(tmaxvec,tminvec,Tlo,Thi):
    gdd = 0.0
    for day in range(tminvec.shape[0]):
        gdd = gdd + calc_agdd_day(tminvec[day],tmaxvec[day],Tlo,Thi)
    return gdd


# In[5]:


# Calculates BGDD for a single day (descending numbers, constrained in both sides)
@jit(nopython=True)
def calc_bgdd_day(Tmin,Tmax,Tlo,Thi):
    res = 0.005 #Resolution (dt, in days) on which to evaluate the T sine curve

    t = np.arange(0,1,res)
    nt = t.shape[0]

    Tamp = (Tmax-Tmin)/2.0
    Tmed = (Tmax+Tmin)/2.0

    T = Tmed + Tamp*np.sin(t*(2.0*np.pi/1))

    gdd = np.sum(np.invert(np.isnan(
        np.where((T>=Tlo)&(T<=Thi),T,np.nan)
    ))*(T-Thi))/nt

    return gdd

@jit(nopython=True)
def calc_bgdd_point(tmaxvec,tminvec,Tlo,Thi):
    gdd = 0.0
    for day in range(tminvec.shape[0]):
        gdd = gdd + calc_bgdd_day(tminvec[day],tmaxvec[day],Tlo,Thi)
    return gdd


# In[6]:


# Calculates CGDD for a single day (ascending numbers, constrained just below)
@jit(nopython=True)
def calc_cgdd_day(Tmin,Tmax,Tlo):
    res = 0.005 #Resolution (dt, in days) on which to evaluate the T sine curve

    t = np.arange(0,1,res)
    nt = t.shape[0]

    Tamp = (Tmax-Tmin)/2.0
    Tmed = (Tmax+Tmin)/2.0

    T = Tmed + Tamp*np.sin(t*(2.0*np.pi/1))

    gdd = np.sum(np.invert(np.isnan(
        np.where((T>=Tlo),T,np.nan)
    ))*(T-Tlo))/nt

    return gdd

@jit(nopython=True)
def calc_cgdd_point(tmaxvec,tminvec,Tlo):
    gdd = 0.0
    for day in range(tminvec.shape[0]):
        gdd = gdd + calc_cgdd_day(tminvec[day],tmaxvec[day],Tlo)
    return gdd


# In[7]:


# Calculates everything for the growing season given numpy arrays. 
# Loops throught the points
# Compiles with Numba parallel
@jit(nopython=True, parallel = True)
# @jit(nopython=True)
def calc_all(planmat,harvmat,tempmat,tmaxmat,tminmat,
             tempmeanmat,tmaxmeanmat,tminmeanmat,
             trngmeanmat,ndaymat,
             agddmat,bgddmat,cgddmat,fgddmat,
             alim,blim,clim,flim):
    for lati in prange(tempmat.shape[1]):
        for lonj in range(tempmat.shape[2]):
            if (np.isnan(planmat[lati,lonj])) or (np.isnan(tempmat[0,lati,lonj])):
                continue
            plan = int(planmat[lati,lonj])
            harv = int(harvmat[lati,lonj])

            tempvec = tempmat[plan:harv,lati,lonj]
            tmaxvec = tmaxmat[plan:harv,lati,lonj]
            tminvec = tminmat[plan:harv,lati,lonj]

            tempmeanmat[lati,lonj] = np.nanmean(tempvec)
            tmaxmeanmat[lati,lonj] = np.nanmean(tmaxvec)
            tminmeanmat[lati,lonj] = np.nanmean(tminvec)
            trngmeanmat[lati,lonj] = np.nanmean(tmaxvec - tminvec)
            ndaymat[lati,lonj] = np.int64(harv-plan)

            agddmat[lati,lonj] = calc_agdd_point(tmaxvec,tminvec,alim,blim)
            bgddmat[lati,lonj] = calc_bgdd_point(tmaxvec,tminvec,blim,clim)
            cgddmat[lati,lonj] = calc_cgdd_point(tmaxvec,tminvec,clim)
            fgddmat[lati,lonj] = calc_agdd_point(tmaxvec,tminvec,clim,flim)


# ### Begin main script

# In[8]:


# Create output folder
if not os.path.exists(outfolder): os.makedirs(outfolder, exist_ok=True)

# Read the mask
mskarr = xr.open_dataarray(mskfname)
    
# crop = crops[0]
for crop in crops:
    print(crop)

    # Open calendar and convert it to two-year based indexes. FIXME: Ignoring leap years here
    caldata = xr.open_dataset(calfolder+crop+calsuf)
    planarr = caldata[planvar]
    harvarr = caldata[harvvar]
    harvarr = xr.where(harvarr < planarr,harvarr + 365,harvarr) - 1 

    # Define DD limits
    alim = alims[crop]
    blim = blims[crop]
    clim = clims[crop]
    flim = flims[crop]


            #FIXME: Loop here
#     hyear = hyears[0]
#     deltat = deltats[0]
    for hyear in tqdm.tqdm(hyears):
        for deltat in deltats:
            print(crop + " deltat = " + str(deltat))
            # Open the climate arrays, concatenating them
            temparr = concat_clim(infolder,masterpref,tempsuf,hyear) + deltat
            tmaxarr = concat_clim(infolder,masterpref,tmaxsuf,hyear) + deltat
            tminarr = concat_clim(infolder,masterpref,tminsuf,hyear) + deltat

            temparr = temparr.where(mskarr == 1)
            tmaxarr = tmaxarr.where(mskarr == 1)
            tminarr = tminarr.where(mskarr == 1)

            # Preallocate the arrays that will be filled
#             lldims = ("latitude","longitude")
            lldims = ("lat","lon")
            coords = [(i,temparr.coords[i].data,temparr.coords[i].attrs) for i in lldims] # Tuples with lat and lon dimension specs

            # 2D arrays
            tempmean = xr.DataArray(coords = coords, name = "tempmean")
            tmaxmean = xr.DataArray(coords = coords, name = "tmaxmean")
            tminmean = xr.DataArray(coords = coords, name = "tminmean")

            trngmean = xr.DataArray(coords = coords, name = "trngmean")

            ndayarr = xr.DataArray(coords = coords, name = "ndays")

            agddarr = xr.DataArray(coords = coords, name = "agdd")
            bgddarr = xr.DataArray(coords = coords, name = "bgdd")
            cgddarr = xr.DataArray(coords = coords, name = "cgdd")
            fgddarr = xr.DataArray(coords = coords, name = "fgdd")

            # This basically creates pointers to the numpy arrays inside the xr.Dataarrays
            # We need those for numba to work. An alternative would be passing the .data in the function call
            planmat = planarr.data
            harvmat = harvarr.data

            tempmat = temparr.data
            tmaxmat = tmaxarr.data
            tminmat = tminarr.data

            tempmeanmat = tempmean.data
            tmaxmeanmat = tmaxmean.data
            tminmeanmat = tminmean.data

            trngmeanmat = trngmean.data

            ndaymat = ndayarr.data

            agddmat = agddarr.data
            bgddmat = bgddarr.data
            cgddmat = cgddarr.data
            fgddmat = fgddarr.data



            # Calculates everything.
            calc_all(planmat,harvmat,tempmat,tmaxmat,tminmat,
                         tempmeanmat,tmaxmeanmat,tminmeanmat,
                         trngmeanmat,ndaymat,
                         agddmat,bgddmat,cgddmat,fgddmat,
                         alim,blim,clim,flim)


            # Merge everything in a single Dataset
            #         outdata = xr.merge([tempmean,tmaxmean,tminmean,trngmean,ndayarr,tempdist,tempgdds])
            outdata = xr.merge([tempmean,tmaxmean,tminmean,trngmean,ndayarr,agddarr,bgddarr,cgddarr,fgddarr])
            outdata.attrs['Crop'] = crop
            outdata.attrs['harvest_year'] = hyear
            outdata.attrs['calendar_path'] = calfolder+crop+calsuf
            outdata.attrs['climdata_path'] = infolder
            #         outdata.attrs['climdata_ex_path'] = infolder+temppref+str(hyear)+".nc"

            # Write output
            outfname = outfolder + crop + ".computed.deltat." + str(deltat) + "." + str(hyear) + ".nc"
            outdata.to_netcdf(outfname,
                      engine = "netcdf4")


# In[9]:


Tmax = 25
Tmin = 15
Tlo = 15
Thi = 25


res = 0.005 #Resolution (dt, in days) on which to evaluate the T sine curve

t = np.arange(0,1,res)
nt = t.shape[0]

Tamp = (Tmax-Tmin)/2.0
Tmed = (Tmax+Tmin)/2.0

T = Tmed + Tamp*np.sin(t*(2.0*np.pi/1))

gdd = np.sum(np.invert(np.isnan(
    np.where((T>=Tlo)&(T<=Thi),T,np.nan)
))*(T-Tlo))/nt

gdd


# In[10]:


# lati = 260
# lonj = 150
# planmat[lati,lonj]


# In[11]:


# calc_all(planmat,harvmat,tempmat,tmaxmat,tminmat,
#              tempmeanmat,tmaxmeanmat,tminmeanmat,
#              trngmeanmat,ndaymat,
#              agddmat,bgddmat,cgddmat,fgddmat,
#              alim,blim,clim,flim)


# In[12]:


# # Calculates everything for the growing season given numpy arrays. 
# # Loops throught the points and calculates both distribution and regular gs means
# # Compiles with Numba parallel
# @jit(nopython=True, parallel = True)
# # @jit(nopython=True)
# def calc_all(planmat,harvmat,tempmat,tmaxmat,tminmat,
#              tempmeanmat,tmaxmeanmat,tminmeanmat,
#              trngmeanmat,ndaymat,
#              agddmat,bgddmat,cgddmat,fgddmat,
#              alim,blim,clim,flim):
#     for lati in prange(tempmat.shape[1]):
#         for lonj in range(tempmat.shape[2]):
#             if (np.isnan(planmat[lati,lonj])) or (np.isnan(tempmat[0,lati,lonj])):
#                 continue
#             plan = int(planmat[lati,lonj])
#             harv = int(harvmat[lati,lonj])

#             tempvec = tempmat[plan:harv,lati,lonj]
#             tmaxvec = tmaxmat[plan:harv,lati,lonj]
#             tminvec = tminmat[plan:harv,lati,lonj]

#             tempmeanmat[lati,lonj] = np.nanmean(tempvec)
#             tmaxmeanmat[lati,lonj] = np.nanmean(tmaxvec)
#             tminmeanmat[lati,lonj] = np.nanmean(tminvec)
#             trngmeanmat[lati,lonj] = np.nanmean(tmaxvec - tminvec)
#             ndaymat[lati,lonj] = np.int64(harv-plan)

#             agddmat[lati,lonj] = calc_agdd_point(tmaxvec,tminvec,alim,blim)
#             bgddmat[lati,lonj] = calc_bgdd_point(tmaxvec,tminvec,blim,clim)
#             cgddmat[lati,lonj] = calc_cgdd_point(tmaxvec,tminvec,clim)
#             fgddmat[lati,lonj] = calc_agdd_point(tmaxvec,tminvec,clim,flim)


# In[13]:


# # Create output folder
# if not os.path.exists(outfolder): os.makedirs(outfolder, exist_ok=True)

# # Read the mask
# mskarr = xr.open_dataarray(mskfname)
    
# crop = crops[0]
# # for crop in crops:
# print(crop)

# # Open calendar and convert it to two-year based indexes. FIXME: Ignoring leap years here
# caldata = xr.open_dataset(calfolder+crop+calsuf)
# planarr = caldata[planvar]
# harvarr = caldata[harvvar]
# harvarr = xr.where(harvarr < planarr,harvarr + 365,harvarr) - 1 

# # Define DD limits
# alim = alims[crop]
# blim = blims[crop]
# clim = clims[crop]
# flim = flims[crop]


#         #FIXME: Loop here
# hyear = hyears[0]
# deltat = deltats[0]
        
# print(crop + " deltat = " + str(deltat))
# # Open the climate arrays, concatenating them
# temparr = concat_clim(infolder,masterpref,tempsuf,hyear) + deltat
# tmaxarr = concat_clim(infolder,masterpref,tmaxsuf,hyear) + deltat
# tminarr = concat_clim(infolder,masterpref,tminsuf,hyear) + deltat

# temparr = temparr.where(mskarr == 1)
# tmaxarr = tmaxarr.where(mskarr == 1)
# tminarr = tminarr.where(mskarr == 1)

# # Generate a vector of lower T bounds to use as metadata and speed up computation
# Tlos = np.arange(Tlo,Thi,Tint)

# # Preallocate the arrays that will be filled
# #             lldims = ("latitude","longitude")
# lldims = ("lat","lon")
# coords = [(i,temparr.coords[i].data,temparr.coords[i].attrs) for i in lldims] # Tuples with lat and lon dimension specs

# # 2D arrays
# tempmean = xr.DataArray(coords = coords, name = "tempmean")
# tmaxmean = xr.DataArray(coords = coords, name = "tmaxmean")
# tminmean = xr.DataArray(coords = coords, name = "tminmean")

# trngmean = xr.DataArray(coords = coords, name = "trngmean")

# ndayarr = xr.DataArray(coords = coords, name = "ndays")

# agddarr = xr.DataArray(coords = coords, name = "agdd")
# bgddarr = xr.DataArray(coords = coords, name = "bgdd")
# cgddarr = xr.DataArray(coords = coords, name = "cgdd")
# fgddarr = xr.DataArray(coords = coords, name = "fgdd")

# # This basically creates pointers to the numpy arrays inside the xr.Dataarrays
# # We need those for numba to work. An alternative would be passing the .data in the function call
# planmat = planarr.data
# harvmat = harvarr.data

# tempmat = temparr.data
# tmaxmat = tmaxarr.data
# tminmat = tminarr.data

# tempmeanmat = tempmean.data
# tmaxmeanmat = tmaxmean.data
# tminmeanmat = tminmean.data

# trngmeanmat = trngmean.data

# ndaymat = ndayarr.data

# agddmat = agddarr.data
# bgddmat = bgddarr.data
# cgddmat = cgddarr.data
# fgddmat = fgddarr.data



# # Calculates everything.
# calc_all(planmat,harvmat,tempmat,tmaxmat,tminmat,
#              tempmeanmat,tmaxmeanmat,tminmeanmat,
#              trngmeanmat,ndaymat,
#              agddmat,bgddmat,cgddmat,fgddmat,
#              alim,blim,clim,flim)


# # Merge everything in a single Dataset
# #         outdata = xr.merge([tempmean,tmaxmean,tminmean,trngmean,ndayarr,tempdist,tempgdds])
# outdata = xr.merge([tempmean,tmaxmean,tminmean,trngmean,ndayarr,agddarr,bgddarr,cgddarr,fgddarr])
# outdata.attrs['Crop'] = crop
# outdata.attrs['harvest_year'] = hyear
# outdata.attrs['calendar_path'] = calfolder+crop+calsuf
# outdata.attrs['climdata_path'] = infolder
# #         outdata.attrs['climdata_ex_path'] = infolder+temppref+str(hyear)+".nc"

# # Write output
# outfname = outfolder + crop + ".computed.deltat." + str(deltat) + "." + str(hyear) + ".nc"
# outdata.to_netcdf(outfname,
#           engine = "netcdf4")


# In[14]:


# outdata["agdd"].plot()


# In[15]:


# calc_gdd_point(tmaxvec,tminvec,Tlos)


# In[16]:


# These plots are good for understanding the distribution of T during the day
# Tl1 = 15
# Tl2 = Tl1 + 1
# l =(np.invert(np.isnan(
#             np.where((T>=Tl1) & (T<=Tl2) ,T,np.nan)
#         )))
# print(np.sum(l)/nt)

# plt.subplot(2, 1, 1)
# plt.plot(t,T)
# plt.axhline(y=Tl1)
# plt.axhline(y=Tl2)
# plt.ylim((10,30))
# plt.subplot(2, 1, 2)
# plt.plot(t,l)
# plt.show()

# Tmin = 15
# Tmax = 25
# Tlo = -5.0
# Thi = 50.0
# Tint = 1.0
# plt.subplot(1, 2, 1)
# plt.plot(t,T,)
# plt.ylim((10,30))
# plt.subplot(1, 2, 2)
# plt.plot(calc_dist_point(Tmin,Tmax,Tlo,Thi,Tint),Tl1s)
# plt.ylim((10,30))
# plt.show()

