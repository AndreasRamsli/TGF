#!/usr/bin/env python
# coding: utf-8

# In[171]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from math import radians, cos, sin, asin, sqrt
import os
import re
import datetime
import csv


# # Prepping of the data

# ## AGILE

# In[264]:


#####  Prepping the AGILE data frame  #######
    
# AGILE csv file gotten from: https://www.ssdc.asi.it/mcal3tgfcat/#start_page

AGILE_data = pd.read_csv("AGILE_plot.csv") 

#Deleting unnecessary columns
AGILE_data = AGILE_data.drop(columns=['Unnamed: 8'], axis = 1)

AGILE_data = AGILE_data.rename(columns={' GeoLon(deg) ':'Lon', ' GeoLat(deg) ': "Lat", ' Satellite Altitude(km) ': "Alt",' Date(UTC) ': "Date(UTC)", ' T0 Local Time(h.dec) ': 'Local Hour'})


# Local time given as a decimal. Rounding the number no nearest hour ####
AGILE_data["Local Hour"] = np.round(AGILE_data["Local Hour"])


# Finding the local date and hour for the TGF phenomena
for index, row in AGILE_data.iterrows():
    
    date = row['Date(UTC)'].split('-') 
    year, month, day = int(date[0]),int(date[1]),int(date[2][0:2])
    time = row['Date(UTC)'].split(':')
    hour, minute, second = int(time[0][-2:]), int(time[1]), int(time[2])
    
    date_UTC = datetime.datetime(year, month, day, hour, minute, second)
    data_add = datetime.timedelta(minutes = 4*row['Lon'])
    
    local_date = date_UTC + data_add
    local_hour = local_date.hour
    
    AGILE_data.at[index,"Local Date"] = local_date.strftime("%d-%m-%Y")
    AGILE_data.at[index,"Local Hour"] = local_hour

#Creating a new column for the duration of the TGF. Given in ms. 
AGILE_data["Duration"] = ''

for index, row in AGILE_data.iterrows():
    T50 = row[' T50']
    AGILE_data.at[index, "Duration"] = T50 * 2 #Doubling because this is the duration of the whole event


AGILE_data.head()


# ## FERMI

# In[283]:


#####  Prepping the FERMI data frame  #######
    
#FERMI csv gotten from: https://fermi.gsfc.nasa.gov/ssc/data/access/gbm/tgf/    

FERMI_data = pd.read_csv("./gbm_tgf_catalog_offline.csv")

#Deleting unnecessary columns
FERMI_data = FERMI_data.drop([' MET', ' file', ' P2', ' Alt', ' LST', ' TRIG_ID'], axis = 1)
FERMI_data = FERMI_data.rename(columns={' Width_ms': 'Duration'})

#Creating a new colum for the all the total counts
FERMI_data[" Total Counts"] = FERMI_data[" BGO_0_N"] + FERMI_data[" BGO_1_N"]  #TGF_data[" NAI_N"]
FERMI_data[' Local Date'] = ''
FERMI_data[' Local Hour'] = ''

#Correcting the longitude to be centered around +- 180 degrees:
for label, row in FERMI_data.iterrows():
    if row[' Lon'] > 180:
        FERMI_data.at[label, " Lon"] = row[' Lon'] - 360
    else:
        pass
    
# Rounding the lon and lat values
FERMI_data[" Lat"] = np.round(FERMI_data[" Lat"],1)
FERMI_data[" Lon"] = np.round(FERMI_data[" Lon"],1)


# Finding the local date and hour for the TGF phenomena
for index, row in FERMI_data.iterrows():
    
    date = row[' Date'].split('-') 
    year, month, day = int(date[0]),int(date[1]),int(date[2])
    
    time = row[' UTC'].split(':')
    hour, minute, second = int(time[0]), int(time[1]), int(time[2][:2])
    
    date_UTC = datetime.datetime(year, month, day, hour, minute, second)
    data_add = datetime.timedelta(minutes = 4*row[' Lon'])
    local_date = date_UTC + data_add
    local_hour = local_date.hour
    FERMI_data.iloc[index, 10] = local_date
    FERMI_data.iloc[index, 11] = local_hour
    
FERMI_data.head()


# ## RHESSI

# In[279]:


#####  Prepping the RHESSI data frame  #######
    
#FERMI csv gotten from: http://scipp.pbsci.ucsc.edu/rhessi/    
    
RHESSI_data = pd.read_csv('RHESSI.csv', names=["BS", "Timestamp", "Lat", "Lon", "Counts", "Local Date", "Local Hour"])
        
        # Assigning new values to each column #
for index, row in RHESSI_data.iterrows():
    if index > 0:
        RHESSI_data.loc[index, 'Timestamp'] = row[0][9:32]
        RHESSI_data.at[index,'Lat'] = row[0][39:47]
        RHESSI_data.at[index, 'Lon'] = row[0][53:62]
        RHESSI_data.at[index, 'Counts'] = row[0][69:72]
        RHESSI_data.at[index, 'Duration'] = row[0][78:]
    else:
        pass

            # Trimming the dataframe #
RHESSI_data = RHESSI_data.drop("BS", axis = 1)
RHESSI_data = RHESSI_data.drop([0], axis = 0)

        # Redefining the lon cordinates #
for label, row in RHESSI_data.iterrows():
    if row['Lon'] > 180:
        RHESSI_data.at[label, "Lon"] = row['Lon'] - 360
    else:
        pass
    
        # Rounding the lon and lat values
RHESSI_data["Lat"] = np.round(RHESSI_data["Lat"],1)
RHESSI_data["Lon"] = np.round(RHESSI_data["Lon"],1)

        # Cleaning up the longitude column
for index, row in RHESSI_data.iterrows():
    lon = row[2]
    if abs(lon) > 180:
        RHESSI_data = RHESSI_data.drop(RHESSI_data.index[index-1])
        

        # Finding the local date and hour for the TGF phenomena
for index, row in RHESSI_data.iterrows():
    
    date = str(row['Timestamp']).split('-')
    year, month, day = int(date[0]), int(date[1]), int(date[2][0:2])
    time = row['Timestamp'].split(':')
    hour, minute, second = int(time[0][-2:]), int(time[1]), int(time[2][0:2])
    
    date_UTC = datetime.datetime(year, month, day, hour, minute, second)
    data_add = datetime.timedelta(minutes = 4*row['Lon'])
    
    local_date = date_UTC + data_add
    local_hour = local_date.hour
    
    RHESSI_data.loc[index,"Local Date"] = local_date.strftime("%Y-%m-%d")
    RHESSI_data.at[index,"Local Hour"] = local_hour
    
RHESSI_data.head()


# # PLOTTING

# In[305]:


#### Creating a template for the plots (2x2) ####

fig, axs = plt.subplots(figsize=(14,10), nrows=2, ncols=2)


#colors = ['black', 'red', 'blue']


                                    ### Longitude ###
axs[0, 0].hist(AGILE_data["Lon"], bins = 'auto', edgecolor= "black", histtype="step", hatch='//', label='AGILE')
axs[0, 0].set_title("Distribution of TGF longitudes")
axs[0, 0].hist(FERMI_data[" Lon"], bins = 'auto', edgecolor = "red", histtype="step", label='FERMI')
axs[0, 0].hist(RHESSI_data["Lon"], bins = 'auto', edgecolor = "blue", histtype='step', linestyle='dashed', label='RHESSI')
axs[0, 0].legend(prop={'size': 10})


                                    ### Latitude ###
axs[0, 1].hist(AGILE_data["Lat"], bins = 'auto', edgecolor= "black", histtype="step", hatch='//', label='AGILE')
axs[0, 1].hist(FERMI_data[" Lat"], bins = 'auto', edgecolor = "red", histtype="step", label='FERMI')
axs[0, 1].hist(RHESSI_data["Lat"], bins = 'auto', edgecolor = "blue", histtype='step', linestyle='dashed', label='RHESSI')
axs[0, 1].set_title("Distribution of TGF latitudes")
axs[0, 1].legend(prop={'size': 10})


                                    ### Local Hour ###
x = np.arange(0,25,1)
axs[1, 0].hist(AGILE_data["Local Hour"], bins = x, edgecolor= "black", histtype="step", hatch='//', label='AGILE')
axs[1, 0].hist(FERMI_data[" Local Hour"], bins = x, edgecolor = "red", histtype="step", label='FERMI')
axs[1, 0].hist(RHESSI_data["Local Hour"], bins = x, edgecolor = "blue", histtype='step', linestyle='dashed', label='RHESSI')
axs[1, 0].set_title("Distribution of TGF local hour")
axs[1, 0].legend(prop={'size': 10})

    
                                    ### Counts ###
y = np.arange(1,1e2)
axs[1, 1].hist(AGILE_data[' ML Counts '], bins = y, edgecolor= "black", histtype="step", hatch='//', label='AGILE')
axs[1, 1].hist(FERMI_data[" Total Counts"], bins = y, edgecolor = "red", histtype="step", label='FERMI')
axs[1, 1].hist(RHESSI_data["Counts"], bins = y, edgecolor = "blue", histtype='step', linestyle='dashed', label='RHESSI')
axs[1, 1].set_title("Distribution of TGF counts")
axs[1, 1].legend(prop={'size': 10})
axs[1,1].set_xscale('log')
    
                                    ### Labels ###
plt.setp(axs[0,0], xlabel='Longitude', ylabel='Number of counts')
plt.setp(axs[0,1], xlabel="Latitude", ylabel='Number of counts')
plt.setp(axs[1,0], xlabel="Local Hour", ylabel='Number of counts')
plt.setp(axs[1,1], xlabel="Counts", ylabel='Number of counts in each bin')


plt.show()


# In[317]:


plt.subplots(figsize=(12,8))
x = np.arange(1e-2,1e2,0.5)
plt.hist(AGILE_data["Duration"], bins = 'auto', edgecolor= "black", histtype="step", hatch='//', label='AGILE')
plt.hist(FERMI_data["Duration"], bins = 'auto', edgecolor = "red", histtype="step", label='FERMI')
plt.hist(RHESSI_data["Duration"], bins = x, edgecolor = "blue", histtype='step', linestyle='dashed', label='RHESSI')
plt.title("Distribution of TGF duration")
plt.legend(prop={'size': 10})
plt.xlabel("Duration in ms")
plt.ylabel("Number of counts")
plt.xscale('log')

plt.show()


# In[ ]:




