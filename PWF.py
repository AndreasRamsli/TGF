#!/usr/bin/env python
# coding: utf-8

# # Program Workflow

# ### A script for finding the closest match between TGF and lightningstroke from WWLLN
# #### First find the closest in distance, then account for the light travel time

# In[225]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from math import radians, cos, sin, asin, sqrt
import os
import re
import datetime
import csv


# ### Prepping the data

# #### Need to have the same format for the TGFID for each satelite. Using AGILE as example

# In[155]:


AGILE_data = pd.read_csv("AGILE_data.csv", sep=",")
AGILE_data = AGILE_data.drop(columns=['Unnamed: 9'])
AGILE_data = AGILE_data.rename(columns={' GeoLon(deg) ':'Lon', ' GeoLat(deg) ': "Lat", 
                                        ' Satellite Altitude(km) ': "Alt",' Date(UTC) ': "Date(UTC)",
                                        'New TGFID ': "New TGFID", ' W Longitude(deg) ': 
                                        'W lon',' W Latitude(deg) ': 'W lat', ' W Dist(km) ': 'W Dist (H)'})
        
'''Removing whitespace before and after the string in the New TGFID column'''
AGILE_data['New TGFID'] = AGILE_data['New TGFID'].str.strip()

'''Creating a list of all the WWLLN filenames'''
lightning_data = os.listdir("./WWLLN_data")


AGILE_data.head()


# # Program Workflow

# In[268]:


######## Program Workflow #########
    
    
# 1 Find the date and time of the TGF
# 2 Search the WWLLN folder. -> Create a list of possible connections that day. Filter by the closest to the satelite.
# 3 Check if the stroke is within a +/- 500 us window
# 4 Return a dataframe of TGFs that have a lightning connection


time_match_dict = {} 
"""  time_match_dict: 
Dictionary of TGFs that might have a lightning assosiation. <key> : TGFID from AGILE 
<value>: <list> filenames in WWLLN directory that have a match in hour,minute and second. Have not adjusted for light 
travel time yet. When the light travel time have been calculated it's necessary to update the dictionary (+/- 500us error) 

"""

TGFIDs = []
possible_matches_dict = {}
WWLLN_candidates = {}         # <WWLLN time> : great circle distance


# A nested dictionary for the matches between TGF and lightning
match = {"TGFID": [], "Lon sat": [], "Lat sat": [], "Alt sat": [],
                                "Lon W": [], "Lat W": [], "W dist (km)": [], 
                                "Propagation time (s)": [], "Delta t (us)": []} 


def TGF_search(satelite_df):
    ### ALL FORMATS ARE CORRECT ###
    
    """ This function goes through each row in the AGILE dataset and updating the time_match_dict for each iteration.
    
    Several variables will then be passed to the WWLLN_filter function. That function will updates the time_match_dict
    """ 
    for index, row in satelite_df.iterrows():
        TGFID = row['New TGFID']
        TGFID_list = row['New TGFID'].split('.')
        year, month, day = "20" + TGFID_list[0][:2], TGFID_list[0][2:4], TGFID_list[0][4:]
        hour, minute, second, microsecond = TGFID_list[1][:2], TGFID_list[1][2:4], TGFID_list[1][-2:], TGFID_list[2]
        date = year + month + day
        time = hour + minute + second
        
        # Passes these variables forward to the next function, which finds the corresponding WWLLN file
        
        WWLLN_filter(lightning_data, date, time, TGFID)
        
        date_UTC_AGILE = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), int(microsecond))
        ### Maybe pass the date_UTC_AGILE object forward? ###
        
    return possible_matches_dict

def WWLLN_filter(filename, date, time, TGFID):
    ### RETURNS THE CORRECT FORMAT ###
    """
    This function updates the time_match_dict. key: TGFID, value: list of WWLLN filenames to be further explored.
    
    filenames <str> : List of filenames in the WWLLN folder. (Defined earlier)
    date : date of the TGF
    time : time of the TGF
    new_TGFID : date and time put together, not including microseconds.
    
    """
    
    for WWLLN_filename in filename:
        if WWLLN_filename == "A" + date + "_" + time + ".txt": ### Have to be the exact time ###
            if date not in possible_matches_dict:
                possible_matches_dict[TGFID] = [WWLLN_filename]
            else:
                possible_matches_dict[TGFID].append(WWLLN_filename)

def haversine(lon_AGILE, lat_AGILE, lon_WWLLN, lat_WWLLN):
    
    ### Correct function ####
    """ Function for calculating the distance between two points on a sphere. I will use this function to find the 
    closest lightning stroke to the TGF. Return distance between TGF and lightning stroke in km """
    
    R = 6.3781370e6 # Earth's equatorial radius in meters
    dLat = radians(lat_WWLLN - lat_AGILE)
    dLon = radians(lon_WWLLN - lon_AGILE)
    lat1 = radians(lat_AGILE)
    lat2 = radians(lat_WWLLN)

    a = sin(dLat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2
    c = 2*asin(sqrt(a))
    S = (R * c) 

    return S

def law_of_cosines(S, alt_sat):
    ### Correct function ###
    
    R = 6.3781370e6 # Earth's equatorial radius in meters
    alt_light = 15000   # Approx. alt of lightning (15km)
    R_sat = R + (alt_sat * 1000)
    R_light = R + alt_light
    #print("R sat:",R_sat,"R light:", R_light)
    
    theta = S / R  ###Returns radians###
    #print(theta, "Theta")
    
    """Is the fuction below correct?"""
    d = sqrt((R_sat)**2 + (R_light)**2 - 2*(R_sat)*(R_light)*cos(theta))
    
    return d

def propagation_time(d):
    
    ### Correct function ###
    "Propagation time for photons travelling from lighning to satelite. t measured in s and c is given in km/s"
    
    """Convert to m/s"""
    
    c = 2.99792e8 # Speed of light in m/s
    t = (d / c)
    
    #t_prop = datetime.timedelta(seconds = t)
    
    return t # Returning a float, not a datetime object

def time_diff(t_sat, t_W, t_prop):
    ### Correct format ###
    """
    t_sat <str>: time of the TGF event
    t_W <str> time of the lightning event
    t_prop <float> propagation time in seconds between the two
    return <float> corrected timediff between TGF event and lightning event given in seconds
    """
    
    Sat_time = datetime.datetime.strptime(t_sat, "%H:%M:%S.%f")
    WWLLN_time = datetime.datetime.strptime(t_W, "%H:%M:%S.%f") # datetime object for WWLLN time
    
    diff = Sat_time -  WWLLN_time
    y = divmod(diff.total_seconds(), 60) # Returns the remainder in a <tuple> : (seconds, microseconds)
    time_diff_sat_WWLLN = y[0] + y[1] #timedifferance between the t_sat and t_WWLLN
    delta_t = (time_diff_sat_WWLLN - t_prop)
    
    return delta_t

def great_filter(filename, lon_sat, lat_sat, alt_sat, TGFID, t_sat):
    ### Correct format ###
    
    """ Checking if the distance between the satelite and WWWLLN match are within 1000km.
    
    If the condition is True, then the calculation of distance (law of cosines) and delta t is executed.
    
    return: tuple or nested dictionary key: TGFID, value_0: SATELITE information, value_1: WWLLN information, 
    value_2: measurements (S,d,delta_t)
    """
    
    df = pd.read_csv("./WWLLN_data/" + filename, names=["Date", "Time", "Lat", "Lon", 'Resid', 'Nstn'])
       
    #Now iterate through the df and filter out the ones that are over 1000 km.
    for index, row in df.iterrows():
        lat_WWLLN = row["Lat"]
        lon_WWLLN = row["Lon"]
        date_WWLLN = [row["Date"], row["Time"]]
        t_WWLLN = row["Time"]
        
        'From this point forward i can do the calculations of S, d and delta t'
        S = haversine(lon_sat, lat_sat, lon_WWLLN, lat_WWLLN) #Great circle distance. Correct format
        if (S / 1000) <= 1000:                          # Checking that the distance is less than or equal to 1000 km.
            d = law_of_cosines(S,alt_sat) #Returns distance in m
            t_prop = propagation_time(d)  # Returns propagation time in seconds as a float
            delta_t = time_diff(t_sat, t_WWLLN, t_prop) # Returns the corrected time differance between
                    # satelite and lightning in seconds
            
            if -500 <= (delta_t * 1e6) <= 500:
                ### Append the relevant values (incl. propagation time) to a dictionary ###
                match["TGFID"].append(TGFID)
                match["Lon sat"].append(lon_sat)
                match["Lat sat"].append(lat_sat)
                match["Alt sat"].append(alt_sat)
                match["Lon W"].append(lon_WWLLN)
                match["Lat W"].append(lat_WWLLN)
                match["W dist (km)"].append(S/1000)
                match["Propagation time (s)"].append(t_prop)
                match["Delta t (us)"].append(delta_t*1e6)
                

                ##### Confirmation between the returned values and the AGILE dataframe #####
                
                
                """Finding out that there are lightning centered around the same time, only a coulpe of kms
                between them """  
        else:
            pass
        
    return match

    
def check_match(dictionary):
    """ There are often several lightning ass. for each TGF. This function will choose the lightning that is
    within 1000km to the satelite. (great circle distance)"""
    match = {}
    
    for TGFID, filename in dictionary.items():
        
        sat_row = AGILE_data.loc[AGILE_data["New TGFID"] == TGFID] #TGFID have whitespace before and after characters
        lon_sat = sat_row["Lon"].iloc[0]
        lat_sat = sat_row["Lat"].iloc[0]
        alt_sat = sat_row["Alt"].iloc[0]
        t_sat = sat_row["Date(UTC)"].iloc[0][12:20] + sat_row["New TGFID"].iloc[0][-5:]
        
        
        ### PASSING IN TGFID AS PARAMETER INTO THE DISTANCE_REDUCTION FUNCTION ###
        
        matches = great_filter(filename[0], lon_sat, lat_sat, alt_sat, TGFID, t_sat) #Return a dict of matches
        
    return matches
                
                
def main():
    first_candidates = TGF_search(AGILE_data) # TGF_search returns a dictionary <TGFID> : [WWLLN files]
    "Returns a dictionary of 113 possible files to be explored"
    
    matches = check_match(first_candidates) #Check_match returns a dictionary <WWLLN time> : S (great circle)
    
    match_df = pd.DataFrame.from_dict(matches)
    
    return match_df

"""Finding out that the same TGF have corrolation to several lightning incidences that are centered around an area.
explaination: lightning from the same thundercloud that occur around the same time
Finding 148 TGF matches from 178 TGFs. On the other hand, many of the matches are from the same thundercloud. This is
what i'm going to improve. If there are several matches from the same thundercload, they need to be filtered, leaving
only one match. """
    
main()


# In[269]:


matches_df = main()


# In[270]:


"Need to find a way around the 'not found' problem. Maybe try, except will work"

i = len(matches_df.index) - 1
print(i)
while i >= 0:
    if matches_df.iloc[i - 1][1] == matches_df.iloc[i][1]:
        "Dropping the row that have the highest delta t"
        print("match")
        if abs(matches_df.iloc[i][8]) > abs(matches_df.iloc[i - 1][8]):
            print("bigger: ", matches_df.iloc[i][8])
            matches_df = matches_df.drop(i)
            i -= 1
        elif abs(matches_df.iloc[i][8]) < abs(matches_df.iloc[i - 1][8]):
            print("bigger: ", matches_df.iloc[i][8])
            matches_df = matches_df.drop(i-1)
            i -= 1
    else:
        i -= 1
        print("not match")


# In[242]:


plt.subplots(figsize=(12,8))
plt.hist(AGILE_data["Lon"], bins = 'auto', edgecolor= "black", histtype="step", label='AGILE')
plt.hist(matches_df['Lon sat'], bins = 'auto', edgecolor= "red", histtype='step', hatch='//', label='AGILE/WWLLN Match')
plt.legend(loc='upper left')
plt.show()


# In[ ]:





# In[228]:


df_A = pd.read_csv("AGILE_match.csv")
df_A.head()


# In[ ]:





# ## Creating a scatter plot of the match_df

# In[ ]:





# In[ ]:




ax = plt.axes(projection=ccrs.PlateCarree())
ax.stock_img()


# In[ ]:





# In[ ]:





# In[ ]:





# ## Double checking that the values from the script are correct

# In[197]:


AGILE_data.head()


# In[ ]:





# In[ ]:





# ## Checking if the functions returns the correct distance and delta_t between AGILE and WWLLN

# In[125]:


def main():
    lon_AGILE = 115.76
    lat_AGILE = 0.62
    alt_sat   = 463.1
    lon_WWLLN = 112.63
    lat_WWLLN = -1.59
    
    S = haversine(lon_AGILE, lat_AGILE, lon_WWLLN, lat_WWLLN)
    d = law_of_cosines(S, alt_sat)
    t_prop = propagation_time(d)
    delta_t = time_diff('11:28:02.7849','11:28:02.782773', t_prop)
    
    return print("Haversine:",S/1000, "Fasit:", 429, "Law of cosines:", d/1000, "delta_t:", delta_t)
    

def haversine(lon_AGILE, lat_AGILE, lon_WWLLN, lat_WWLLN):
    """ Function for calculating the distance between two points on a sphere. I will use this function to find the 
    closest lightning stroke to the TGF. Return distance between TGF and lightning stroke in km """
    
    R = 6.3781370e6 # Earth's equatorial radius in meters
    dLat = radians(lat_WWLLN - lat_AGILE)
    dLon = radians(lon_WWLLN - lon_AGILE)
    lat1 = radians(lat_AGILE)
    lat2 = radians(lat_WWLLN)

    a = sin(dLat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2
    c = 2*asin(sqrt(a))
    S = (R * c) 

    return S

def law_of_cosines(S, alt_sat):
    
    R = 6.3781370e6 # Earth's equatorial radius in meters
    alt_light = 15000   # Approx. alt of lightning (15km)
    R_sat = R + (alt_sat * 1000)
    R_light = R + alt_light
    #print("R sat:",R_sat,"R light:", R_light)
    
    theta = S / R  ###Returns radians###
    #print(theta, "Theta")
    
    """Is the fuction below correct?"""
    d = sqrt((R_sat)**2 + (R_light)**2 - 2*(R_sat)*(R_light)*cos(theta))
    
    return d

def propagation_time(d):
    "Propagation time for photons travelling from lighning to satelite. t measured in s and c is given in km/s"
    
    """Convert to m/s"""
    
    c = 2.99792e8 # Speed of light in m/s
    t = (d / c)
    print('t prop:', t)
    
    #t_prop = datetime.timedelta(seconds = t)
    
    return t # Returning a float, not a datetime object

def time_diff(t_sat, t_W, t_prop):
    """
    t_sat <str>: time of the TGF event
    t_W <str> time of the lightning event
    t_prop <float> propagation time in seconds between the two
    return <float> corrected timediff between TGF event and lightning event given in seconds
    """
    
    Sat_time = datetime.datetime.strptime(t_sat, "%H:%M:%S.%f")
    WWLLN_time = datetime.datetime.strptime(t_W, "%H:%M:%S.%f") # datetime object for WWLLN time
    
    
    #print("WWLLN_time:", WWLLN_time,"Sat_time:", Sat_time,"t_prop:", t_prop)
    #print("Sat_time:",Sat_time,"corrected_sat_time:", corrected_sat_time)
    
    diff = Sat_time -  WWLLN_time
    y = divmod(diff.total_seconds(), 60) # Returns the remainder in a <tuple> : (seconds, microseconds)
    time_diff_sat_WWLLN = y[0] + y[1] #timedifferance between the t_sat and t_WWLLN
    delta_t = (time_diff_sat_WWLLN - t_prop)

    
    """date_TGF = [] # A list that stores the date and time of TGF
    
    t_sat = t_sat.strip()
    t_sat = t_sat.split(".")

    date_TGF.extend(["20" + t_sat[0][:2], t_sat[0][2:4], t_sat[0][4:6], #date
                        t_sat[1][:2], t_sat[1][2:4], t_sat[1][4:], #time
                        t_sat[2]]) #microseconds
                                  
    sat_date = datetime.datetime(int(date_TGF[0]), int(date_TGF[1]), int(date_TGF[2]), int(date_TGF[3]), int(date_TGF[4]), int(date_TGF[5]), int(date_TGF[6]))
    
    WWLLN_date = datetime.datetime(int(t_W[0][:4]),int(t_W[0][5:7]),int(t_W[0][8:]),int(t_W[1][:2]),int(t_W[1][3:5]),int(t_W[1][6:8]), int(t_W[1][9:]))
    
    TGF_corrected = sat_date - t_prop
    delta_t = TGF_corrected - WWLLN_date # This is the delta t between lightning and TGF. Needs to be +- 500 us"""
    
    return delta_t





main()


# ## Old functions

# In[ ]:


def haversine(lon_AGILE, lat_AGILE, lon_WWLLN, lat_WWLLN):
    """ Function for calculating the distance between two points on a sphere. I will use this function to find the 
    closest lightning stroke to the TGF. Return distance between TGF and lightning stroke in km """
    
    R = 6.371009e6 # Earth radius in meters
    dLat = radians(lat_WWLLN - lat_AGILE)
    dLon = radians(lon_WWLLN - lon_AGILE)
    lat1 = radians(lat_AGILE)
    lat2 = radians(lat_WWLLN)

    a = sin(dLat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2
    c = 2*asin(sqrt(a))
    S = (R * c)

    return np.round(S)

def law_of_cosines(S, alt_sat):
    """ WWLLN_candidates <dict>: Key is the WWLLN time, value is the great circle distance (S) between TGF and lightning
    
    Parameters for calculating the propagation time between lightning and TGF (calculating in km and seconds): 
    
    1) S = Great circle distance
    2) R = 6.371009e6 # Earth's avrage radius in m
    3) c = 2.99792e8 # Speed of light in m/s
    4) alt_sat =  x (given in km -> need to convert to m)
    5) alt_light = 15 # Altitude of lightning phenomena in km
    
    """   
    ### Need to convert from degrees to radians ###
    ### Is it correct to add the radius of earth to the two variables ### 
    R = 6.371009e6 # Radius of earth (km)
    alt_light = 15000   # Approx. alt of lightning (km)
    alt_sat = alt_sat * 1000  # Converting from km to m
    
    theta = S / R  ###Returns radians### 
    
    """Is the fuction below correct?"""
    d = sqrt((R + alt_sat)**2 + (R + alt_light)**2 - 2*(R + alt_sat)*(R + alt_light)*cos(theta))
    
    return (np.round(d))

def propagation_time(d):
    "Propagation time for photons travelling from lighning to satelite. t measured in s and c is given in km/s"
    
    """Convert to m/s"""
    
    c = 2.99792e8 # Speed of light in m/s
    t = (d / c)
    
    t_prop = datetime.timedelta(seconds = t)
    
    return t_prop

def time_diff(t_sat, t_W, t_prop):
    """t_sat <str>: full date of the TGF event
    t_W <list>: a list of the full date of the lightning stroke
    t_prop <object>: propagation time from lightning location to satelite. Given as datetime object (see function:
    propagation_time). """
    
    date_TGF = [] # A list that stores the date and time of TGF
    
    t_sat = t_sat.strip()
    t_sat = t_sat.split(".")

    date_TGF.extend(["20" + t_sat[0][:2], t_sat[0][2:4], t_sat[0][4:6], #date
                        t_sat[1][:2], t_sat[1][2:4], t_sat[1][4:], #time
                        t_sat[2]]) #microseconds
                                  
    sat_date = datetime.datetime(int(date_TGF[0]), int(date_TGF[1]), int(date_TGF[2]), int(date_TGF[3]), int(date_TGF[4]), int(date_TGF[5]), int(date_TGF[6]))
    
    WWLLN_date = datetime.datetime(int(t_W[0][:4]),int(t_W[0][5:7]),int(t_W[0][8:]),int(t_W[1][:2]),int(t_W[1][3:5]),int(t_W[1][6:8]), int(t_W[1][9:]))
    
    TGF_corrected = sat_date - t_prop
    delta_t = TGF_corrected - WWLLN_date # This is the delta t between lightning and TGF. Needs to be +- 500 us
    
    return delta_t


# In[235]:


""" Need to figure out how to check that the TGF are within +- 500 us. """

t_0 = '180311.114515.931552'

t_W = ['2018/02/04', '05:47:28.016891'] # Need to fix this into the right format. Assuming it's passed in as it is

t_prop = 1957.5634326049794

def time_diff(t_sat, t_W, t_prop):
    """t_sat <str>: full date of the TGF event
    t_W <list>: a list of the full date of the lightning stroke
    t_propagation <int>: propagation time from lightning location to satelite. Given in microseconds (see function:
    propagation_time). """
    
    date_TGF = [] # A list that stores the date and time of TGF
    
    """Maybe split the string instead"""
    
    t_sat = t_sat.split(".")
    print(t_sat)
        
    date_TGF.extend(["20" + t_sat[0][:2], t_sat[0][2:4], t_sat[0][4:6], #date
                        t_sat[1][:2], t_sat[1][2:4], t_sat[1][4:], #time
                        t_sat[2]]) #microseconds
                                  
    sat_date = datetime.datetime(int(date_TGF[0]), int(date_TGF[1]), int(date_TGF[2]), int(date_TGF[3]),
                                 int(date_TGF[4]), int(date_TGF[5]), int(date_TGF[6]))
    
    WWLLN_date = datetime.datetime(int(t_W[0][:4]),int(t_W[0][5:7]),int(t_W[0][8:]),int(t_W[1][:2]),
                                int(t_W[1][3:5]),int(t_W[1][6:8]), int(t_W[1][9:]))
    
    t_propagation = datetime.timedelta(microseconds = t_prop) # propagation time
    
    TGF_corrected = sat_date - t_propagation
    
    delta_t = TGF_corrected - WWLLN_date # This is the delta t between lightning and TGF. Needs to be +- 500 us
    
    return delta_t

print(time_diff(t_0, t_W, t_prop))


# In[104]:





# In[ ]:




