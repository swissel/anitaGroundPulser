# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 16:12:59 2016

@author: Caroline
"""

import matplotlib.pyplot as pyp
import numpy as np
import math as m
import numpy as np

# function to return distance using latitude, longitude, and altitude                                                                                                                                                
def ArcDistCalc(alt_Pulser, DLat_Pulser, DLon_Pulser, alt_ANITA, DLat_ANITA, DLon_ANITA):
    try:
        # lat, lon, alt equaling zero or "None" means that no data was received (they are not usable coordinates)                                                                                                    
        if DLat_Pulser == 0 or DLon_Pulser == 0 or alt_Pulser or alt_ANITA == 'grnd' or DLat_ANITA == '' or DLon_ANITA == '' :
            print("Please use better input values.")
            arcDist = None
        else:
            # this function assumes that lat and lon are in degrees (DDD.ddddddd) ... Convert to radians                                                                                                             
            lat1 = (np.pi/180)*DLat_Pulser
            lat2 = (np.pi/180)*DLat_ANITA
            lon1 = (np.pi/180)*DLon_Pulser
            lon2 = (np.pi/180)*DLon_ANITA

            # use the "great circle distance" formula with reduction in short distance errors found at:                                                                                                              
                # "http://williams.best.vwh.net/avform.htm#Dist"                                                                                                                                                     
            Dshort = 2.*np.arcsin(np.sqrt((np.sin((lat1-lat2)/2.))**2 ))
            Dshort = Dshort + np.cos(lat1)*np.cos(lat2)*(np.sin((lon1-lon2)/2.))**2

            # convert the short approx of distance from radians to km ... assume Earth radius = 6371 km                                                                                                              
            dist = 6371.0 * Dshort
            # pythagorean theorem to estimate arc distance (mi) using dist_mi and alt_mi                                                                                                                             
            arcDist = np.sqrt(dist**2 + (alt_ANITA - alt_Pulser)**2)
        return arcDist
    except:
        print("ERROR")
        pass
    
    
def coords_to_ECEF(latitude, longitude, alt):
    
    lat = m.radians(latitude)   # convert from degrees into radians
    lon = m.radians(longitude)
    
    a = 6378137. #mean Earth radius in m
    f = 1./298.257223563 #flattening coefficient


    p1 = f*(2-f)
    p2 = m.sin(lat)**2
    d = m.sqrt(1-(p2*p1))
    rn = a/d
     
    x = (rn + alt)*m.cos(lat)*m.cos(lon)    # meters
    y = (rn + alt)*m.cos(lat)*m.sin(lon)
    z = ((rn*(1-f*(2.-f)))+alt)*m.sin(lat)
    
    
    
    return x,y,z
    


#####%%%%%%%%%%%                                                                                                                                                                                                     
if __name__ == '__main__':
        alt_p = 0.
        DLat_p = -79 - 27.937/60.
        DLon_p = -112. - 6.744/60.

        alt_a = 37.00
        DLat_a = -75 -  27.937/60.
        DLon_a = -112. - 6.744/60.

        print (DLat_p, DLon_p, DLat_a, DLon_a)
        
        
        [x,y,z] = coords_to_ECEF(-60,60,500)
        print x,y,z
        




