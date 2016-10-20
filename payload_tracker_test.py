# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 16:12:59 2016

@author: Caroline
"""

import matplotlib.pyplot as pyp
import numpy as np
import math as m
import numpy as np

# function to return distance using latitude, longitude, and altitude (in meters)                                                                                                                                                
def ArcDistCalc(alt_Pulser, DLat_Pulser, DLon_Pulser, alt_ANITA, DLat_ANITA, DLon_ANITA):
    try:
        # lat, lon, alt equaling zero or "None" means that no data was received (they are not usable coordinates)                                                                                                    
        if DLat_Pulser == 0 or DLon_Pulser == 0 or alt_ANITA == 'grnd' or DLat_ANITA == '' or DLon_ANITA == '' :
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

            # convert the short approx of distance from radians to meters
                # use the same value of R_earth as next function ... assume Earth radius = 66378137 m                                                                                                        
            dist = 66378137 * Dshort
            
            # find elevation angle --> inverse tangent of alt/dist --> convert from rad to degrees
            elev_Angle = m.degrees(m.atan((alt_ANITA - alt_Pulser)/dist))
            
            # pythagorean theorem to estimate arc distance (m)                                                                                                                              
            arcDist = np.sqrt(dist**2 + (alt_ANITA - alt_Pulser)**2)
        return (arcDist, elev_Angle)
    except:
        print("ERROR")
        pass
    
# function to convert degrees gps coords to (x,y,z)ECEF ... altitude in meters
def coord_degToECEF(latitude, longitude, alt):
    
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
 
    return (x,y,z)
   
    
# function for azimuth angle using ECEF x y z coordinates
def azimuthAngle(x, y, z, xG, yG, zG):
    # sine and cosine of azimuth angle ... (xG, yG, zG --> "ground")
    (dx, dy, dz) = (x-xG, y-yG, z-zG)
    print "delta coords = ", (dx, dy, dz)
    # north "level" vector (projection of north vector)
    (xN, yN, zN) = (-zG*xG, -zG*yG, xG**2+yG**2)
    print "north pointing vector = ", (xN, yN, zN)
    # east "level" vector (projection of east vector)
    (xE, yE, zE) = (-yG, xG, 0)
    print "east pointing vector = ", (xE, yE, zE)

    cosAz = (xN*dx+yN*dy+zN*dz)/m.sqrt(zN*(x**2+y**2+z**2)*(dx**2+dy**2+dz**2))
    print "cosAz = ", cosAz
    sinAz = (xE*dx+yE*dy)/m.sqrt(zN*(dx**2+dy**2+dz**2))
    print "sinAz = ", sinAz

    # return azimuth as inverse tangent of sine and cosine
    AzAngle = m.degrees(m.atan(sinAz/cosAz))

    return AzAngle


#####%%%%%%%%%%%
if __name__ == '__main__':
        
        # convert coordinates to degrees (DDD.ddddd) ... altitude in m
        # alt_p = 0.
        # DLat_p = -79 - 27.937/60.
        # DLon_p = -112. - 6.744/60.

        # alt_a = 37.00
        # DLat_a = -75 -  27.937/60.
        # DLon_a = -112. - 6.744/60.

        # TESTER WITH EXAMPLE NUMBERS FROM "http://gis.stackexchange.com/questions/58923/calculate-view-angle"
        alt_p = 4000
        DLat_p = 39.
        DLon_p = -75.

        alt_a = 12000
        DLat_a = 39.
        DLon_a = -76.

        print "coordinates in degrees = ", (DLat_p, DLon_p, DLat_a, DLon_a)

        [dist_arc, angle_elevation] = ArcDistCalc(alt_p, DLat_p, DLon_p, alt_a, DLat_a, DLon_a)
        print "\n elevation angle (degrees) = ", angle_elevation

        # test input coordinates [-60, 60, 500]
        [x_a, y_a, z_a] = coord_degToECEF(DLat_a, DLon_a, alt_a)
        # print "\n ECEF coordinates of anita payload (meters) = ", (x_a,y_a,z_a)
        print "\n coords of jet = ", (x_a,y_a,z_a)

        [x_p, y_p, z_p] = coord_degToECEF(DLat_p, DLon_p, alt_p)
        print "coords of plane = ", (x_p, y_p, z_p)

        AzAngle = azimuthAngle(x_a, y_a, z_a, x_p, y_p, z_p)

        print "\n azimuth angle (degrees) = ", AzAngle

        print "positive value of azimuth angle (degrees)", 360 + AzAngle
        
        
        
