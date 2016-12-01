# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 16:12:59 2016

@author: Caroline
"""

import matplotlib.pyplot as pyp
import numpy as np
import math as m
import numpy as np
import argparse

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

def test():
     # convert coordinates to degrees (DDD.ddddd) ... altitude in m
     # alt_p = 0.
     # DLat_p = -79 - 27.937/60.
     # DLon_p = -112. - 6.744/60.

     # alt_a = 37.00
     # DLat_a = -75 -  27.937/60.
     # DLon_a = -112. - 6.744/60.

     # TESTER WITH EXAMPLE NUMBERS FROM "http://gis.stackexchange.com/questions/58923/calculate-view-angle"

     #For accurate calculations, convert (lat, lon, elevation) directly to earth-centered (x,y,z). (If you don't do this, you need to retain additional information about the local normal ["up"] directions in order to compute angles accurately at nonzero elevations.)

     #Elevation

     #Given two points (x,y,z) and (x',y',z') in an earth-centered coordinate system, the vector from the first to the second is (dx,dy,dz) = (x'-x, y'-y, z'-z), whence the cosine of the angle made to the normal at (x,y,z) is the inner product of the unit length versions of those vectors:

     #Cos(elevation) = (x*dx + y*dy + z*dz) / Sqrt((x^2+y^2+z^2)*(dx^2+dy^2+dz^2))
     #Obtain its principal inverse cosine. Subtract this from 90 degrees if you want the angle of view relative to a nominal horizon. This is the "elevation."

     #Azimuth

     #A similar calculation obtains the local direction of view ("azimuth"). We need a level vector (u,v,w) that points due north. One such vector at the location (x,y,z) is (-z*x, -z*y, x^2+y^2). (The inner product of these two vectors is zero, proving it is level. Its projection onto the Equatorial plane is proportional to (-x,-y) which points directly inward, making it the projection of a north-pointing vector. These two calculations confirm that this is indeed the desired vector). Therefore

     #Cos(azimuth) = (-z*x*dx - z*y*dy + (x^2+y^2)*dz) / Sqrt((x^2+y^2)(x^2+y^2+z^2)(dx^2+dy^2+dz^2))
     #We also need the sine of the azimuth, which is similarly obtained once we know a vector pointing due East (locally). Such a vector is (-y, x, 0), because it clearly is perpendicular to (x,y,z) (the up direction) and the northern direction. Therefore

     #Sin(azimuth) = (-y*dx + x*dy) / Sqrt((x^2+y^2)(dx^2+dy^2+dz^2))
     #These values enable us to recover the azimuth as the inverse tangent of the cosine and sine.

     # EXAMPLE
     # A pilot in an airplane flying flying west at 4000 meters, located at (lat, lon) = (39, -75), sees a jet far ahead located at (39, -76) and flying at 12000 meters. What is are the angles of view (relative to the level direction at the pilot's location)?

     #The XYZ coordinates of the airplanes are (x,y,z) = (1285410, -4797210, 3994830) and (x',y',z') = (1202990, -4824940, 3999870), respectively (in the ITRF00 datum, which uses the GRS80 ellipsoid). The pilot's view vector therefore is (dx,dy,dz) = (-82404.5, -27735.3, 5034.56). Applying the formula gives the cosine of the view angle as 0.0850706. Its inverse cosine is 85.1199 degrees, whence the elevation is 4.88009 degrees: the pilot is looking up by that much.

     #A north-pointing level vector is (-5.13499, 19.1641, 24.6655) (times 10^12) and an east-pointing level vector is (4.79721, 1.28541, 0) (times 10^6). Therefore, applying the last two formulas, the cosine of the azimuth equals 0.00575112 and its sine equals -0.996358. The ArcTangent function tells us the angle for the direction (0.00575112, -0.996358) is 270.331 degrees: almost due west. (It's not exactly west because the two planes lie on the same circle of latitude, which is curving toward the North pole: see Why is the 'straight line' path across continent so curved? for an extended explanation.)

    #By the way, this example confirms we got the orientation correct for the azimuth calculation: although it was clear the east-pointing vector was orthogonal to the other two vectors, until now it was not plain that it truly points east and not west.
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


#####%%%%%%%%%%%
if __name__ == '__main__':
        
     opts = argparse.ArgumentParser(description='process commandline options')
	
     # default pulser location is LDB near McMurdo Station
     # you can change the default as you wish.
     opts.add_argument('pulser_lat', help='Latitude of pulser in +/- decimal degrees E.', nargs='?', type=float, default=-77.853836167)
     opts.add_argument('pulser_lon', help='Longitude of pulser in +/- decimal degrees N.',nargs='?', type=float, default=167.202818)
     opts.add_argument('pulser_alt', help='Altitude of pulser in m above sea level.', nargs='?', type=float, default=0.019812)
     opts.add_argument('lat', help='Latitude of ANITA payload in +/- decimal degrees E.', nargs=1, type=float)
     opts.add_argument('lon', help='Longitude of ANITA payload in +/- decimal degress N.', nargs=1, type=float)
     opts.add_argument('alt', help='Altitude of ANITA payload in m above sea level.', nargs=1, type=float)
        
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


   
