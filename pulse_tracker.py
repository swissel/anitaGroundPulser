
# pulse_tracker.py
# script to calculate elevation and azimuth angles from two GPS coordinates


import matplotlib.pyplot as pyp
import numpy as np
import math as m



# function to convert degrees gps coords to (x,y,z)ECEF ... altitude in meters
def coord_degToECEF(lat, lon, alt):
    # constant in meters
    a = 6378137
    # flattening coefficient 
    f = 1/298.257223563
    # convert to radians
    lat = m.radians(lat)
    lon = m.radians(lon)

    r = a/m.sqrt(1-(m.sin(lat)**2)*(f*(2-f)))
    x = (r+alt)*m.cos(lat)*m.cos(lon)
    y = (r+alt)*m.cos(lat)*m.sin(lon)
    z = ((r*(1-f*(2-f)))+alt)*m.sin(lat)

    return x,y,z

# function for elevation angle using ECEF x y z coordinates
def elevationAngle(x, y, z, xG, yG, zG):
    # sine and cosine of azimuth angle ... (xG, yG, zG --> "ground")
    (dx, dy, dz) = (x-xG, y-yG, z-zG)

    cosEl = (x*dx+y*dy+z*dz)/np.sqrt((x**2+y**2+z**2)*(dx**2+dy**2+dz**2))
    # elevation from the horizon
    horiz_elev = 90-m.acos(cosEl)

    return horiz_elev

# function for azimuth angle using ECEF x y z coordinates
def azimuthAngle(x, y, z, xG, yG, zG):
    # sine and cosine of azimuth angle ... (xG, yG, zG --> "ground")
    (dx, dy, dz) = (x-xG, y-yG, z-zG)
    # print "delta coords = ", (dx, dy, dz)

    # north "level" vector (projection of north vector)
    (xN, yN, zN) = (-zG*xG, -zG*yG, xG**2+yG**2)
    # print "north pointing vector = ", (xN, yN, zN)
    
    # east "level" vector (projection of east vector)
    (xE, yE, zE) = (-yG, xG, 0)
    # print "east pointing vector = ", (xE, yE, zE)

    cosAz = (xN*dx+yN*dy+zN*dz)/m.sqrt(zN*(x**2+y**2+z**2)*(dx**2+dy**2+dz**2))
    # print "cosAz = ", cosAz
    sinAz = (xE*dx+yE*dy)/m.sqrt(zN*(dx**2+dy**2+dz**2))
    # print "sinAz = ", sinAz

    # return azimuth as inverse tangent of sine and cosine
    AzAngle = m.degrees(m.atan(sinAz/cosAz))

    return AzAngle


#####%%%%%%%%%%%
if __name__ == '__main__':
        
        # TESTER WITH EXAMPLE NUMBERS FROM "http://gis.stackexchange.com/questions/58923/calculate-view-angle"
        # convert coordinates to degrees (DDD.ddddd) ... altitude in m
        # alt_p = 0.
        # DLat_p = -79 - 27.937/60.
        # DLon_p = -112. - 6.744/60.
        # alt_a = 37.00
        # DLat_a = -75 -  27.937/60.
        # DLon_a = -112. - 6.744/60.



        alt_p = 0.019812
        print 'Pulser altitude hard-coded as ', alt_p

        DLat_p = -77.853836167
        print 'Pulser latitude hard-coded as ', DLat_p

        DLon_p = 167.202818
        print 'Pulser longitude hard-coded as ', DLon_p

        print 'ENTER LAT AND LON AS DECIMAL DEGREES'

        alt_a = input('Enter the altitude of ANITA: ')

        DLat_a = input('Enter the latitude of ANITA: ')

        DLon_a = input('Enter the longitude of ANITA: ')



        print "coordinates in degrees (pulser lat, pulser lon, ANITA lat, ANITA lon)= ", (DLat_p, DLon_p, DLat_a, DLon_a)

        # test input coordinates [-60, 60, 500]
        [x_a, y_a, z_a] = coord_degToECEF(DLat_a, DLon_a, alt_a)
        print "ECEF coordinates of ANITA payload (meters) = ", (x_a,y_a,z_a)
        [x_p, y_p, z_p] = coord_degToECEF(DLat_p, DLon_p, alt_p)
        print "ECEF coordinates of pulser (meters) = ", (x_p, y_p, z_p)

        AzAngle = azimuthAngle(x_a, y_a, z_a, x_p, y_p, z_p)
        angle_elevation = elevationAngle(x_a, y_a, z_a, x_p, y_p, z_p)

        print "\n elevation angle (degrees) = ", angle_elevation
        if AzAngle < 0:
            print "azimuth angle (degrees) = ", abs(AzAngle), " S .... (%s)" %AzAngle 
        else:
            print "azimuth angle (degrees) = %s N" %AzAngle
