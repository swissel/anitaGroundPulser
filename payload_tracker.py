
# pulse_tracker.py
# script to calculate elevation and azimuth angles from two GPS coordinates


import matplotlib.pyplot as pyp
import numpy as np
import math as m
import argparse

# function to convert degrees gps coords to (x,y,z)ECEF ... altitude in meters
def coord_degToECEF(lat, lon, alt):
    # constant in meters
    a = 6378137.
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

def distance(x, y, z, xG, yG, zG):
    (dx, dy, dz) = (x-xG, y-yG, z-zG)	
    return np.sqrt( dx**2 + dy**2 + dz**2)

# function for elevation angle using ECEF x y z coordinates
def elevationAngle(x, y, z, xG, yG, zG):
    # elevation angle from the ground to the payload

    # dv is vector from the pulser vG = (xG, yG, zG) to the payload v = (x, y, z)
    # v and vG are in ECEF (earth-centered coordinates)
    (dx, dy, dz) = (x-xG, y-yG, z-zG)
    #print "elevationAngle:\n\t (dx, dy, dz): ", (dx,dy,dz)
    
    # cos zenith angle is the inner product of 
    # unit vectors of vG and dv 
    cosZe = (xG*dx + yG*dy + zG*dz) / np.sqrt( (xG**2+yG**2+zG**2)*(dx**2+dy**2+dz**2) )

    # elevation from the horizon
    horiz_elev = 90. - m.acos(cosZe)*180./np.pi
    #print "\televation angle ", horiz_elev
    return horiz_elev

# function for azimuth angle using ECEF x y z coordinates
def azimuthAngle(x, y, z, xG, yG, zG):
    # sine and cosine of azimuth angle ... (xG, yG, zG --> "ground")
    (dx, dy, dz) = (x-xG, y-yG, z-zG)
    #print "delta coords = ", (dx, dy, dz)

    # north "level" vector (projection of north vector)
    (xN, yN, zN) = (-zG*xG, -zG*yG, xG**2+yG**2)
    #print "north pointing vector = ", (xN, yN, zN)
    
    # east "level" vector (projection of east vector)
    (xE, yE, zE) = (-yG, xG, 0)
    #print "east pointing vector = ", (xE, yE, zE)

    cosAz = (xN*dx+yN*dy+zN*dz)/m.sqrt(zN*(x**2+y**2+z**2)*(dx**2+dy**2+dz**2))
    # print "cosAz = ", cosAz
    sinAz = (xE*dx+yE*dy)/m.sqrt(zN*(dx**2+dy**2+dz**2))
    #print "sinAz = ", sinAz

    # return azimuth as inverse tangent of sine and cosine
    AzAngle = m.degrees(m.atan(sinAz/cosAz))

    return AzAngle


#####%%%%%%%%%%%
if __name__ == '__main__':
        
        # TESTER WITH EXAMPLE NUMBERS FROM "http://gis.stackexchange.com/questions/58923/calculate-view-angle"
        # convert coordinates to degrees (DDD.ddddd) ... altitude in m
        # alt_p = 4000.
        # DLat_p = 39 
        # DLon_p = -75
        # alt_a = 12000
        # DLat_a = 39
        # DLon_a = -76
	# --> Gives elevation angle from horizontal : 4.88 degrees
	# 	    azimuthal angle in degrees clockwise from from N : 270.331 / -89.67

	opts = argparse.ArgumentParser(description='process commandline options')
	
     	# default pulser location is LDB near McMurdo Station
     	# you can change the default as you wish.
     	opts.add_argument('anita_lat', help='Latitude of ANITA payload in +/- decimal degrees E.', type=float)
     	opts.add_argument('anita_lon', help='Longitude of ANITA payload in +/- decimal degress N.', type=float)
     	opts.add_argument('anita_alt', help='Altitude of ANITA payload in m above sea level.', type=float)
        opts.add_argument('pulser_lat', help='Latitude of pulser in +/- decimal degrees E.', nargs='?', type=float, default=-77.853836167)
     	opts.add_argument('pulser_lon', help='Longitude of pulser in +/- decimal degrees N.',nargs='?', type=float, default=167.202818)
     	opts.add_argument('pulser_alt', help='Altitude of pulser in m above sea level.', nargs='?', type=float, default=0.019812)

     	args = opts.parse_args()

	alt_p = args.pulser_alt
	DLat_p = args.pulser_lat
	DLon_p = args.pulser_lon

	alt_a = args.anita_alt
	DLat_a = args.anita_lat
	DLon_a = args.anita_lon

	# Print out the geometery
	print "-------------------------------------------------"
	print "Pulser location:"
	print "-------------------------------------------------"
	print "\tLatitude : ", DLat_p, " degrees"
	print "\tLongitude : ", DLon_p, " degrees"
	print "\tAltitude : ", alt_p, " m"
	print "-------------------------------------------------"
	print "\tANITA location:"
	print "\tLatitude : ", DLat_a, " degrees"
	print "\tLongitude : ", DLon_a, " degrees"
	print "\tAltitude : ", alt_a, " m"
	print "-------------------------------------------------"

        #print "coordinates in degrees (pulser lat, pulser lon, ANITA lat, ANITA lon)= ", (DLat_p, DLon_p, DLat_a, DLon_a)

        [x_a, y_a, z_a] = coord_degToECEF(DLat_a, DLon_a, alt_a)
        #print "ECEF coordinates of ANITA payload (meters) = ", (x_a,y_a,z_a)
        [x_p, y_p, z_p] = coord_degToECEF(DLat_p, DLon_p, alt_p)
        #print "ECEF coordinates of pulser (meters) = ", (x_p, y_p, z_p)
	#print "(dx, dy, dz): ", (x_a-x_p, y_a-y_p, z_a-z_p)
	
	print "Distance from ANITA to Pulser : ", distance(x_a, y_a, z_a, x_p, y_p, z_p)/1000., " km."
	print "-------------------------------------------------"

        AzAngle = azimuthAngle(x_a, y_a, z_a, x_p, y_p, z_p)
        angle_elevation = elevationAngle(x_a, y_a, z_a, x_p, y_p, z_p)
	
	print "Point pulser antenna in the following direction:"
	print "-------------------------------------------------"
        print "\tElevation angle:\t%2.2f degrees from horizontal" % (angle_elevation)
        if AzAngle < 0:
            print "\tAzimuth angle:\t%3.2f / %3.2f  degrees CW from N" % (360. + AzAngle, AzAngle)
        else:
            print "\tAzimuth angle:\t%3.2f degrees CW from N" % (AzAngle)
	print "-------------------------------------------------"
