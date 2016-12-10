
# pulse_tracker.py
# script to calculate elevation and azimuth angles from two GPS coordinates


#import matplotlib.pyplot as pyp
import numpy as np
import math as m
import argparse

# class to handle degreesDminutesSsecondsSN/S/E/W or decimal degrees
class HandleDegrees(argparse.Action):
       def __call__(self, parser, namespace, value, option_string=None):
	   try:
		# must be decimal degrees already
		setattr(namespace, self.dest, float(value)) 
	   except:
		   dind = value.find('D', 0, -1) 
		   if(  dind == -1):
			try:
				# must be decimal degrees already
				setattr(namespace, self.dest, float(value)) 
			except:
				print "degrees must be in either decimal degrees (20.5) or degDminD[secS]N/S/E/W (20D30M00SN or 20D30MN):", value
		   else:
			mind = value.find('M', 0, -1)
			sind = value.find('S', 0, -1)
			direction = value[-1]

			degrees = float(value[0:dind])
			minutes = float(value[dind+1:mind])
			seconds = 0.
			if( sind != -1 ):
				seconds = float(value[mind+1:sind])
			
			decdeg = degrees + minutes/60. + seconds/3600.
			if( direction == "S" or direction == "W"):
				decdeg *= -1.
			setattr(namespace, self.dest, decdeg) 


		

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

def greatcircle_distance(lat1, lon1, lat2, lon2):
    # using the haversine formula
    
    phi1    = np.deg2rad(lat1)
    phi2    = np.deg2rad(lat2)
    dphi    = np.deg2rad(lat2 - lat1)
    dlambda = np.deg2rad(lon2 - lon1)
    
    #R = 6378137.
    R = 6371e3

    a = pow(np.sin(dphi/2.),2.) + np.cos(phi1) * np.cos(phi2) * pow(np.sin(dlambda/2.),2.)
    c = 2*np.arctan2(np.sqrt(a), np.sqrt(1.-a))
    d = R * c
    return d

def distance(x, y, z, xG, yG, zG):
    # this gives the wrong distance compared with the haversine formula
    # do not understand why?
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
	###############################################################################
	# Check with Mt. Erebus from LDB should give answer that is due North and <10deg in elevation
	#
	# bash-3.2$ python payload_tracker.py -77.516666667 167.15 3794
	# Mt. Erebus lat  -77.516666667 lon  167.15 altitude  3794.0 m
	# Gives elevation angle of 5.66 deg. from horizontal / 358.06 / -1.94  degrees CW from N
	###############################################################################

	opts = argparse.ArgumentParser(description='process commandline options')
	
     	# default pulser location is LDB near McMurdo Station
     	# you can change the default as you wish.
     	opts.add_argument('anita_lat', help='Latitude of ANITA payload in +/- decimal degrees N. or DegDMinM[SecS]E/W', action=HandleDegrees, metavar='deg[DminM[secS]]')
     	opts.add_argument('anita_lon', help='Latitude of ANITA payload in +/- decimal degrees N. or DegDMinM[SecS]E/W', action=HandleDegrees, metavar='deg[DminM[secS]]')
     	opts.add_argument('anita_alt', help='Altitude of ANITA payload in ft above sea level.', type=float)
        # LDB
	#opts.add_argument('pulser_lat', help='Latitude of ANITA payload in +/- decimal degrees N. or DegDMinM[SecS]E/W', action=HandleDegrees, metavar='deg[DminM[secS]]', nargs='?', default=-77.853836167)
	#opts.add_argument('pulser_lat', help='Latitude of ANITA payload in +/- decimal degrees N. or DegDMinM[SecS]E/W', action=HandleDegrees, metavar='deg[DminM[secS]]', nargs='?',  default=167.202818)
	#opts.add_argument('pulser_lat', help='Latitude of ANITA payload in +/- decimal degrees N. or DegDMinM[SecS]E/W', action=HandleDegrees, metavar='deg[DminM[secS]]', nargs='?', type=float, default=0.019812)
	# WAIS
	opts.add_argument('pulser_lat', help='Latitude of ANITA payload in +/- decimal degrees N. or DegDMinM[SecS]E/W', action=HandleDegrees, metavar='deg[DminM[secS]]', nargs='?', default=-79.465616667)
     	opts.add_argument('pulser_lon', help='Latitude of ANITA payload in +/- decimal degrees N. or DegDMinM[SecS]E/W', action=HandleDegrees, metavar='deg[DminM[secS]]', nargs='?', default=-112.1124)
     	opts.add_argument('pulser_alt', help='Altitude of pulser in m above sea level.', nargs='?', type=float, default=1775.68)
	# Siple
	#opts.add_argument('pulser_lat', help='Latitude of ANITA payload in +/- decimal degrees N. or DegDMinM[SecS]E/W', action=HandleDegrees, metavar='deg[DminM[secS]]', nargs='?', default=-81.652316667)
	#opts.add_argument('pulser_lat', help='Latitude of ANITA payload in +/- decimal degrees N. or DegDMinM[SecS]E/W', action=HandleDegrees, metavar='deg[DminM[secS]]', nargs='?', default=-149.000166667)
	#opts.add_argument('pulser_lat', help='Latitude of ANITA payload in +/- decimal degrees N. or DegDMinM[SecS]E/W', action=HandleDegrees, metavar='deg[DminM[secS]]', nargs='?', default=0.019812)


     	args = opts.parse_args()

	alt_p = args.pulser_alt
	DLat_p = args.pulser_lat
	DLon_p = args.pulser_lon

	alt_a = args.anita_alt * 0.3048
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
	
	#print "Distance from ANITA to Pulser : ", distance(x_a, y_a, z_a, x_p, y_p, z_p)/1000., " km."
	d_a = distance(x_a, y_a, z_a, x_p, y_p, z_p)
	d_cf = greatcircle_distance(DLat_p, DLon_p, DLat_a, DLon_a)
	dt_cf = d_cf/2.99e8 * 1000 # ms
	dt_a = d_a/2.99e8 * 1000 # ms
	print "Distance from ANITA to Pulser : as-crow-flies %4.2f / with altitude %4.2f km." %( d_cf/1000.,  d_a/1000.)
	print "Arrival time at ANITA : as-crow-flies %2.2f / with altitude %2.2f ms"%(dt_cf, dt_a)
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

