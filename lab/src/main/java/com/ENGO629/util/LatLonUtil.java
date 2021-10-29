package com.ENGO629.util;

public class LatLonUtil {
	// All are WGS-84 params
	// Semi-major axis or Equatorial radius
	private static final double a = 6378137;
	// flattening
	private static final double f = 1 / 298.257223563;
	// Semi-minor axis or Polar radius
	private static final double b = 6356752.314245;
	private static final double e = Math.sqrt((Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2));
	private static final double e2 = Math.sqrt((Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(b, 2));

	public static double[] ecef2lla(double[] ECEF) {

		double x = ECEF[0];
		double y = ECEF[1];
		double z = ECEF[2];
		double[] lla = { 0, 0, 0 };
		double lat = 0, lon, height, N, theta, p;

		p = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2));

		theta = Math.atan((z * a) / (p * b));

		lon = Math.atan(y / x);

		if (x < 0) {
			if (y > 0) {
				lon = Math.PI + lon;
			} else {
				lon = -Math.PI + lon;
			}
		}
		for (int i = 0; i < 3; i++) {
			lat = Math.atan(((z + Math.pow(e2, 2) * b * Math.pow(Math.sin(theta), 3))
					/ ((p - Math.pow(e, 2) * a * Math.pow(Math.cos(theta), 3)))));
			theta = Math.atan((Math.tan(lat) * b) / (a));

		}

		N = a / (Math.sqrt(1 - (Math.pow(e, 2) * Math.pow(Math.sin(lat), 2))));

		height = (p * Math.cos(lat)) + ((z + (Math.pow(e, 2) * N * Math.sin(lat))) * Math.sin(lat)) - N;

		lon = lon * 180 / Math.PI;

		lat = lat * 180 / Math.PI;
		lla[0] = lat;
		lla[1] = lon;
		lla[2] = height;
		return lla;
	}

	public static double getHaversineDistance(double[] LatLon1, double[] LatLon2) {

		double lat1 = LatLon1[0];
		double lon1 = LatLon1[1];
		double lat2 = LatLon2[0];
		double lon2 = LatLon2[1];
		// The math module contains a function named toRadians which converts from
		// degrees to radians.
		lon1 = Math.toRadians(lon1);
		lon2 = Math.toRadians(lon2);
		lat1 = Math.toRadians(lat1);
		lat2 = Math.toRadians(lat2);

		// Haversine formula
		double dlon = lon2 - lon1;
		double dlat = lat2 - lat1;
		double a = Math.pow(Math.sin(dlat / 2), 2) + Math.cos(lat1) * Math.cos(lat2) * Math.pow(Math.sin(dlon / 2), 2);

		double c = 2 * Math.asin(Math.sqrt(a));

		// Used Mean Earth Radius
		double r = 6371000;

		// calculate the result
		return (c * r);
	}

}
