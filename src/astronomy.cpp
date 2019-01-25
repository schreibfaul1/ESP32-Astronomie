/*
 * astronomy.cpp
 *
 *  Created on: 19.12.2017
 *  Updated on: 25.01.2019
 *  Version 1.1.0
 *      Author: Wolle
 */

#include <Arduino.h>
#include "astronomy.h"


Astronomy::Astronomy(as_geo geoa, int8_t deltaT){
	m_Lat     = geoa.latitude;
	m_Lon     = geoa.longitude;
	m_Zone    = geoa.timezone;
	m_DeltaT  = deltaT; // time lag to Universal Time Coordinated [UTC] seconds
}

Astronomy::~Astronomy(){

}
// Calculate Julian date: valid only from 1.3.1901 to 28.2.2100
double Astronomy::CalcJD(int day, int month, int year){
	double jd = 2415020.5 - 64; // 1.1.1900 - correction of algorithm
	if (month <= 2) { year--; month += 12; }
	jd += Int(((year - 1900)) * 365.25);
	jd += Int(30.6001 * ((1 + month)));
	return jd + day;
}
// Julian Date to Greenwich Mean Sidereal Time
double Astronomy::CalcGMST(double JD){
	double UT = frac(JD - 0.5) * 24.0; // UT in hours
	JD = floor(JD - 0.5) + 0.5;   // JD at 0 hours UT
	double T = (JD - 2451545.0) / 36525.0;
	double T0 = 6.697374558 + T * (2400.051336 + T * 0.000025862);
	return(Mod(T0 + UT * 1.002737909, 24.0));
}
// Local Mean Sidereal Time, geographical longitude in radians, East is positive
double Astronomy::GMST2LMST(double gmst, double lon){
	double res=RAD * lon / 15;
	return Mod((gmst + res), 24.0);
}
// Convert Greenwich mean sidereal time to UT
double Astronomy::GMST2UT(double JD, double gmst){
	JD = floor(JD - 0.5) + 0.5;   // JD at 0 hours UT
	double T = (JD - 2451545.0) / 36525.0;
	double T0 = Mod(6.697374558 + T * (2400.051336 + T * 0.000025862), 24.0);
	return 0.9972695663 * ((gmst - T0));
}
// Find GMST of rise/set of object from the two calculates
// (start)points (day 1 and 2) and at midnight UT(0)
double Astronomy::InterpolateGMST(double gmst0, double gmst1, double gmst2, double timefactor)
{
	return ((timefactor * 24.07 * gmst1 - gmst0 * (gmst2 - gmst1)) / (timefactor * 24.07 + gmst1 - gmst2));
}

// Calculate observers cartesian equatorial coordinates (x,y,z in celestial frame)
// from geodetic coordinates (longitude, latitude, height above WGS84 ellipsoid)
// Currently only used to calculate distance of a body from the observer
Astronomy::coor Astronomy::Observer2EquCart(double lon, double lat, double height, double gmst)
{
	double flat = 298.257223563;        // WGS84 flatening of earth
	double aearth = 6378.137;           // GRS80/WGS84 semi major axis of earth ellipsoid
	coor xyz;
	// Calculate geocentric latitude from geodetic latitude
	double co = cos(lat);
	double si = sin(lat);
	double fl = 1.0 - 1.0 / flat;
	fl = fl * fl;
	si = si * si;
	double u = 1.0 / sqrt(co * co + fl * si);
	double a = aearth * u + height;
	double b = aearth * fl * u + height;
	double radius = sqrt(a * a * co * co + b * b * si); // geocentric distance from earth center
	xyz.y = acos(a * co / radius); // geocentric latitude, rad
	xyz.x = lon; // longitude stays the same
	if (lat < 0.0) { xyz.y = -xyz.y; } // adjust sign
	xyz = EquPolar2Cart(xyz.x, xyz.y, radius); // convert from geocentric polar to geocentric cartesian, with regard to Greenwich
	// rotate around earth's polar axis to align coordinate system from Greenwich to vernal equinox
	double x = xyz.x;
	double y = xyz.y;
	double rotangle = gmst / 24.0 * 2.0 * M_PI; // sideral time gmst given in hours. Convert to radians
	xyz.x = x * cos(rotangle) - y * sin(rotangle);
	xyz.y = x * sin(rotangle) + y * cos(rotangle);
	xyz.r = radius;
	xyz.lon = lon;
	xyz.lat = lat;
	return xyz;
}

Astronomy::SIGN Astronomy::Sign(double lon){
	//char* signs[] = { "Widder", "Stier", "Zwillinge", "Krebs", "L�we", "Jungfrau", "Waage", "Skorpion", "Sch�tze", "Steinbock", "Wassermann", "Fische" };
	return (Astronomy::SIGN)((int)floor(lon * RAD / 30.0));
}



// Calculate cartesian from polar coordinates
Astronomy::coor Astronomy::EquPolar2Cart(double lon, double lat, double distance){
	coor xyz;
	double rcd = cos(lat) * distance;
	xyz.x = rcd * cos(lon);
	xyz.y = rcd * sin(lon);
	xyz.z = distance * sin(lat);
	return xyz;
}

// Calculate coordinates for Sun
// Coordinates are accurate to about 10s (right ascension)
// and a few minutes of arc (declination)
Astronomy::coor Astronomy::SunPosition(double TDT, double geolat, double lmst){

	double D = TDT - 2447891.5;

	double eg = 279.403303 * DEG;
	double wg = 282.768422 * DEG;
	double e = 0.016713;
	double a = 149598500; // km
	double diameter0 = 0.533128 * DEG; // angular diameter of Moon at a distance
	double MSun = 360 * DEG / 365.242191 * D + eg - wg;
	double nu = MSun + 360.0 * DEG / M_PI * e * sin(MSun);

	Astronomy::coor sunCoor;
	sunCoor.lon = Mod2Pi(nu + wg);
	sunCoor.lat= 0;
	sunCoor.anomalyMean = MSun;
	sunCoor.distance = (1 - e*e) / (1 + e * cos(nu)); // distance in astronomical units
	sunCoor.diameter = diameter0 / sunCoor.distance; // angular diameter in radians
	sunCoor.distance *=  a;         // distance in km
	sunCoor.parallax = 6378.137 / sunCoor.distance;  // horizonal parallax
	sunCoor = Ecl2Equ(sunCoor, TDT);

	// Calculate horizonal coordinates of sun, if geographic positions is given
	if (!isnan(geolat) && !isnan(lmst))
	{
		sunCoor = Equ2Altaz(sunCoor, TDT, geolat, lmst);
	}
	return sunCoor;
}

// Transform ecliptical coordinates (lon/lat) to equatorial coordinates (RA/dec)
Astronomy::coor Astronomy::Ecl2Equ(Astronomy::coor co, double TDT){
	double T = (TDT - 2451545.0) / 36525.0; // Epoch 2000 January 1.5
	double eps = (23.0 + (26 + 21.45 / 60.0) / 60.0 + T * (-46.815 + T * (-0.0006 + T * 0.00181)) / 3600.0) * DEG;
	double coseps = cos(eps);
	double sineps = sin(eps);
	double sinlon = sin(co.lon);
	co.ra = Mod2Pi(atan2((sinlon * coseps - tan(co.lat) * sineps), cos(co.lon)));
	co.dec = asin(sin(co.lat) * coseps + cos(co.lat) * sineps * sinlon);

	return co;
}





// Transform equatorial coordinates (RA/Dec) to horizonal coordinates (azimuth/altitude)
// Refraction is ignored
Astronomy::coor Astronomy::Equ2Altaz(Astronomy:: coor co, double TDT, double geolat, double lmst){
	double cosdec = cos(co.dec);
	double sindec = sin(co.dec);
	double lha = lmst - co.ra;
	double coslha = cos(lha);
	double sinlha = sin(lha);
	double coslat = cos(geolat);
	double sinlat = sin(geolat);

	double N = -cosdec * sinlha;
	double D = sindec * coslat - cosdec * coslha * sinlat;
	co.az = Mod2Pi(atan2(N, D));
	co.alt = asin(sindec * sinlat + cosdec * coslha * coslat);

	return co;
}

// Calculate data and coordinates for the Moon
// Coordinates are accurate to about 1/5 degree (in ecliptic coordinates)
Astronomy::coor Astronomy::MoonPosition(Astronomy::coor sunCoor, double TDT, Astronomy::coor observer, double lmst){
	double D = TDT - 2447891.5;

	// Mean Moon orbit elements as of 1990.0
	double l0 = 318.351648 * DEG;
	double P0 = 36.340410 * DEG;
	double N0 = 318.510107 * DEG;
	double i = 5.145396 * DEG;
	double e = 0.054900;
	double a = 384401; // km
	double diameter0 = 0.5181 * DEG; // angular diameter of Moon at a distance
	double parallax0 = 0.9507 * DEG; // parallax at distance a

	double l = 13.1763966 * DEG * D + l0;
	double MMoon = l - 0.1114041 * DEG * D - P0; // Moon's mean anomaly M
	double N = N0 - 0.0529539 * DEG * D;       // Moon's mean ascending node longitude
	double C = l - sunCoor.lon;
	double Ev = 1.2739 * DEG * sin(2 * C - MMoon);
	double Ae = 0.1858 * DEG * sin(sunCoor.anomalyMean);
	double A3 = 0.37 * DEG * sin(sunCoor.anomalyMean);
	double MMoon2 = MMoon + Ev - Ae - A3;  // corrected Moon anomaly
	double Ec = 6.2886 * DEG * sin(MMoon2);  // equation of centre
	double A4 = 0.214 * DEG * sin(2 * MMoon2);
	double l2 = l + Ev + Ec - Ae + A4; // corrected Moon's longitude
	double V = 0.6583 * DEG * sin(2 * (l2 - sunCoor.lon));
	double l3 = l2 + V; // true orbital longitude;

	double N2 = N - 0.16 * DEG * sin(sunCoor.anomalyMean);

	Astronomy::coor moonCoor;
	moonCoor.lon = Mod2Pi(N2 + atan2(sin(l3 - N2) * cos(i), cos(l3 - N2)));
	moonCoor.lat = asin(sin(l3 - N2) * sin(i));
	moonCoor.orbitLon = l3;

	moonCoor = Ecl2Equ(moonCoor, TDT);
	// relative distance to semi mayor axis of lunar oribt
	moonCoor.distance = (1 - e*e) / (1 + e * cos(MMoon2 + Ec));
	moonCoor.diameter = diameter0 / moonCoor.distance; // angular diameter in radians
	moonCoor.parallax = parallax0 / moonCoor.distance; // horizontal parallax in radians
	moonCoor.distance = moonCoor.distance * a; // distance in km

	// Calculate horizonal coordinates of sun, if geographic positions is given
	if ((observer.lat+observer.r) > 0 && !isnan(lmst))
	{
		// transform geocentric coordinates into topocentric (==observer based) coordinates
		moonCoor = GeoEqu2TopoEqu(moonCoor, observer, lmst);
		moonCoor.raGeocentric = moonCoor.ra; // backup geocentric coordinates
		moonCoor.decGeocentric = moonCoor.dec;
		moonCoor.ra = moonCoor.raTopocentric;
		moonCoor.dec = moonCoor.decTopocentric;
		moonCoor = Equ2Altaz(moonCoor, TDT, observer.lat, lmst); // now ra and dec are topocentric
	}

	// Age of Moon in radians since New Moon (0) - Full Moon (pi)
	moonCoor.moonAge = Mod2Pi(l3 - sunCoor.lon);
	moonCoor.phase = 0.5 * (1 - cos(moonCoor.moonAge)); // Moon phase, 0-1

	double mainPhase = 1.0 / 29.53 * 360 * DEG; // show 'Newmoon, 'Quarter' for +/-1 day arond the actual event
	double p = Mod(moonCoor.moonAge, 90.0 * DEG);
	if (p < mainPhase || p > 90 * DEG - mainPhase) p = 2 * roundl(moonCoor.moonAge / (90.0 * DEG));
	else p = 2 * floor(moonCoor.moonAge / (90.0 * DEG)) + 1;
	moonCoor.moonPhase = (int)p;

	return moonCoor;
}

// Transform geocentric equatorial coordinates (RA/Dec) to topocentric equatorial coordinates
Astronomy::coor Astronomy::GeoEqu2TopoEqu(Astronomy::coor co, Astronomy::coor observer, double lmst){
	double cosdec = cos(co.dec);
	double sindec = sin(co.dec);
	double coslst = cos(lmst);
	double sinlst = sin(lmst);
	double coslat = cos(observer.lat); // we should use geocentric latitude, not geodetic latitude
	double sinlat = sin(observer.lat);
	double rho = observer.r; // observer-geocenter in Kilometer

	double x = co.distance * cosdec * cos(co.ra) - rho * coslat * coslst;
	double y = co.distance * cosdec * sin(co.ra) - rho * coslat * sinlst;
	double z = co.distance * sindec - rho * sinlat;

	co.distanceTopocentric = sqrt(x * x + y * y + z * z);
	co.decTopocentric = asin(z / co.distanceTopocentric);
	co.raTopocentric = Mod2Pi(atan2(y, x));

	return co;
}
Astronomy::timespan Astronomy::TimeSpan(double tdiff){
	Astronomy::timespan ts;
	char buf[10];String hms;
	m_hh=0; m_mm=0; m_ss=0; m_dv=tdiff;
	if (tdiff == 0.0 || isnan(tdiff)) return ts;
	double m = (tdiff - floor(tdiff)) * 60.0;
	m_hh = Int(tdiff);
	double s = (m - floor(m)) * 60.0;
	m_mm = Int(m);
	if (s >= 59.5) { m_mm++; s -= 60.0; }
	if (m_mm >= 60) { m_hh++; m_mm -= 60; }
	m_ss = (int)roundl(s);
	ts.Hour=m_hh;
	ts.Minute=m_mm;
	ts.Second=m_ss;
	sprintf(buf, "%02d:%02d:%02d", m_hh, m_mm, m_ss);
	ts.HHMMSS= String(buf);
	ts.TotalHour=((float)m_hh + ((float)m_mm + (float)m_ss / 60.0f) / 60.0f);
	ts.TotalMinute=((float)(m_hh * 60 + m_mm) + (float)m_ss / 60.0f);
	ts.TotalSecond=((float)((m_hh * 60 + m_mm) * 60 + m_ss));
	return ts;
}
// Rough refraction formula using standard atmosphere: 1015 mbar and 10°C
// Input true altitude in radians, Output: increase in altitude in degrees
double Astronomy::Refraction(double alt){
	double altdeg = alt * RAD;
	if (altdeg < -2 || altdeg >= 90) return 0.0;

	double pressure = 1015;
	double temperature = 10;
	if (altdeg > 15) return (0.00452 * pressure / ((273 + temperature) * tan(alt)));

	double y = alt;
	double D = 0.0;
	double P = (pressure - 80.0) / 930.0;
	double Q = 0.0048 * (temperature - 10.0);
	double y0 = y;
	double D0 = D;

	for (int i = 0; i < 3; i++)
	{
		double N = y + (7.31 / (y + 4.4));
		N = 1.0 / tan(N * DEG);
		D = N * P / (60.0 + Q * (N + 39.0));
		N = y - y0;
		y0 = D - D0 - N;
		if ((N != 0.0) && (y0 != 0.0)) { N = y - N * (alt + D - y) / y0; }
		else { N = alt + D; }
		y0 = y;
		D0 = D;
		y = N;
	}
	return D; // Hebung durch Refraktion in radians
}
// returns Greenwich sidereal time (hours) of time of rise
// and set of object with coordinates coor.ra/coor.dec
// at geographic position lon/lat (all values in radians)
// Correction for refraction and semi-diameter/parallax of body is taken care of in function RiseSet
// h is used to calculate the twilights. It gives the required elevation of the disk center of the sun
Astronomy::coor Astronomy::GMSTRiseSet(Astronomy::coor co, double lon, double lat, double hn){
	double h = isnan(hn) ? 0.0: hn; // set default value
	Astronomy::coor riseset;
	//  double tagbogen = std::acos(-std::tan(lat)*std::tan(coor["dec"])); // simple formula if twilight is not required
	double tagbogen = acos((sin(h) - sin(lat) * sin(co.dec)) / (cos(lat) * cos(co.dec)));

	riseset.transit = RAD / 15 * (+co.ra - lon);
	riseset.rise = 24.0 + RAD / 15 * (-tagbogen + co.ra - lon); // calculate GMST of rise of object
	riseset.set = RAD / 15 * (+tagbogen + co.ra - lon); // calculate GMST of set of object

	// using the modulo function Mod, the day number goes missing. This may get a problem for the moon
	riseset.transit = Mod(riseset.transit, 24);
	riseset.rise = Mod(riseset.rise, 24);
	riseset.set = Mod(riseset.set, 24);

	return riseset;
}
// Find GMST of rise/set of object from the two calculates
// (start)points (day 1 and 2) and at midnight UT(0)
double InterpolateGMST(double gmst0, double gmst1, double gmst2, double timefactor){
	return ((timefactor * 24.07 * gmst1 - gmst0 * (gmst2 - gmst1)) / (timefactor * 24.07 + gmst1 - gmst2));
}
// JD is the Julian Date of 0h UTC time (midnight)
Astronomy::coor Astronomy::RiseSet(double jd0UT, Astronomy::coor  coor1, Astronomy::coor  coor2, double lon, double lat, double timeinterval, double naltitude)
{
	// altitude of sun center: semi-diameter, horizontal parallax and (standard) refraction of 34'
	double alt = 0.0; // calculate
	double altitude = isnan(naltitude) ? 0.0 : naltitude; // set default value

	// true height of sun center for sunrise and set calculation. Is kept 0 for twilight (ie. altitude given):
	if (altitude == 0.0) alt = 0.5 * coor1.diameter - coor1.parallax + 34.0 / 60 * DEG;

	Astronomy::coor rise1 = GMSTRiseSet(coor1, lon, lat, altitude);
	Astronomy::coor rise2 = GMSTRiseSet(coor2, lon, lat, altitude);

	Astronomy::coor rise;

	// unwrap GMST in case we move across 24h -> 0h
	if (rise1.transit > rise2.transit && abs(rise1.transit - rise2.transit) > 18) rise2.transit += 24.0;
	if (rise1.rise > rise2.rise && abs(rise1.rise - rise2.rise) > 18) rise2.rise += 24.0;
	if (rise1.set > rise2.set && abs(rise1.set - rise2.set) > 18) rise2.set += 24.0;
	double T0 = CalcGMST(jd0UT);
	//  var T02 = T0-zone*1.002738; // Greenwich sidereal time at 0h time zone (zone: hours)

	// Greenwich sidereal time for 0h at selected longitude
	double T02 = T0 - lon * RAD / 15 * 1.002738;
	if (T02 < 0) T02 += 24.0;

	if (rise1.transit < T02) { rise1.transit += 24.0; rise2.transit += 24.0; }
	if (rise1.rise < T02) { rise1.rise += 24.0; rise2.rise += 24.0; }
	if (rise1.set < T02) { rise1.set += 24.0; rise2.set += 24.0; }

	// Refraction and Parallax correction
	double decMean = 0.5 * (coor1.dec + coor2.dec);
	double psi = acos(sin(lat) / cos(decMean));
	double y = asin(sin(alt) / sin(psi));
	double dt = 240 * RAD * y / cos(decMean) / 3600; // time correction due to refraction, parallax
	rise.transit = GMST2UT(jd0UT, InterpolateGMST(T0, rise1.transit, rise2.transit, timeinterval));
	rise.rise = GMST2UT(jd0UT, InterpolateGMST(T0, rise1.rise, rise2.rise, timeinterval) - dt);
	rise.set = GMST2UT(jd0UT, InterpolateGMST(T0, rise1.set, rise2.set, timeinterval) + dt);

	return (rise);
}
// Find local time of moonrise and moonset
// JD is the Julian Date of 0h local time (midnight)
// Accurate to about 5 minutes or better
// recursive: 1 - calculate rise/set in UTC
// recursive: 0 - find rise/set on the current local day (set could also be first)
// returns '' for moonrise/set does not occur on selected day
Astronomy::coor Astronomy::CalcMoonRise(double JD, double deltaT, double lon, double lat, int zone, bool recursive){
	double timeinterval = 0.5;
	double jd0UT = floor(JD - 0.5) + 0.5;   // JD at 0 hours UT
	Astronomy::coor suncoor1 = SunPosition(jd0UT + deltaT / 24.0 / 3600.0);
	Astronomy::coor coor1 = MoonPosition(suncoor1, jd0UT + deltaT / 24.0 / 3600.0);
	Astronomy::coor suncoor2 = SunPosition(jd0UT + timeinterval + deltaT / 24.0 / 3600.0); // calculations for noon
	// calculations for next day's midnight
	Astronomy::coor coor2 = MoonPosition(suncoor2, jd0UT + timeinterval + deltaT / 24.0 / 3600.0);

	Astronomy::coor risetemp;
	// rise/set time in UTC, time zone corrected later.
	// Taking into account refraction, semi-diameter and parallax
	Astronomy::coor rise = RiseSet(jd0UT, coor1, coor2, lon, lat, timeinterval);

	if (!recursive)
	{ // check and adjust to have rise/set time on local calendar day
		if (zone > 0)
		{
			// recursive call to MoonRise returns events in UTC
			risetemp = CalcMoonRise(JD - 1.0, deltaT, lon, lat, zone, true);
			if (rise.transit >= 24.0 - zone || rise.transit < -zone)
			{ // transit time is tomorrow local time
				if (risetemp.transit < 24.0 - zone) rise.transit = NAN_DOUBLE; // there is no moontransit today
				else rise.transit = risetemp.transit;
			}

			if (rise.rise >= 24.0 - zone || rise.rise < -zone)
			{ // rise time is tomorrow local time
				if (risetemp.rise < 24.0 - zone) rise.rise = NAN_DOUBLE; // there is no moontransit today
				else rise.rise = risetemp.rise;
			}

			if (rise.set >= 24.0 - zone || rise.set < -zone)
			{ // set time is tomorrow local time
				if (risetemp.set < 24.0 - zone) rise.set = NAN_DOUBLE; // there is no moontransit today
				else rise.set = risetemp.set;
			}

		}
		else if (zone < 0)
		{
			// rise/set time was tomorrow local time -> calculate rise time for former UTC day
			if (rise.rise < -zone || rise.set < -zone || rise.transit < -zone)
			{
				risetemp = CalcMoonRise(JD + 1.0, deltaT, lon, lat, zone, true);

				if (rise.rise < -zone)
				{
					if (risetemp.rise > -zone) rise.rise = NAN_DOUBLE; // there is no moonrise today
					else rise.rise = risetemp.rise;
				}

				if (rise.transit < -zone)
				{
					if (risetemp.transit > -zone) rise.transit = NAN_DOUBLE; // there is no moonset today
					else rise.transit = risetemp.transit;
				}

				if (rise.set < -zone)
				{
					if (risetemp.set > -zone) rise.set = NAN_DOUBLE; // there is no moonset today
					else rise.set = risetemp.set;
				}

			}
		}

		if (rise.rise != NAN_DOUBLE) rise.rise = Mod(rise.rise + zone, 24.0);    // correct for time zone, if time is valid
		if (rise.transit != NAN_DOUBLE) rise.transit = Mod(rise.transit + zone, 24.0); // correct for time zone, if time is valid
		if (rise.set != NAN_DOUBLE) rise.set = Mod(rise.set + zone, 24.0);    // correct for time zone, if time is valid
	}
	return rise;
}



// Find (local) time of sunrise and sunset, and twilights
// JD is the Julian Date of 0h local time (midnight)
// Accurate to about 1-2 minutes
// recursive: 1 - calculate rise/set in UTC in a second run
// recursive: 0 - find rise/set on the current local day. This is set when doing the first call to this function
Astronomy::coor Astronomy::CalcSunRise(double JD, double deltaT, double lon, double lat, int zone, bool recursive){
	double jd0UT = floor(JD - 0.5) + 0.5;   // JD at 0 hours UT
	Astronomy::coor coor1 = SunPosition(jd0UT + deltaT / 24.0 / 3600.0);
	Astronomy::coor coor2 = SunPosition(jd0UT + 1.0 + deltaT / 24.0 / 3600.0); // calculations for next day's UTC midnight

	Astronomy::coor risetemp;
	// rise/set time in UTC.
	Astronomy::coor rise = RiseSet(jd0UT, coor1, coor2, lon, lat, 1);
	if (!recursive)
	{ // check and adjust to have rise/set time on local calendar day
		if (zone > 0)
		{
			// rise time was yesterday local time -> calculate rise time for next UTC day
			if (rise.rise >= 24 - zone || rise.transit >= 24 - zone || rise.set >= 24 - zone)
			{
				risetemp = CalcSunRise(JD + 1, deltaT, lon, lat, zone, true);
				if (rise.rise >= 24 - zone) rise.rise = risetemp.rise;
				if (rise.transit >= 24 - zone) rise.transit = risetemp.transit;
				if (rise.set >= 24 - zone) rise.set = risetemp.set;
			}
		}
		else if (zone < 0)
		{
			// rise time was yesterday local time -> calculate rise time for next UTC day
			if (rise.rise < -zone || rise.transit < -zone || rise.set < -zone)
			{
				risetemp = CalcSunRise(JD - 1, deltaT, lon, lat, zone, true);
				if (rise.rise < -zone) rise.rise = risetemp.rise;
				if (rise.transit < -zone) rise.transit = risetemp.transit;
				if (rise.set < -zone) rise.set = risetemp.set;
			}
		}

		rise.transit = Mod(rise.transit + zone, 24.0);
		rise.rise = Mod(rise.rise + zone, 24.0);
		rise.set = Mod(rise.set + zone, 24.0);

		// Twilight calculation
		// civil twilight time in UTC.
		risetemp = RiseSet(jd0UT, coor1, coor2, lon, lat, 1, -6.0 * DEG);
		rise.cicilTwilightMorning = Mod(risetemp.rise + zone, 24.0);
		rise.cicilTwilightEvening = Mod(risetemp.set + zone, 24.0);

		// nautical twilight time in UTC.
		risetemp = RiseSet(jd0UT, coor1, coor2, lon, lat, 1, -12.0 * DEG);
		rise.nauticalTwilightMorning = Mod(risetemp.rise + zone, 24.0);
		rise.nauticalTwilightEvening = Mod(risetemp.set + zone, 24.0);

		// astronomical twilight time in UTC.
		risetemp = RiseSet(jd0UT, coor1, coor2, lon, lat, 1, -18.0 * DEG);
		rise.astronomicalTwilightMorning = Mod(risetemp.rise + zone, 24.0);
		rise.astronomicalTwilightEvening = Mod(risetemp.set + zone, 24.0);
	}
	return rise;
}


void Astronomy::setInput(as_date d, as_time t){
	String res="";  char buf[20];

	double JD0 = CalcJD(d.day, d.month, d.year);
	double jd = JD0 + (t.hour - m_Zone + t.minute / 60.0 + t.second / 3600.0) / 24.0;
	double TDT = jd + m_DeltaT / 24.0 / 3600.0;
	double lat = m_Lat * DEG; // geodetic latitude of observer on WGS84
	double lon = m_Lon * DEG; // latitude of observer
	double height = 0 * 0.001; // altiude of observer in meters above WGS84 ellipsoid (and converted to kilometers)
	double gmst = CalcGMST(jd);
	double lmst = GMST2LMST(gmst, lon);
	observerCart = Observer2EquCart(lon, lat, height, gmst); // geocentric cartesian coordinates of observer
	sunCoor = SunPosition(TDT, lat, lmst * 15.0 * DEG);   // Calculate data for the Sun at given time
	moonCoor = MoonPosition(sunCoor, TDT, observerCart, lmst * 15.0 * DEG);    // Calculate data for the Moon at given time

	m_JD = round100000(jd);
	m_GMST = TimeSpan(gmst);
	m_LMST = TimeSpan(lmst);

	m_SunLon = round1000(sunCoor.lon * RAD);
	m_SunRA = TimeSpan(sunCoor.ra * RAD / 15);
	m_SunDec = round1000(sunCoor.dec * RAD);
	m_SunAz = round100(sunCoor.az * RAD);
	m_SunAlt = round10(sunCoor.alt * RAD + Refraction(sunCoor.alt));  // including refraction

	m_SunSign = Sign(sunCoor.lon);
	m_SunDiameter = round100(sunCoor.diameter * RAD * 60.0); // angular diameter in arc seconds
	m_SunDistance = round10(sunCoor.distance);

	// Calculate distance from the observer (on the surface of earth) to the center of the sun
	sunCart = EquPolar2Cart(sunCoor.ra, sunCoor.dec, sunCoor.distance);
	double sunCardxSqr=(sunCart.x - observerCart.x) * (sunCart.x - observerCart.x);
	double sunCardySqr=(sunCart.y - observerCart.y) * (sunCart.y - observerCart.y);
	double sunCartzSqr=(sunCart.z - observerCart.z) * (sunCart.z - observerCart.z);
	m_SunDistanceObserver = round10(sqrt(sunCardxSqr  + sunCardySqr  + sunCartzSqr));

	sunRise = CalcSunRise(JD0, m_DeltaT, lon, lat, m_Zone, false);
	m_SunTransit = TimeSpan(sunRise.transit);
	m_SunRise = TimeSpan(sunRise.rise);
	m_SunSet = TimeSpan(sunRise.set);
	m_SunCivilTwilightMorning = TimeSpan(sunRise.cicilTwilightMorning);
	m_SunCivilTwilightEvening = TimeSpan(sunRise.cicilTwilightEvening);
	m_SunNauticalTwilightMorning = TimeSpan(sunRise.nauticalTwilightMorning);
	m_SunNauticalTwilightEvening = TimeSpan(sunRise.nauticalTwilightEvening);
	m_SunAstronomicalTwilightMorning = TimeSpan(sunRise.astronomicalTwilightMorning);
	m_SunAstronomicalTwilightEvening = TimeSpan(sunRise.astronomicalTwilightEvening);

	m_MoonLon = round1000(moonCoor.lon * RAD);
	m_MoonLat = round1000(moonCoor.lat * RAD);
	m_MoonRA = TimeSpan(moonCoor.ra * RAD / 15.0);
	m_MoonDec = round1000(moonCoor.dec * RAD);
	m_MoonAz = round100(moonCoor.az * RAD);
	m_MoonAlt = round10(moonCoor.alt * RAD + Refraction(moonCoor.alt));  // including refraction
	m_MoonAge = round1000(moonCoor.moonAge * RAD);
	m_MoonPhaseNumber = round1000(moonCoor.phase);

	int phase = (int)moonCoor.moonPhase;
	if (phase == 8) phase = 0;
	m_MoonPhase = (LUNARPHASE)phase;

	m_MoonSign = Sign(moonCoor.lon);
	m_MoonDistance = round10(moonCoor.distance);
	m_MoonDiameter = round100(moonCoor.diameter * RAD * 60.0); // angular diameter in arc seconds

	// Calculate distance from the observer (on the surface of earth) to the center of the moon
	moonCart = EquPolar2Cart(moonCoor.raGeocentric, moonCoor.decGeocentric, moonCoor.distance);
	double moonCardxSqr=(moonCart.x - observerCart.x) * (moonCart.x - observerCart.x);
	double moonCardySqr=(moonCart.y - observerCart.y) * (moonCart.y - observerCart.y);
	double moonCartzSqr=(moonCart.z - observerCart.z) * (moonCart.z - observerCart.z);
	m_MoonDistanceObserver = round10(sqrt(moonCardxSqr + moonCardySqr + moonCartzSqr));

	moonRise = CalcMoonRise(JD0, m_DeltaT, lon, lat, m_Zone, false);

	m_MoonTransit = TimeSpan(moonRise.transit);
	m_MoonRise = TimeSpan(moonRise.rise);
	m_MoonSet = TimeSpan(moonRise.set);

    sprintf(buf,"%02d:%02d:%02d",t.hour, t.minute, t.second); m_Time=String(buf);
    sprintf(buf,"%02d.%02d.%04d",d.day, d.month, d.year); m_Date=String(buf);
    m_os = "\n";
    m_os+= "Eastern geographic length: .............. " + String(m_Lon) + " deg\n" ;
    m_os+= "Geographic latitude: .................... " + String(m_Lat) + " deg\n";
    m_os+= "Date/Time: .............................. " + m_Date + " " + m_Time + "\n";
    m_os+= "Time difference to world time: .......... " + String(m_Zone) + " h\n";
    m_os+= "deltaT: ................................. " + String(m_DeltaT) +  " sec\n";
    m_os+= "Julian Date: ............................ " + String(m_JD) + " days\n";
    m_os+= "Greenwich sidereal  GMST: ............... " + m_GMST.HHMMSS + " = " + String(m_GMST.TotalHour, 4) + " h\n";
    m_os+= "Local sidereal time LMST: ............... " + m_LMST.HHMMSS + " = " +String(m_LMST.TotalHour, 4) + " h\n";
    m_os+= "Distance to the sun (Earth's center): ... " + String(m_SunDistance) + " km\n";
    m_os+= "Distance to the sun (from the observer):. " + String(m_SunDistanceObserver) + " km\n";
	m_os+= "Ecliptic length of the sun: ............. " + String(m_SunLon) + " deg\n";
	m_os+= "Right ascension of the sun: ............. " + m_SunRA.HHMMSS + " = " + String(m_SunRA.TotalHour, 4) + " h\n";
	m_os+= "Declination of the sun: ................. " + String(m_SunDec) + " deg\n";
	m_os+= "Azimuth of the sun: ..................... " + String(m_SunAz) + " deg\n";
	m_os+= "Height of the sun above the horizon: .... " + String(m_SunAlt) + " deg\n";
	m_os+= "Diameter of the sun: .................... " + String(m_SunDiameter) + " '\n";
	m_os+= "Astronomical dawn: ...................... " + m_SunAstronomicalTwilightMorning.HHMMSS + " = " + String(m_SunAstronomicalTwilightMorning.TotalHour) + " h\n";
	m_os+= "Nautical dawn: .......................... " + m_SunNauticalTwilightMorning.HHMMSS + " = " + String(m_SunNauticalTwilightMorning.TotalHour) + " h\n";
	m_os+= "Civil dawn: ............................. " + m_SunCivilTwilightMorning.HHMMSS + " = " + String(m_SunCivilTwilightMorning.TotalHour) + " h\n";
	m_os+= "Sunrise: ................................ " + m_SunRise.HHMMSS + " = " + String(m_SunRise.TotalHour) + " h\n";
	m_os+= "Culmination of the sun: ................. " + m_SunTransit.HHMMSS + " = " + String(m_SunTransit.TotalHour) + " h\n";
	m_os+= "Sunset: ................................. " + m_SunSet.HHMMSS + " = " + String(m_SunSet.TotalHour) + " h\n";
	m_os+= "Civil dusk: ............................. " + m_SunCivilTwilightEvening.HHMMSS + " = " + String(m_SunCivilTwilightEvening.TotalHour) + " h\n";
	m_os+= "Nautical dusk: .......................... " + m_SunNauticalTwilightEvening.HHMMSS + " = " + String(m_SunNauticalTwilightEvening.TotalHour) + " h\n";
	m_os+= "Astronomical dusk: ...................... " + m_SunAstronomicalTwilightEvening.HHMMSS + " = " + String(m_SunAstronomicalTwilightEvening.TotalHour) + " h\n";
	m_os+= "Signs of the zodiac: .................... " + TKZ[(int)m_SunSign] + "\n";
	m_os+= "Distance to the moon (Earth's center): .. " + String(m_MoonDistance) + " km\n";
	m_os+= "Distance to the moon (from the observer): " + String(m_MoonDistanceObserver) + " km\n";
	m_os+= "Ecliptic longitude of the moon: ......... " + String(m_MoonLon) + " deg\n";
	m_os+= "Ecliptic latitude of the moon: .......... " + String(m_MoonLat) + " deg\n";
	m_os+= "Right ascension of the moon: ............ " + m_MoonRA.HHMMSS + " = " + String(m_MoonRA.TotalHour) + " h\n";
	m_os+= "Declination of the moon: ................ " + String(m_MoonDec) + " deg\n";
	m_os+= "Azimuth of the moon: .................... " + String(m_MoonAz) + " deg\n";
	m_os+= "Height of the moon above the horizon: ... " + String(m_MoonAlt) + " deg\n";
	m_os+= "Diameter of the moon: ................... " + String(m_MoonDiameter) + " '\n";
	m_os+= "Moonrise: ............................... " + m_MoonRise.HHMMSS + " = " + String(m_MoonRise.TotalHour) + " h\n";
	m_os+= "Culmination of the moon: ................ " + m_MoonTransit.HHMMSS + " = " + String(m_MoonTransit.TotalHour) + " h\n";
	m_os+= "Moonset: ................................ " + m_MoonSet.HHMMSS + " = " + String(m_MoonSet.TotalHour) + " h\n";
	m_os+= "Moon phase number: ...................... " + String(m_MoonPhaseNumber)+ "\n";
	m_os+= "Mon age: ................................ " + String(m_MoonAge) + " deg\n";
	m_os+= "Moon phase: ............................. " + phases[(int)m_MoonPhase] + "\n";
	m_os+= "Mon sign: ............................... " + TKZ[(int)m_MoonSign] + "\n";
}

