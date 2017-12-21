/*
 * astronomy.h
 *
 *  Created on: 19.12.2017
 *      Author: Wolle
 *
 *   Entfernen Sie folgende Informationen auf keinen Fall: / Do not remove following text:
 *   Source code by Arnold Barmettler, www.astronomie.info / www.CalSky.com
 *   based on algorithms by Peter Duffett-Smith's great and easy book
 *   'Practical Astronomy with your Calculator'. ESP32 adaption by schreibfaul1
 *
 */

#ifndef ASTRONOMY_H_
#define ASTRONOMY_H_

#include<Arduino.h>
#include <limits>

struct as_geo{
	double longitude;
	double latitude;
	int timezone;
};

struct as_time{
	uint8_t hour;
	uint8_t minute;
	uint8_t second;
};

struct as_date{
	uint8_t  day;
	uint8_t  month;
	uint16_t year;
};

class Astronomy{

#define NAN_DOUBLE std::numeric_limits<double>::quiet_NaN()
#define NAN_INT std::numeric_limits<int>::quiet_NaN()

private:
	const double DEG=(M_PI/180.0);
	const double RAD=(180.0/M_PI);


	struct coor{
		double x;
		double y;
		double z;
		double r;
		double az;
		double ra;
		double alt;
		double dec;
		double lon;
		double lat;
		double set;
		double phase;
		double rise;
		double moonAge;
		double moonPhase;
		double anomalyMean;
		double distance;
		double diameter;
		double parallax;
		double orbitLon;
		double transit;
		double raGeocentric;
		double decGeocentric;
		double raTopocentric;
		double decTopocentric;
		double distanceTopocentric;
		double cicilTwilightMorning;
		double cicilTwilightEvening;
		double nauticalTwilightMorning;
		double nauticalTwilightEvening;
		double astronomicalTwilightMorning;
		double astronomicalTwilightEvening;
	};
	coor observerCart, sunCoor, moonCoor, sunCart, sunRise, moonCart, moonRise;

	struct timespan{
		uint32_t Hour;
		uint32_t Minute;
		uint32_t Second;
		String HHMMSS;
		double TotalHour;
		double TotalMinute;
		double TotalSecond;
	};

	enum SIGN
	{
		SIGN_ARIES,			//!< Widder
		SIGN_TAURUS,		//!< Stier
		SIGN_GEMINI,		//!< Zwillinge
		SIGN_CANCER,		//!< Krebs
		SIGN_LEO,			//!< Löwe
		SIGN_VIRGO,			//!< Jungfrau
		SIGN_LIBRA,			//!< Waage
		SIGN_SCORPIO,		//!< Skorpion
		SIGN_SAGITTARIUS,	//!< Schütze
		SIGN_CAPRICORNUS,	//!< Steinbock
		SIGN_AQUARIUS,		//!< Wassermann
		SIGN_PISCES			//!< Fische
	};
	String TKZ[12] {"Widder", "Stier", "Zwillinge", "Krebs", "Loewe", "Jungfrau",
					"Waage", "Skorpion", "Schuetze", "Steinbock", "Wassermann", "Fische" };

	String phases[8] {  "Neumond", "Zunehmende Sichel", "Erstes Viertel", "Zunehmender Mond",
						"Vollmond", "Abnehmender Mond", "Letztes Viertel", "Abnehmende Sichel" };

	enum LUNARPHASE
	{
		LP_NEW_MOON,                //!< Neumond
		LP_WAXING_CRESCENT_MOON,    //!< Zunehmende Sichel
		LP_FIRST_QUARTER_MOON,      //!< Erstes Viertel
		LP_WAXING_GIBBOUS_MOON,     //!< Zunehmender Mond
		LP_FULL_MOON,               //!< Vollmond
		LP_WANING_GIBBOUS_MOON,     //!< Abnehmender Mond
		LP_LAST_QUARTER_MOON,       //!< Letztes Viertel
		LP_WANING_CRESCENT_MOON,    //!< Abnehmende Sichel
	};

	double m_Lat=0;
	double m_Lon=0;
	double m_Zone=0;
	double m_DeltaT=0;
	double m_JD=0;
	double m_SunLon=0;
	double m_SunDistance=0;
	double m_SunDec=0;
	double m_SunAz=0;
	double m_SunAlt=0;
	double m_SunDiameter=0;
	double m_SunDistanceObserver=0;
	double m_MoonDistance=0;
	double m_MoonDistanceObserver=0;
	double m_MoonLon=0;
	double m_MoonLat=0;
	double m_MoonDec=0;
	double m_MoonAz=0;
	double m_MoonAlt=0;
	double m_MoonDiameter=0;
	double m_MoonPhaseNumber=0;
	double m_MoonAge=0;
	double m_dv=0;
	uint32_t m_hh=0;
	uint32_t m_mm=0;
	uint32_t m_ss=0;
	String m_os;
	String m_Date="";
	String m_Time="";
	LUNARPHASE m_MoonPhase=LP_NEW_MOON;
	SIGN m_MoonSign=SIGN_ARIES;
	SIGN m_SunSign=SIGN_ARIES;
	timespan m_GMST;
	timespan m_LMST;
	timespan m_SunRA;
	timespan m_SunTransit;
	timespan m_SunRise, m_SunSet;
	timespan m_SunCivilTwilightMorning, m_SunCivilTwilightEvening;
	timespan m_SunNauticalTwilightMorning, m_SunNauticalTwilightEvening;
	timespan m_SunAstronomicalTwilightMorning, m_SunAstronomicalTwilightEvening;
	timespan m_MoonRA;
	timespan m_MoonRise;
	timespan m_MoonTransit;
	timespan m_MoonSet;


public:

	Astronomy(as_geo, int8_t deltaT=65);
	~Astronomy();
	String setInput(as_date, as_time);
	inline String GetAll() {return m_os;}
	inline double GetLat() {return m_Lat;}
	inline double GetLon() {return m_Lon;}
	inline String GetDate() {return m_Date;}
	inline String GetTime() {return m_Time;}
	inline double GetJD() {return m_JD;}
	inline double GetZone() {return m_Zone;}
	inline double GetDeltaT() {return m_DeltaT;}
	inline double GetSunDistance() {return m_SunDistance;}
	inline double GetSunDistanceObserver() {return m_SunDistanceObserver;}
	inline double GetSunLon() {return m_SunLon;}
	inline double GetSunDec() {return m_SunDec;}
	inline double GetSunAz() {return m_SunAz;}
	inline double GetSunAlt() {return m_SunAlt;}
	inline double GetSunDiameter() {return m_SunDiameter;}
	inline double GetMoonDistance() {return m_MoonDistance; }
	inline double GetMoonDistanceObserver() {return m_MoonDistanceObserver;}
	inline double GetMoonLon() {return m_MoonLon;}
	inline double GetMoonLat() {return m_MoonLat;}
	inline double GetMoonDec() {return m_MoonDec;}
	inline double GetMoonAz() {return m_MoonAz;}
	inline double GetMoonAlt() {return m_MoonAlt;}
	inline double GetMoonDiameter() {return m_MoonDiameter;}
	inline double GetMoonPhaseNumber() {return m_MoonPhaseNumber;}
	inline double GetMoonAge() {return m_MoonAge;}
	inline double GetSunAstronomicalTwilightMorning() { return m_SunAstronomicalTwilightMorning.TotalHour;}
	inline double GetSunNauticalTwilightMorning() { return m_SunNauticalTwilightMorning.TotalHour;}
	inline double GetSunCivilTwilightMorning() { return m_SunCivilTwilightMorning.TotalHour;}
	inline double GetSunRise() { return m_SunRise.TotalHour;}
	inline double GetSunTransit() { return m_SunTransit.TotalHour;}
	inline double GetSunSet() { return m_SunSet.TotalHour;}
	inline double GetSunCivilTwilightEvening() { return m_SunCivilTwilightEvening.TotalHour;}
	inline double GetSunNauticalTwilightEvening() { return m_SunNauticalTwilightEvening.TotalHour;}
	inline double GetSunAstronomicalTwilightEvening() { return m_SunAstronomicalTwilightEvening.TotalHour;}
	inline double GetMoonRA() { return m_MoonRA.TotalHour;}
	inline double GetMoonRise() { return m_MoonRise.TotalHour;}
	inline double GetMoonTransit() { return m_MoonTransit.TotalHour;}
	inline double GetMoonSet() { return m_MoonSet.TotalHour;}
	inline String GetSunAstronomicalTwilightMorning_s() {return m_SunAstronomicalTwilightMorning.HHMMSS;}
	inline String GetSunNauticalTwilightMorning_s() {return m_SunNauticalTwilightMorning.HHMMSS;}
	inline String GetSunCivilTwilightMorning_s() {return m_SunCivilTwilightMorning.HHMMSS;}
	inline String GetSunRise_s() {return m_SunRise.HHMMSS;}
	inline String GetSunTransit_s() {return m_SunTransit.HHMMSS;}
	inline String GetSunSet_s() {return m_SunSet.HHMMSS;}
	inline String GetSunCivilTwilightEvening_s() {return m_SunCivilTwilightEvening.HHMMSS;}
	inline String GetSunNauticalTwilightEvening_s() {return m_SunNauticalTwilightEvening.HHMMSS;}
	inline String GetSunAstronomicalTwilightEvening_s() {return m_SunAstronomicalTwilightEvening.HHMMSS;}
	inline String GetMoonRA_s() {return m_MoonRA.HHMMSS;}
	inline String GetMoonRise_s() {return m_MoonRise.HHMMSS;}
	inline String GetMoonTransit_s() {return m_MoonTransit.HHMMSS;}
	inline String GetMoonSet_s() {return m_MoonSet.HHMMSS;}
	inline String GetMoonPhase() {return phases[(int)m_MoonPhase];}
	inline String GetMoonSign() {return TKZ[(int)m_MoonSign];}
	inline String GetSunSign() {return TKZ[(int)m_SunSign];}

private:


protected:
	double CalcJD(int day, int month, int year); // Calculate Julian date: valid only from 1.3.1901 to 28.2.2100
	double CalcGMST(double JD);
	double GMST2LMST(double gmst, double lon);
	double Refraction(double alt);
	double GMST2UT(double JD, double gmst);
	double InterpolateGMST(double gmst0, double gmst1, double gmst2, double timefactor);
	coor EquPolar2Cart(double lon, double lat, double distance);
	coor Observer2EquCart(double lon, double lat, double height, double gmst);
	coor SunPosition(double TDT, double geolat = NAN_DOUBLE, double lmst = NAN_DOUBLE);
	coor Equ2Altaz(coor co, double TDT, double geolat, double lmst);
	coor Ecl2Equ(coor co, double TDT);
	coor MoonPosition(coor sunCoor, double TDT, coor observer=coor(), double lmst=NAN_DOUBLE);
	coor GeoEqu2TopoEqu(coor co, coor observer, double lmst);
	coor RiseSet(double jd0UT, coor coor1, coor coor2, double lon, double lat, double timeinterval, double naltitude = NAN_DOUBLE);
	coor GMSTRiseSet(coor co, double lon, double lat, double hn = NAN_DOUBLE);
	coor CalcSunRise(double JD, double deltaT, double lon, double lat, int zone, bool recursive);
	coor CalcMoonRise(double JD, double deltaT, double lon, double lat, int zone, bool recursive);
	SIGN Sign(double lon);
	inline int Int(double x) {return (x < 0) ? (int)ceil(x) : (int)floor(x);}
	inline double frac(double x) {return (x - floor(x));}
	inline double Mod(double a, double b) {return (a - floor(a / b) * b);}
	inline double Mod2Pi(double x) {return Mod(x, 2.0 * M_PI);}
	inline double round10(double x) {return (round(10.0 * x) / 10.0);}
	inline double round100(double x) {return (roundl(100.0 * x) / 100.0);}
	inline double round1000(double x) {return (roundl(1000.0 * x) / 1000.0);}
	inline double round10000(double x) {return (roundl(10000.0 * x) / 10000.0);}
	inline double round100000(double x) {return (roundl(100000.0 * x) / 100000.0);}
	timespan TimeSpan(double tdiff);


};























#endif /* ASTRONOMY_H_ */
