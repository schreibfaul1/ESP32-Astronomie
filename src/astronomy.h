/*
 * astronomy.h
 *
 *  Created on: 19.12.2017
 *  Updated on: 25.01.2019
 *  Version 1.1.0
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
//	String TKZ[12] {"Widder", "Stier", "Zwillinge", "Krebs", "Loewe", "Jungfrau",
//					"Waage", "Skorpion", "Schuetze", "Steinbock", "Wassermann", "Fische" };

	String TKZ[12] {"Aries", "Taurus",  "Gemini",      "Cancer",    "Leo",      "Virgo",
	                "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces"};

//	String phases[8] {  "Neumond", "Zunehmende Sichel", "Erstes Viertel", "Zunehmender Mond",
//						"Vollmond", "Abnehmender Mond", "Letztes Viertel", "Abnehmende Sichel" };

	String phases[8] { "New Moon", "Waxing Crescent Moon", "First Quarter Moon", "Waxing Gibbous Moon,"
	                   "Full Moon","Waning Gibbous Moon",  "Third Quarter Moon", "Waning Crescent Moon"};



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
	void setInput(as_date, as_time);
	String GetAll() {return m_os;}
	double GetLat() {return m_Lat;}
	double GetLon() {return m_Lon;}
	String GetDate() {return m_Date;}
	String GetTime() {return m_Time;}
	double GetJD() {return m_JD;}
	double GetZone() {return m_Zone;}
	double GetDeltaT() {return m_DeltaT;}
	double GetSunDistance() {return m_SunDistance;}
	double GetSunDistanceObserver() {return m_SunDistanceObserver;}
	double GetSunLon() {return m_SunLon;}
	double GetSunDec() {return m_SunDec;}
	double GetSunAz() {return m_SunAz;}
	double GetSunAlt() {return m_SunAlt;}
	double GetSunDiameter() {return m_SunDiameter;}
	double GetMoonDistance() {return m_MoonDistance; }
	double GetMoonDistanceObserver() {return m_MoonDistanceObserver;}
	double GetMoonLon() {return m_MoonLon;}
	double GetMoonLat() {return m_MoonLat;}
	double GetMoonDec() {return m_MoonDec;}
	double GetMoonAz() {return m_MoonAz;}
	double GetMoonAlt() {return m_MoonAlt;}
	double GetMoonDiameter() {return m_MoonDiameter;}
	double GetMoonPhaseNumber() {return m_MoonPhaseNumber;}
	double GetMoonAge() {return m_MoonAge;}
	double GetSunAstronomicalTwilightMorning() { return m_SunAstronomicalTwilightMorning.TotalHour;}
	double GetSunNauticalTwilightMorning() { return m_SunNauticalTwilightMorning.TotalHour;}
	double GetSunCivilTwilightMorning() { return m_SunCivilTwilightMorning.TotalHour;}
	double GetSunRise() { return m_SunRise.TotalHour;}
	double GetSunTransit() { return m_SunTransit.TotalHour;}
	double GetSunSet() { return m_SunSet.TotalHour;}
	double GetSunCivilTwilightEvening() { return m_SunCivilTwilightEvening.TotalHour;}
	double GetSunNauticalTwilightEvening() { return m_SunNauticalTwilightEvening.TotalHour;}
	double GetSunAstronomicalTwilightEvening() { return m_SunAstronomicalTwilightEvening.TotalHour;}
	double GetMoonRA() { return m_MoonRA.TotalHour;}
	double GetMoonRise() { return m_MoonRise.TotalHour;}
	double GetMoonTransit() { return m_MoonTransit.TotalHour;}
	double GetMoonSet() { return m_MoonSet.TotalHour;}
	String GetSunAstronomicalTwilightMorning_s() {return m_SunAstronomicalTwilightMorning.HHMMSS;}
	String GetSunNauticalTwilightMorning_s() {return m_SunNauticalTwilightMorning.HHMMSS;}
	String GetSunCivilTwilightMorning_s() {return m_SunCivilTwilightMorning.HHMMSS;}
	String GetSunRise_s() {return m_SunRise.HHMMSS;}
	String GetSunTransit_s() {return m_SunTransit.HHMMSS;}
	String GetSunSet_s() {return m_SunSet.HHMMSS;}
	String GetSunCivilTwilightEvening_s() {return m_SunCivilTwilightEvening.HHMMSS;}
	String GetSunNauticalTwilightEvening_s() {return m_SunNauticalTwilightEvening.HHMMSS;}
	String GetSunAstronomicalTwilightEvening_s() {return m_SunAstronomicalTwilightEvening.HHMMSS;}
	String GetMoonRA_s() {return m_MoonRA.HHMMSS;}
	String GetMoonRise_s() {return m_MoonRise.HHMMSS;}
	String GetMoonTransit_s() {return m_MoonTransit.HHMMSS;}
	String GetMoonSet_s() {return m_MoonSet.HHMMSS;}
	String GetMoonPhase() {return phases[(int)m_MoonPhase];}
	String GetMoonSign() {return TKZ[(int)m_MoonSign];}
	String GetSunSign() {return TKZ[(int)m_SunSign];}

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
