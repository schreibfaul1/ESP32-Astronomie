# ESP32-Astronomie
Compute astronomical data from given geoaraphic coordinates, date and time
Library zum Bau einer astronomischen Uhr mit ESP32. Berechnet für einen geografischen Ort für die angegebene Zeit den Sonnenaufgang, Sonnenuntergang, die Mondphase, das Tierkreiszeichen u.v.m.
Das Datum muss zwischen 1.3.1901 und 28.2.2100 liegen. Eine ausführliche Beschreibung der verwendeten Formeln ist hier zu finden: http://lexikon.astronomie.info/java/sunmoon/
Beispielprogramm: 

````c++
#include <Arduino.h>
#include "astronomy.h"

as_geo Berlin ={ 13.40, 52.52,  1}; // longitude, latitude, timezone
as_geo NewYork={-74.00, 40.71, -5};
as_geo Moscow ={ 37.61, 55.75,  3};

Astronomy astro(Berlin);

char buf[100];

void setup(){
	Serial.begin(115200);

	//Set the Date/Time (UTC)
    as_date d={  1, 12, 2017}; // day, month, year
	as_time t={ 10, 00,   00}; // hour, minute, second
	astro.setInput(d, t);
	Serial.println(astro.GetAll());

}

void loop()
{
//Add your repeated code here
}
````

Ausgabe:
Oestl. geografische Laenge: ........... 13.40 Grad
Geografische Breite: .................. 52.52 Grad
Datum/Zeit: ........................... 01.12.2017 10:00:00
Zeitdifferenz zu Weltzeit: ............ 1.00 h
deltaT: ............................... 65.00 sek
Julianisches Datum: ................... 2458088.88 Tage
Greenwich Sternzeit  GMST: ............ 13:41:39 = 13.6942 h
Lokale Sternzeit LMST: ................ 14:35:15 = 14.5875 h
Entfernung der Sonne (Erdmittelpunkt):  147498071.20 km
Entfernung der Sonne (vom Beobachter):  147496787.00 km
Eklipt. Laenge der Sonne: ............. 249.40 Grad
Rektaszension der Sonne: .............. 16:30:52 = 16.5144 h
Deklination der Sonne: ................ -21.86 Grad
Azimut der Sonne: ..................... 152.76 Grad
Hoehe der Sonne ueber Horizont: ....... 11.90 Grad
Durchmesser der Sonne: ................ 32.44 '
Astronomische Morgendaemmerung: ....... 05:49:54 = 5.83 h
Nautische Morgendaemmerung: ........... 06:31:06 = 6.52 h
Buergerliche Morgendaemmerung: ........ 07:14:30 = 7.24 h
Sonnenaufgang: ........................ 07:54:44 = 7.91 h
Sonnenkulmination: .................... 11:55:39 = 11.93 h
Sonnenuntergang: ...................... 15:56:12 = 15.94 h
Buergerliche Abenddaemmerung: ......... 16:36:26 = 16.61 h
Nautische Abenddaemmerung: ............ 17:19:49 = 17.33 h
Astronomische Abenddaemmerung: ........ 18:00:59 = 18.02 h
Tierkreiszeichen: ..................... Schuetze
Entfernung des Mondes (Erdmittelpunkt): 371229.90 km
Entfernung des Mondes (vom Beobachter): 374301.30 km
Eklipt. Laenge des Mondes: ............ 37.45 Grad
Eklipt. Breite des Mondes: ............ -5.05 Grad
Rektaszension des Mondes: ............. 02:27:07 = 2.45 h
Deklination des Mondes: ............... 8.36 Grad
Azimut des Mondes: .................... 2.30 Grad
Hoehe des Mondes ueber Horizont: ...... -29.10 Grad
Durchmesser des Mondes: ............... 32.19 '
Mondaufgang: .......................... 15:09:53 = 15.16 h
Mondkulmination: ...................... 22:18:21 = 22.31 h
Monduntergang: ........................ 04:21:56 = 4.37 h
Mondphase: ............................ 0.92
Mondalter: ............................ 148.10 Grad
Mondphase: ............................ Zunehmender Mond
Mondzeichen: .......................... Stier

Weiteres Beispiel, Berechnung des Sonnenaufgangs und Sonnenuntergangs für Berlin vom 01.Dez-11.Det. 2017;
````c++
#include <Arduino.h>
#include "astronomy.h"

as_geo Berlin ={ 13.40, 52.52,  1}; // longitude, latitude, timezone
as_geo NewYork={-74.00, 40.71, -5};
as_geo Moscow ={ 37.61, 55.75,  3};

Astronomy astro(Berlin);

char buf[100];

void setup(){
	Serial.begin(115200);

	//Set the Date/Time (UTC)
    as_date d={ 1, 12, 2017}; // day, month, year
	as_time t;

	Serial.println("\n Berlin, Dezember 2017");
	Serial.println(" Tag Sonnenaufgang Sonnenuntergang");
	for(int i=1; i<12; i++){
		d.day=i;
		astro.setInput(d, t);
		sprintf(buf, " %02d    %s        %s\n", i, astro.GetSunRise_s().c_str(), astro.GetSunSet_s().c_str());
		Serial.print(buf);
	}
}

void loop()
{
//Add your repeated code here
}
````




