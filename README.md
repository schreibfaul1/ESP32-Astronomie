# ESP32-Astronomie
Compute astronomical data from given geoaraphic coordinates, date and time.<br>
<br>
Library zum Bau einer astronomischen Uhr mit ESP32. Berechnet für einen geografischen Ort für die angegebene Zeit den Sonnenaufgang, Sonnenuntergang, die Mondphase, das Tierkreiszeichen u.v.m.
Das Datum muss zwischen 1.3.1901 und 28.2.2100 liegen. Eine ausführliche Beschreibung der verwendeten Formeln ist hier zu finden: http://lexikon.astronomie.info/java/sunmoon/ <br>
<br>
Library for the construction of an astronomical clock with ESP32. Calculates the sunrise, sunset, moon phase, zodiac sign, etc. for a given geographical location for the specified time. The date must be between 1.3.1901 and 28.2.2100. A detailed description of the formulas used can be found here: http://lexikon.astronomie.info/java/sunmoon/ <br>


### Example 1: 

````c++
// preparation: avoid insufficient stack
//
// Guru Meditation Error: Core  1 panic'ed (Unhandled debug exception)
// Debug exception reason: Stack canary watchpoint triggered (loopTask)
//
// double the stack in main.cpp      ...packages\esp32\hardware\esp32\version\cores\esp32\main.cpp
// set: xTaskCreatePinnedToCore(loopTask, "loopTask", 8192*2, NULL, 1, &loopTaskHandle, ARDUINO_RUNNING_CORE);
//
#include <Arduino.h>
#include "astronomy.h"

as_geo Berlin ={ 13.40, 52.52,  1}; // longitude, latitude, timezone
as_geo NewYork={-74.00, 40.71, -5};
as_geo Moscow ={ 37.61, 55.75,  3};

Astronomy astro(Berlin);

char buf[100];

void setup(){
    Serial.begin(115200);

    //Set the Date/Time
    as_date d={ 25, 01, 2019}; // day, month, year
    as_time t={ 10, 00,   00}; // hour, minute, second
    astro.setInput(d, t);
    Serial.println("Berlin, 25.Jan.2019 at 10 O'clock");
    Serial.println("=================================");
    Serial.println(astro.GetAll());

}

void loop()
{
//Add your repeated code here
}
````
Output:: <br>
![Ausgabe](https://github.com/schreibfaul1/ESP32-Astronomie/blob/master/images/astronomy1.jpg)
 <br>
### Example 2:
````c++
// double the stack in main.cpp
// set: xTaskCreatePinnedToCore(loopTask, "loopTask", 8192*2, NULL, 1, &loopTaskHandle, ARDUINO_RUNNING_CORE);
//
#include <Arduino.h>
#include "astronomy.h"

as_geo Berlin ={ 13.40, 52.52,  1}; // longitude, latitude, timezone
as_geo NewYork={-74.00, 40.71, -5};
as_geo Moscow ={ 37.61, 55.75,  3};

Astronomy astro(Berlin);

char buf[100];

void setup(){
    Serial.begin(115200);

    //Set the Date/Time
    as_date d={ 1, 1, 2018}; // day, month, year
    as_time t;

    Serial.println("\n Berlin, January 2018");
    Serial.println(" Day   Sunrise         Sunset  ");
    Serial.println(" ===   ========        ========");
    for(int i=1; i<11; i++){
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
Output: <br>
![Ausgabe](https://github.com/schreibfaul1/ESP32-Astronomie/blob/master/images/astronomy2.jpg)
<br>



