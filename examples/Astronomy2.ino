// preparation: avoid insufficient stack if You can see
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

