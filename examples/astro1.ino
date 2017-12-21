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
