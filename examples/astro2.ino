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
    as_date d={  1, 12, 2017}; // day, month, year
	as_time t={ 10, 00,   00}; // hour, minute, second
	astro.setInput(d, t);
	Serial.println(astro.GetAll());

}

void loop()
{
//Add your repeated code here
}
