//
// timescan.C - a utility function (C++) to convert a free-format date/time
//              string into unix epoch time.
//
// author: richard.t.jones at uconn.edu
// version: november 8, 2014
//
// usage:
//    long epochtime = timescan(tstring);
//
// where any of the following examples for tstring would work.
//     "Jan 4, 2012 07:54:18"
//     "07:54:18 Fri Jan 4 2012"
//     "07:54:18 Fri Jan 4 2012"
//     "Fri Jan  4, 2012 at 07:54:18"
//     "January 4, 2012 at 07:54:18"
// If the algorithm is unable to find the required information in the string
// then the return value is -1, otherwise it returns the unix epoch time.

#include <stdio.h>
#include <string.h>
#include <time.h>

#include <iostream>
#include <string>
#include <map>

std::map< std::string, int > monthshort;

std::map< std::string, int > monthlong;

const char delims[] = " \t,;";

void initialize() {
   monthshort["Jan"] = 0;
   monthshort["Feb"] = 1;
   monthshort["Mar"] = 2;
   monthshort["Apr"] = 3;
   monthshort["May"] = 4;
   monthshort["Jun"] = 5;
   monthshort["Jul"] = 6;
   monthshort["Aug"] = 7;
   monthshort["Sep"] = 8;
   monthshort["Oct"] = 9;
   monthshort["Nov"] = 10;
   monthshort["Dec"] = 11;
   monthlong["January"] = 0;
   monthlong["February"] = 1;
   monthlong["March"] = 2;
   monthlong["April"] = 3;
   monthlong["May"] = 4;
   monthlong["June"] = 5;
   monthlong["July"] = 6;
   monthlong["August"] = 7;
   monthlong["September"] = 8;
   monthlong["October"] = 9;
   monthlong["November"] = 10;
   monthlong["December"] = 11;
}

time_t timescan(const std::string tstring)
{
   char *tstringcopy = new char[tstring.size()];
   char *tstringtoken;
   initialize();

   struct tm tstruct;

   // extract the hour:minute:second
   strcpy(tstringcopy, tstring.c_str());
   for (tstringtoken = strtok(tstringcopy, delims);
        tstringtoken != 0; tstringtoken = strtok(0, delims)) 
   {
      if (sscanf(tstringtoken, "%d:%d:%d", &tstruct.tm_hour,
                 &tstruct.tm_min, &tstruct.tm_sec) == 3)
         break;
   }
   if (tstringtoken == 0) {
      delete [] tstringcopy;
      std::cerr << "hour:minute:second not found!" << std::endl;
      return -1;
   }

   // extract the month day
   strcpy(tstringcopy, tstring.c_str());
   for (tstringtoken = strtok(tstringcopy, delims);
        tstringtoken != 0; tstringtoken = strtok(0, delims)) 
   {
      char monstr[10];
      if (strlen(tstringtoken) < 10 &&
          sscanf(tstringtoken, "%s", monstr) == 1)
      {
         if (monthlong.find(monstr) != monthlong.end()) {
            tstruct.tm_mon = monthlong[monstr];
            tstringtoken = strtok(0, delims);
            if (sscanf(tstringtoken, "%d", &tstruct.tm_mday) != 0)
               break;
         }
         else if (monthshort.find(monstr) != monthshort.end()) {
            tstruct.tm_mon = monthshort[monstr];
            tstringtoken = strtok(0, delims);
            if (sscanf(tstringtoken, "%d", &tstruct.tm_mday) != 0)
               break;
         }
      }
   }
   if (tstringtoken == 0) {
      delete [] tstringcopy;
      std::cerr << "month mday found!" << std::endl;
      return -1;
   }

   // extract the year
   strcpy(tstringcopy, tstring.c_str());
   for (tstringtoken = strtok(tstringcopy, delims);
        tstringtoken != 0; tstringtoken = strtok(0, delims)) 
   {
      if (sscanf(tstringtoken, "%d", &tstruct.tm_year) == 1 &&
          tstruct.tm_year > 2000 && tstruct.tm_year < 2100)
            break;
   }
   if (tstringtoken == 0) {
      delete [] tstringcopy;
      std::cerr << "year not found!" << std::endl;
      return -1;
   }
   tstruct.tm_year -= 1900;
   tstruct.tm_isdst = -1;

   delete [] tstringcopy;
   return mktime(&tstruct);
}

#ifndef COMPILE_TEST_MAIN
int main()
{
   std::string tstring("Time:     Fri Nov  7 03:08:09 2014");
   time_t timecode = timescan(tstring);
   struct tm *tstruct = localtime(&timecode);
   std::cout << "\"" << tstring << "\" returns "  << timecode
             << " which translates back to local time as "
             << ((tstruct->tm_isdst == 1)? "DST " : "EST ") 
             << asctime(tstruct)
             << std::endl;
   return 0;
}
#endif
