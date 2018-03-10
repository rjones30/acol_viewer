#include <iostream>
#include <acol_viewer.h>
#include <TSystem.h>
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>

#include <stdio.h>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>

#include <stdlib.h> 
#include <math.h>
#include <time.h>

const char statedir[] = "/gluonfs1/home/hdops/acol_viewer/state/"; 

float sum[8] = {0,0,0,0,0,0,0,0};
float pedestal[8] = {0,0,0,0,0,0,0,0};
float gain[8] = {0,0,0,0,0,0,0,0};
float noise[8] = {0,0,0,0,0,0,0,0};

const int IXP = acol_viewer::kIXplus;
const int IXM = acol_viewer::kIXminus;
const int IYP = acol_viewer::kIYplus;
const int IYM = acol_viewer::kIYminus;
const int OXP = acol_viewer::kOXplus;
const int OXM = acol_viewer::kOXminus;
const int OYP = acol_viewer::kOYplus;
const int OYM = acol_viewer::kOYminus;

void play_acol_viewer(int verbosity=0) {

   // initialize the viewer
   acol_viewer *viewer = new acol_viewer;
   acol_viewer::record record;
   double current_nA = 99.;

   record.wIXplus = 0;
   record.wIXminus = 0;
   record.wIYplus = 0;
   record.wIYminus = 0;
   record.wOXplus = 0;
   record.wOXminus = 0;
   record.wOYplus = 0;
   record.wOYminus = 0;
   viewer->setZero(record);

   int ok = 1;
   while (ok) {
      std::string cmd("./epics_pv.sh > ");
      std::string payload(statedir);
      payload += "latest_epics_payload.txt";
      cmd += payload;
      system(cmd.c_str());
      std::ifstream infile(payload.c_str());
      int Vlevel[8] = {0,0,0,0,0,0,0,0};
      std::string datestamp;
      std::string timestamp;
      if (infile.good()) {
         std::string line;
         while (getline(infile, line)) {
            std::stringstream ss(line);
            std::string varname;
            double value;
            ss >> varname >> datestamp >> timestamp >> value;
            if (varname.compare("IBCAD00CRCUR6") == 0)
               current_nA = value;
            else if (varname.compare("IOCHDCOL:VMICADC1_1") == 0)
               Vlevel[0] = value;
            else if (varname.compare("IOCHDCOL:VMICADC2_1") == 0)
               Vlevel[1] = value;
            else if (varname.compare("IOCHDCOL:VMICADC3_1") == 0)
               Vlevel[2] = value;
            else if (varname.compare("IOCHDCOL:VMICADC4_1") == 0)
               Vlevel[3] = value;
            else if (varname.compare("IOCHDCOL:VMICADC1_2") == 0)
               Vlevel[4] = value;
            else if (varname.compare("IOCHDCOL:VMICADC2_2") == 0)
               Vlevel[5] = value;
            else if (varname.compare("IOCHDCOL:VMICADC3_2") == 0)
               Vlevel[6] = value;
            else if (varname.compare("IOCHDCOL:VMICADC4_2") == 0)
               Vlevel[7] = value;
            else
               std::cerr << "warning: unknown epics variable "
                         << varname << " found in epics record."
                         << std::endl;
         }
      }
      infile.close();

      struct tm timest;
      float tsec;
      sscanf(datestamp.c_str(), "%i-%i-%i", &timest.tm_year, &timest.tm_mon,
                                                             &timest.tm_mday);
      sscanf(timestamp.c_str(), "%i:%i:%f", &timest.tm_hour, &timest.tm_min,
                                                             &tsec);
      timest.tm_sec = int(tsec);
      record.tnsec = (tsec - timest.tm_sec) * 1000000000;
      timest.tm_mon -= 1;
      timest.tm_year -= 1900;
      timest.tm_isdst = -1;
      record.tsec = mktime(&timest);

      payload = statedir;
      payload += "calibration.txt";
      infile.open(payload.c_str());
      if (infile.good()) {
         std::string line;
         getline(infile, line);
         for (int n=0; n < 8; ++n) {
            getline(infile, line);
            std::stringstream ss(line);
            double ped;
            ss >> gain[n] >> ped >> noise[n];
            pedestal[n] = (pedestal[n] == 0)? ped : pedestal[n];
         }
         infile.close();
      }

      record.wIXplus = noise[0];
      record.wIXminus = noise[1];
      record.wIYplus = noise[2];
      record.wIYminus = noise[3];
      record.wOXplus = noise[4];
      record.wOXminus = noise[5];
      record.wOYplus = noise[6];
      record.wOYminus = noise[7];
      viewer->setNoise(record);

      record.wIXplus = gain[0] * Vlevel[0] * (0.005/16) - pedestal[0];
      record.wIXminus = gain[1] * Vlevel[1] * (0.005/16) - pedestal[1];
      record.wIYplus = gain[2] * Vlevel[2] * (0.005/16) - pedestal[2];
      record.wIYminus = gain[3] * Vlevel[3] * (0.005/16) - pedestal[3];
      record.wOXplus = gain[4] * Vlevel[4] * (0.005/16) - pedestal[4];
      record.wOXminus = gain[5] * Vlevel[5] * (0.005/16) - pedestal[5];
      record.wOYplus = gain[6] * Vlevel[6] * (0.005/16) - pedestal[6];
      record.wOYminus = gain[7] * Vlevel[7] * (0.005/16) - pedestal[7];
      viewer->setFullScale(record);
      viewer->display(record);
      if (verbosity > 0)
         std::cout << record.wIXplus << ", "
                   << record.wIXminus << ", "
                   << record.wIYplus << ", "
                   << record.wIYminus << ", "
                   << record.wOXplus << ", "
                   << record.wOXminus << ", "
                   << record.wOYplus << ", "
                   << record.wOYminus << std::endl;

      if (current_nA == 0) {
         float *Vlevel = &record.wIXplus;
         for (int n=0; n < 8; ++n) {
            pedestal[n] += Vlevel[n] / 20;
         }
      }
      sleep(1);
   }
}

