#include <iostream>
#include <fstream>
#include <string>
#include <string.h>

#include <TFile.h>
#include <TChain.h>

TChain *enchain(const char* treename, 
                const char* treetitle,
                const char* listfile=0)
{
   TChain *ch = new TChain(treename, treetitle);

   ifstream iflist;
   if (listfile != 0)
      iflist.open(listfile);
   else
      iflist.open("enchain.list");

   while (iflist.good()) {
      std::string ifline;
      getline(iflist,ifline);
      std::string rfile = ifline.substr(ifline.find_last_of(' ') + 1);
      if (rfile.size() > 0) {
         rfile = "rawdata/" + rfile;
         ch->AddFile(rfile.c_str());
      }
   }
   iflist.close();
   return ch;
}
