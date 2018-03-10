#include <iostream>
#include <acol_image.h>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>

acol_image *img=0;

void paint_scan()
{
   TTree *tscantree = (TTree*)gROOT->FindObject("tscan");
   if (tscantree == 0) {
      std::cerr << "Error - unable to find input tree tscan." << std::endl;
      return;
   }

   if (img == 0)
      img = new acol_image();

   struct tscan_t {
      double tscan_s;
      double Iscan_nA;
      double Xscan_bpu_mm;
      double Yscan_bpu_mm;
      double Vscan_ixp;
      double Vscan_ixm;
      double Vscan_iyp;
      double Vscan_iym;
      double Vscan_oxp;
      double Vscan_oxm;
      double Vscan_oyp;
      double Vscan_oym;
      double ATrate_Hz;
   } tscan;
   tscantree->SetBranchAddress("tscan", &tscan.tscan_s);

   int nentries = tscantree->GetEntries();
   double x[nentries];
   double y[nentries];
   for (int entry=0; entry < nentries; ++entry) {
      tscantree->GetEntry(entry);
      x[entry] = tscan.Xscan_bpu_mm + 1;
      y[entry] = tscan.Yscan_bpu_mm + 4;
   }
   img->PaintStripe(nentries, x, y);
}

void paint_scans()
{
   TString rfile[10];
   rfile[0] = "beam_0026.root";
   rfile[1] = "beam_0027.root";
   rfile[2] = "beam_0028.root";
   rfile[3] = "beam_0029.root";
   rfile[4] = "beam_0030.root";
   rfile[5] = "beam_0031.root";
   rfile[6] = "beam_0032.root";
   rfile[7] = "beam_0033.root";
   rfile[8] = "beam_0035.root";
   rfile[9] = "beam_0036.root";
   for (int i=0; i < 10; ++i) {
      TFile *fi = new TFile(rfile[i]);
      paint_scan();
      delete fi;
   }
}
