//
// show_acoltree.C - root script to open a chain of root files generated
//                   by Yi Qiang containing 8 separate trees, one for each
//                   sector of the active collimator, and condense these
//                   into a single chain with one voltage value for each
//                   sector computed as the average of the 8192 samples
//                   in the original tree.
//
// author: richard.t.jones at uconn.edu
// version: november 7, 2014
//
// Usage:
// The following sequence of commands will create a chain of the files
// listed in standard "ls -l" format in enchain.list and condense them
// into an acoltree tree, and then generate a time series plot of the
// voltage level from one of the wedges, identified in the string argument
// to the show_acoltree function.
//   root > .L enchain.C+g
//   root > .L acoltree.C+g
//   root > .L show_acoltree.C+g
//   root > show_acoltree("oyp");
// A side effect of the above command is to create a new root file called
// acoltree.root in the working directory and save the condensed tree in
// it, called tsum. If acoltree.root already exists in the working directory
// then the tsum tree is read from the file instead of being regenerated
// by a scan over the input chain.

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TChain.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>

static TCanvas *c1 = 0;

#include <acoltree.h>

TChain *enchain(const char* treename, 
                const char* treetitle,
                const char* listfile=0);
static bool initialized = true;

TH1 *show_acoltree(const char *wedge, double t0 = 1416163221, 
                                      double t1 = 1416163821)
{
   if (! initialized) {
      gSystem->CompileMacro("enchain.C");
      gSystem->CompileMacro("acoltree.C");
      initialized = true;
   }

   TFile tsumfile("acoltree.root", "update");

   struct tsum_t {
      double t_s;
      double I_nA;
      double Vixp;
      double Vixm;
      double Viyp;
      double Viym;
      double Voxp;
      double Voxm;
      double Voyp;
      double Voym;
   } tsum;

   TTree *tsumtree = (TTree*)gROOT->FindObject("tsum");
   if (tsumtree == 0) {
      TChain *rawch[8];
      rawch[0] = enchain("N9:raw_XP", "inner X+ wedge");
      rawch[1] = enchain("N9:raw_XM", "inner X- wedge");
      rawch[2] = enchain("N9:raw_YP", "inner Y+ wedge");
      rawch[3] = enchain("N9:raw_YM", "inner Y- wedge");
      rawch[4] = enchain("N10:raw_XP", "outer X+ wedge");
      rawch[5] = enchain("N10:raw_XM", "outer X- wedge");
      rawch[6] = enchain("N10:raw_YP", "outer Y+ wedge");
      rawch[7] = enchain("N10:raw_YM", "outer Y- wedge");
      TTree *sumtr[8];
      for (int n=0; n < 8; ++n) {
         std::string tname("h_");
         tname += rawch[n]->GetName();
         sumtr[n] = (TTree*)gROOT->FindObject(tname.c_str());
         while (sumtr[n] == 0) {
            acoltree rawtr(rawch[n]);
            rawtr.Loop();
            sumtr[n] = (TTree*)gROOT->FindObject(tname.c_str());
         }
      }

      struct acoltree_t {
         double t_s;
         double I_nA;
         double V_average;
      } sums[8];
      for (int n=0; n < 8; ++n) {
         sumtr[n]->SetBranchAddress("def", &sums[n].t_s);
         sums[n].t_s = 0;
      }
      tsumtree = new TTree("tsum", "active collimator summary");
      tsumtree->Branch("tsum", &tsum.t_s, "t_s/D:I_nA/D"
                                          ":Vixp/D:Vixm/D:Viyp/D:Viym/D"
                                          ":Voxp/D:Voxm/D:Voyp/D:Voym/D",
                       32768);

      Long64_t jentry[8];
      for (int n=1; n < 8; ++n)
         jentry[n] = 0;
      Long64_t nentries = sumtr[0]->GetEntries();
      for (jentry[0]=0; jentry[0] < nentries; ++jentry[0]) {
         if (sumtr[0]->GetEntry(jentry[0]) > 0) {
            int n;
            for (n=1; n < 8; ++n) {
               while (sums[n].t_s < sums[0].t_s)
                  if (sumtr[n]->GetEntry(++jentry[n]) <= 0) {
                     std::cout << "synchronization failed with tree " << n
                               << " at event " << jentry[0] << std::endl;
                     break;
                  }
               if (fabs(sums[0].t_s - sums[n].t_s) > 2.) {
                  std::cout << "time mismatch at event " << jentry[0]
                            << " in tree " << n << std::endl;
                  std::cout << " " << std::setprecision(10) << sums[0].t_s;
                  for (int nn=1; nn <= n; ++nn)
                     std::cout << ", " << std::setprecision(10) 
                               << sums[nn].t_s;
                  std::cout << std::endl;
                  break;
               }
            }
            if (n == 8) {
               tsum.t_s = sums[0].t_s;
               tsum.I_nA = sums[0].I_nA;
               tsum.Vixp = sums[0].V_average;
               tsum.Vixm = sums[1].V_average;
               tsum.Viyp = sums[2].V_average;
               tsum.Viym = sums[3].V_average;
               tsum.Voxp = sums[4].V_average;
               tsum.Voxm = sums[5].V_average;
               tsum.Voyp = sums[6].V_average;
               tsum.Voym = sums[7].V_average;
               tsumtree->Fill();
            }
         }
      }
      tsumtree->Write();
   }
   else {
      tsumtree->SetBranchAddress("tsum", &tsum.t_s);
   }
   
   double tlim[2] = {1.e20,0.};
   double Vlim[2] = {1.e20,-1.e20};
   double *Vsamp = &tsum.Vixp;
   Long64_t nentries = tsumtree->GetEntries();
   for (Long64_t n=0; n < nentries; ++n) {
      tsumtree->GetEntry(n);
      tlim[0] = (tsum.t_s < tlim[0])? tsum.t_s : tlim[0];
      tlim[1] = (tsum.t_s > tlim[1])? tsum.t_s : tlim[1];
      for (int i=0; i < 8; ++i) {
         Vlim[0] = (Vsamp[i] < Vlim[0])? Vsamp[i] : Vlim[0];
         Vlim[1] = (Vsamp[i] > Vlim[1])? Vsamp[i] : Vlim[1];
      }
   }
   tlim[0] = (t0 == 0)? tlim[0]/3600 : 0;
   tlim[1] = (t1 == 0)? tlim[1]/3600 : (t1-t0)/3600.;
   Vlim[0] = int(Vlim[0] * 10 - 9.99) / 10.;
   Vlim[1] = int(Vlim[1] * 10 + 9.99) / 10.;

   gStyle->SetPadLeftMargin(0.15);
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 550, 500);

   std::ostringstream varexpr;
   if (strncmp(wedge, "ixp", 3) == 0 || strncmp(wedge, "IXP", 3) == 0)
      varexpr << std::setprecision(12) << "Vixp:(t_s-" << t0 << ")/3600"
              << ">>hwedge(100," << tlim[0] << "," << tlim[1]
              << ",100," << Vlim[0] << "," << Vlim[1] << ")";
   else if (strncmp(wedge, "ixm", 3) == 0 || strncmp(wedge, "IXM", 3) == 0)
      varexpr << std::setprecision(12) << "Vixm:(t_s-" << t0 << ")/3600"
              << ">>hwedge(100," << tlim[0] << "," << tlim[1]
              << ",100," << Vlim[0] << "," << Vlim[1] << ")";
   else if (strncmp(wedge, "iyp", 3) == 0 || strncmp(wedge, "IYP", 3) == 0)
      varexpr << std::setprecision(12) << "Viyp:(t_s-" << t0 << ")/3600"
              << ">>hwedge(100," << tlim[0] << "," << tlim[1]
              << ",100," << Vlim[0] << "," << Vlim[1] << ")";
   else if (strncmp(wedge, "iym", 3) == 0 || strncmp(wedge, "IYM", 3) == 0)
      varexpr << std::setprecision(12) << "Viym:(t_s-" << t0 << ")/3600"
              << ">>hwedge(100," << tlim[0] << "," << tlim[1]
              << ",100," << Vlim[0] << "," << Vlim[1] << ")";
   else if (strncmp(wedge, "oxp", 3) == 0 || strncmp(wedge, "OXP", 3) == 0)
      varexpr << std::setprecision(12) << "Voxp:(t_s-" << t0 << ")/3600"
              << ">>hwedge(100," << tlim[0] << "," << tlim[1]
              << ",100," << Vlim[0] << "," << Vlim[1] << ")";
   else if (strncmp(wedge, "oxm", 3) == 0 || strncmp(wedge, "OXM", 3) == 0)
      varexpr << std::setprecision(12) << "Voxm:(t_s-" << t0 << ")/3600"
              << ">>hwedge(100," << tlim[0] << "," << tlim[1]
              << ",100," << Vlim[0] << "," << Vlim[1] << ")";
   else if (strncmp(wedge, "oyp", 3) == 0 || strncmp(wedge, "OYP", 3) == 0)
      varexpr << std::setprecision(12) << "Voyp:(t_s-" << t0 << ")/3600"
              << ">>hwedge(100," << tlim[0] << "," << tlim[1]
              << ",100," << Vlim[0] << "," << Vlim[1] << ")";
   else if (strncmp(wedge, "oym", 3) == 0 || strncmp(wedge, "OYM", 3) == 0)
      varexpr << std::setprecision(12) << "Voym:(t_s-" << t0 << ")/3600"
              << ">>hwedge(100," << tlim[0] << "," << tlim[1]
              << ",100," << Vlim[0] << "," << Vlim[1] << ")";
   else
      return 0;
   std::ostringstream condexpr;
   condexpr << std::setprecision(12) << "t_s > " << t0 << " && t_s < " << t1;
std::cout << "varexpr is " << varexpr.str() << std::endl;
std::cout << "condexpr is " << condexpr.str() << std::endl;
   tsumtree->Draw(varexpr.str().c_str(), condexpr.str().c_str());

   TH1 *hwedge = (TH1*)c1->GetPrimitive("hwedge");
   if (hwedge) {
      std::string name(wedge);
      name = "h" + name;
      hwedge->SetName(name.c_str());
      std::string title(wedge);
      title = "wedge " + title;
      hwedge->SetTitle(title.c_str());
      hwedge->GetYaxis()->SetTitle("1 second signal average (V)");
      hwedge->GetYaxis()->SetTitleOffset(2.0);
      hwedge->GetXaxis()->SetTitle("time (hr)");
      c1->Update();
   }

   return hwedge;
}

void show_all()
{
   show_acoltree("ixp");
   c1->Print("hixp.png");
   show_acoltree("ixm");
   c1->Print("hixm.png");
   show_acoltree("iyp");
   c1->Print("hiyp.png");
   show_acoltree("iym");
   c1->Print("hiym.png");
   show_acoltree("oxp");
   c1->Print("hoxp.png");
   show_acoltree("oxm");
   c1->Print("hoxm.png");
   show_acoltree("oyp");
   c1->Print("hoyp.png");
   show_acoltree("oym");
   c1->Print("hoym.png");
}
