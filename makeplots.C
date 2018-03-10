#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>

void Vinde(const char* varname, double Vmin=-0.2, double Vmax=+0.2, 
           double tmin=1414930823)
{
   TTree *tsum = (TTree*)gROOT->FindObject("tsum");
   if (tsum == 0) {
      std::cout << "unable to load tsum tree, "
                   "please cd to the correct directory and try again."
                << std::endl;
      return;
   }

   gStyle->SetPadLeftMargin(0.15);
   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 550, 500);

   std::ostringstream varexpr;
   varexpr << varname << ":(t_s-" << tmin << ")/3600.";
   std::ostringstream cutexpr;
   cutexpr << varname << " > " << Vmin << " && " << varname << " < " << Vmax
           << " && t_s > " << tmin;
   tsum->Draw(varexpr.str().c_str(), cutexpr.str().c_str());

   TH1 *htemp = (TH1*)gROOT->FindObject("htemp");
   if (htemp != 0) {
      htemp->GetXaxis()->SetTitle("t (hr)");
      htemp->SetTitle("");
   }
}

void set_profile_errors(TProfile *prof)
{
   int nbins = prof->GetNbinsX();
   for (int i=1; i <= nbins; ++i) {
      double cont = prof->GetBinContent(i); 
      double err = prof->GetBinError(i);
      if (cont != 0 && err >= 0)
         prof->SetBinError(i, sqrt(1e-4 + cont * cont));
   }
}

void pairxscan(const char* varname1, const char* varname2)
{
   TFile *tscanfile = 0;
   TTree *tscantree = (TTree*)gROOT->FindObject("tscan");
   if (tscantree == 0) {
      tscanfile= new TFile("xbeamscan2-11-07_rad_2e-5.root");
      tscantree = (TTree*)gROOT->FindObject("tscan");
      if (tscantree == 0) {
         std::cout << "unable to load tscan tree, "
                      "please cd to the correct directory and try again."
                   << std::endl;
         return;
      }
   }

   double xlow = tscantree->GetMinimum("Xscan_bpu_mm");
   double xhigh = tscantree->GetMaximum("Xscan_bpu_mm");
   xlow = int((xlow + 100.) / 5) * 5 - 100.;
   xhigh = int((xhigh + 100.) / 5 + 1) * 5 - 100.;
   TProfile *pleft = new TProfile("pleft", tscantree->GetTitle(),
                                  100, xlow, xhigh, -11., 11.);
   TProfile *pright = new TProfile("pright", tscantree->GetTitle(),
                                  100, xlow, xhigh, -11., 11.);
   pleft->BuildOptions(-11., 11., "s");
   pright->BuildOptions(-11., 11., "s");

   gStyle->SetPadLeftMargin(0.15);
   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 550, 500);

   std::string varexpr1(varname1);
   varexpr1 += ":Xscan_bpu_mm";
   varexpr1 += ">>pleft";
   std::string varexpr2(varname2);
   varexpr2 += ":Xscan_bpu_mm";
   varexpr2 += ">>pright";
   std::string cutexpr("t_s<tscan_s+5");
   tscantree->Draw(varexpr1.c_str(), cutexpr.c_str());
   tscantree->Draw(varexpr2.c_str(), cutexpr.c_str());

   set_profile_errors(pleft);
   set_profile_errors(pright);

   pleft->GetXaxis()->SetTitle("x_bpu_mean (mm)");
   pleft->GetYaxis()->SetTitle("inner x wedges (V)");
   pleft->GetYaxis()->SetTitleOffset(1.5);
   pleft->SetStats(0);
   pleft->SetLineColor(kRed);
   pright->SetLineColor(kBlue);
   pleft->Draw(); 
   pright->Draw("same"); 
}

void pairyscan(const char* varname1, const char* varname2)
{
   TFile *tscanfile = 0;
   TTree *tscantree = (TTree*)gROOT->FindObject("tscan");
   if (tscantree == 0) {
      tscanfile= new TFile("ybeamscan2-11-07_rad_2e-5.root");
      tscantree = (TTree*)gROOT->FindObject("tscan");
      if (tscantree == 0) {
         std::cout << "unable to load tscan tree, "
                      "please cd to the correct directory and try again."
                   << std::endl;
         return;
      }
   }

   double xlow = tscantree->GetMinimum("Xscan_bpu_mm");
   double xhigh = tscantree->GetMaximum("Xscan_bpu_mm");
   xlow = int((xlow + 100.) / 5) * 5 - 100.;
   xhigh = int((xhigh + 100.) / 5 + 1) * 5 - 100.;
   TProfile *pleft = new TProfile("pleft", tscantree->GetTitle(),
                                  150, xlow, xhigh, -11., 11.);
   TProfile *pright = new TProfile("pright", tscantree->GetTitle(),
                                  150, xlow, xhigh, -11., 11.);
   pleft->BuildOptions(-11., 11., "s");
   pright->BuildOptions(-11., 11., "s");

   gStyle->SetPadLeftMargin(0.15);
   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 550, 500);

   pleft->GetXaxis()->SetTitle("y_bpu_mean (mm)");
   std::string varexpr1(varname1);
   varexpr1 += ":Yscan_bpu_mm";
   varexpr1 += ">>pleft";
   std::string varexpr2(varname2);
   varexpr2 += ":Yscan_bpu_mm";
   varexpr2 += ">>pright";
   std::string cutexpr("t_s<tscan_s+5");
   tscantree->Draw(varexpr1.c_str(), cutexpr.c_str());
   tscantree->Draw(varexpr2.c_str(), cutexpr.c_str());

   pleft->GetYaxis()->SetTitle("inner y wedges (V)");
   pleft->GetYaxis()->SetTitleOffset(1.5);
   pleft->SetStats(0);
   pleft->Draw(); 
   pright->Draw("same"); 
}

void Xasym(const char* varname, const char* cutstr=0)
{
   TFile *tscanfile = 0;
   TTree *tscantree = (TTree*)gROOT->FindObject("tscan");
   if (tscantree == 0) {
      tscanfile= new TFile("xbeamscan2-11-07_rad_2e-5.root");
      tscantree = (TTree*)gROOT->FindObject("tscan");
      if (tscantree == 0) {
         std::cout << "unable to load tscan tree, "
                      "please cd to the correct directory and try again."
                   << std::endl;
         return;
      }
   }

   TProfile *xasym = new TProfile("xasym", tscantree->GetTitle(),
                                  200, -1.0, 1.0, -30., 30.);
   xasym->BuildOptions(-30., 30., "s");

   gStyle->SetPadLeftMargin(0.15);
   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 550, 500);

   std::string varexpr(varname);
   varexpr = "Xscan_bpu_mm:" + varexpr;
   varexpr += ">>xasym";
   std::string cutexpr("t_s<tscan_s+5");
   if (cutstr != 0) {
      cutexpr += "&&";
      cutexpr += cutstr;
   }
   tscantree->Draw(varexpr.c_str(), cutexpr.c_str());

   xasym->GetYaxis()->SetTitle("x_bpu_mean (mm)");
   xasym->GetXaxis()->SetTitle("V asymmetry");
   xasym->GetYaxis()->SetTitleOffset(1.5);
   xasym->SetStats(0);
   xasym->Draw(); 
}

void Yasym(const char* varname, const char* cutstr=0)
{
   TFile *tscanfile = 0;
   TTree *tscantree = (TTree*)gROOT->FindObject("tscan");
   if (tscantree == 0) {
      tscanfile= new TFile("ybeamscan2-11-07_rad_2e-5.root");
      tscantree = (TTree*)gROOT->FindObject("tscan");
      if (tscantree == 0) {
         std::cout << "unable to load tscan tree, "
                      "please cd to the correct directory and try again."
                   << std::endl;
         return;
      }
   }

   TProfile *yasym = new TProfile("yasym", tscantree->GetTitle(),
                                  200, -1.0, 1.0, -15., 15.);
   yasym->BuildOptions(-15., 15., "s");

   gStyle->SetPadLeftMargin(0.15);
   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 550, 500);

   std::string varexpr(varname);
   varexpr = "Yscan_bpu_mm:" + varexpr;
   varexpr += ">>yasym";
   std::string cutexpr("t_s<tscan_s+5");
   if (cutstr != 0) {
      cutexpr += "&&";
      cutexpr += cutstr;
   }
   tscantree->Draw(varexpr.c_str(), cutexpr.c_str());

   yasym->GetYaxis()->SetTitle("y_bpu_mean (mm)");
   yasym->GetXaxis()->SetTitle("V asymmetry");
   yasym->GetYaxis()->SetTitleOffset(1.5);
   yasym->SetStats(0);
   yasym->Draw(); 
}

void yinfitxscan()
{
   double fitpars[6] = {-2.61809,
                        8.37568,
                        -1.30637,
                        1.47959,
                        0.615187,
                        31.0245};

   TFile *tscanfile = 0;
   TTree *tscantree = (TTree*)gROOT->FindObject("tscan");
   if (tscantree == 0) {
      tscanfile= new TFile("beam_0018.root");
      tscantree = (TTree*)gROOT->FindObject("tscan");
      if (tscantree == 0) {
         std::cout << "unable to load the tscan tree, "
                      "please cd to the correct directory and try again."
                   << std::endl;
         return;
      }
   }

   TProfile *yinfit = new TProfile("yinfit", tscantree->GetTitle(),
                                  50, -15., 15., -20., 20.);
   yinfit->BuildOptions(-20., 20., "s");

   gStyle->SetPadLeftMargin(0.15);
   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 550, 500);

   std::string yasym("(Viyp-Viym)/(Viyp+Viym)");
   std::stringstream varexpr;
   int n;
   for (n=0; n < 5; ++n)
      varexpr << "(";
   varexpr << fitpars[n--];
   for (; n >= 0; )
      varexpr << ")*" << yasym << "+" << fitpars[n--];
   varexpr << ":Xscan_bpu_mm>>yinfit";
   tscantree->Draw(varexpr.str().c_str());

   yinfit->GetXaxis()->SetTitle("111.7 - x_col_motor (mm)");
   yinfit->GetYaxis()->SetTitle("y from acol inner y-wedges (mm)");
   yinfit->GetYaxis()->SetTitleOffset(1.3);
   yinfit->SetStats(0);
   yinfit->Draw(); 
}

void youtfitxscan()
{
   double fitpars[6] = {-2.68086,
                        13.9037,
                        -0.712466,
                        15.3233,
                        -0.758043,
                        4.94836};

   TFile *tscanfile = 0;
   TTree *tscantree = (TTree*)gROOT->FindObject("tscan");
   if (tscantree == 0) {
      tscanfile= new TFile("beam_0018.root");
      tscantree = (TTree*)gROOT->FindObject("tscan");
      if (tscantree == 0) {
         std::cout << "unable to load the tscan tree, "
                      "please cd to the correct directory and try again."
                   << std::endl;
         return;
      }
   }

   TProfile *youtfit = new TProfile("youtfit", tscantree->GetTitle(),
                                  50, -15., 15., -20., 20.);
   youtfit->BuildOptions(-20., 20., "s");

   gStyle->SetPadLeftMargin(0.15);
   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 550, 500);

   std::string yasym("(Voyp-Voym)/(Voyp+Voym)");
   std::stringstream varexpr;
   int n;
   for (n=0; n < 5; ++n)
      varexpr << "(";
   varexpr << fitpars[n--];
   for (; n >= 0; )
      varexpr << ")*" << yasym << "+" << fitpars[n--];
   varexpr << ":Xscan_bpu_mm>>youtfit";
   tscantree->Draw(varexpr.str().c_str());

   youtfit->GetXaxis()->SetTitle("111.7 - x_col_motor (mm)");
   youtfit->GetYaxis()->SetTitle("y from acol outer y-wedges (mm)");
   youtfit->GetYaxis()->SetTitleOffset(1.3);
   youtfit->SetStats(0);
   youtfit->Draw(); 
}

void plotratescan(const char* filename, int ndatum)
{
   std::string filenames(filename);
   filenames = "rawdata/RateScan/" + filenames;
   std::ifstream ifs(filenames.c_str());
   if (! ifs.good())
      return;

   double datum[999];
   int npoints;
   for (npoints=0; ifs.good(); ++npoints) {
      double data[5];
      ifs >> data[0] >> data[1] >> data[2] >> data[3] >> data[4];
      datum[npoints] = data[ndatum];
   }
   ifs.close();
   --npoints;

   TProfile *pscan = new TProfile("pscan", filename, npoints, 0., npoints);
   for (int point=0; point < npoints; ++point)
      pscan->Fill(point, datum[point]);

   set_profile_errors(pscan);

   gStyle->SetPadLeftMargin(0.15);
   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 550, 500);
   pscan->SetStats(0);
   pscan->Draw();
}
