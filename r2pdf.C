//
// r2pdf.C - generates the PDF for the radial distribution of beam
//           photons incident on the face of the active collimator.
//
// author: richard.t.jones at uconn.edu
// version: november 25, 2014
//
// usage:
//  root > .L r2pdf.C+O
//  root > TF1 *f1 = new TF1("hr2", "r2pdf(x/100)", 0, 1e4);
//  root> f1->Draw("c");
//
// notes:
// 1) These functions assume length units in cm.
// 2) The currently implemented model fits well the distribution of beam
//    photons at the Hall D collimator that were generated from an amorphous
//    radiator at the nominal goniometer position.
//

#include <iostream>
#include <sstream>
#include <math.h>

#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TTree.h>
#include <TCanvas.h>

double r2pdf(double r2)
{
   // Return the unnormalized PDF, assuming r2 in cm^2
   // The normalization is set such that r2pdf(0) = 1.
 
   double z = r2 / 0.25;
   return  1 / (1 + 3.8 * z + 1.11 * pow(z, 1.9));
}

TH1D *r2dist(TH1D *htemp)
{
   // Generate dP/d[ln(z)] such that it can be overlaid
   // on the output from r2sim(), used for refining the model.

   TH1D *hfit = (TH1D*)htemp->Clone("hfit");
   Long64_t nbins = hfit->GetNbinsX();
   for (int n=1; n <= nbins; ++n) {
      double z = exp(hfit->GetXaxis()->GetBinCenter(n));
      double dPdlnz = z / (1 + 3.8 * z + 1.11 * pow(z, 1.9));
      hfit->SetBinContent(n, dPdlnz);
   }
   double height1 = htemp->GetMaximum();
   double height0 = hfit->GetMaximum();
   hfit->Scale(height1 / height0);
   hfit->SetLineColor(kRed);

   gStyle->SetPadLeftMargin(0.15);
   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 600, 500);
   hfit->Draw("same");
   return hfit;
}

TH1D *r2sim(TH1D* hsim=0)
{
   // Generate a projection onto the r2 axis of a tree containing
   // Monte Carlo generated beam impacts on the active collimator.
   // The Monte Carlo tree can be made by running hdgeant with the
   // BACKGROUND_PROFILING switch enabled in gustep, and then running
   //  $ h2root geant.hbook geant.root 1 1 65536
   // on the output file from a beam-generator simulation run.

   TTree *h10 = (TTree*)gROOT->FindObject("h10");
   if (h10 == 0) {
      std::cerr << "Error in r2sim: unable to find input tree h10,"
                << " please open geant.root and try again."
                << std::endl;
      return 0;
   }
   TString name("log(x[0]*x[0]+x[1]*x[1])-2*log(0.50)");
   if (hsim == 0)
      TH1D *hsim = new TH1D("hsim", "r{^2} distribution from simulation",
                            120, -10., 10.);

   gStyle->SetPadLeftMargin(0.15);
   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 600, 500);
   h10->Draw(name + ">>hsim", "x[6]>0.1");
   return hsim;
}
