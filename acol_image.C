/*
 * acol_image - a modeling and visualization class for representing
 *              the internal geometry of the GlueX active collimator.
 *              Methods are provided for visualizing scans of the beam
 *              across the collimator face, and for predicting what
 *              the currents should be as a sum of contributions from
 *              the basic geometric elements.
 *
 * author: richard.t.jones at uconn.edu
 * version: october 24, 2014
 */

#include <acol_image.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include <TVirtualFFT.h>
#include <TMatrixT.h>
#include <TVectorT.h>
#include <TStyle.h>
#include <TF1.h>

acol_image::acol_image()
 : fXscanOffset(-1.2), fYscanOffset(-2.5), fChiSquare(0)
{
   /// Makes a new canvas for drawing the active collimator viewer
   /// and draws a blank image showing the outlines of the display.

   fCanvasWidth = 500;
   fCanvasHeight = 500;
   fCanvas = (TCanvas*)gROOT->FindObject("acol_image_canvas");
   if (fCanvas)
      delete fCanvas;
   fCanvas = new TCanvas("acol_image_canvas",
                         "Hall D Active Collimator", 
                         50, 20, fCanvasWidth, fCanvasHeight);
   fCanvas->Range(-70.0,-70.0,70.0,70.0);
   fPline = new TPolyLine();
   fArc = new TArc();

   // Draw collimator outline with wedge outlines
   fCollimatorDiam = 0.197 * 25.4;
   fCoverPlateDiam = 5.4 * 25.4;
   fWedgeHalfAngle = 30.0;
   fWedgeAngleNseg = 30;
   fInnerWedgeR[0] = 0.0984 * 25.4;
   fInnerWedgeR[1] = 0.9842 * 25.4;
   fInnerWedgeX2 = 0.3200 * 25.4;
   fInnerHoleR = 1.6 / 2 * 25.4;
   fInnerConnR = 2.0 / 2 * 25.4;
   fOuterWedgeR[0] = 2.3604 / 2 * 25.4;
   fOuterWedgeR[1] = 4.7244 / 2 * 25.4;
   fOuterWedgeX2 = 1.4500 * 25.4;
   fOuterHoleR = 4.2 / 2 * 25.4;
   fOuterConnR = 4.6 / 2 * 25.4;
   fWedgeHoleDiam = 0.125 * 25.4;
   fWedgeConnDiam = 0.25 * 25.4;
   fWedgeHoleAngle = 15.;

   for (int phi=0; phi < 360; phi += 90) {
      draw_wedge(fInnerWedgeR, fInnerWedgeX2, fWedgeHalfAngle, phi);
      draw_wedge(fOuterWedgeR, fOuterWedgeX2, fWedgeHalfAngle, phi);
   }

   load_couplings();

   gStyle->SetPadRightMargin(0.20);
   fCanvas1 = (TCanvas*)gROOT->FindObject("c1");
   if (fCanvas1 == 0)
      fCanvas1 = new TCanvas("c1", "c1", 30, 150, 600, 500);
}

acol_image::~acol_image()
{
   delete fCanvas;
   delete fPline;
   delete fArc;
   delete fCanvas1;
}

void acol_image::draw_wedge(double rlim[2], double x2, 
                            double hangle, double phi0)
{
   // Draws an image of one wedge on the current canvas (in mm).

   double x[2*fWedgeAngleNseg];
   double y[2*fWedgeAngleNseg];
   double phi = phi0 - hangle;
   double dphi = 2 * hangle / fWedgeAngleNseg;

   // pins region
   int seg;
   double r1 = rlim[0] + 0.25;
   for (seg=0; seg <= fWedgeAngleNseg; ++seg) {
      x[seg] = r1 * cos((phi + seg * dphi) * 3.1415926/180.);
      y[seg] = r1 * sin((phi + seg * dphi) * 3.1415926/180.);
   }
   --seg;
   double r2 = x2 / cos(hangle * 3.1415926/180.);
   x[seg+1] = r2 * cos((phi + seg * dphi) * 3.1415926/180.);
   y[seg+1] = r2 * sin((phi + seg * dphi) * 3.1415926/180.);
   x[seg+2] = r2 * cos(phi * 3.1415926/180.);
   y[seg+2] = r2 * sin(phi * 3.1415926/180.);
   x[seg+3] = r1 * cos(phi * 3.1415926/180.);
   y[seg+3] = r1 * sin(phi * 3.1415926/180.);
   fPline->SetFillColor(41);
   fPline->DrawPolyLine(fWedgeAngleNseg + 4, x, y, "f");

   // inner arc
   for (int nseg=0; nseg <= fWedgeAngleNseg; ++nseg) {
      x[nseg] = rlim[0] * cos((phi + nseg * dphi) * 3.1415926/180.);
      y[nseg] = rlim[0] * sin((phi + nseg * dphi) * 3.1415926/180.);
   }
   fPline->DrawPolyLine(fWedgeAngleNseg + 1, x, y);

   // outer arc
   for (int nseg=0; nseg <= fWedgeAngleNseg; ++nseg) {
      x[nseg] *= rlim[1] / rlim[0];
      y[nseg] *= rlim[1] / rlim[0];
   }
   fPline->DrawPolyLine(fWedgeAngleNseg + 1, x, y);

   // 2 radial edges
   x[0] = rlim[0] * cos(phi * 3.1415926/180);
   y[0] = rlim[0] * sin(phi * 3.1415926/180);
   x[1] = x[0] * rlim[1] / rlim[0];
   y[1] = y[0] * rlim[1] / rlim[0];
   fPline->DrawPolyLine(2, x, y);
   phi += 2 * hangle;
   x[0] = rlim[0] * cos(phi * 3.1415926/180);
   y[0] = rlim[0] * sin(phi * 3.1415926/180);
   x[1] = x[0] * rlim[1] / rlim[0];
   y[1] = y[0] * rlim[1] / rlim[0];
   fPline->DrawPolyLine(2, x, y);

   // add the mounting holes
   double rhole = (rlim[0] < fInnerHoleR)? fInnerHoleR : fOuterHoleR;
   fArc->DrawArc(rhole * cos((phi0 + fWedgeHoleAngle) * 3.1415926/180.),
                 rhole * sin((phi0 + fWedgeHoleAngle) * 3.1415926/180.),
                 fWedgeHoleDiam / 2.);
   fArc->DrawArc(rhole * cos((phi0 - fWedgeHoleAngle) * 3.1415926/180.),
                 rhole * sin((phi0 - fWedgeHoleAngle) * 3.1415926/180.),
                 fWedgeHoleDiam / 2.);

   fCanvas->Update();
}

void acol_image::PaintScan(TChain *tscantree)
{
   // Read 2D scan data from an input tree and plot the trajectory
   // of the scan on top of the active collimator image.

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

   if (tscantree == 0) {
      std::cerr << "Error in PaintScan - called without a valid"
                << " input tree, cannot continue." << std::endl;
      return;
   }
   tscantree->SetBranchAddress("tscan", &tscan.tscan_s);
   Long64_t nentries = tscantree->GetEntries();
   double x[nentries], y[nentries];
   double &Xscan = tscan.Xscan_bpu_mm;
   double &Yscan = tscan.Yscan_bpu_mm;

   TFile *seg = 0;
   int npoints = 0;
   for (Long64_t i=0; i < nentries; ++i) {
      tscantree->GetEntry(i);
      if (tscantree->GetCurrentFile() != seg) {
         if (npoints > 0)
            paint_stripe(npoints, x, y);
         seg = tscantree->GetCurrentFile();
         npoints = 0;
      }
      x[npoints] = Xscan - fXscanOffset;
      y[npoints] = Yscan - fYscanOffset;
      ++npoints;
   }
   if (npoints > 0)
      paint_stripe(npoints, x, y);
}

void acol_image::paint_stripe(int npoints, double *x, double *y)
{
   // Generates a visualization of a collimator scan in the form
   // of a red stripe painted over the collimator image.

   fCanvas->cd();
   fPline->SetLineColor(2);
   fPline->SetLineWidth(3);
   fPline->DrawPolyLine(npoints, x, y);
   fCanvas1->cd();
}

TH2D *acol_image::BeamConvolution(TH2D *hin)
{
   // Clones input histogram hin and fills the clone with the original
   // image convoluted with the transverse beam intensity profile,
   // as represetned in the Beam_dIdr2() method. This implementation
   // uses the convolution theorem and the TVirtualFFT package, so
   // there can be wrap-around artifacts at the edges that come from
   // the assumption of infinite periodicity of the 2D distributions.

   int ipoint;
   int nbins[2];
   nbins[0] = hin->GetNbinsX();
   nbins[1] = hin->GetNbinsY();
   double xmin = hin->GetXaxis()->GetBinLowEdge(1);
   double xmax = hin->GetXaxis()->GetBinUpEdge(hin->GetNbinsX());
   double ymin = hin->GetYaxis()->GetBinLowEdge(1);
   double ymax = hin->GetYaxis()->GetBinUpEdge(hin->GetNbinsY());
   double dx = (xmax - xmin) / nbins[0];
   double dy = (ymax - ymin) / nbins[1];

   // set up the FFTs of the input histogram and the beam profiles
   TVirtualFFT *fft_hin = TVirtualFFT::FFT(2, nbins, "R2C M K");
   TVirtualFFT *fft_beam = TVirtualFFT::FFT(2, nbins, "R2C M K");
   TVirtualFFT *fft_spot = TVirtualFFT::FFT(2, nbins, "R2C M K");
   ipoint = 0;
   for (int ix=1; ix <= nbins[0]; ++ix) {
      double x = hin->GetXaxis()->GetBinCenter(ix);
      for (int iy=1; iy <= nbins[1]; ++iy) {
         double y = hin->GetYaxis()->GetBinCenter(iy);
         double r2 = x*x + y*y;
         fft_hin->SetPoint(ipoint, hin->GetBinContent(ix, iy));
         fft_beam->SetPoint(ipoint, Beam_dIdr2(r2));
         if (ix - 1 == (nbins[0] + 1)/2 && iy - 1 == (nbins[1] + 1)/2)
            fft_spot->SetPoint(ipoint, 12);
         else
            fft_spot->SetPoint(ipoint, 0.1 / sqrt(0.01 + r2) +
                                       0.030 * exp(-0.5 * r2 / 35));
         ++ipoint;
      }
   }
   fft_hin->Transform();
   fft_beam->Transform();
   fft_spot->Transform();

   // take the complex product of the two FTs
   int Nyqvist_limit = (nbins[0] + 1) * (nbins[1] + 1) / 2 + 1;
   double areal[Nyqvist_limit], aimag[Nyqvist_limit];
   double breal[Nyqvist_limit], bimag[Nyqvist_limit];
   double creal[Nyqvist_limit], cimag[Nyqvist_limit];
   double dreal[Nyqvist_limit], dimag[Nyqvist_limit];
   fft_hin->GetPointsComplex(areal, aimag);
   fft_beam->GetPointsComplex(breal, bimag);
   fft_spot->GetPointsComplex(creal, cimag);
   TVirtualFFT *fft_inv = TVirtualFFT::FFT(2, nbins, "C2R M K");
   for (ipoint = 0; ipoint < Nyqvist_limit; ++ipoint) {
      double re = areal[ipoint]*breal[ipoint] - aimag[ipoint]*bimag[ipoint];
      double im = areal[ipoint]*bimag[ipoint] + aimag[ipoint]*breal[ipoint];
      dreal[ipoint] = re * creal[ipoint] - im * cimag[ipoint];
      dimag[ipoint] = re * cimag[ipoint] + im * creal[ipoint];
      dreal[ipoint] /= nbins[0] * nbins[1];
      dimag[ipoint] /= nbins[0] * nbins[1];
   }
   fft_inv->SetPointsComplex(dreal, dimag);

   // take the inverse transform and save the result in hout
   fft_inv->Transform();
   TH2D *hout = (TH2D*)TH2D::TransformHisto(fft_inv, 0, "RE");

   // clean up the FFT objects
   delete fft_hin;
   delete fft_beam;
   delete fft_spot;
   delete fft_inv;

   std::string name = hin->GetName();
   hout->SetName((name + "_conv").c_str());
   std::string title = hin->GetTitle();
   hout->SetTitle((title + " convolved with beam").c_str());
   hout->GetXaxis()->SetLimits(xmin, xmax);
   hout->GetYaxis()->SetLimits(ymin, ymax);
   hout->SetContour(100);
   hout->SetStats(0);
   return hout;
}

TH2D *acol_image::BeamConvolution2(TH2D *hin)
{
   // Reimplementation of BeamConvolution, but using an explicit
   // convolution sum algorithm instead of using FFTs and the convolution
   // theorem. It is much slower than BeamConvolution but it does not
   // suffer from the wrap-around effects that are inherent in the
   // application of the convolution theorem.

   std::string name = hin->GetName();
   TH2D *hout = (TH2D*)hin->Clone((name + "_conv").c_str());
   hout->Reset();
   int nbinx = hout->GetNbinsX();
   int nbiny = hout->GetNbinsY();
   for (int ix=1; ix <= nbinx; ++ix) {
      std::cout << ix << '\r' << std::flush;
      double x = hout->GetXaxis()->GetBinCenter(ix);
      for (int iy=1; iy <= nbiny; ++iy) {
         double y = hout->GetYaxis()->GetBinCenter(iy);
         for (int jx=1; jx <= nbinx; ++jx) {
            double dx = hout->GetXaxis()->GetBinCenter(jx) - x;
            for (int jy=1; jy <= nbiny; ++jy) {
               double dy = hout->GetYaxis()->GetBinCenter(jy) - y;
               double r2 = dx*dx + dy*dy;
               hout->Fill(x, y, hin->GetBinContent(jx,jy) * Beam_dIdr2(r2));
            }
         }
      }
   }
   hout->SetContour(100);
   return hout;
}

TH2D *acol_image::MapInnerXplus(TH2D *htemp)
{
   // Generate a 2D histogram populated with the response 
   // of the inner X+ wedge to a point-like beam directed
   // into each pixel in the map.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("map_ixp")))
         delete  htemp;
      htemp = new TH2D("map_ixp", "inner X+ wedge response map",
                       281, -70, 70, 281, -70, 70);
   }

   int nbinx = htemp->GetNbinsX();
   int nbiny = htemp->GetNbinsY();
   for (int ix=1; ix < nbinx; ++ix) {
      for (int iy=1; iy < nbiny; ++iy) {
         double x = htemp->GetXaxis()->GetBinCenter(ix);
         double y = htemp->GetYaxis()->GetBinCenter(iy);
         double resp = fCoupling[0] * hits_inner_wedge_pins(x, y, 0) +
                       fCoupling[1] * hits_inner_wedge_flat(x, y, 0) +
                       fCoupling[2] * hits_inner_wedge_hole(x, y, 0) +
                       fCoupling[3] * hits_inner_wedge_conn(x, y, 0) +
                       fCoupling[4] * hits_inner_wedge_pins(x, y, 90) +
                       fCoupling[5] * hits_inner_wedge_flat(x, y, 90) +
                       fCoupling[6] * hits_inner_wedge_hole(x, y, 90) +
                       fCoupling[7] * hits_inner_wedge_conn(x, y, 90) +
                       fCoupling[8] * hits_inner_wedge_pins(x, y, 180) +
                       fCoupling[9] * hits_inner_wedge_flat(x, y, 180) +
                       fCoupling[10] * hits_inner_wedge_hole(x, y, 180) +
                       fCoupling[11] * hits_inner_wedge_conn(x, y, 180) +
                       fCoupling[4] * hits_inner_wedge_pins(x, y, 270) +
                       fCoupling[5] * hits_inner_wedge_flat(x, y, 270) +
                       fCoupling[6] * hits_inner_wedge_hole(x, y, 270) +
                       fCoupling[7] * hits_inner_wedge_conn(x, y, 270) +
                       fCoupling[12] * hits_outer_wedge_pins(x, y, 0) +
                       fCoupling[13] * hits_outer_wedge_flat(x, y, 0) +
                       fCoupling[14] * hits_outer_wedge_hole(x, y, 0) +
                       fCoupling[15] * hits_outer_wedge_conn(x, y, 0) +
                       fCoupling[16] * hits_outer_wedge_pins(x, y, 90) +
                       fCoupling[17] * hits_outer_wedge_flat(x, y, 90) +
                       fCoupling[18] * hits_outer_wedge_hole(x, y, 90) +
                       fCoupling[19] * hits_outer_wedge_conn(x, y, 90) +
                       fCoupling[20] * hits_outer_wedge_pins(x, y, 180) +
                       fCoupling[21] * hits_outer_wedge_flat(x, y, 180) +
                       fCoupling[22] * hits_outer_wedge_hole(x, y, 180) +
                       fCoupling[23] * hits_outer_wedge_conn(x, y, 180) +
                       fCoupling[16] * hits_outer_wedge_pins(x, y, 270) +
                       fCoupling[17] * hits_outer_wedge_flat(x, y, 270) +
                       fCoupling[18] * hits_outer_wedge_hole(x, y, 270) +
                       fCoupling[19] * hits_outer_wedge_conn(x, y, 270) +
                       fCoupling[36] * hits_cover_plate(x, y);
         htemp->Fill(x,y,resp);
      }
   }
   htemp->SetStats(0);
   htemp->SetContour(100);
   return htemp;
}

TH2D *acol_image::MapInnerXminus(TH2D *htemp)
{
   // Generate a 2D histogram populated with the response 
   // of the inner X- wedge to a point-like beam directed
   // into each pixel in the map.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("map_ixm")))
         delete  htemp;
      htemp = new TH2D("map_ixm", "inner X- wedge response map",
                       281, -70, 70, 281, -70, 70);
   }

   htemp = MapInnerXplus(htemp);
   TH2D *hcopy = (TH2D*)htemp->Clone("map_ixm_copy");
   int nbinx = hcopy->GetNbinsX();
   int nbiny = hcopy->GetNbinsY();
   for (int ix=1; ix < nbinx; ++ix)
      for (int iy=1; iy < nbinx; ++iy)
         htemp->SetBinContent(nbinx - ix + 1, 
                              nbiny - iy + 1, hcopy->GetBinContent(ix, iy));
   delete hcopy;
   return htemp;
}

TH2D *acol_image::MapInnerYplus(TH2D *htemp)
{
   // Generate a 2D histogram populated with the response 
   // of the inner Y+ wedge to a point-like beam directed
   // into each pixel in the map.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("map_iyp")))
         delete  htemp;
      htemp = new TH2D("map_iyp", "inner Y+ wedge response map",
                       281, -70, 70, 281, -70, 70);
   }

   htemp = MapInnerXplus(htemp);
   TH2D *hcopy = (TH2D*)htemp->Clone("map_iyp_copy");
   int nbinx = hcopy->GetNbinsX();
   int nbiny = hcopy->GetNbinsY();
   for (int ix=1; ix < nbinx; ++ix)
      for (int iy=1; iy < nbinx; ++iy)
         htemp->SetBinContent(nbiny - iy + 1, ix, hcopy->GetBinContent(ix, iy));
   delete hcopy;
   return htemp;
}

TH2D *acol_image::MapInnerYminus(TH2D *htemp)
{
   // Generate a 2D histogram populated with the response 
   // of the inner Y- wedge to a point-like beam directed
   // into each pixel in the map.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("map_iym")))
         delete  htemp;
      htemp = new TH2D("map_iym", "inner Y- wedge response map",
                       281, -70, 70, 281, -70, 70);
   }

   htemp = MapInnerXplus(htemp);
   TH2D *hcopy = (TH2D*)htemp->Clone("map_iym_copy");
   int nbinx = hcopy->GetNbinsX();
   int nbiny = hcopy->GetNbinsY();
   for (int ix=1; ix < nbinx; ++ix)
      for (int iy=1; iy < nbinx; ++iy)
         htemp->SetBinContent(iy, nbinx - ix + 1, hcopy->GetBinContent(ix, iy));
   delete hcopy;
   return htemp;
}

TH2D *acol_image::MapOuterXplus(TH2D *htemp)
{
   // Generate a 2D histogram populated with the response 
   // of the outer X+ wedge to a point-like beam directed
   // into each pixel in the map.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("map_oxp")))
         delete  htemp;
      htemp = new TH2D("map_oxp", "outer X+ wedge response map",
                       281, -70, 70, 281, -70, 70);
   }

   int nbinx = htemp->GetNbinsX();
   int nbiny = htemp->GetNbinsY();
   for (int ix=1; ix < nbinx; ++ix) {
      for (int iy=1; iy < nbiny; ++iy) {
         double x = htemp->GetXaxis()->GetBinCenter(ix);
         double y = htemp->GetYaxis()->GetBinCenter(iy);
         double resp = fCoupling[24] * hits_outer_wedge_pins(x, y, 0) +
                       fCoupling[25] * hits_outer_wedge_flat(x, y, 0) +
                       fCoupling[26] * hits_outer_wedge_hole(x, y, 0) +
                       fCoupling[27] * hits_outer_wedge_conn(x, y, 0) +
                       fCoupling[28] * hits_outer_wedge_pins(x, y, 90) +
                       fCoupling[29] * hits_outer_wedge_flat(x, y, 90) +
                       fCoupling[30] * hits_outer_wedge_hole(x, y, 90) +
                       fCoupling[31] * hits_outer_wedge_conn(x, y, 90) +
                       fCoupling[32] * hits_outer_wedge_pins(x, y, 180) +
                       fCoupling[33] * hits_outer_wedge_flat(x, y, 180) +
                       fCoupling[34] * hits_outer_wedge_hole(x, y, 180) +
                       fCoupling[35] * hits_outer_wedge_conn(x, y, 180) +
                       fCoupling[28] * hits_outer_wedge_pins(x, y, 270) +
                       fCoupling[29] * hits_outer_wedge_flat(x, y, 270) +
                       fCoupling[30] * hits_outer_wedge_hole(x, y, 270) +
                       fCoupling[31] * hits_outer_wedge_conn(x, y, 270) +
                       fCoupling[12] * hits_inner_wedge_pins(x, y, 0) +
                       fCoupling[13] * hits_inner_wedge_flat(x, y, 0) +
                       fCoupling[14] * hits_inner_wedge_hole(x, y, 0) +
                       fCoupling[15] * hits_inner_wedge_conn(x, y, 0) +
                       fCoupling[16] * hits_inner_wedge_pins(x, y, 90) +
                       fCoupling[17] * hits_inner_wedge_flat(x, y, 90) +
                       fCoupling[18] * hits_inner_wedge_hole(x, y, 90) +
                       fCoupling[19] * hits_inner_wedge_conn(x, y, 90) +
                       fCoupling[20] * hits_inner_wedge_pins(x, y, 180) +
                       fCoupling[21] * hits_inner_wedge_flat(x, y, 180) +
                       fCoupling[22] * hits_inner_wedge_hole(x, y, 180) +
                       fCoupling[23] * hits_inner_wedge_conn(x, y, 180) +
                       fCoupling[16] * hits_inner_wedge_pins(x, y, 270) +
                       fCoupling[17] * hits_inner_wedge_flat(x, y, 270) +
                       fCoupling[18] * hits_inner_wedge_hole(x, y, 270) +
                       fCoupling[19] * hits_inner_wedge_conn(x, y, 270) +
                       fCoupling[37] * hits_cover_plate(x, y);
         htemp->Fill(x,y,resp);
      }
   }
   htemp->SetStats(0);
   htemp->SetContour(100);
   return htemp;
}

TH2D *acol_image::MapOuterXminus(TH2D *htemp)
{
   // Generate a 2D histogram populated with the response 
   // of the outer X- wedge to a point-like beam directed
   // into each pixel in the map.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("map_oxm")))
         delete  htemp;
      htemp = new TH2D("map_oxm", "outer X- wedge response map",
                       281, -70, 70, 281, -70, 70);
   }

   htemp = MapOuterXplus(htemp);
   TH2D *hcopy = (TH2D*)htemp->Clone("map_oxm_copy");
   int nbinx = hcopy->GetNbinsX();
   int nbiny = hcopy->GetNbinsY();
   for (int ix=1; ix < nbinx; ++ix)
      for (int iy=1; iy < nbinx; ++iy)
         htemp->SetBinContent(nbinx - ix + 1,
                              nbiny - iy + 1, hcopy->GetBinContent(ix, iy));
   delete hcopy;
   return htemp;
}

TH2D *acol_image::MapOuterYplus(TH2D *htemp)
{
   // Generate a 2D histogram populated with the response 
   // of the outer Y+ wedge to a point-like beam directed
   // into each pixel in the map.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("map_oyp")))
         delete  htemp;
      htemp = new TH2D("map_oyp", "outer Y+ wedge response map",
                       281, -70, 70, 281, -70, 70);
   }

   htemp = MapOuterXplus(htemp);
   TH2D *hcopy = (TH2D*)htemp->Clone("map_oyp_copy");
   int nbinx = hcopy->GetNbinsX();
   int nbiny = hcopy->GetNbinsY();
   for (int ix=1; ix < nbinx; ++ix)
      for (int iy=1; iy < nbinx; ++iy)
         htemp->SetBinContent(nbiny - iy + 1, ix, hcopy->GetBinContent(ix, iy));
   delete hcopy;
   return htemp;
}

TH2D *acol_image::MapOuterYminus(TH2D *htemp)
{
   // Generate a 2D histogram populated with the response 
   // of the outer Y- wedge to a point-like beam directed
   // into each pixel in the map.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("map_oym")))
         delete  htemp;
      htemp = new TH2D("map_oym", "inner Y- wedge response map",
                       281, -70, 70, 281, -70, 70);
   }

   htemp = MapOuterXplus(htemp);
   TH2D *hcopy = (TH2D*)htemp->Clone("map_oym_copy");
   int nbinx = hcopy->GetNbinsX();
   int nbiny = hcopy->GetNbinsY();
   for (int ix=1; ix < nbinx; ++ix)
      for (int iy=1; iy < nbinx; ++iy)
         htemp->SetBinContent(iy, nbinx - ix + 1, hcopy->GetBinContent(ix, iy));
   delete hcopy;
   return htemp;
}

int acol_image::FitScanData(TChain *tscantree)
{
   // Read 2D scan data from an input tree and do a linear least-squares
   // fit of the 8 measured currents to the acol_image model, using the
   // fCoupling factors as adjustable parameters. Upon completion, the
   // fCoupling data members are loaded with the best-fit values. Return
   // value is 0 for success, or non-zero for failure to find a solution.

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

   if (tscantree == 0) {
      std::cerr << "Error in FitScanData - called without a valid"
                << " input tree, cannot continue." << std::endl;
      return -1;
   }
   tscantree->SetBranchAddress("tsum", &tsum.t_s);
   tscantree->SetBranchAddress("tscan", &tscan.tscan_s);
   Long64_t nentries = tscantree->GetEntries();
   if (nentries < 10) {
      std::cerr << "Error in FitScanData - called without sufficient"
                << " data points to perform a fit, cannot continue."
                << std::endl;
      return -1;
   }
   double *Vmeas = &tsum.Vixp;
   double &Xscan = tscan.Xscan_bpu_mm;
   double &Yscan = tscan.Yscan_bpu_mm;

   // Create an array of 2D maps Z_kl where k=0..7 is the index of the
   // sensing wedge and l=0..24 is the index of the source element. The
   // order of the sensing wedge is ix-, ix+, iy+, iy-, ox-, ox+, oy+, oy-.
   // The order of the source elements is specified by the order of the
   // fCoupling data members as listed in the header file. The model
   // function for the current on wedge k for the scan point i is
   //
   //     y_ki = sum_l{ a_l Z_kl(x_i) }
   //
   // where x_i is a 2D vector indicating where the beam was located at
   // point i during the scan, and must be within the range of map Z_kl.

   int nbins = 281;
   double xmax = 70.;
   TH2D *Zmap[8][fNpars];
   for (int l=0; l < fNpars; l++) {
      for (int ll=0; ll < fNpars; ll++) {
         fCoupling[ll] = (l == ll)? 1 : 0;
      }
      TH2D *hmap[8];
      for (int k=0; k < 8; ++k) {
         std::stringstream name;
         name << "Z_" << k << "_" << l;
         hmap[k] = (TH2D*)gROOT->FindObject(name.str().c_str());
         if (hmap[k] != 0)
            delete hmap[k];
         hmap[k] = new TH2D(name.str().c_str(), "",
                            nbins, -xmax, xmax, nbins, -xmax, xmax);
      }
      Zmap[0][l] = BeamConvolution(MapInnerXminus(hmap[0]));
      Zmap[1][l] = BeamConvolution(MapInnerXplus(hmap[1]));
      Zmap[2][l] = BeamConvolution(MapInnerYplus(hmap[2]));
      Zmap[3][l] = BeamConvolution(MapInnerYminus(hmap[3]));
      Zmap[4][l] = BeamConvolution(MapOuterXminus(hmap[4]));
      Zmap[5][l] = BeamConvolution(MapOuterXplus(hmap[5]));
      Zmap[6][l] = BeamConvolution(MapOuterYplus(hmap[6]));
      Zmap[7][l] = BeamConvolution(MapOuterYminus(hmap[7]));
      std::cout << l << '\r' << std::flush;
   }

   // Form the second derivative matrix of the chisquare. It is
   // real symmetric, so only the upper half-matrix needs to be
   // filled. There is one matrix for the inner wedges and another
   // for the outer.

   TMatrixTSym<double> SD(fNpars);
   TVectorT<double> yZ(fNpars);
   for (int l=0; l < fNpars; ++l) {
      for (int ll=0; ll < fNpars; ++ll) {
         SD[l][ll] = 0;
      }
      yZ[l] = 0;
   }

   for (Long64_t i=0; i < nentries; ++i) {
      tscantree->GetEntry(i);
      double sigma2 = 0.01;
      for (int k=0; k < 8; ++k) {
         double Zki[fNpars];
         sigma2 = 0.02 * Vmeas[k];
         for (int l=0; l < fNpars; ++l) {
            double x = Xscan - fXscanOffset;
            double y = Yscan - fYscanOffset;
            Zki[l] = Zmap[k][l]->Interpolate(x, y);
            yZ[l] += Vmeas[k] * Zki[l] / sigma2;
         }
         for (int l=0; l < fNpars; ++l) {
            for (int ll=0; ll < fNpars; ++ll) {
              SD[l][ll] += Zki[l] * Zki[ll] / sigma2;
            }
         }
      }
   }

   // Now solve for the least-squares coefficients by linear methods

   TMatrixTSym<double> SDinv(SD);
   double detSD;
   SDinv.Invert(&detSD);
   if (detSD < 1e-16) {
      std::cerr << "FitScanData error - SD matrix is singular,"
                << " determinant is " << detSD
                << ", cannot continue." << std::endl;
      return -1;
   }

   TVectorT<double> fitvec(yZ);
   fitvec *= SDinv;
   for (int l=0; l < fNpars; ++l) {
      fCoupling[l] = fitvec[l];
   }

   // Compute the least-squares value

   double ysqr = 0;
   double chisqr = 0;
   for (Long64_t i=0; i < nentries; ++i) {
      tscantree->GetEntry(i);
      double sigma2 = 0.01;
      double x = Xscan - fXscanOffset;
      double y = Yscan - fYscanOffset;
      for (int k=0; k < 8; ++k) {
         sigma2 = 0.02 * Vmeas[k];
         ysqr += Vmeas[k] * Vmeas[k] / sigma2;
         double adotZ = 0;
         for (int l=0; l < fNpars; ++l)
            adotZ += fitvec[l] * Zmap[k][l]->Interpolate(x,y);
         chisqr += pow(Vmeas[k] - adotZ, 2) / sigma2;
      }
   }
 
   double asqr = 0;
   for (int l=0; l < fNpars; ++l)
      for (int ll=0; ll < fNpars; ++ll)
         asqr += fitvec[l] * SD[l][ll] * fitvec[ll];
 
   // Print a summary message and return

   std::cout << "FitScanData succeeds, best-fit parameters follow below"
             << std::endl;
   for (int i=0; i < 10; ++i)
      std::cout << "  coupling  " << i << " = " << fCoupling[i]
                << " +/- " << sqrt(SDinv[i][i]) << std::endl;
   for (int i=10; i < fNpars; ++i)
      std::cout << "  coupling " << i << " = " << fCoupling[i]
                << " +/- " << sqrt(SDinv[i][i]) << std::endl;
   std::cout << "least-square sum = " << ysqr << " - " << asqr
             << " = " << chisqr << std::endl;
   fChiSquare = chisqr;

   return 0;
}

TH1D *acol_image::PlotScanData(TChain *tscantree, int wedge)
{
   // Read 2D scan data from an input tree and plot the model response
   // as a function of the x displacement of the beam. The measured
   // response is superimposed as a set of data points with errors.

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

   if (tscantree == 0) {
      std::cerr << "Error in PlotScanData - called without a valid"
                << " input tree, cannot continue." << std::endl;
      return 0;
   }
   tscantree->SetBranchAddress("tsum", &tsum.t_s);
   tscantree->SetBranchAddress("tscan", &tscan.tscan_s);
   Long64_t nentries = tscantree->GetEntries();
   double *Vmeas = &tsum.Vixp;
   double &Xscan = tscan.Xscan_bpu_mm;
   double &Yscan = tscan.Yscan_bpu_mm;
   double xmin = tscantree->GetMinimum("Xscan_bpu_mm");
   double xmax = tscantree->GetMaximum("Xscan_bpu_mm");
   xmin -= (xmax - xmin) * 0.10;
   xmax += (xmax - xmin) * 0.09;

   TH1D *hmodel = (TH1D*)gROOT->FindObject("hscan_model");
   if (hmodel != 0)
      delete hmodel;
   hmodel = new TH1D("hscan_model", "active collimator scan", 100, xmin, xmax);
   TH1D *hdata = (TH1D*)gROOT->FindObject("hscan_data");
   if (hdata != 0)
      delete hdata;
   hdata = (TH1D*)hmodel->Clone("hscan_data");
   TH2D *hmap;
   if (wedge == 0)
      hmap = BeamConvolution(MapInnerXminus());
   else if (wedge == 1)
      hmap = BeamConvolution(MapInnerXplus());
   else if (wedge == 2)
      hmap = BeamConvolution(MapInnerYplus());
   else if (wedge == 3)
      hmap = BeamConvolution(MapInnerYminus());
   else if (wedge == 4)
      hmap = BeamConvolution(MapOuterXminus());
   else if (wedge == 5)
      hmap = BeamConvolution(MapOuterXplus());
   else if (wedge == 6)
      hmap = BeamConvolution(MapOuterYplus());
   else if (wedge == 7)
      hmap = BeamConvolution(MapOuterYminus());
   else {
      std::cerr << "Error in PlotScanData - called with invalid"
                << " wedge index " << wedge << ", cannot continue."
                << std::endl;
      return 0;
   }

   for (Long64_t i=0; i < nentries; ++i) {
      tscantree->GetEntry(i);
      double xscan = Xscan - fXscanOffset;
      double yscan = Yscan - fYscanOffset;
      hdata->Fill(xscan, Vmeas[wedge]);
      hmodel->Fill(xscan, yscan);
   }
   double yscan = 999;
   for (int bin=1; bin <= hmodel->GetNbinsX(); ++bin) {
      double x = hmodel->GetXaxis()->GetBinCenter(bin);
      double y = hmodel->GetBinContent(bin);
      if (y != 0)
         yscan = y;
      if (yscan != 999)
         hmodel->SetBinContent(bin, hmap->Interpolate(x, yscan));
   }

   int nbins = hdata->GetNbinsX();
   double err[nbins];
   for (int i=1; i <= nbins; ++i)
      err[i] = (hdata->GetBinContent(i) > 0)? 0.02 : 0;
   hdata->SetError(err);

   double min_model = hmodel->GetMinimum();
   double max_model = hmodel->GetMaximum();
   double min_data = hdata->GetMinimum();
   double max_data = hdata->GetMaximum();
   hmodel->SetMinimum(((min_model < min_data)? min_model : min_data) * 1.1); 
   hmodel->SetMaximum(((max_model > max_data)? max_model : max_data) * 1.1); 
   hmodel->Draw("c");
   hdata->Draw("esame");

   delete hmap;
   return hmodel;
}

TH2D *acol_image::PlotBeamSpot(TH2D *htemp)
{
   // Generate a 2D histogram of the intensity map of the photon beam
   // that the model uses to generate currents in the wedge elements
   // as a function of the beam offset from the collimator axis. If
   // the user provides an input histogram htemp then fill that, 
   // otherwise create a new one with nominal dimensions and resolution.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("beam_spot")))
         delete  htemp;
      htemp = new TH2D("beam_spot", "model collimator beam spot",
                       281, -70, 70, 281, -70, 70);
   }

   int nbinx = htemp->GetNbinsX();
   int nbiny = htemp->GetNbinsY();
   htemp->Reset();
   htemp->SetBinContent((nbinx + 1) / 2, (nbiny + 1) / 2, 1);
   return BeamConvolution(htemp);
}

TH1D *acol_image::PlotBeamSpotProfile(TH2D *htemp, double xoff)
{
   // Generate a photon beam intensity map with PlotBeamSpot and then
   // take a slice through it at offset xoff from the center of the beam
   // and return it as a 1D histogram.

   TH2D *h2 = PlotBeamSpot(htemp);
   TH1D *h1 = (TH1D*)gROOT->FindObject("beam_spot_profile");
   if (h1)
      delete h1;

   int nbins = h2->GetNbinsX();
   int xmin = h2->GetXaxis()->GetBinLowEdge(1);
   int xmax = h2->GetXaxis()->GetBinUpEdge(nbins);
   h1 = new TH1D("beam_spot_profile",
                 "model collimator beam spot profile", nbins, xmin, xmax);
   int iy = h2->GetYaxis()->FindBin(xoff);
   for (int ix=1; ix < nbins; ++ix)
      h1->SetBinContent(ix, h2->GetBinContent(ix, iy));
   h1->Draw();
   return h1;
}

TH2D *acol_image::PlotInnerXasym(TH2D *htemp)
{
   // Generate the asymmetry for the inner X wedges as a function of
   // beam displacement x,y from the collimator axis, and return
   // the result as a 2D histogram.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("innerXasym")))
         delete  htemp;
      htemp = new TH2D("innerXasym", "inner X wedges asymmetry",
                       181, -20, 20, 181, -20, 20);
   }
   TH2D *hminus = BeamConvolution(MapInnerXminus());
   TH2D* hplus = BeamConvolution(MapInnerXplus());
   TH2D *hwork = (TH2D*)hplus->Clone("hwork");
   hwork->Add(hminus, -1);
   hplus->Add(hminus, 1);
   hwork->Divide(hplus);
   int nbinx = htemp->GetNbinsX();
   int nbiny = htemp->GetNbinsY();
   for (int ix=1; ix <= nbinx; ++ix) {
      double x = htemp->GetXaxis()->GetBinCenter(ix);
      for (int iy=1; iy <= nbiny; ++iy) {
         double y = htemp->GetYaxis()->GetBinCenter(iy);
         htemp->SetBinContent(ix, iy, hwork->Interpolate(x, y));
      }
   }
   htemp->SetStats(0);
   htemp->Draw("colz");
   delete hplus;
   delete hminus;
   delete hwork;
   return htemp;
}

TH2D *acol_image::PlotInnerYasym(TH2D *htemp)
{
   // Generate the asymmetry for the inner Y wedges as a function of
   // beam displacement x,y from the collimator axis, and return
   // the result as a 2D histogram.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("innerYasym")))
         delete  htemp;
      htemp = new TH2D("innerYasym", "inner Y wedges asymmetry",
                       181, -20, 20, 181, -20, 20);
   }
   TH2D *hminus = BeamConvolution(MapInnerYminus());
   TH2D* hplus = BeamConvolution(MapInnerYplus());
   TH2D *hwork = (TH2D*)hplus->Clone("hwork");
   hwork->Add(hminus, -1);
   hplus->Add(hminus, 1);
   hwork->Divide(hplus);
   int nbinx = htemp->GetNbinsX();
   int nbiny = htemp->GetNbinsY();
   for (int ix=1; ix <= nbinx; ++ix) {
      double x = htemp->GetXaxis()->GetBinCenter(ix);
      for (int iy=1; iy <= nbiny; ++iy) {
         double y = htemp->GetYaxis()->GetBinCenter(iy);
         htemp->SetBinContent(ix, iy, hwork->Interpolate(x, y));
      }
   }
   htemp->SetStats(0);
   htemp->Draw("colz");
   delete hplus;
   delete hminus;
   delete hwork;
   return htemp;
}

TH2D *acol_image::PlotOuterXasym(TH2D *htemp)
{
   // Generate the asymmetry for the outer X wedges as a function of
   // beam displacement x,y from the collimator axis, and return
   // the result as a 2D histogram.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("outerXasym")))
         delete  htemp;
      htemp = new TH2D("outerXasym", "outer X wedges asymmetry",
                       281, -30, 30, 281, -30, 30);
   }
   TH2D *hminus = BeamConvolution(MapOuterXminus());
   TH2D* hplus = BeamConvolution(MapOuterXplus());
   TH2D *hwork = (TH2D*)hplus->Clone("hwork");
   hwork->Add(hminus, -1);
   hplus->Add(hminus, 1);
   hwork->Divide(hplus);
   int nbinx = htemp->GetNbinsX();
   int nbiny = htemp->GetNbinsY();
   for (int ix=1; ix <= nbinx; ++ix) {
      double x = htemp->GetXaxis()->GetBinCenter(ix);
      for (int iy=1; iy <= nbiny; ++iy) {
         double y = htemp->GetYaxis()->GetBinCenter(iy);
         htemp->SetBinContent(ix, iy, hwork->Interpolate(x, y));
      }
   }
   htemp->SetStats(0);
   htemp->Draw("colz");
   delete hplus;
   delete hminus;
   delete hwork;
   return htemp;
}

TH2D *acol_image::PlotOuterYasym(TH2D *htemp)
{
   // Generate the asymmetry for the outer Y wedges as a function of
   // beam displacement x,y from the collimator axis, and return
   // the result as a 2D histogram.

   if (htemp == 0) {
      if ((htemp = (TH2D*)gROOT->FindObject("outerYasym")))
         delete  htemp;
      htemp = new TH2D("outerYasym", "outer Y wedges asymmetry",
                       281, -30, 30, 281, -30, 30);
   }
   TH2D *hminus = BeamConvolution(MapOuterYminus());
   TH2D* hplus = BeamConvolution(MapOuterYplus());
   TH2D *hwork = (TH2D*)hplus->Clone("hwork");
   hwork->Add(hminus, -1);
   hplus->Add(hminus, 1);
   hwork->Divide(hplus);
   int nbinx = htemp->GetNbinsX();
   int nbiny = htemp->GetNbinsY();
   for (int ix=1; ix <= nbinx; ++ix) {
      double x = htemp->GetXaxis()->GetBinCenter(ix);
      for (int iy=1; iy <= nbiny; ++iy) {
         double y = htemp->GetYaxis()->GetBinCenter(iy);
         htemp->SetBinContent(ix, iy, hwork->Interpolate(x, y));
      }
   }
   htemp->SetStats(0);
   htemp->Draw("colz");
   delete hplus;
   delete hminus;
   delete hwork;
   return htemp;
}

TH2D *acol_image::PlotInnerSolution(double xasym, double yasym)
{
   // Plot the solution for the beam centroid given a set of
   // measured inner wedge asymmetries xasym and yasym.

   TH2D *ixa = PlotInnerXasym();
   TH2D *iya = PlotInnerYasym();
   TH2D *sol = (TH2D*)gROOT->FindObject("insol");
   if (sol != 0)
      delete sol;
   sol = (TH2D*)ixa->Clone("insol");
   int nbinx = sol->GetNbinsX();
   int nbiny = sol->GetNbinsY();
   for (int ix=1; ix <= nbinx; ++ix) {
      for (int iy=1; iy <= nbiny; ++iy) {
         double dax = ixa->GetBinContent(ix, iy) - xasym;
         double day = iya->GetBinContent(ix, iy) - yasym;
         sol->SetBinContent(ix, iy, exp(-0.5 * (dax*dax + day*day) / 1e-4));
      }
   }
   sol->SetStats(0);
   sol->Draw("colz");

   fXsolution = sol->GetMean(1);
   fYsolution = sol->GetMean(2);
   std::cout << "solution found at x,y = " << fXsolution << "," << fYsolution
             << " mm" << std::endl
             << "with corresponding asymmetries "
             << ixa->Interpolate(fXsolution, fYsolution) << ", "
             << iya->Interpolate(fXsolution, fYsolution) << std::endl;
   delete ixa;
   delete iya;
   return sol;
}

TH2D *acol_image::PlotOuterSolution(double xasym, double yasym)
{
   // Plot the solution for the beam centroid given a set of
   // measured outer wedge asymmetries xasym and yasym.

   TH2D *oxa = PlotOuterXasym();
   TH2D *oya = PlotOuterYasym();
   TH2D *sol = (TH2D*)gROOT->FindObject("outsol");
   if (sol != 0)
      delete sol;
   sol = (TH2D*)oxa->Clone("outsol");
   int nbinx = sol->GetNbinsX();
   int nbiny = sol->GetNbinsY();
   for (int ix=1; ix <= nbinx; ++ix) {
      for (int iy=1; iy <= nbiny; ++iy) {
         double dax = oxa->GetBinContent(ix, iy) - xasym;
         double day = oya->GetBinContent(ix, iy) - yasym;
         sol->SetBinContent(ix, iy, exp(-0.5 * (dax*dax + day*day) / 1e-4));
      }
   }
   sol->SetStats(0);
   sol->Draw("colz");

   fXsolution = sol->GetMean(1);
   fYsolution = sol->GetMean(2);
   std::cout << "solution found at x,y = " << fXsolution << "," << fYsolution
             << " mm" << std::endl
             << "with corresponding asymmetries "
             << oxa->Interpolate(fXsolution, fYsolution) << ", "
             << oya->Interpolate(fXsolution, fYsolution) << std::endl;
   delete oxa;
   delete oya;
   return sol;
}

void acol_image::CheckScanData(TChain *tscantree)
{
   // Read 2D scan data from an input tree and compare the coordinates
   // recorded in the tree with those computed based on the measured
   // active collimator asymmetries.

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

   if (tscantree == 0) {
      std::cerr << "Error in CheckScanData - called without a valid"
                << " input tree, cannot continue." << std::endl;
      return;
   }
   tscantree->SetBranchAddress("tscan", &tscan.tscan_s);
   tscantree->SetBranchAddress("tsum", &tsum.t_s);
   Long64_t nentries = tscantree->GetEntries();
   double *Vmeas = &tsum.Vixp;
   double &Xscan = tscan.Xscan_bpu_mm;
   double &Yscan = tscan.Yscan_bpu_mm;

   TH2D *ixa = PlotInnerXasym();
   TH2D *iya = PlotInnerYasym();
   TH2D *oxa = PlotOuterXasym();
   TH2D *oya = PlotOuterYasym();
   ixa = (TH2D*)ixa->Clone("my_private_ixa");
   iya = (TH2D*)iya->Clone("my_private_iya");
   oxa = (TH2D*)oxa->Clone("my_private_oxa");
   oya = (TH2D*)oya->Clone("my_private_oya");

   double ixmin = ixa->GetXaxis()->GetBinLowEdge(1);
   double ixmax = ixa->GetXaxis()->GetBinUpEdge(ixa->GetNbinsX());
   double iymin = ixa->GetYaxis()->GetBinLowEdge(1);
   double iymax = ixa->GetYaxis()->GetBinUpEdge(ixa->GetNbinsY());
   double oxmin = oxa->GetXaxis()->GetBinLowEdge(1);
   double oxmax = oxa->GetXaxis()->GetBinUpEdge(oxa->GetNbinsX());
   double oymin = oxa->GetYaxis()->GetBinLowEdge(1);
   double oymax = oxa->GetYaxis()->GetBinUpEdge(oxa->GetNbinsY());
   TH1D *hix = new TH1D("hix","x inner model vs x scan", 100, ixmin, ixmax);
   TH1D *hiy = new TH1D("hiy","y inner model vs y scan", 100, iymin, iymax);
   TH1D *hox = new TH1D("hox","x outer model vs x scan", 100, oxmin, oxmax);
   TH1D *hoy = new TH1D("hoy","y outer model vs y scan", 100, oymin, oymax);

   for (Long64_t entry=0; entry < nentries; ++entry) {
      tscantree->GetEntry(entry);
      double ixasym = (Vmeas[1] - Vmeas[0]) / (Vmeas[1] + Vmeas[0]);
      double iyasym = (Vmeas[2] - Vmeas[3]) / (Vmeas[2] + Vmeas[3]);
      double oxasym = (Vmeas[5] - Vmeas[4]) / (Vmeas[5] + Vmeas[4]);
      double oyasym = (Vmeas[6] - Vmeas[7]) / (Vmeas[6] + Vmeas[7]);
      double xscan = Xscan - fXscanOffset;
      double yscan = Yscan - fYscanOffset;
      std::cout << std::endl
                << "* entry " << entry << " has recorded x,y = "
                << xscan << ", " << yscan << " mm" << std::endl;
      std::cout << "  with measured inner wedge asymmetries x,y = "
                << ixasym << "," << iyasym << std::endl;
      if (xscan > ixmin && xscan < ixmax && yscan > iymin && yscan < iymax) {
         std::cout << "  and model inner wedge asymmetries x,y = "
                   << ixa->Interpolate(xscan, yscan) << ", "
                   << iya->Interpolate(xscan, yscan) << std::endl
                   << "    inner wedge ";
         PlotInnerSolution(ixasym, iyasym);
         hix->Fill(xscan, fXsolution);
         hiy->Fill(xscan, fYsolution);
      }
      else {
         std::cout << "  *** outside limits of inner wedge current model ***"
                   << std::endl;
      }
      std::cout << "  with measured outer wedge asymmetries x,y = "
                << oxasym << "," << oyasym << std::endl;
      if (xscan > oxmin && xscan < oxmax && yscan > oymin && yscan < oymax) {
         std::cout << "  and model outer wedge asymmetries x,y = "
                   << oxa->Interpolate(xscan, yscan) << ", "
                   << oya->Interpolate(xscan, yscan) << std::endl
                   << "    outer wedge ";
         PlotOuterSolution(oxasym, oyasym);
         hox->Fill(xscan, fXsolution);
         hoy->Fill(xscan, fYsolution);
      }
      else {
         std::cout << "  *** outside limits of outer wedge current model ***"
                   << std::endl;
      }
   }
   delete ixa;
   delete iya;
   delete oxa;
   delete oya;

   if (hix) {
      int nbins = hix->GetNbinsX();
      double err[nbins];
      for (int i=1; i <= nbins; ++i)
         err[i] = (hix->GetBinContent(i) != 0)? 0.2 : 0;
      hix->SetError(err);
   }
   if (hiy) {
      int nbins = hiy->GetNbinsX();
      double err[nbins];
      for (int i=1; i <= nbins; ++i)
         err[i] = (hiy->GetBinContent(i) != 0)? 0.2 : 0;
      hiy->SetError(err);
   }
   if (hox) {
      int nbins = hox->GetNbinsX();
      double err[nbins];
      for (int i=1; i <= nbins; ++i)
         err[i] = (hox->GetBinContent(i) != 0)? 0.2 : 0;
      hox->SetError(err);
   }
   if (hoy) {
      int nbins = hoy->GetNbinsX();
      double err[nbins];
      for (int i=1; i <= nbins; ++i)
         err[i] = (hoy->GetBinContent(i) != 0)? 0.2 : 0;
      hoy->SetError(err);
   }
   hix->SetStats(0);
   hiy->SetStats(0);
   hox->SetStats(0);
   hoy->SetStats(0);
   new TF1("feq", "x", -50, 50);
}

void acol_image::load_couplings(char *filename)
{
   // Loads the values of the model coupling parameters from an input
   // file, in the format of one floating point parameter per line.
   // If no input file is specified, default values are loaded.

   if (filename == 0) {
      fCoupling[0] = 1;
      fCoupling[1] = 0.35;
      fCoupling[2] = 0.05;
      fCoupling[3] = 0.50;
      fCoupling[4] = 0;
      fCoupling[5] = 0;
      fCoupling[6] = 0;
      fCoupling[7] = 0;
      fCoupling[8] = 0;
      fCoupling[9] = 0;
      fCoupling[10] = 0;
      fCoupling[11] = 0;
      fCoupling[12] = 0;
      fCoupling[13] = 0;
      fCoupling[14] = 0;
      fCoupling[15] = 0;
      fCoupling[16] = 0;
      fCoupling[17] = 0;
      fCoupling[18] = 0;
      fCoupling[19] = 0;
      fCoupling[20] = 0;
      fCoupling[21] = 0;
      fCoupling[22] = 0;
      fCoupling[23] = 0;
      fCoupling[24] = 1;
      fCoupling[25] = 0.35;
      fCoupling[26] = 0.05;
      fCoupling[27] = 0.50;
      fCoupling[28] = 0;
      fCoupling[29] = 0;
      fCoupling[30] = 0;
      fCoupling[31] = 0;
      fCoupling[32] = 0;
      fCoupling[33] = 0;
      fCoupling[34] = 0;
      fCoupling[35] = 0;
      fCoupling[36] = 0;
      fCoupling[37] = 0;
      return;
   }

   std::ifstream parfile(filename);
   if (! parfile.good()) {
      std::cerr << "acol_image::load_couplings error - unable to open"
                << " input file " << filename << std::endl;
      return;
   }
   for (int i=0; i < fNpars; ++i)
      parfile >> fCoupling[i];
   parfile.close();
}

void acol_image::save_couplings(char *filename)
{
   // Saves the values of the model coupling parameters to an output
   // file, in the format of one floating point parameter per line.
   // If the file already exists, it is overwritten.

   std::ofstream parfile(filename);
   for (int i=0; i < fNpars; ++i)
      parfile << fCoupling[i] << std::endl;
   parfile.close();
}

void acol_image::Print()
{
   std::cout << "acol_image object property      value" << std::endl
             << "     fCollimatorDiam          = " 
             << fCollimatorDiam << " mm" << std::endl
             << "     fCoverPlateDiam          = "
             << fCoverPlateDiam << " mm" << std::endl
             << "     fInnerWedgeR[2]          = "
             << fInnerWedgeR[0] << ", " << fInnerWedgeR[1]
             << " mm" << std::endl
             << "     fOuterWedgeR[2]          = "
             << fOuterWedgeR[0] << ", " << fOuterWedgeR[1]
             << " mm" << std::endl
             << "     fInnerWedgeX2            = "
             << fInnerWedgeX2 << " mm" << std::endl
             << "     fOuterWedgeX2            = "
             << fOuterWedgeX2 << " mm" << std::endl
             << "     fWedgeHalfAngle          = "
             << fWedgeHalfAngle << " deg" << std::endl
             << "     fInnerHoleR              = "
             << fInnerHoleR << " mm" << std::endl
             << "     fInnerConnR              = "
             << fInnerConnR << " mm" << std::endl
             << "     fOuterHoleR              = "
             << fOuterHoleR << " mm" << std::endl
             << "     fOuterConnR              = "
             << fOuterConnR << " mm" << std::endl
             << "     fWedgeHoleDiam           = "
             << fWedgeHoleDiam << " mm" << std::endl
             << "     fWedgeHoleAngle          = "
             << fWedgeHoleAngle << " deg" << std::endl
             << "     fWedgeConnDiam           = "
             << fWedgeConnDiam << " mm" << std::endl
             << "     fWedgeAngleNseg          = "
             << fWedgeAngleNseg << std::endl
             << "     fXscanOffset             = "
             << fXscanOffset << " mm" << std::endl
             << "     fYscanOffset             = "
             << fYscanOffset << " mm" << std::endl
             << std::endl;

   std::cout << "Model parameters have the following values:" << std::endl
             << "   fCoupling[0..3]     = "
             << std::setw(12) << fCoupling[0] << ", " 
             << std::setw(12) << fCoupling[1] << ", "
             << std::setw(12) << fCoupling[2] << ", "
             << std::setw(12) << fCoupling[3] << std::endl
             << "   fCoupling[4..7]     = "
             << std::setw(12) << fCoupling[4] << ", " 
             << std::setw(12) << fCoupling[5] << ", "
             << std::setw(12) << fCoupling[6] << ", "
             << std::setw(12) << fCoupling[7] << std::endl
             << "   fCoupling[8..11]    = "
             << std::setw(12) << fCoupling[8] << ", " 
             << std::setw(12) << fCoupling[9] << ", "
             << std::setw(12) << fCoupling[10] << ", "
             << std::setw(12) << fCoupling[11] << std::endl
             << "   fCoupling[12..15]   = "
             << std::setw(12) << fCoupling[12] << ", " 
             << std::setw(12) << fCoupling[13] << ", "
             << std::setw(12) << fCoupling[14] << ", "
             << std::setw(12) << fCoupling[15] << std::endl
             << "   fCoupling[16..19]   = "
             << std::setw(12) << fCoupling[16] << ", " 
             << std::setw(12) << fCoupling[17] << ", "
             << std::setw(12) << fCoupling[18] << ", "
             << std::setw(12) << fCoupling[19] << std::endl
             << "   fCoupling[20..23]   = "
             << std::setw(12) << fCoupling[20] << ", " 
             << std::setw(12) << fCoupling[21] << ", "
             << std::setw(12) << fCoupling[22] << ", "
             << std::setw(12) << fCoupling[23] << std::endl
             << "   fCoupling[24..27]   = "
             << std::setw(12) << fCoupling[24] << ", " 
             << std::setw(12) << fCoupling[25] << ", "
             << std::setw(12) << fCoupling[25] << ", "
             << std::setw(12) << fCoupling[27] << std::endl
             << "   fCoupling[28..31]   = "
             << std::setw(12) << fCoupling[28] << ", " 
             << std::setw(12) << fCoupling[29] << ", "
             << std::setw(12) << fCoupling[30] << ", "
             << std::setw(12) << fCoupling[31] << std::endl
             << "   fCoupling[32..35]   = "
             << std::setw(12) << fCoupling[32] << ", " 
             << std::setw(12) << fCoupling[33] << ", "
             << std::setw(12) << fCoupling[34] << ", "
             << std::setw(12) << fCoupling[35] << std::endl
             << "   fCoupling[36..37]   = "
             << std::setw(12) << fCoupling[36] << ", "
             << std::setw(12) << fCoupling[37] << std::endl;

   if (fChiSquare > 0)
      std::cout << "chi-square from latest fit was " << fChiSquare
                << std::endl;
}

int acol_image::hits_inner_wedge_pins(double x, double y, double phi0)
{
   // Logical function indicating whether or not point x,y lies
   // inside the region of the pins on inner wedge centered at phi0.

   double phi0rad = phi0 * M_PI/180;
   double phi1rad = (phi0 - 360) * M_PI/180;
   double phi = atan2(y,x);
   return ((fabs(phi - phi0rad) > fWedgeHalfAngle * M_PI/180 &&
            fabs(phi - phi1rad) > fWedgeHalfAngle * M_PI/180) ||
           x*x + y*y < fInnerWedgeR[0] * fInnerWedgeR[0] ||
           x*cos(phi0rad) + y*sin(phi0rad) > fInnerWedgeX2)? 0 : 1;
}

int acol_image::hits_inner_wedge_flat(double x, double y, double phi0)
{
   // Logical function indicating whether or not point x,y lies
   // inside the region of the flat on inner wedge centered at phi0.

   double phi0rad = phi0 * M_PI/180;
   double phi1rad = (phi0 - 360) * M_PI/180;
   double phi = atan2(y,x);
   double r2 = x*x + y*y;
   return ((fabs(phi - phi0rad) > fWedgeHalfAngle * M_PI/180 &&
            fabs(phi - phi1rad) > fWedgeHalfAngle * M_PI/180) ||
           r2 < fInnerWedgeR[0] * fInnerWedgeR[0] ||
           r2 > fInnerWedgeR[1] * fInnerWedgeR[1])? 0 : 1;
}

int acol_image::hits_inner_wedge_hole(double x, double y, double phi0)
{
   // Logical function indicating whether or not point x,y lies
   // inside the region of the hole on inner wedge centered at phi0.

   double phi0rad = (phi0 - fWedgeHoleAngle) * M_PI/180;
   double phi1rad = (phi0 + fWedgeHoleAngle) * M_PI/180;
   double dx0 = fInnerHoleR * cos(phi0rad) - x;
   double dy0 = fInnerHoleR * sin(phi0rad) - y;
   double dx1 = fInnerHoleR * cos(phi1rad) - x;
   double dy1 = fInnerHoleR * sin(phi1rad) - y;
   return (sqrt(dx0 * dx0 + dy0 * dy0) > fWedgeHoleDiam / 2 &&
           sqrt(dx1 * dx1 + dy1 * dy1) > fWedgeHoleDiam / 2)? 0 : 1;
}

int acol_image::hits_inner_wedge_conn(double x, double y, double phi0)
{
   // Logical function indicating whether or not point x,y lies
   // inside the region of the connector on inner wedge centered at phi0.

   double phi0rad = (phi0 + fWedgeHoleAngle) * M_PI/180;
   double dx = fInnerConnR * cos(phi0rad) - x;
   double dy = fInnerConnR * sin(phi0rad) - y;
   return (sqrt(dx*dx + dy*dy) > fWedgeConnDiam / 2)? 0 : 1;
}

int acol_image::hits_outer_wedge_pins(double x, double y, double phi0)
{
   // Logical function indicating whether or not point x,y lies
   // inside the region of the pins on outer wedge centered at phi0.

   double phi0rad = phi0 * M_PI/180;
   double phi1rad = (phi0 - 360) * M_PI/180;
   double phi = atan2(y,x);
   return ((fabs(phi - phi0rad) > fWedgeHalfAngle * M_PI/180 &&
            fabs(phi - phi1rad) > fWedgeHalfAngle * M_PI/180) ||
           x*x + y*y < fOuterWedgeR[0] * fOuterWedgeR[0] ||
           x*cos(phi0rad) + y*sin(phi0rad) > fOuterWedgeX2)? 0 : 1;
}

int acol_image::hits_outer_wedge_flat(double x, double y, double phi0)
{
   // Logical function indicating whether or not point x,y lies
   // inside the region of the flat on outer wedge centered at phi0.

   double phi0rad = phi0 * M_PI/180;
   double phi1rad = (phi0 - 360) * M_PI/180;
   double phi = atan2(y,x);
   double r2 = x*x + y*y;
   return ((fabs(phi - phi0rad) > fWedgeHalfAngle * M_PI/180 &&
            fabs(phi - phi1rad) > fWedgeHalfAngle * M_PI/180) ||
           r2 < fOuterWedgeR[0] * fOuterWedgeR[0] ||
           r2 > fOuterWedgeR[1] * fOuterWedgeR[1])? 0 : 1;
}

int acol_image::hits_outer_wedge_hole(double x, double y, double phi0)
{
   // Logical function indicating whether or not point x,y lies
   // inside the region of the hole on outer wedge centered at phi0.

   double phi0rad = (phi0 - fWedgeHoleAngle) * M_PI/180;
   double phi1rad = (phi0 + fWedgeHoleAngle) * M_PI/180;
   double dx0 = fOuterHoleR * cos(phi0rad) - x;
   double dy0 = fOuterHoleR * sin(phi0rad) - y;
   double dx1 = fOuterHoleR * cos(phi1rad) - x;
   double dy1 = fOuterHoleR * sin(phi1rad) - y;
   return (sqrt(dx0 * dx0 + dy0 * dy0) > fWedgeHoleDiam / 2 &&
           sqrt(dx1 * dx1 + dy1 * dy1) > fWedgeHoleDiam / 2)? 0 : 1;
}

int acol_image::hits_outer_wedge_conn(double x, double y, double phi0)
{
   // Logical function indicating whether or not point x,y lies
   // inside the region of the connector on outer wedge centered at phi0.

   double phi0rad = (phi0 + fWedgeHoleAngle) * M_PI/180;
   double dx = fOuterConnR * cos(phi0rad) - x;
   double dy = fOuterConnR * sin(phi0rad) - y;
   return (sqrt(dx*dx + dy*dy) > fWedgeConnDiam / 2)? 0 : 1;
}

int acol_image::hits_cover_plate(double x, double y)
{
   // Logical function indicating whether or not point x,y lies
   // inside the region of the aluminum cover plate.

   double r = sqrt(x*x + y*y);
   return (r < fCollimatorDiam / 2 || r > fCoverPlateDiam / 2)? 0 : 1;
}
