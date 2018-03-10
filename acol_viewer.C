/*
 * acol_viewer - a visualization class for analyzing the currents
 *               in the 8 wedges of the GlueX photon beam active
 *               collimator and displaying the photon beam spot
 *               centroid when the currents meet the expected
 *               conditions.
 *
 * author: richard.t.jones at uconn.edu
 * version: october 24, 2014
 */

#include <acol_viewer.h>

#include <TBox.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTimeStamp.h>
#include <TPaveText.h>
#include <TText.h>

#include <iostream>
#include <math.h>

acol_viewer::acol_viewer() {
   /// Makes a new canvas for drawing the active collimator viewer
   /// and draws a blank image showing the outlines of the display.

   fInitialized = 0;

   fCanvasWidth = 500;
   fCanvasHeight = 500;
   fCanvas = (TCanvas*)gROOT->FindObject("acol_viewer_canvas");
   if (fCanvas)
      delete fCanvas;
   fCanvas = new TCanvas("acol_viewer_canvas",
                         "Hall D Active Collimator Viewer",
                         fCanvasWidth, fCanvasHeight);
   fCanvas->Range(-0.5,-0.6,0.6,0.5);

   // Draw collimator aperture with cross-hairs
   fApertureRadius = 0.25;
   fAperture = new TArc(0.,0.,fApertureRadius);
   fAperture->Draw();
   fXhair = new TGaxis(-fApertureRadius, 0., fApertureRadius, 0.,
                       -fApertureRadius, fApertureRadius, 505, "U+-");
   fXhair->Draw();
   fYhair = new TGaxis(0., -fApertureRadius, 0., fApertureRadius,
                       -fApertureRadius, fApertureRadius, 505, "U+-");
   fYhair->Draw();

   // No beam spot yet
   fBeamSpot = 0;
   fBeamCircle = 0;
   fBeamSpotSize = 0.02;
   fBeamCircleSize = 0.05;

   // No timestamp yet
   fTimestamp = 0;
   fScaleMax = 0;
   fCentInner = 0;
   fCentOuter = 0;

   // Draw images of wedges
   fSegmentsPerArc = 10;
   fInnerWedgeRadius[0] = 0.28;
   fInnerWedgeRadius[1] = 0.32;
   fOuterWedgeRadius[0] = 0.35;
   fOuterWedgeRadius[1] = 0.39;
   for (int i=0; i < 8; ++i) {
      fWedgeBow[i] = makeWedge(i, kRed - 10);
      fWedgeBow[i]->Draw("f");
      fWedgeBow[i]->Draw();
   }

   // Chose the color palette for wedge current display
   fNumberColors = 50;
   gStyle->SetPalette(1);
   double xbar[2] = {0.52, 0.58};
   double ybar[2] = {-0.5, 0.45};
   double dybar = (ybar[1] - ybar[0]) / fNumberColors;
   TBox *bar = new TBox();
   for (int i=0; i < fNumberColors; ++i) {
      bar->SetFillColor(gStyle->GetColorPalette(i));
      bar->DrawBox(xbar[0], ybar[0] + i * dybar, 
                  xbar[1], ybar[0] + (i+1) * dybar);
   }
   fCanvas->Update();

   // Add text labels for the wedges
   TText *text1 = new TText( 0.40, -0.02, "+X");
   TText *text2 = new TText(-0.03,  0.40, "+Y");
   TText *text3 = new TText(-0.45, -0.02, "-X");
   TText *text4 = new TText(-0.03, -0.45, "-Y");
   text1->Draw();
   text2->Draw();
   text3->Draw();
   text4->Draw();
   fCanvas->Update();

   // Default constants
   for (int i=0; i < 8; ++i) {
      fGainFactor[i] = 1;
   }
   fBaseline = 0;
   fFullScaleInner = 10;
   fFullScaleOuter = 10;
   fEffRadiusInner = 6.3;
   fEffRadiusOuter = 11.3;
   fMaximumRadius = 2.5;
}

acol_viewer::~acol_viewer() {
   delete fCanvas;
   delete fAperture;
   delete fXhair;
   delete fYhair;
   if (fBeamSpot)
      delete fBeamSpot;
   if (fBeamCircle)
      delete fBeamCircle;
   if (fTimestamp)
      delete fTimestamp;
   if (fScaleMax)
      delete fScaleMax;
   if (fCentInner)
      delete fCentInner;
   if (fCentOuter)
      delete fCentOuter;
   for (int i=0; i < 8; ++i)
      delete fWedgeBow[i];
}

int acol_viewer::setZero(record &data) {
   // The data provided are used to set pedestal currents in
   // the 8 wedges. This must be called before the viewer
   // can display position information. It is assumed to
   // represent the state of the inputs without beam present.

   fPedestal[kIXplus] = data.wIXplus;
   fPedestal[kIXminus] = data.wIXminus;
   fPedestal[kIYplus] = data.wIYplus;
   fPedestal[kIYminus] = data.wIYminus;
   fPedestal[kOXplus] = data.wOXplus;
   fPedestal[kOXminus] = data.wOXminus;
   fPedestal[kOYplus] = data.wOYplus;
   fPedestal[kOYminus] = data.wOYminus;
   fInitialized = fInitialized | 1;
   return 0;
}

int acol_viewer::setNoise(record &data) {
   // The data provided are used to set the rms noise level on 
   // the currents in the 8 wedges. This must be called before
   // the viewer can display position information. The rms noise
   // is used to form a threshold below which a sector is
   // considered to be seeing no beam flux.

   fRMSnoise[kIXplus] = data.wIXplus;
   fRMSnoise[kIXminus] = data.wIXminus;
   fRMSnoise[kIYplus] = data.wIYplus;
   fRMSnoise[kIYminus] = data.wIYminus;
   fRMSnoise[kOXplus] = data.wOXplus;
   fRMSnoise[kOXminus] = data.wOXminus;
   fRMSnoise[kOYplus] = data.wOYplus;
   fRMSnoise[kOYminus] = data.wOYminus;
   fInitialized = fInitialized | 2;
   return 0;
}

int acol_viewer::setGains(record &data) {
   // The data provided are used to set the gain factors on 
   // the currents in the 8 wedges. If not called, all of the
   // gain factors are assumed to be unity.

   fGainFactor[kIXplus] = data.wIXplus;
   fGainFactor[kIXminus] = data.wIXminus;
   fGainFactor[kIYplus] = data.wIYplus;
   fGainFactor[kIYminus] = data.wIYminus;
   fGainFactor[kOXplus] = data.wOXplus;
   fGainFactor[kOXminus] = data.wOXminus;
   fGainFactor[kOYplus] = data.wOYplus;
   fGainFactor[kOYminus] = data.wOYminus;
   return 0;
}

double acol_viewer::setBaseline(record &data) {
   // Use the least two of the outer wedges to estimate
   // the common-mode reverse current on the device.
   // Return value is the value for the inner wedges;

   double Imin[2] = {0,0};
   float *wouter = &data.wIXplus;
   for (int i=kOXplus; i <= kOYminus; ++i) {
      double Iouter = (wouter[i] - fPedestal[i]) * fGainFactor[i];
      if (Iouter < Imin[0]) {
         Imin[1] = Imin[0];
         Imin[0] = Iouter;
      }
      else if (Iouter < Imin[1]) {
         Imin[1] = Iouter;
      }
   }
   return fBaseline = (Imin[0] + Imin[1]) / 2;
}

double acol_viewer::setFullScale(record &data) {
   // Use the maximum currents on the wedges to set
   // the full-scale current for the display.

   double Imax = 0;
   float *wedges = &data.wIXplus;
   for (int i=kIXplus; i <= kIYminus; ++i) {
      double Iwedge = (wedges[i] - fPedestal[i]) * fGainFactor[i];
      Imax = (Iwedge > Imax)? Iwedge : Imax;
   }
   fFullScaleInner = (Imax < 1.0)? 1.0 :
                     (Imax < 3.0)? 3.0 :
                     (Imax < 5.0)? 5.0 : 10.;
   Imax = 0;
   for (int i=kOXplus; i <= kOYminus; ++i) {
      double Iwedge = (wedges[i] - fPedestal[i]) * fGainFactor[i];
      Imax = (Iwedge > Imax)? Iwedge : Imax;
   }
   fFullScaleOuter = (Imax < 1.0)? 1.0 :
                     (Imax < 3.0)? 3.0 :
                     (Imax < 5.0)? 5.0 : 10.;
   return fFullScaleInner;
}

int acol_viewer::display(record &data) {
   // Update the display with the contents of the data provided,
   // assuming that zero() has already been called.

   if (fInitialized != 3) {
      return 1;
   }

   setBaseline(data);
   if (fBeamSpot) {
      delete fBeamSpot;
      fBeamSpot = 0;
   }
   if (fBeamCircle) {
      delete fBeamCircle;
      fBeamCircle = 0;
   }

   // Draw the beam centroid: middle dot
   double xIn, yIn;
   if (getBeamCentroidInner(data, xIn, yIn) == 0) {
      fBeamSpot = new TArc(xIn/10, yIn/10, fBeamSpotSize);
      fBeamSpot->SetFillColor(kRed);
      fBeamSpot->Draw("f");
      fBeamSpot->Draw();
   }
 
   // Draw the beam centroid: outer circle
   double xOut, yOut;
   if (getBeamCentroidOuter(data, xOut, yOut) == 0) {
      fBeamCircle = new TArc(xOut/10, yOut/10, fBeamCircleSize);
      fBeamCircle->SetLineColor(kMagenta + 2);
      fBeamCircle->SetLineWidth(3);
      fBeamCircle->SetFillStyle(0);
      fBeamCircle->Draw();
   }

   // Color the wedges
   float *wedges = &data.wIXplus;
   for (int i=kIXplus; i <= kOYminus; ++i) {
      int color = 0;
      double Iwedge = (wedges[i] - fPedestal[i]) * fGainFactor[i];
      if (Iwedge > fBaseline + 5 * fRMSnoise[i]) {
         double fullscale = (i < kOXplus)? fFullScaleInner : fFullScaleOuter;
         int level = (Iwedge / (fullscale + 0.1)) * fNumberColors;
         color = gStyle->GetColorPalette(level);
      }
      fWedgeBow[i] = makeWedge(i, color);
      fWedgeBow[i]->Draw("f");
      fWedgeBow[i]->Draw();
   }

   // Update the timestamp 
   TTimeStamp tstamp((time_t)data.tsec, (int)data.tnsec);
   if (fTimestamp == 0)
      fTimestamp = new TPaveText(-0.5, -0.55, 0.5, -0.5, "NB");
   fTimestamp->Clear();
   fTimestamp->AddText(tstamp.AsString());
   fTimestamp->Draw();

   // Update the full scale
   char scalestr[99];
   sprintf(scalestr, "%2.0f V", fFullScaleInner);
   if (fScaleMax == 0)
      fScaleMax = new TPaveText(0.50, 0.45, 0.60, 0.50, "NB");
   fScaleMax->Clear();
   fScaleMax->AddText(scalestr);
   fScaleMax->Draw();
 
   // Update the numerical centroid displays
   char xInner[99], xOuter[99];
   sprintf(xInner, "inner: %6.3f,%6.3f mm", xIn, yIn);
   sprintf(xOuter, "outer: %6.3f,%6.3f mm", xOut, yOut);
   if (fCentInner == 0)
      fCentInner = new TPaveText(-0.50, 0.45, -0.18, 0.50, "");
   fCentInner->Clear();
   fCentInner->AddText(xInner);
   fCentInner->Draw();
   if (fCentOuter == 0)
      fCentOuter = new TPaveText(-0.50, 0.38, -0.18, 0.43, "");
   fCentOuter->Clear();
   fCentOuter->AddText(xOuter);
   fCentOuter->Draw();
 
   fCanvas->Update();
   return 0;
}

int acol_viewer::getBeamCentroidInner(record &data, double &x_mm, double &y_mm)
{
   // Analyzes the data record and tries to compute the centroid
   // position of the beam. If the currents are too small or outside
   // the ranges allowed by the algorithm, it returns a non-zero value.
 
   float *wedges = &data.wIXplus;
   double Inorm[4], Inorm_sigma[4];
   for (int i=kIXplus; i <= kIYminus; ++i) {
      Inorm[i] = (wedges[i] - fPedestal[i]) * fGainFactor[i];
      Inorm_sigma[i] = fRMSnoise[i] * fGainFactor[i];
      Inorm[i] -= fBaseline;
   }
   if (Inorm[kIXplus] < 3 * Inorm_sigma[kIXplus] ||
       Inorm[kIYplus] < 3 * Inorm_sigma[kIYplus])
   {
      return -1;
   }

   x_mm = (Inorm[kIXplus] - Inorm[kIXminus]) /
          (Inorm[kIXplus] + Inorm[kIXminus] + 1e-10) *
          fEffRadiusInner;
   y_mm = (Inorm[kIYplus] - Inorm[kIYminus]) /
          (Inorm[kIYplus] + Inorm[kIYminus] + 1e-10) *
          fEffRadiusInner;
   return (sqrt(x_mm * x_mm + y_mm * y_mm + 1e-10) > fMaximumRadius)? 1 : 0;
}

int acol_viewer::getBeamCentroidOuter(record &data, double &x_mm, double &y_mm)
{
   // Analyzes the data record and tries to compute the centroid
   // position of the beam. If the currents are too small or outside
   // the ranges allowed by the algorithm, it returns a non-zero value.
 
   float *wedges = &data.wIXplus;
   double Inorm[8], Inorm_sigma[8];
   for (int i=kOXplus; i <= kOYminus; ++i) {
      Inorm[i] = (wedges[i] - fPedestal[i]) * fGainFactor[i];
      Inorm_sigma[i] = fRMSnoise[i] * fGainFactor[i];
      Inorm[i] -= fBaseline;
   }
   if (Inorm[kOXplus] < 3 * Inorm_sigma[kOXplus] ||
       Inorm[kOYplus] < 3 * Inorm_sigma[kOYplus])
   {
      return -1;
   }

   x_mm = (Inorm[kOXplus] - Inorm[kOXminus]) /
          (Inorm[kOXplus] + Inorm[kOXminus] + 1e-10) *
          fEffRadiusOuter;
   y_mm = (Inorm[kOYplus] - Inorm[kOYminus]) /
          (Inorm[kOYplus] + Inorm[kOYminus] + 1e-10) *
          fEffRadiusOuter;
   return (sqrt(x_mm * x_mm + y_mm * y_mm + 1e-10) > fMaximumRadius)? 1 : 0;
}

TPolyLine *acol_viewer::makeWedge(int wedge, int color) {
   // Makes a filled arc on the active collimator viewer
   // outside the aperture to represent one of the 8 wedges,
   // using color to represent the current in the wedge.

   int sgn;
   double *r;
   if (wedge == kIXplus || wedge == kIYplus) {
      r = fInnerWedgeRadius;
      sgn = 1;
   }
   else if (wedge == kIXminus || wedge == kIYminus) {
      r = fInnerWedgeRadius;
      sgn = -1;
   }
   else if (wedge == kOXplus || wedge == kOYplus) {
      r = fOuterWedgeRadius;
      sgn = 1;
   }
   else if (wedge == kOXminus || wedge == kOYminus) {
      r = fOuterWedgeRadius;
      sgn = -1;
   }
   else {
      return 0;
   }
  
   int nseg = 2*fSegmentsPerArc+3;
   double x[nseg];
   double y[nseg];
   double phi = TMath::Pi() / 6.;
   double dphi = 2 * phi / fSegmentsPerArc;
   int vert = 0;
   for (int i = 0; i <= fSegmentsPerArc; ++i, ++vert) {
      x[vert] = sgn * r[0] * cos(phi);
      y[vert] = sgn * r[0] * sin(phi);
      phi -= dphi;
   }
   for (int i = 0; i <= fSegmentsPerArc; ++i, ++vert) {
      phi += dphi;
      x[vert] = sgn * r[1] * cos(phi);
      y[vert] = sgn * r[1] * sin(phi);
   }
   x[vert] = x[0];
   y[vert] = y[0];

   TPolyLine *pline;
   if (wedge == kIXplus || wedge == kIXminus ||
       wedge == kOXplus || wedge == kOXminus)
   {
      pline = new TPolyLine(nseg, x, y);
   }
   else {
      pline = new TPolyLine(nseg, y, x);
   }
   pline->SetFillColor(color);
   pline->SetLineColor(1);
   pline->SetLineWidth(1);
   return pline;
}
