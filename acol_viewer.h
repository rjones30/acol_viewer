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

#ifndef acol_viewer_h
#define acol_viewer_h

#include <TROOT.h>
#include <TCanvas.h>
#include <TArc.h>
#include <TGaxis.h>
#include <TColor.h>
#include <TPolyLine.h>
#include <TPaveText.h>

class acol_viewer {
 public:
   static const int kIXplus = 0;
   static const int kIXminus = 1;
   static const int kIYplus = 2;
   static const int kIYminus = 3;
   static const int kOXplus = 4;
   static const int kOXminus = 5;
   static const int kOYplus = 6;
   static const int kOYminus = 7;

   double fPedestal[8];      // adc counts
   double fRMSnoise[8];      // adc counts
   double fGainFactor[8];
   double fBaseline;         // normalized counts
   double fFullScaleInner;   // normalized counts
   double fFullScaleOuter;   // normalized counts
   double fEffRadiusInner;   // mm
   double fEffRadiusOuter;   // mm
   double fMaximumRadius;    // mm
   int fInitialized;

   TCanvas *fCanvas;
   int fCanvasWidth;
   int fCanvasHeight;
   TArc *fAperture;
   TArc *fBeamSpot;
   TArc *fBeamCircle;
   TGaxis *fXhair;
   TGaxis *fYhair;
   TPaveText *fTimestamp;
   TPaveText *fScaleMax;
   TPaveText *fCentInner;
   TPaveText *fCentOuter;
   int fSegmentsPerArc;
   double fApertureRadius;
   double fInnerWedgeRadius[2];
   double fOuterWedgeRadius[2];
   double fBeamSpotSize;
   double fBeamCircleSize;
   TPolyLine *fWedgeBow[8];
   int fNumberColors;

   class record {
    public:
      Long64_t tsec;
      Long64_t tnsec;
      float wIXplus;   // inner X+
      float wIXminus;  // inner X-
      float wIYplus;   // inner Y+
      float wIYminus;  // inner Y-
      float wOXplus;   // outer X+
      float wOXminus;  // outer X-
      float wOYplus;   // outer Y+
      float wOYminus;  // outer Y-
   };

   acol_viewer();
   ~acol_viewer();
   int setZero(record &data);
   int setNoise(record &data);
   int setGains(record &data);
   double setBaseline(record &data);
   double setFullScale(record &data);
   int display(record &data);
   int getBeamCentroidInner(record &data, double &x_mm, double &y_mm);
   int getBeamCentroidOuter(record &data, double &x_mm, double &y_mm);

 protected:
   TPolyLine *makeWedge(int wedge, int color);
};

#endif
