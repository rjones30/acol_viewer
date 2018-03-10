/*
 * acol_image - a modeling and visualization class for representing
 *              the internal geometry of the GlueX active collimator.
 *
 * author: richard.t.jones at uconn.edu
 * version: november 24, 2014
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
#include <TH2D.h>
#include <TChain.h>

class acol_image {
 public:
   TCanvas *fCanvas;
   int fCanvasWidth;
   int fCanvasHeight;
   TArc *fArc;
   TPolyLine *fPline;
   double fCollimatorDiam;
   double fCoverPlateDiam;
   double fInnerWedgeR[2];
   double fOuterWedgeR[2];
   double fInnerWedgeX2;
   double fOuterWedgeX2;
   double fWedgeHalfAngle;
   double fInnerHoleR;
   double fInnerConnR;
   double fOuterHoleR;
   double fOuterConnR;
   double fWedgeHoleDiam;
   double fWedgeHoleAngle;
   double fWedgeConnDiam;
   int fWedgeAngleNseg;

   double fXscanOffset;
   double fYscanOffset;
   TCanvas *fCanvas1;

   // Coupling factor matrix between different wedge elements,
   // where the coupling elements to an inner wedge are numbered as follows.
   //   0 - inner wedge pins region, self
   //   1 - inner wedge flat region, self
   //   2 - inner wedge hole region, self
   //   3 - inner wedge hole with connector, self
   //   4 - inner wedge pins region, adjacent
   //   5 - inner wedge flat region, adjacent
   //   6 - inner wedge hole region, adjacent
   //   7 - inner wedge hole with connector, adjacent
   //   8 - inner wedge pins region, opposite
   //   9 - inner wedge flat region, opposite
   //  10 - inner wedge hole region, opposite
   //  11 - inner wedge hole with connector, opposite
   //  12 - inout wedge pins region, same sector
   //  13 - inout wedge flat region, same sector
   //  14 - inout wedge hole region, same sector
   //  15 - inout wedge hole with connector, same sector
   //  16 - inout wedge pins region, opposite
   //  17 - inout wedge flat region, opposite
   //  18 - inout wedge hole region, opposite
   //  19 - inout wedge hole with connector, opposite
   //  20 - inout wedge pins region, opposite
   //  21 - inout wedge flat region, opposite
   //  22 - inout wedge hole region, opposite
   //  23 - inout wedge hole with connector, opposite
   //  24 - outer wedge pins region, same sector
   //  25 - outer wedge flat region, same sector
   //  26 - outer wedge hole region, same sector
   //  27 - outer wedge hole with connector, same sector
   //  28 - outer wedge pins region, opposite
   //  29 - outer wedge flat region, opposite
   //  30 - outer wedge hole region, opposite
   //  31 - outer wedge hole with connector, opposite
   //  32 - outer wedge pins region, opposite
   //  33 - outer wedge flat region, opposite
   //  34 - outer wedge hole region, opposite
   //  35 - outer wedge hole with connector, opposite
   //  36 - cover plate to inner wedge
   //  37 - cover plate to outer wedge
   static const int fNpars = 38;
   double fCoupling[fNpars];
   double fChiSquare;
   double fXsolution;
   double fYsolution;

   acol_image();
   ~acol_image();
   void PaintScan(TChain *ch);
   double Beam_dIdr2(double r2);
   TH2D *BeamConvolution(TH2D *hin);
   TH2D *BeamConvolution2(TH2D *hin);
   TH2D *MapInnerXplus(TH2D *htemp=0);
   TH2D *MapInnerXminus(TH2D *htemp=0);
   TH2D *MapInnerYplus(TH2D *htemp=0);
   TH2D *MapInnerYminus(TH2D *htemp=0);
   TH2D *MapOuterXplus(TH2D *htemp=0);
   TH2D *MapOuterXminus(TH2D *htemp=0);
   TH2D *MapOuterYplus(TH2D *htemp=0);
   TH2D *MapOuterYminus(TH2D *htemp=0);

   int FitScanData(TChain *tscantree);
   TH1D *PlotScanData(TChain *tscantree, int wedge);
   TH2D *PlotBeamSpot(TH2D *htemp=0);
   TH1D *PlotBeamSpotProfile(TH2D *htemp=0, double xoff=0);
   TH2D *PlotInnerXasym(TH2D *htemp=0);
   TH2D *PlotInnerYasym(TH2D *htemp=0);
   TH2D *PlotOuterXasym(TH2D *htemp=0);
   TH2D *PlotOuterYasym(TH2D *htemp=0);
   TH2D *PlotInnerSolution(double xasym, double yasym);
   TH2D *PlotOuterSolution(double xasym, double yasym);
   void CheckScanData(TChain *tscantree);
   void Print();

   void resetFactor(int index1, int index2=0);
   void setFactor(int index, double value);
   double getFactor(int index);
   void setScanOffsets(double xoffset, double yoffset);
   void save_couplings(char *filename);
   void load_couplings(char *filename=0);

 protected:
   void draw_wedge(double rlim[2], double x2, double hangle, double phi0);
   void paint_stripe(int npoints, double *x, double *y);
   int hits_inner_wedge_pins(double x, double y, double phi0);
   int hits_inner_wedge_flat(double x, double y, double phi0);
   int hits_inner_wedge_hole(double x, double y, double phi0);
   int hits_inner_wedge_conn(double x, double y, double phi0);
   int hits_outer_wedge_pins(double x, double y, double phi0);
   int hits_outer_wedge_flat(double x, double y, double phi0);
   int hits_outer_wedge_hole(double x, double y, double phi0);
   int hits_outer_wedge_conn(double x, double y, double phi0);
   int hits_cover_plate(double x, double y);
};

inline double acol_image::Beam_dIdr2(double r2)
{
   // Intensity profile of the photon beam as a function of radius
   // squared r2 (mm^2) from the beam center. This shape contains
   // both the intrinsic bremsstrahlung characteristic spread and
   // also the effect of the electron beam emittance and multiple
   // scattering in the radiator target. Intensity is normalized
   // to unity as r2 -> 0.

   double z = r2 / 25.;
   double dIdr2 = 1 / (1 + 3.8 * z + 1.11 * pow(z, 1.9));
   return dIdr2;
}

inline void acol_image::resetFactor(int index1, int index2)
{
   for (int i = index1; (i == index1 || i <= index2) && i < fNpars; ++i)
      fCoupling[i] = 0;
}

inline void acol_image::setFactor(int index, double value)
{
   if (index >= 0 && index < fNpars)
      fCoupling[index] = value;
}

inline double acol_image::getFactor(int index)
{
   if (index >= 0 && index < fNpars)
      return fCoupling[index];
   return 0;
}

inline void acol_image::setScanOffsets(double xoffset, double yoffset)
{
   fXscanOffset = xoffset;
   fYscanOffset = yoffset;
}

#endif
