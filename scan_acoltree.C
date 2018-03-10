//
// scan_acoltree.C - root script to read a text file containing the log
//                   of a "burt" scan of the beam over the face of the
//                   active collimator, and correlate the times in the
//                   log with records in the acoltree to form a new tree
//                   that merges beam position information with acoltree
//                   info.
//
// author: richard.t.jones at uconn.edu
// version: november 8, 2014
//
// Usage:
// The following sequence of commands will create a chain of the files
// listed in standard "ls -l" format in enchain.list and condense them
// into an acoltree tree, and then generate a time series plot of the
// voltage level from one of the wedges, identified in the string argument
// to the show_acoltree function.
//   root > .L enchain.C+g
//   root > .L acoltree.C+g
//   root > .L scan_acoltree.C+g
//   root > scan_acoltree("rawdata/ybeamscan2-11-06_rad_1e-4.txt");
// A side effect of the above command is to create a new root file called
// ybeamscan2-11-06_rad_1e-4.root in the working directory and save the 
// merged tree in it, called tscan. If the output root file already exists
// in the working directory, it is overwritten.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <math.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TChain.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>

static TCanvas *c1 = 0;

#include <acoltree.h>
#include <timescan.C>

static bool initialized = true;

// offsets set for 11/07 at 3:40 am
//double offset[8] = {-0.080, +0.045, -0.095, -0.060,
//                    -0.055, +0.065, +0.150, +0.040}; 

// offsets set for 11/07 at 10:00 am
//double offset[8] = {-0.080, +0.045, -0.095, -0.050,
//                    -0.025, +0.115, +0.195, +0.062}; 

// offsets set for 11/16 at 12:45 am
double offset[8] = {-0.120, +0.005, -0.095, -0.105,
                    -0.195, -0.133, -0.240, -0.195}; 

// gains set for 11/07 at 3:40 am
//double gain[8] = {1.12, 1.00, 1.10, 1.00,
//                  1.25, 1.00, 0.91, 1.00};

// gains set for 11/07 at 10:00 am
//double gain[8] = {1.02, 1.00, 1.10, 1.00,
//                  1.15, 1.00, 0.96, 1.00};

// gains set for 11/16 at 13:40 am
double gain[8] = {1.04, 1.00, 1.10, 1.00,
                  1.15, 1.00, 0.96, 1.00};

#define VOLTS_PER_ADC_COUNT (0.005)
#define CENTER_OF_5MM_COLLIMATOR (-111.70)

TH1 *scan_acoltree(const char *logfile)
{
   if (! initialized) {
      gSystem->CompileMacro("enchain.C");
      gSystem->CompileMacro("acoltree.C");
      initialized = true;
   }

   double pedestal[8];
   for (int i=0; i < 8; ++i)
      pedestal[i] = offset[i] * gain[i];

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

   std::string rootfile(logfile);
   size_t dot_txt = rootfile.rfind(".txt");
   if (dot_txt != rootfile.npos)
      rootfile.replace(dot_txt, 4, ".root");
   else
      rootfile += ".root";
   TFile tscanfile(rootfile.c_str(), "recreate");
   std::string logfilestr(logfile);
   logfilestr = "burt scan " + logfilestr;
   TTree *tscantree = new TTree("tscan", logfilestr.c_str());
   tscantree->Branch("tsum", &tsum.t_s, "t_s/D:I_nA/D"
                                        ":Vixp/D:Vixm/D:Viyp/D:Viym/D"
                                        ":Voxp/D:Voxm/D:Voyp/D:Voym/D",
                     32768);
   tscantree->Branch("tscan", &tscan.tscan_s, "tscan_s/D:Iscan_nA/D"
                                        ":Xscan_bpu_mm/D:Yscan_bpu_mm/D"
                                        ":Vscan_ixp/D:Vscan_ixm/D"
                                        ":Vscan_iyp/D:Vscan_iym/D"
                                        ":Vscan_oxp/D:Vscan_oxm/D"
                                        ":Vscan_oyp/D:Vscan_oym/D"
                                        ":ATrate_Hz/D",
                     32768);

   TFile tsumfile("acoltree.root");
   TTree *tsumtree = (TTree*)tsumfile.Get("tsum");
   if (tsumtree == 0) {
      std::cerr << "scan_acoltree error - unable to open acoltree.root"
                << std::endl;
      return 0;
   }
   tsumtree->SetBranchAddress("tsum", &tsum.t_s);

   logfilestr = std::string(logfile);
   logfilestr = "rawdata/" + logfilestr;
   std::ifstream inlog(logfilestr.c_str());
   if (! inlog.good()) {
      std::cerr << "scan_acoltree error - unable to open input log file "
                << logfilestr << std::endl;
      return 0;
   }
   Long64_t jentry = 0;
   tscan.tscan_s = 1e20;
   double tnext_s = 1e20;
   char readbuf[999];
   while (inlog.good()) {
      inlog.getline(readbuf, 999);
      std::string line(readbuf);
      if (line.find("Time: ") == 0) {
         tscan.tscan_s = tnext_s;
         tnext_s = timescan(line);
         if (tnext_s == -1)
            std::cerr << "timescan error in line \"" << line << "\""
                      << std::endl;
      }
      else if (line.find("IBCAD00CRCUR6 ") == 0)
         sscanf(readbuf, "%*s %*i %le", &tscan.Iscan_nA);
      else if (line.find("bpu_mean_x ") == 0)
         sscanf(readbuf, "%*s %*i %le", &tscan.Xscan_bpu_mm);
      else if (line.find("bpu_mean_y ") == 0)
         sscanf(readbuf, "%*s %*i %le", &tscan.Yscan_bpu_mm);
      else if (line.find("activeTarget:rate ") == 0)
         sscanf(readbuf, "%*s %*i %le", &tscan.ATrate_Hz);
      else if (line.find("IOCHDCOL:VMICADC1_1 ") == 0)
         sscanf(readbuf, "%*s %*i %le", &tscan.Vscan_ixp);
      else if (line.find("IOCHDCOL:VMICADC2_1 ") == 0)
         sscanf(readbuf, "%*s %*i %le", &tscan.Vscan_ixm);
      else if (line.find("IOCHDCOL:VMICADC3_1 ") == 0)
         sscanf(readbuf, "%*s %*i %le", &tscan.Vscan_iyp);
      else if (line.find("IOCHDCOL:VMICADC4_1 ") == 0)
         sscanf(readbuf, "%*s %*i %le", &tscan.Vscan_iym);
      else if (line.find("IOCHDCOL:VMICADC1_2 ") == 0)
         sscanf(readbuf, "%*s %*i %le", &tscan.Vscan_oxp);
      else if (line.find("IOCHDCOL:VMICADC2_2 ") == 0)
         sscanf(readbuf, "%*s %*i %le", &tscan.Vscan_oxm);
      else if (line.find("IOCHDCOL:VMICADC3_2 ") == 0)
         sscanf(readbuf, "%*s %*i %le", &tscan.Vscan_oyp);
      else if (line.find("IOCHDCOL:VMICADC4_2 ") == 0)
         sscanf(readbuf, "%*s %*i %le", &tscan.Vscan_oym);
      while (tnext_s > tscan.tscan_s) {
         if (jentry < tsumtree->GetEntries())
            tsumtree->GetEntry(jentry++);
         else
            break;
         if (jentry/100 * 100 == jentry) 
            std::cout << jentry << ": " << tsum.t_s
                      << std::setprecision(10)
                      << " (" << tscan.tscan_s << ", "
                      << tnext_s << ")\n";
         if (tsum.t_s > tscan.tscan_s && tsum.t_s < tnext_s) {
            double *V12 = &tsum.Vixp;
            double *V16 = &tscan.Vscan_ixp;
            for (int n=0; n < 8; ++n) {
               V12[n] *= gain[n];
               V12[n] -= pedestal[n];
               V16[n] /= 16;
               V16[n] *= VOLTS_PER_ADC_COUNT * gain[n];
               V16[n] -= pedestal[n];
            }
            tscantree->Fill();
            std::cout << "filled with " << tscan.tscan_s << " < " << tsum.t_s
                      << " < " << tnext_s << std::endl;
            std::cout << "  Iscan_nA = " << tscan.Iscan_nA << std::endl
                      << "  Xscan_bpu_mm = " << tscan.Xscan_bpu_mm << std::endl
                      << "  Yscan_bpu_mm = " << tscan.Yscan_bpu_mm << std::endl
                      << "  ATrate_Hz = " << tscan.ATrate_Hz << std::endl;
         }   
         else if (tsum.t_s > tnext_s) {
            tscan.tscan_s = tnext_s;
         }
      }
   }
   tscanfile.cd();
   tscantree->Write();

   double xlim[2] = {1.e20,0.};
   double ylim[2] = {1.e20,0.};
   long int nentries = tscantree->GetEntries();
   for (long int entry=0; entry < nentries; ++entry) {
      tscantree->GetEntry(entry);
      xlim[0] = (tscan.Xscan_bpu_mm < xlim[0])? tscan.Xscan_bpu_mm : xlim[0];
      xlim[1] = (tscan.Xscan_bpu_mm > xlim[1])? tscan.Xscan_bpu_mm : xlim[1];
      ylim[0] = (tscan.Yscan_bpu_mm < ylim[0])? tscan.Yscan_bpu_mm : ylim[0];
      ylim[1] = (tscan.Yscan_bpu_mm > ylim[1])? tscan.Yscan_bpu_mm : ylim[1];
   }

   gStyle->SetPadLeftMargin(0.15);
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 550, 500);

   xlim[0] = int(xlim[0] - 0.999);
   xlim[1] = int(xlim[1] + 0.999);
   ylim[0] = int(ylim[0] - 0.999);
   ylim[1] = int(ylim[1] + 0.999);
   ostringstream varexpr;
   varexpr << "Yscan_bpu_mm:Xscan_bpu_mm>>xyscan"
           << "(10," << xlim[0] << "," << xlim[1]
           << ",10," << ylim[0] << "," << ylim[1] << ")";
   tscantree->Draw(varexpr.str().c_str());

   TH1 *xyscan = (TH1*)c1->GetPrimitive("xyscan");
   if (xyscan) {
      std::string title(logfile);
      title = "scan coordinates for " + title;
      xyscan->SetTitle(title.c_str());
      xyscan->GetYaxis()->SetTitle("y (mm)");
      xyscan->GetYaxis()->SetTitleOffset(2.0);
      xyscan->GetXaxis()->SetTitle("x (mm)");
      c1->Update();
   }

   return xyscan;
}

TH1 *scan_thetree(const char *rootfile)
{
   if (! initialized) {
      gSystem->CompileMacro("enchain.C");
      gSystem->CompileMacro("acoltree.C");
      initialized = true;
   }

   double pedestal[8];
   for (int i=0; i < 8; ++i)
      pedestal[i] = offset[i] * gain[i];

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

   std::string inrootfile(rootfile);
   inrootfile = "rawdata/sscanData/" + inrootfile;
   std::string outrootfile(rootfile);
   size_t dot_txt = outrootfile.rfind(".asc");
   if (dot_txt != outrootfile.npos)
      outrootfile.replace(dot_txt, 9, ".root");
   else
      outrootfile += ".root";
   TFile tscanfile(outrootfile.c_str(), "recreate");
   std::string rootfilestr(rootfile);
   rootfilestr = "collimator scan " + rootfilestr;
   TTree *tscantree = new TTree("tscan", rootfilestr.c_str());
   tscantree->Branch("tsum", &tsum.t_s, "t_s/D:I_nA/D"
                                        ":Vixp/D:Vixm/D:Viyp/D:Viym/D"
                                        ":Voxp/D:Voxm/D:Voyp/D:Voym/D",
                     32768);
   tscantree->Branch("tscan", &tscan.tscan_s, "tscan_s/D:Iscan_nA/D"
                                        ":Xscan_bpu_mm/D:Yscan_bpu_mm/D"
                                        ":Vscan_ixp/D:Vscan_ixm/D"
                                        ":Vscan_iyp/D:Vscan_iym/D"
                                        ":Vscan_oxp/D:Vscan_oxm/D"
                                        ":Vscan_oyp/D:Vscan_oym/D"
                                        ":ATrate_Hz/D",
                     32768);

   struct thetree_t {
      Float_t         Entry;
      Float_t         US1_2_BOT_m1_VAL;
      Float_t         US1_2_BOT_m1_RBV;
      Float_t         IBCAD00CRCUR6;
      Float_t         PSC_disc_scaler_1;
      Float_t         PSC_disc_scaler_2;
      Float_t         PSC_disc_scaler_3;
      Float_t         PSC_disc_scaler_4;
      Float_t         PSC_disc_scaler_5;
      Float_t         PSC_disc_scaler_6;
      Float_t         PSC_disc_scaler_7;
      Float_t         PSC_disc_scaler_8;
      Float_t         PSC_disc_scaler_9;
      Float_t         PSC_disc_scaler_10;
      Float_t         PSC_disc_scaler_11;
      Float_t         PSC_disc_scaler_12;
      Float_t         PSC_disc_scaler_13;
      Float_t         PSC_disc_scaler_14;
      Float_t         PSC_disc_scaler_15;
      Float_t         PSC_disc_scaler_16;
      Float_t         HALO_g_col_l_rate;
      Float_t         HALO_g_col_t_rate;
      Float_t         HALO_g_col_r_rate;
      Float_t         HALO_g_col_b_rate;
      Float_t         HALO_g_tgt_l_rate;
      Float_t         HALO_g_tgt_t_rate;
      Float_t         HALO_g_tgt_r_rate;
      Float_t         HALO_g_tgt_b_rate;
      Float_t         IOCHDCOL_VMICADC1_1;
      Float_t         IOCHDCOL_VMICADC2_1;
      Float_t         IOCHDCOL_VMICADC3_1;
      Float_t         IOCHDCOL_VMICADC4_1;
      Float_t         IOCHDCOL_VMICADC1_2;
      Float_t         IOCHDCOL_VMICADC2_2;
      Float_t         IOCHDCOL_VMICADC3_2;
      Float_t         IOCHDCOL_VMICADC4_2;
      Float_t         ST_disc_scaler_1;
      Float_t         ST_disc_scaler_2;
      Float_t         ST_disc_scaler_3;
      Float_t         ST_disc_scaler_4;
      Float_t         ST_disc_scaler_5;
      Float_t         ST_disc_scaler_6;
      Float_t         ST_disc_scaler_7;
      Float_t         ST_disc_scaler_8;
      Float_t         ST_disc_scaler_9;
      Float_t         ST_disc_scaler_10;
      Float_t         ST_disc_scaler_11;
      Float_t         ST_disc_scaler_12;
      Float_t         ST_disc_scaler_13;
      Float_t         ST_disc_scaler_14;
      Float_t         ST_disc_scaler_15;
      Float_t         ST_disc_scaler_16;
      Float_t         ST_disc_scaler_17;
      Float_t         ST_disc_scaler_18;
      Float_t         ST_disc_scaler_19;
      Float_t         ST_disc_scaler_20;
      Float_t         ST_disc_scaler_21;
      Float_t         ST_disc_scaler_22;
      Float_t         ST_disc_scaler_23;
      Float_t         ST_disc_scaler_24;
      Float_t         ST_disc_scaler_25;
      Float_t         ST_disc_scaler_26;
      Float_t         ST_disc_scaler_27;
      Float_t         ST_disc_scaler_28;
      Float_t         ST_disc_scaler_29;
      Float_t         ST_disc_scaler_30;
      Float_t         activeTarget_rate;
      Float_t         IPMAD00_XPOS;
      Float_t         IPMAD00_YPOS;
      Float_t         bpu_mean_x;
      Float_t         bpu_mean_y;
   } thetree;

   TFile thetreefile(inrootfile.c_str());
   TTree *thetreetree = (TTree*)thetreefile.Get("theTree");
   if (thetreetree == 0) {
      std::cerr << "scan_thetree error - unable to open " << inrootfile
                << std::endl;
      return 0;
   }
   thetreetree->SetBranchAddress("Entry", &thetree.Entry);
   thetreetree->SetBranchAddress("US1_2_BOT_m1.VAL", &thetree.US1_2_BOT_m1_VAL);
   thetreetree->SetBranchAddress("US1_2_BOT_m1.RBV", &thetree.US1_2_BOT_m1_RBV);
   thetreetree->SetBranchAddress("IBCAD00CRCUR6", &thetree.IBCAD00CRCUR6);
   thetreetree->SetBranchAddress("PSC_disc_scaler_1", &thetree.PSC_disc_scaler_1);
   thetreetree->SetBranchAddress("PSC_disc_scaler_2", &thetree.PSC_disc_scaler_2);
   thetreetree->SetBranchAddress("PSC_disc_scaler_3", &thetree.PSC_disc_scaler_3);
   thetreetree->SetBranchAddress("PSC_disc_scaler_4", &thetree.PSC_disc_scaler_4);
   thetreetree->SetBranchAddress("PSC_disc_scaler_5", &thetree.PSC_disc_scaler_5);
   thetreetree->SetBranchAddress("PSC_disc_scaler_6", &thetree.PSC_disc_scaler_6);
   thetreetree->SetBranchAddress("PSC_disc_scaler_7", &thetree.PSC_disc_scaler_7);
   thetreetree->SetBranchAddress("PSC_disc_scaler_8", &thetree.PSC_disc_scaler_8);
   thetreetree->SetBranchAddress("PSC_disc_scaler_9", &thetree.PSC_disc_scaler_9);
   thetreetree->SetBranchAddress("PSC_disc_scaler_10", &thetree.PSC_disc_scaler_10);
   thetreetree->SetBranchAddress("PSC_disc_scaler_11", &thetree.PSC_disc_scaler_11);
   thetreetree->SetBranchAddress("PSC_disc_scaler_12", &thetree.PSC_disc_scaler_12);
   thetreetree->SetBranchAddress("PSC_disc_scaler_13", &thetree.PSC_disc_scaler_13);
   thetreetree->SetBranchAddress("PSC_disc_scaler_14", &thetree.PSC_disc_scaler_14);
   thetreetree->SetBranchAddress("PSC_disc_scaler_15", &thetree.PSC_disc_scaler_15);
   thetreetree->SetBranchAddress("PSC_disc_scaler_16", &thetree.PSC_disc_scaler_16);
   thetreetree->SetBranchAddress("HALO_g_col_l_rate", &thetree.HALO_g_col_l_rate);
   thetreetree->SetBranchAddress("HALO_g_col_t_rate", &thetree.HALO_g_col_t_rate);
   thetreetree->SetBranchAddress("HALO_g_col_r_rate", &thetree.HALO_g_col_r_rate);
   thetreetree->SetBranchAddress("HALO_g_col_b_rate", &thetree.HALO_g_col_b_rate);
   thetreetree->SetBranchAddress("HALO_g_tgt_l_rate", &thetree.HALO_g_tgt_l_rate);
   thetreetree->SetBranchAddress("HALO_g_tgt_t_rate", &thetree.HALO_g_tgt_t_rate);
   thetreetree->SetBranchAddress("HALO_g_tgt_r_rate", &thetree.HALO_g_tgt_r_rate);
   thetreetree->SetBranchAddress("HALO_g_tgt_b_rate", &thetree.HALO_g_tgt_b_rate);
   thetreetree->SetBranchAddress("IOCHDCOL_VMICADC1_1", &thetree.IOCHDCOL_VMICADC1_1);
   thetreetree->SetBranchAddress("IOCHDCOL_VMICADC2_1", &thetree.IOCHDCOL_VMICADC2_1);
   thetreetree->SetBranchAddress("IOCHDCOL_VMICADC3_1", &thetree.IOCHDCOL_VMICADC3_1);
   thetreetree->SetBranchAddress("IOCHDCOL_VMICADC4_1", &thetree.IOCHDCOL_VMICADC4_1);
   thetreetree->SetBranchAddress("IOCHDCOL_VMICADC1_2", &thetree.IOCHDCOL_VMICADC1_2);
   thetreetree->SetBranchAddress("IOCHDCOL_VMICADC2_2", &thetree.IOCHDCOL_VMICADC2_2);
   thetreetree->SetBranchAddress("IOCHDCOL_VMICADC3_2", &thetree.IOCHDCOL_VMICADC3_2);
   thetreetree->SetBranchAddress("IOCHDCOL_VMICADC4_2", &thetree.IOCHDCOL_VMICADC4_2);
   thetreetree->SetBranchAddress("ST_disc_scaler_1", &thetree.ST_disc_scaler_1);
   thetreetree->SetBranchAddress("ST_disc_scaler_2", &thetree.ST_disc_scaler_2);
   thetreetree->SetBranchAddress("ST_disc_scaler_3", &thetree.ST_disc_scaler_3);
   thetreetree->SetBranchAddress("ST_disc_scaler_4", &thetree.ST_disc_scaler_4);
   thetreetree->SetBranchAddress("ST_disc_scaler_5", &thetree.ST_disc_scaler_5);
   thetreetree->SetBranchAddress("ST_disc_scaler_6", &thetree.ST_disc_scaler_6);
   thetreetree->SetBranchAddress("ST_disc_scaler_7", &thetree.ST_disc_scaler_7);
   thetreetree->SetBranchAddress("ST_disc_scaler_8", &thetree.ST_disc_scaler_8);
   thetreetree->SetBranchAddress("ST_disc_scaler_9", &thetree.ST_disc_scaler_9);
   thetreetree->SetBranchAddress("ST_disc_scaler_10", &thetree.ST_disc_scaler_10);
   thetreetree->SetBranchAddress("ST_disc_scaler_11", &thetree.ST_disc_scaler_11);
   thetreetree->SetBranchAddress("ST_disc_scaler_12", &thetree.ST_disc_scaler_12);
   thetreetree->SetBranchAddress("ST_disc_scaler_13", &thetree.ST_disc_scaler_13);
   thetreetree->SetBranchAddress("ST_disc_scaler_14", &thetree.ST_disc_scaler_14);
   thetreetree->SetBranchAddress("ST_disc_scaler_15", &thetree.ST_disc_scaler_15);
   thetreetree->SetBranchAddress("ST_disc_scaler_16", &thetree.ST_disc_scaler_16);
   thetreetree->SetBranchAddress("ST_disc_scaler_17", &thetree.ST_disc_scaler_17);
   thetreetree->SetBranchAddress("ST_disc_scaler_18", &thetree.ST_disc_scaler_18);
   thetreetree->SetBranchAddress("ST_disc_scaler_19", &thetree.ST_disc_scaler_19);
   thetreetree->SetBranchAddress("ST_disc_scaler_20", &thetree.ST_disc_scaler_20);
   thetreetree->SetBranchAddress("ST_disc_scaler_21", &thetree.ST_disc_scaler_21);
   thetreetree->SetBranchAddress("ST_disc_scaler_22", &thetree.ST_disc_scaler_22);
   thetreetree->SetBranchAddress("ST_disc_scaler_23", &thetree.ST_disc_scaler_23);
   thetreetree->SetBranchAddress("ST_disc_scaler_24", &thetree.ST_disc_scaler_24);
   thetreetree->SetBranchAddress("ST_disc_scaler_25", &thetree.ST_disc_scaler_25);
   thetreetree->SetBranchAddress("ST_disc_scaler_26", &thetree.ST_disc_scaler_26);
   thetreetree->SetBranchAddress("ST_disc_scaler_27", &thetree.ST_disc_scaler_27);
   thetreetree->SetBranchAddress("ST_disc_scaler_28", &thetree.ST_disc_scaler_28);
   thetreetree->SetBranchAddress("ST_disc_scaler_29", &thetree.ST_disc_scaler_29);
   thetreetree->SetBranchAddress("ST_disc_scaler_30", &thetree.ST_disc_scaler_30);
   thetreetree->SetBranchAddress("activeTarget_rate", &thetree.activeTarget_rate);
   thetreetree->SetBranchAddress("IPMAD00.XPOS", &thetree.IPMAD00_XPOS);
   thetreetree->SetBranchAddress("IPMAD00.YPOS", &thetree.IPMAD00_YPOS);
   thetreetree->SetBranchAddress("bpu_mean_x", &thetree.bpu_mean_x);
   thetreetree->SetBranchAddress("bpu_mean_y", &thetree.bpu_mean_y);

   long int nentries = thetreetree->GetEntries();
   for (long int n=0; n < nentries; ++n) {
      thetreetree->GetEntry(n);
      std::cout << "Just got entry " << thetree.Entry << " from theTree"
                << " with xcollimator = " << thetree.US1_2_BOT_m1_VAL
                << std::endl;
      float *adcvalue = &thetree.IOCHDCOL_VMICADC1_1;
      double *Vraw = &tscan.Vscan_ixp;
      double *Vtransf = &tsum.Vixp;
      for (int i=0; i < 8; ++i) {
         Vtransf[i] = Vraw[i] = adcvalue[i];
         Vtransf[i] /= 16;
         Vtransf[i] *= VOLTS_PER_ADC_COUNT * gain[i];
         Vtransf[i] -= pedestal[i];
      }
      tsum.t_s = tscan.tscan_s = 0;
      tsum.I_nA = tscan.Iscan_nA = thetree.IBCAD00CRCUR6;
      tscan.Xscan_bpu_mm = thetree.bpu_mean_x;
      tscan.Yscan_bpu_mm = thetree.bpu_mean_y;
      tsum.I_nA = tscan.Xscan_bpu_mm;
      tscan.Xscan_bpu_mm -= thetree.US1_2_BOT_m1_VAL - CENTER_OF_5MM_COLLIMATOR;
      tscan.ATrate_Hz = thetree.activeTarget_rate;
      if (fabs(tscan.Xscan_bpu_mm) < 50)
         tscantree->Fill();
   }
   tscanfile.cd();
   tscantree->Write();

   double xlim[2] = {1.e20,0.};
   double ylim[2] = {1.e20,0.};
   nentries = tscantree->GetEntries();
   for (long int entry=0; entry < nentries; ++entry) {
      tscantree->GetEntry(entry);
      xlim[0] = (tscan.Xscan_bpu_mm < xlim[0])? tscan.Xscan_bpu_mm : xlim[0];
      xlim[1] = (tscan.Xscan_bpu_mm > xlim[1])? tscan.Xscan_bpu_mm : xlim[1];
      ylim[0] = (tscan.Yscan_bpu_mm < ylim[0])? tscan.Yscan_bpu_mm : ylim[0];
      ylim[1] = (tscan.Yscan_bpu_mm > ylim[1])? tscan.Yscan_bpu_mm : ylim[1];
   }

   gStyle->SetPadLeftMargin(0.15);
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1", 40, 20, 550, 500);

   xlim[0] = int(xlim[0] - 0.999);
   xlim[1] = int(xlim[1] + 0.999);
   ylim[0] = int(ylim[0] - 0.999);
   ylim[1] = int(ylim[1] + 0.999);
   ostringstream varexpr;
   varexpr << "Yscan_bpu_mm:Xscan_bpu_mm>>xyscan"
           << "(10," << xlim[0] << "," << xlim[1]
           << ",10," << ylim[0] << "," << ylim[1] << ")";
   tscantree->Draw(varexpr.str().c_str());

   TH1 *xyscan = (TH1*)c1->GetPrimitive("xyscan");
   if (xyscan) {
      std::string title(rootfile);
      title = "scan coordinates for " + title;
      xyscan->SetTitle(title.c_str());
      xyscan->GetYaxis()->SetTitle("y (mm)");
      xyscan->GetYaxis()->SetTitleOffset(2.0);
      xyscan->GetXaxis()->SetTitle("x (mm)");
      c1->Update();
   }

   return xyscan;
}
