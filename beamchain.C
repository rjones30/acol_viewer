#include <TFile.h>
#include <TChain.h>
#include <TString.h>

#include <math.h>

TChain **beamchain()
{
   TChain **ch = new TChain*[11];
   TString fi[10];
   fi[0] = "beam_0026.root";
   fi[1] = "beam_0027.root";
   fi[2] = "beam_0028.root";
   fi[3] = "beam_0029.root";
   fi[4] = "beam_0030.root";
   fi[5] = "beam_0031.root";
   fi[6] = "beam_0032.root";
   fi[7] = "beam_0033.root";
   fi[8] = "beam_0035.root";
   fi[9] = "beam_0036.root";
   ch[0] = new TChain("tscan");
   for (int i=1; i <= 10; ++i) {
      ch[i] = new TChain("tscan");
      ch[i]->Add(fi[i-1]);
      ch[0]->Add(fi[i-1]);
   }
   ch[0]->Add("beam_const.root");
   return ch;
}

void make_beam_const()
{
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

   TFile *tscanfile = new TFile("beam_const.root", "recreate");
   TTree *tscantree = new TTree("tscan", "constraints");
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

   // add some fake scan entries to act as constraints

   double rhole[2] = {1.6 / 2 * 25.4, 4.2 / 2 * 25.4};
   double rconn[2] = {2.2 / 2 * 25.4, 4.5 / 2 * 25.4};
   for (double phi=0; phi < 360; phi += 90) {
      tsum.t_s = tscan.tscan_s = 0;
      tsum.I_nA = tscan.Iscan_nA = 0;
      double *Vmeas = &tsum.Vixp;
      double *Vraw = &tscan.Vscan_ixp;
      for (int i=0; i < 8; i++)
         Vraw[i] = 0;
       if (fabs(phi - 0) < 1e-10) {
          Vmeas[0] = 0.10;  Vmeas[4] = 0.05;
          Vmeas[1] = 0.20;  Vmeas[5] = 2.00;
          Vmeas[2] = 0.15;  Vmeas[6] = 0.10;
          Vmeas[3] = 0.15;  Vmeas[7] = 0.10;
       }
       else if (fabs(phi - 90) < 1e-10) {
          Vmeas[0] = 0.15;  Vmeas[4] = 0.10;
          Vmeas[1] = 0.15;  Vmeas[5] = 0.10;
          Vmeas[2] = 0.20;  Vmeas[6] = 2.00;
          Vmeas[3] = 0.10;  Vmeas[7] = 0.05;
      }
      else if (fabs(phi - 180) < 1e-10) {
          Vmeas[0] = 0.20;  Vmeas[4] = 2.00;
          Vmeas[1] = 0.10;  Vmeas[5] = 0.05;
          Vmeas[2] = 0.15;  Vmeas[6] = 0.10;
          Vmeas[3] = 0.15;  Vmeas[7] = 0.10;
       }
       else if (fabs(phi - 270) < 1e-10) {
          Vmeas[0] = 0.15;  Vmeas[4] = 0.10;
          Vmeas[1] = 0.15;  Vmeas[5] = 0.10;
          Vmeas[2] = 0.10;  Vmeas[6] = 0.05;
          Vmeas[3] = 0.20;  Vmeas[7] = 2.00;
      }
      tscan.ATrate_Hz = 0;
      double phirad = (phi - 15) * M_PI/180.;
      tscan.Xscan_bpu_mm = rhole[1] * cos(phirad);
      tscan.Yscan_bpu_mm = rhole[1] * sin(phirad);
      tscantree->Fill();
      phirad = (phi + 15) * M_PI/180.;
      tscan.Xscan_bpu_mm = rhole[1] * cos(phirad);
      tscan.Yscan_bpu_mm = rhole[1] * sin(phirad);
      tscantree->Fill();
   }
   tscantree->Write();
   tscanfile->Close();
}
