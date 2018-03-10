#define acoltree_cxx
#include "acoltree.h"
#include <iostream>
#include <string>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#define VOLTS_PER_ADC_COUNT 0.005

void acoltree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L acoltree.C
//      Root > acoltree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   if (fChain == 0) return;

   struct summary_record {
      double t_s;
      double I_nA;
      double V_average;
   } summary;

   if (fSummary == 0) {
      std::string sname(fChain->GetName());
      sname = "h_" + sname;
      std::string stitle(fChain->GetTitle());
      stitle += " summary";
      fSummary = new TTree(sname.c_str(), stitle.c_str());
      fSummary->Branch("def", &summary.t_s, 
                       "t_s/D:I_nA/D:V_average/D", 32768);
   }

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      if (nb == 0) {
         std::cerr << "Error in acoltree::Loop() - "
                   << "GetEntry(" << jentry << ") returns 0, "
                   << "cannot continue." << std::endl;
         return;
      }
      nbytes += nb;
      summary.t_s = record.tsec + record.tnsec / 1.e9;
      summary.I_nA = record.current;
      double ave = 0;
      int n;
      for (n=0; n < 8192; ++n)
         ave += record.data[n];
      ave /= n;
      summary.V_average =  ave * VOLTS_PER_ADC_COUNT;
      fSummary->Fill();
   }
   fSummary->Write();
   std::cout << "acoltree::Loop() - created new tree " << fSummary->GetName()
             << " with " << fSummary->GetEntries() << " rows." << std::endl;
}
