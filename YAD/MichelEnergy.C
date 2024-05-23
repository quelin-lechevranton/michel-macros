#include "YADtools.h"

const size_t n_file=1;
const vector<string> filelist = yad::ReadFileList(n_file,"jeremy.list");

void MichelEnergy() {

    TH1D* hDep = new TH1D("hDep",";Energy (GeV?);count",10,0,0.5);
    TH1D* hCal = new TH1D("hCal",";dEdx (GeV/cm?);count",40,0,40);

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());
        yad::Reco R(filename.c_str());

        size_t n_evt = R.GetEntries();
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {
        // size_t i_evt=45; {

            T.GetEntry(i_evt);
            R.GetEntry(i_evt);
            
            for (size_t i_dep=0; i_dep < T.NDep; i_dep++) {
                // if (T.DepPdg->at(i_dep)!=11) {continue;}

                // hDep->Fill(T.DepE->at(i_dep));

            } //end deposit loop

            // for (size_t i_trk=0; i_trk < R.NTrk; i_trk++) {
            size_t i_trk=0; {

                for (size_t i_cal=0; i_cal < R.TrkCalNPt->at(i_trk); i_cal++) {
                    hCal->Fill((*R.TrkCaldEdx)[i_trk][i_cal]);
                } //end calorimetry loop
            } //end track loop
 
        } //end event loop
    } //end file loop

    TCanvas* c1 = new TCanvas("c1","MichelEnergy");
    c1->Divide(2,1);
    c1->cd(1);
    // hDep->Draw("hist");
    c1->cd(2);
    hCal->Draw("hist");
}