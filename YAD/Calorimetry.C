#include "YAD_tools.h"

const size_t n_file=3;
const vector<string> filelist = yad::ReadFileList(n_file,"ijclab.list");

void Calorimetry() {

    TH2D* hMuon = new TH2D("hMuon",";Residual Range (cm);dE/dx (MeV/cm)",150,0,300,50,0,6);

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());
        yad::Reco R(filename.c_str());

        size_t n_evt = R.GetEntries();
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            // cout << "Event#" << i_evt << endl;

            T.GetEntry(i_evt);
            R.GetEntry(i_evt);

            for (size_t i_trk=0; i_trk < R.NTrk; i_trk++) {

                bool is_inside = true;
                for (size_t i_tpt=0; i_tpt < R.TrkNPt->at(i_trk); i_tpt++) {              

                    double X = (*R.TrkPtX)[i_trk][i_tpt];     
                    double Y = (*R.TrkPtY)[i_trk][i_tpt];     
                    double Z = (*R.TrkPtZ)[i_trk][i_tpt];     
                
                    bool x_inside = -320 <= X && X <= 350;
                    bool y_inside = -317 <= Y && Y <= 317;
                    bool z_inside = 20 <= Z && Z <= 280;

                    is_inside = is_inside && x_inside && y_inside && z_inside;

                    if (!is_inside) {break;}
                } //end of spt loop 
		if(!is_inside) continue;

                for (size_t i_cal=0; i_cal < R.TrkCalNPt->at(i_trk); i_cal++) {
                    hMuon->Fill((*R.TrkCalResRange)[i_trk][i_cal],(*R.TrkCaldEdx)[i_trk][i_cal]);
                }
            }

         
        } //end event loop
    } //end file loop

    TCanvas* c1 = new TCanvas("c1","dEdx");
    c1->cd();
    c1->cd();
    hMuon->Draw("colZ");
}
