#include "YAD_tools.h"

const size_t n_file=3;
const vector<string> filelist = yad::readFileList(n_file,"list/ijclab.list");

void Calorimetry() {
    clock_t start_time=clock();

    TH2D* hMuon = new TH2D("hMuon",";Residual Range (cm);dE/dx (MeV/cm)",150,0,300,50,0,6);

    size_t N_trk=0;
    double avg_tpt=0;

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());
        yad::Reco R(filename.c_str());

        size_t n_evt = R.GetEntries();
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            cout << "Event#" << i_evt+1 << "/" << n_evt << "\r" << flush;

            T.GetEntry(i_evt);
            R.GetEntry(i_evt);

            for (size_t i_trk=0; i_trk < R.NTrk; i_trk++) {

                if(!yad::isInside(
                    R.TrkPtX->at(i_trk),
                    R.TrkPtY->at(i_trk),
                    R.TrkPtZ->at(i_trk),
                    -320, 350,
                    -317, 317,
                    20, 280
                )) continue;

                N_trk++;
                avg_tpt+=R.TrkNPt->at(i_trk);

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
    
    avg_tpt /= N_trk;
    cout << N_trk << " tracks with a average of " << avg_tpt << " points per track " << endl;
    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}
