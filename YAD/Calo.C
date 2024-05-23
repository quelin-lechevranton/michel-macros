#include "YADtools.h"

const size_t n_file=1;
const vector<string> filelist = yad::ReadFileList(n_file,"ijclab.list");

void Calo() {

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());
        yad::Reco R(filename.c_str());

        size_t n_evt = R.GetEntries();
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            cout << "Event#" << i_evt << endl;

            T.GetEntry(i_evt);
            R.GetEntry(i_evt);

            for (size_t i_trk=0; i_trk < R.NTrk; i_trk++) {

                cout << "\tTrack#" << i_trk << "\tCalPlane" << R.TrkCalPlane->at(i_trk) << "\tCalNPt=" << R.TrkCalNPt->at(i_trk) << endl;

                // for (size_t i_tpt=0; i_tpt < R.TrkNPt->at(i_trk); i_tpt++) {

                // } //end trackpoint loop

                // for (size_t i_cal=0; i_cal < R.TrkCalNPt->at(i_trk); i_cal++) {

                // } //end calorimetry loop
            } //end track loop
 
            // for (size_t i_shw=0; i_shw < R.NSHw; i_shw++) {

            // } //end shower loop
        } //end event loop
    } //end file loop

    // TCanvas* c1 = new TCanvas("c1","Calo");
    // c1->cd();
}