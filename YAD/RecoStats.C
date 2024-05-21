#include "YADtools.h"

const size_t n_file=2;
const vector<string> filelist = yad::ReadFileList(n_file,"ijclab.list");

void RecoStats() {

    size_t  N_evt=0,      N_trk=0,
            N_trkless=0,  N_trkone=0,   N_trkful=0,
            N_trklong=0,  N_trkshort=0, N_trkmuless=0,
            N_muless=0,   N_muone=0,    N_muful=0;


    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Reco R(filename.c_str());
        yad::Truth T(filename.c_str());

        size_t n_evt = R.GetEntries();
        N_evt += n_evt;
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            R.GetEntry(i_evt);
            T.GetEntry(i_evt);

            size_t n_trklong=0, n_trk=R.NTrk;
            N_trk += n_trk;

            switch (n_trk) {
                case 0: N_trkless++; break;
                case 1: N_trkone++;  break;
                default: 
                    N_trkful++; 
                    // for (size_t i_trk=0; i_trk < R.NTrk; i_trk++) {
                    for (size_t i_pfp=0; i_pfp < R.NPfp; i_pfp++) {

                        int i_trk = R.PfpTrkID->at(i_pfp);
                        if (i_trk < 0) {continue;}

                        if (R.TrkLength->at(i_trk) > 25) {
                            N_trklong++;
                            n_trklong++;
                        }
                        else {
                            N_trkshort++;
                        }

                        if (R.PfpPdg->at(i_pfp)!=13) {
                            N_trkmuless++;
                        }
                    }
            }
            if (n_trklong > 1) { cout << "Multiple long tracks ! at Event #" << i_evt << endl;}
            
            size_t n_mu = count(R.PfpPdg->begin(), R.PfpPdg->end(), 13);
            switch (n_mu) {
                case 0:  N_muless++; break;
                case 1:  N_muone++;  break;
                default: N_muful++;
            }

            // cout << "Event#" << i_evt+1 << ": #pfp=" << R.NPfp << "\t#muon=" << count(R.PfpPdg->begin(), R.PfpPdg->end(), 13) << "\t#genmu=" << count(T.PrtPdg->begin(), T.PrtPdg->end(), 13) << endl;

        } //end event loop
    } //end file loop


    cout    << "#Event:     " << N_evt << endl
            << "#TrkLess:   " << N_trkless << " (" << 100.*N_trkless/N_evt << "%)" << endl
            << "#TrkOne:    " << N_trkone << " (" << 100.*N_trkone/N_evt << "%)" << endl
            << "#TrkFul:    " << N_trkful << " (" << 100.*N_trkful/N_evt << "%)" << endl
            << "#TrkLong:   " << N_trklong << " (" << 100.*N_trklong/N_trkful << "%)" << endl
            << "#TrkShort:  " << N_trkshort << " (" << 100.*N_trkshort/N_trkful << "%)" << endl
            << "#TrkMuLess: " << N_trkmuless << " (" << 100.*N_trkmuless/N_trk << "%)" << endl
            << "#MuLess:    " << N_muless << " (" << 100.*N_muless/N_evt << "%)" << endl
            << "#MuOne:     " << N_muone << " (" << 100.*N_muone/N_evt << "%)" << endl
            << "#MuFul:     " << N_muful << " (" << 100.*N_muful/N_evt << "%)" << endl;

}