#include "OmegaLight_tools.h"

const size_t n_file=1;
const vector<string> filelist = omega::ReadFileList(n_file,"list/pdvd_cosmics.list");

void RecoStats() {

    size_t  N_evt=0,      N_trk=0,
            N_trkless=0,  N_trkone=0,   N_trkful=0,
            N_trklong=0,  N_trkshort=0, N_trkmuless=0,
            N_muless=0,   N_muone=0,    N_muful=0;


    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        omega::Reco R(filename.c_str());
        //omega::Truth T(filename.c_str());

        size_t n_evt = R.N;
        N_evt += n_evt;
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {
            
            R.GetEvtPfp(i_evt);
            //T.GetEvt(i_evt);

            size_t n_trklong=0, n_trk=R.Trk.N;
            N_trk += n_trk;

            switch (n_trk) {
                case 0: N_trkless++; break;
                case 1: N_trkone++;  break;
                default: 
                    N_trkful++; 
                    // for (size_t i_trk=0; i_trk < R.NTrk; i_trk++) {
                    
                    for (size_t i_pfp=0; i_pfp < R.Pfp.N; i_pfp++) {

                        if (!R.Pfp.isTrk[i_pfp]) continue;

                        R.GetPfpTrk(i_pfp);

                        if (R.Trk.Length > 25) {
                            N_trklong++;
                            n_trklong++;
                        }
                        else {
                            N_trkshort++;
                        }
                    }
            }
            if (n_trklong > 1) { cout << "Multiple long tracks ! at Event #" << i_evt << endl;}
            
            // cout << "Event#" << i_evt+1 << ": #pfp=" << R.NPfp << "\t#muon=" << count(R.PfpPdg->begin(), R.PfpPdg->end(), 13) << "\t#genmu=" << count(T.PrtPdg->begin(), T.PrtPdg->end(), 13) << endl;

        } //end event loop
    } //end file loop


    cout    << "#Event:     " << N_evt << endl
            << "#TrkLess:   " << N_trkless << " (" << 100.*N_trkless/N_evt << "%)" << endl
            << "#TrkOne:    " << N_trkone << " (" << 100.*N_trkone/N_evt << "%)" << endl
            << "#TrkFul:    " << N_trkful << " (" << 100.*N_trkful/N_evt << "%)" << endl
            << "#TrkLong:   " << N_trklong << " (" << 100.*N_trklong/N_trkful << "%)" << endl
            << "#TrkShort:  " << N_trkshort << " (" << 100.*N_trkshort/N_trkful << "%)" << endl;

}
