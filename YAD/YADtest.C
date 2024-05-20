#include "YADtools.h"

void YADtest() {

    yad::Reco R("~/Code/out/pdvd_1k_mu_1GeV_YAdumped.root");

    for (int i=0; i<10; i++) {
        R.GetEntry(i);
        cout << "Event#" << i << endl;

        for(int i_trk=0; i_trk<R.NTrk; i_trk++){ 
            cout << "\t TrackID=" << R.TrkID->at(i_trk) << endl;
        }
    }
}