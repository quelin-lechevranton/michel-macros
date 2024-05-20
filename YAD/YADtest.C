#include "YADtools.h"

void YADtest() {

    YAD::Reco R("~/Code/out/pdvd_1k_mu_1GeV_YAdumped.root");

    cout << R.GetEntries() << endl;
    for (int i=0; i<10; i++) {
        R.GetEntry(i);
        cout << R.TrkLength->at(0) << endl;
    }

}