#include <vector>
#include <string>
#include <iostream>

using namespace std;

const vector<string> filelist = {
    "~/Code/out/pdvd_1k_mu_1GeV_Jdumped.root",
};

void YADref() {

    vector<double>  *TrkLen=nullptr;


    for (int i_file=0; i_file < filelist.size(); i_file++) {
        string filename = filelist[i_file];

        cout << "\e[3mOpening file #" << i_file+1 << ": " << filename << "\e[0m" << endl;

        TFile file(filename.c_str());
        TTree *Reco = file.Get<TTree>("JDumper/Reco");
        int n_evt=Reco->GetEntries();

        Reco->SetBranchAddress("fTrackLength",  &TrkLen);

        for (int i_evt=0; i_evt < 10; i_evt++) {

            Reco->GetEntry(i_evt);
            cout << TrkLen->at(0) << endl;

        file.Close();

    } //end of file loop
}