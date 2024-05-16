#include <vector>
#include <string>

using namespace std;

// const vector<string> filelist = {
//     "~/Code/out/protodunevd_10_muon_500MeV_Jdumped.root",
//     "~/Code/out/protodunevd_100_muon_800MeV_Jdumped.root",
//     "~/Code/out/pdvd_1k_muon_2GeV_allangles_Jdumped.root",
//     "~/Code/out/pdvd_1k_muon_1500MeV_Jdumped.root",
//     "~/Code/out/pdvd_100_muon_1GeV_Jdumped.root"
// };

const vector<string> filelist = {
    "/silver/DUNE/quelin-lechevranton/out/protodunevd_10_muon_500MeV_Jdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/protodunevd_100_muon_800MeV_Jdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/pdvd_1k_muon_2GeV_allangles_Jdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/pdvd_1k_muon_1500MeV_Jdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/pdvd_100_muon_1GeV_Jdumped.root"
};

void Macro() {

    vector<double>  *X=nullptr;

    // Hist ==================================================
    vector<TH1D*> h(3);
    int n_bin=30, x_min=0, x_max=700;

    h[0] = new TH1D("h0","",n_bin,x_min,x_max);

    h[1] = new TH1D("h1","",n_bin,x_min,x_max);

    h[2] = new TH1D("h2","",n_bin,x_min,x_max);

    THStack* hsTrkLen = new THStack("hs",";X;Y");
    for (int i=0; i<h.size(); i++) { hs->Add(h[i]);}

    int nEvent=0;

    for (int i_file=0; i_file < filelist.size(); i_file++) {
        string filename = filelist[i_file];

        cout << "\e[3mOpening file #" << i_file+1 << ": " << filename << "\e[0m" << endl;

        TFile file(filename.c_str());
        TTree *Reco = file.Get<TTree>("JDumper/Reco");
        TTree *Truth = file.Get<TTree>("JDumper/Truth");
        int n_ev=Reco->GetEntries();
        nEvent+=n_ev;

        Reco->SetBranchAddress("", &);

        Truth->SetBranchAddress("", &);

        for (int i_ev=0; i_ev < n_ev; i_ev++) {

            Reco->GetEntry(i_ev);
            Truth->GetEntry(i_ev);

            for(int i_trk=0; i_trk < nTracks ; i_trk++) {

            }
        
        } //end of event loop

        file.Close();
    
    } //end of file loop

    cout << "nEvent=" << nEvent << endl;

    TCanvas* c1 = new TCanvas("c1","");
    c1->cd();
}