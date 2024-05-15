#include <vector>
#include <string>

#include <ROOT/RDataFrame.hxx>
#include <TChain.h>

using namespace std;

const vector<string> filelist = {
    "/silver/DUNE/quelin-lechevranton/out/protodunevd_10_muon_500MeV_Jdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/protodunevd_100_muon_800MeV_Jdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/pdvd_1k_muon_2GeV_allangles_Jdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/pdvd_1k_muon_1500MeV_Jdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/pdvd_100_muon_1GeV_Jdumped.root"
};

// void Zenith();
void TrackLength();

void DataFrame() {

    // Zenith();
    TrackLength();

}

void TrackLength() {

    TChain CReco("JDumper/Reco");
    TChain CTruth("JDumper/Truth");

    for (const string filename : filelist) {
        CReco.Add(filename.c_str()); 
        CTruth.Add(filename.c_str()); 
    }

    ROOT::RDataFrame DFReco(CReco);
    ROOT::RDataFrame DFTruth(CTruth);

    int nEvent = *DFReco.Count();
    int nTrackless = *DFReco.Filter("fNTracks==0").Count();
    auto g1 = DFReco.Filter("fNTracks==1").Histo1D("fTrackLength");

    double max (vector<double> v);

    auto DFMaxLength = DFReco.Define("fTrackMaxLength", max, {"fTrackLength"});
    auto gm = DFMaxLength.Filter("fNTracks>1").Histo1D("fTrackMaxLength");

    cout << "nEvent=" << nEvent << endl;
    cout << "nTrackless=" << nTrackless << endl;

    g1->DrawClone();
    gm->DrawClone();

    // auto c1 = ROOT::RCanvas::Create("c1");
    // c1->Divide(2,1);
    // c1->cd(1);
    // c1->Draw(g1);
    // c1->cd(2);
    // c1->Draw(gm);

}

double max (vector<double> v) {
    return *max_element(v.begin(),v.end());
}


// void Zenith() {

//     TChain CReco("JDumper/Reco");
//     TChain CTruth("JDumper/Truth");

//     for (const string filename : filelist) {
//         CReco.Add(filename.c_str()); 
//         CTruth.Add(filename.c_str()); 
//     }

//     ROOT::RDataFrame DFReco(CReco);
//     ROOT::RDataFrame DFTruth(CTruth);

//     DFReco.Define("fTrackStartZenithX","TMath::RadToDeg()*TMath::ACos(-fTrackStartPx/fTrackStartP)");
//     auto g1 = DFReco1.Graph("fTrackStartP","fTrackLength");
//     auto g2 = DFReco1.Graph("fTrackStartZenithX","fTrackLength");

//     TCanvas* c1 = new TCanvas("c1","");
//     c1->Divide(2,1);
//     c1->cd(1);
//     g1->Draw("AP");
//     c1->cd(2);
//     g2->Draw("AP");
// }