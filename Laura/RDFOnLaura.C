vector<string> filelist = {
    "/afs/cern.ch/work/j/jquelinl/out/protodunevd_10_muon_500MeV_dumped.root"
    // "/afs/cern.ch/work/j/jquelinl/out/pdvd_100_muon_1GeV_dumped.root",
    // "/afs/cern.ch/work/j/jquelinl/out/protodunevd_100_muon_800MeV_dumped.root",
    // "/afs/cern.ch/work/j/jquelinl/out/pdvd_1k_muon_1500MeV_dumped.root",
    // "/afs/cern.ch/work/j/jquelinl/out/pdvd_1k_muon_2GeV_allangles_dumped.root"
};

void TrackLength();

void RDFOnLaura() {
    TrackLength();
}

void TrackLength() {
    TFile file(filelist[0]);
    TChain chain("Reco");
    TChain.Add(filelist);
    TTree *TReco = file.Get<TTree>("LauraPDumper/Reco");
    TTree *TTruth = file.Get<TTree>("LauraPDumper/Truth");
    // ROOT::RDataFrame RReco("LauraPDumper/Reco",filelist);
    // ROOT::RDataFrame RTruth("LauraPDumper/Truth",filelist);
    // TTree* TReco("LauraPDumper/Reco",filelist);
    ROOT::RDataFrame R(TReco.Merge(TTruth));


    auto graph R.Filter("pfpTrackID.size()==1").Graph("genTrueMomentum","pfpTrackLength");

    TCanvas* canvas = new TCanvas("c",""); 
    canvas->cd();
    graph->Draw("AP");
}