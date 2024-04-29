// #include "../includes.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

// string path = "/silver/DUNE/quelin-lechevranton/out/";
// string file = "PDVD_100_muon_800MeV_LauraP_dumped.root";

int test(int);

void MacroOnLaura() {
    test(1);
}

int test(int i_event) {
    // TFile file(path+file);
    // TFile* file=TFile::Open("/silver/DUNE/quelin-lechevranton/out/PDVD_10_muon_500MeV_LauraP_dumped.root");
    TFile file("/eos/user/t/thoudy/pdvd/sims/out/PDVD_10_muon_500MeV_LauraP_dumped.root");
    // cout << "bonjour" << endl;
    TTree* Reco=(TTree*) file.Get("LauraPDumper/Reco");
    cout << "bonjour" << endl;

    int n_event = Reco->GetEntries();
    cout << n_event << endl;
    if (i_event<0 || i_event>n_event) {
        cout << "event index out of bound" << endl; 
        file.Close();
        return 1;
    }

    int j_total=0;

    // // vector<vector<double>*> TrackStart;
    // // vector<vector<double>*> TrackEnd;

    // // Reco->SetBranchAddress("pfpTrackStartX",&(TrackStart[0]));
    // // Reco->SetBranchAddress("pfpTrackStartY",&(TrackStart[1]));
    // // Reco->SetBranchAddress("pfpTrackStartZ",&(TrackStart[2]));
    // // Reco->SetBranchAddress("pfpTrackEndX",&(TrackEnd[0]));
    // // Reco->SetBranchAddress("pfpTrackEndY",&(TrackEnd[1]));
    // // Reco->SetBranchAddress("pfpTrackEndZ",&(TrackEnd[2]));

    vector<double> *TrackStartX, *TrackStartY, *TrackStartZ, *TrackEndX, *TrackEndY, *TrackEndZ;
    
    Reco->SetBranchAddress("pfpTrackStartX",&TrackStartX);
    Reco->SetBranchAddress("pfpTrackStartY",&TrackStartY);
    Reco->SetBranchAddress("pfpTrackStartZ",&TrackStartZ);
    Reco->SetBranchAddress("pfpTrackEndX",&TrackEndX);
    Reco->SetBranchAddress("pfpTrackEndY",&TrackEndY);
    Reco->SetBranchAddress("pfpTrackEndZ",&TrackEndZ);

    Reco->GetEntry(i_event);

    // cout << Reco->GetEntries() << endl;
    // Int_t kiki=(Int_t) i_event;
    // Reco->GetEntry(kiki);
    // for (Int_t iev=0; iev < Reco->GetEntries(); ++iev) { //Loop over the events
    //     Reco->GetEntry(iev);
    //     cout<<"Treating event number: "<<iev<<endl;
    // }
    // file->Close();

    vector<TGraph2D*> graph(2);
    graph[0] = new TGraph2D();
    graph[0]->SetName("Track Start Positions");
    graph[0]->SetMarkerColor(kBlue);

    graph[1] = new TGraph2D();
    graph[1]->SetName("Track End Positions");
    graph[1]->SetMarkerColor(kRed);

    graph[0]->GetXaxis()->SetTitle("X (cm)");
    graph[0]->GetYaxis()->SetTitle("Y (cm)");
    graph[0]->GetZaxis()->SetTitle("Z (cm)");

    for (int j=0; j<-TrackStartX->size(); j++) {        
        graph[0]->SetPoint(
            j_total++,
            TrackStartX->at(j),
            TrackStartY->at(j),
            TrackStartZ->at(j)
        );
        graph[1]->SetPoint(
            j_total++,
            TrackEndX->at(j),
            TrackEndY->at(j),
            TrackEndZ->at(j)
        );
    }

    stringstream title;
    title << "Tracks Ends of Event #" << i_event;
    TCanvas* canvas = new TCanvas("c",title.str().c_str());
    canvas->cd();
    graph[0]->Draw("AP");
    graph[1]->Draw("P");     

    file.Close();
    return 0;
}
