// #include "../includes.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

// const char filename[] = "/silver/DUNE/quelin-lechevranton/out/PDVD_10_muon_500MeV_LauraP_dumped.root"
const char filename[] = "/eos/user/t/thoudy/pdvd/sims/out/PDVD_10_muon_500MeV_LauraP_dumped.root";

int TrackEnds(void);
int Clusters(void);

void MacroOnLaura() {
    // TrackEnds();
    Clusters();

}

int Clusters() {
    TFile file(filename);
    TTree *Reco = (TTree*) file.Get("LauraPDumper/Reco");
    Int_t n_event=Reco->GetEntries();

    //POURQUOI PAS UN POINTER ?
    unsigned int nParticles=0; 

    vector<int> *nClusters=nullptr,
        *TrackID=nullptr,
        *PdgCode=nullptr;

    vector<double> *TrackStartDirectionX=nullptr,
        *TrackStartDirectionY=nullptr,
        *TrackStartDirectionZ=nullptr,
        *TrackVertexDirectionX=nullptr,
        *TrackVertexDirectionY=nullptr,
        *TrackVertexDirectionZ=nullptr;

    vector<vector<double>> *CluPlane=nullptr,
        *CluView=nullptr,
        *CluNHits=nullptr,
        *CluSummedADC=nullptr,
        *CluIntegral=nullptr,
        *CluWidth=nullptr;

    Reco->SetBranchAddress("pfpTrackStartDirectionX",   &TrackStartDirectionX);

    // TBranch* B_nClusters = Reco->GetBranch("pfpNClusters");
    // B_nClusters->SetAddress(&nClusters);
    // B_nClusters->SetAutoDelete(true);

    Reco->SetBranchAddress("pfpNClusters",              &nClusters);
    Reco->SetBranchAddress("nPFParticles",              &nParticles);
    Reco->SetBranchAddress("pfpTrackID",                &TrackID);
    Reco->SetBranchAddress("pfpPdgCode",                &PdgCode);

    Reco->SetBranchAddress("pfpCluPlane",               &CluPlane);
    Reco->SetBranchAddress("pfpCluView",                &CluView);
    Reco->SetBranchAddress("pfpCluNHits",               &CluNHits);
    Reco->SetBranchAddress("pfpCluSummedADC",           &CluSummedADC);
    Reco->SetBranchAddress("pfpCluIntegral",            &CluIntegral);
    Reco->SetBranchAddress("pfpCluWidth",               &CluWidth);


    int n_bin=100, x_min=0, x_max=1000;
    vector<TH1F*> histo();
    histo[0] = new TH1F("hSum","CluSummedADC",n_bin,x_min,x_max);
    histo[0]->SetLineColor(kRed+1);
    // histo[0]->SetName("muon CluSummedADC on collection");
    // histo[0]->GetXaxis()->SetTitle("SummedADC");
    // histo[0]->GetXaxis()->SetMaximum(1000);
    // histo[0]->GetYaxis()->SetTitle("count");

    histo[1] = new TH1F("hInt","CluIntegral",n_bin,x_min,x_max);
    histo[0]->SetLineColor(kBlue-3);
    // histo[1]->SetName("muon CluIntegral on collection");
    // histo[0]->GetXaxis()->SetTitle("Integral");
    // histo[0]->GetXaxis()->SetMaximum(1000);


    for (Int_t i_event=0; i_event < n_event; i_event++) {
        cout << "Event #" << i_event << ": ";

        Reco->GetEntry(i_event);
        cout << "\tn_track=" << TrackID->size();
        cout << "\tn_particle=" << nParticles;
        cout << "\tn_particle=" << nClusters->size();

        for(Int_t i_part=0; i_part < PdgCode->size(); i_part++) {
            if (PdgCode[i_part]!=13) continue;
            for (Int_t i_clu=0; i_clu < CluPlane->at(i_part).size(); i_clu++) {
                if(CluPlane->at(i_part)[i_clu]!=0) continue;
                histo[0]->Fill(CluSummedADC->at(i_part)[i_clu]);
                histo[1]->Fill(CluIntegral->at(i_part)[i_clu]);
            }
        }

        cout << endl;
    } 


    // Reco->GetEntry(9);
    // cout << "\tn_track=" << TrackID->size();
    // cout << "\tn_particle=" << nParticles;
    // cout << "\tn_particle=" << nClusters->size();

    // for (int i=0; i<TrackID->size(); i++) {
    //     cout << TrackID->at(i) << " ";
    // }
    // cout << endl;


    auto canvas = new TCanvas("c1","muon dE/dx on collection");
    canvas->cd();
    histo[0]->Draw();
    histo[1]->Draw();
    // Reco->Draw("pfpCluSummedADC/pfpCluWidth","pfpPdgCode==13 && pfpCluPlane==0");

    return 0;



}




















int TrackEnds() {
    TFile* file = TFile::Open(filename);
    TTree* Reco=(TTree*) file->Get("LauraPDumper/Reco");

    // if (i_event<0 || i_event>n_event) {
    //     cout << "event index out of bound" << endl; 
    //     file.Close();
    //     return 1;
    // }

    Int_t j_total=0;

    // // vector<vector<double>*> TrackStart;
    // // vector<vector<double>*> TrackEnd;

    // // Reco->SetBranchAddress("pfpTrackStartX",&(TrackStart[0]));
    // // Reco->SetBranchAddress("pfpTrackStartY",&(TrackStart[1]));
    // // Reco->SetBranchAddress("pfpTrackStartZ",&(TrackStart[2]));
    // // Reco->SetBranchAddress("pfpTrackEndX",&(TrackEnd[0]));
    // // Reco->SetBranchAddress("pfpTrackEndY",&(TrackEnd[1]));
    // // Reco->SetBranchAddress("pfpTrackEndZ",&(TrackEnd[2]));

    vector<double> *TrackStartX=0, *TrackStartY=0, *TrackStartZ=0, *TrackEndX=0, *TrackEndY=0, *TrackEndZ=0;
    
    Reco->SetBranchAddress("pfpTrackStartX",&TrackStartX);
    Reco->SetBranchAddress("pfpTrackStartY",&TrackStartY);
    Reco->SetBranchAddress("pfpTrackStartZ",&TrackStartZ);
    Reco->SetBranchAddress("pfpTrackEndX",&TrackEndX);
    Reco->SetBranchAddress("pfpTrackEndY",&TrackEndY);
    Reco->SetBranchAddress("pfpTrackEndZ",&TrackEndZ);

    // vector<int> *TrackStartX=0;
    // Reco->SetBranchAddress("pfpTrackEndZ",&TrackStartX);
    // cout << (*TrackStartX).size() << endl;

    // for (Int_t iev=0; iev < Reco->GetEntries(); ++iev) { //Loop over the events
    //     Reco->GetEntry(iev);
    //     cout<<"Treating event number: "<<iev<<endl;
    // }
    // file->Close();

    vector<TGraph2D*> graph(2);
    graph[0] = new TGraph2D();
    graph[0]->SetName("Track Start Positions");
    graph[0]->SetMarkerColor(kBlue);
    graph[0]->SetMarkerStyle(20);
    graph[0]->SetMarkerSize(2);

    graph[1] = new TGraph2D();
    graph[1]->SetName("Track End Positions");
    graph[1]->SetMarkerColor(kRed);
    graph[1]->SetMarkerStyle(20);
    graph[1]->SetMarkerSize(2);

    graph[0]->GetXaxis()->SetTitle("X (cm)");
    graph[0]->GetYaxis()->SetTitle("Y (cm)");
    graph[0]->GetZaxis()->SetTitle("Z (cm)");


    Int_t n_event = Reco->GetEntries();
    for (Int_t i_event=0; i_event < n_event; i_event++) {
        Reco->GetEntry(i_event);

        for (size_t j=0; j< TrackStartX->size(); j++) {        
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
    }

    stringstream title;
    // title << "Tracks Ends of Event #" << i_event;
    title << "Tracks Ends";
    TCanvas* canvas = new TCanvas("c",title.str().c_str());
    canvas->cd();
    graph[0]->Draw("AP");
    graph[1]->Draw("P");     

    file->Close();
    return 0;
}
