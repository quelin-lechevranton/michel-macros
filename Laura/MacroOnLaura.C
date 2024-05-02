// #include "../includes.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

vector<string> filelist = {
    // "/silver/DUNE/quelin-lechevranton/out/PDVD_10_muon_500MeV_LauraP_dumped.root",
    "/afs/cern.ch/work/j/jquelinl/out/protodunevd_10_muon_500MeV_dumped.root",
    "/afs/cern.ch/work/j/jquelinl/out/pdvd_100_muon_1GeV_dumped.root",
    "/afs/cern.ch/work/j/jquelinl/out/protodunevd_100_muon_800MeV_dumped.root",
    "/afs/cern.ch/work/j/jquelinl/out/pdvd_1k_muon_1500MeV_dumped.root"
};

void TrackEnds();
void Clusters();

void MacroOnLaura() {
    TrackEnds();
    // Clusters();
}

void Clusters() {

    int n_bin=50, x_min=0, x_max=2000;
    vector<TH1D*> histo(2);
    histo[0] = new TH1D("hSum",";SummedADC/Width;count",n_bin,x_min,x_max);
    histo[0]->SetLineColor(kRed+1);
    histo[0]->SetLineWidth(2);
    // histo[0]->SetName("muon CluSummedADC on collection");
    // histo[0]->GetXaxis()->SetTitle("SummedADC");
    // histo[0]->GetXaxis()->SetMaximum(1000);
    // histo[0]->GetYaxis()->SetTitle("count");

    histo[1] = new TH1D("hInt",";Integral/Width;count",n_bin,x_min,x_max);
    histo[1]->SetLineColor(kBlue-3);
    histo[1]->SetLineWidth(2);
    // histo[1]->SetName("muon CluIntegral on collection");
    // histo[0]->GetXaxis()->SetTitle("Integral");
    // histo[0]->GetXaxis()->SetMaximum(1000);
    

    //POURQUOI PAS UN POINTER ?
    unsigned int nParticles=0; 

    vector<int> *nClusters=nullptr,
        *TrackID=nullptr,
        *PdgCode=nullptr;

    // vector<double> *TrackStartDirectionX=nullptr,
    //     *TrackStartDirectionY=nullptr,
    //     *TrackStartDirectionZ=nullptr,
    //     *TrackVertexDirectionX=nullptr,
    //     *TrackVertexDirectionY=nullptr,
    //     *TrackVertexDirectionZ=nullptr;

    vector<vector<double>*> TrackEnd = {nullptr,nullptr,nullptr};

    vector<vector<double>> *CluPlane=nullptr,
        *CluView=nullptr,
        *CluNHits=nullptr,
        *CluSummedADC=nullptr,
        *CluIntegral=nullptr,
        *CluWidth=nullptr;    


    for (string filename : filelist) {

        TFile file(filename.c_str());
        TTree *Reco = (TTree*) file.Get("LauraPDumper/Reco");
        Int_t n_event=Reco->GetEntries();

        cout << "Opening: " << filename << "==============" << endl;


        Reco->SetBranchAddress("pfpTrackEndX",   &(TrackEnd[0]));
        Reco->SetBranchAddress("pfpTrackEndY",   &(TrackEnd[1]));
        Reco->SetBranchAddress("pfpTrackEndZ",   &(TrackEnd[2]));

        Reco->SetBranchAddress("pfpNClusters",   &nClusters);
        Reco->SetBranchAddress("nPFParticles",   &nParticles);
        Reco->SetBranchAddress("pfpTrackID",     &TrackID);
        Reco->SetBranchAddress("pfpPdgCode",     &PdgCode);

        Reco->SetBranchAddress("pfpCluPlane",    &CluPlane);
        Reco->SetBranchAddress("pfpCluView",     &CluView);
        Reco->SetBranchAddress("pfpCluNHits",    &CluNHits);
        Reco->SetBranchAddress("pfpCluSummedADC",&CluSummedADC);
        Reco->SetBranchAddress("pfpCluIntegral", &CluIntegral);
        Reco->SetBranchAddress("pfpCluWidth",    &CluWidth);



        double Sum,Int,Width;


        for (Int_t i_event=0; i_event < n_event; i_event++) {
            // cout << "Event #" << i_event << ": ";

            cout << "Event #" << i_event << endl;

            Reco->GetEntry(i_event);
            // cout << "\tn_track=" << TrackID->size();
            // cout << "\tn_particle=" << nParticles;
            // cout << "\tn_particle=" << nClusters->size();

            for (i_track) {

            }

            for(Int_t i_part=0; i_part < nParticles; i_part++) {
                
                // cout << "\tPart #" << i_part << endl;

                if (PdgCode->at(i_part)!=13) continue;

                cout << "\tMuon #" << i_part << endl;

                for (Int_t i_clu=0; i_clu < CluPlane->at(i_part).size(); i_clu++) {

                    // cout << "\t\tClu #" << i_clu;

                    if(CluPlane->at(i_part)[i_clu]!=0) continue;

                    cout << "\t\tCollec.Clu #" << i_clu;

                    Sum=CluSummedADC->at(i_part)[i_clu];
                    Int=CluIntegral->at(i_part)[i_clu];
                    Width=CluWidth->at(i_part)[i_clu];

                    cout << "\t|\tSum=" << Sum << "\tInt=" << Int << "\tWidth=" << Width << endl;

                    histo[0]->Fill(Sum/Width);
                    histo[1]->Fill(Int/Width);
                    // cout << histo[0]->Fill(Sum) << " " << histo[1]->Fill(Int) << endl;
                } //end of cluster loop
            } //end of particule loop

            for(int i_track : *TrackID) {


            } //end of track loop
            // cout << endl;
        } //end of event loop

        file.Close();

    } //end of file loop

    // Reco->GetEntry(9);
    // cout << "\tn_track=" << TrackID->size();
    // cout << "\tn_particle=" << nParticles;
    // cout << "\tn_particle=" << nClusters->size();

    // for (int i=0; i<TrackID->size(); i++) {
    //     cout << TrackID->at(i) << " ";
    // }
    // cout << endl;

    cout << histo[0]->GetEntries() << endl;
    cout << histo[1]->GetEntries() << endl;

    auto canvas = new TCanvas("c1","muon dE/dx on collection");
    canvas->cd();
    histo[1]->Draw("hist");
    histo[0]->Draw("samehist");
    // Reco->Draw("pfpCluSummedADC/pfpCluWidth","pfpPdgCode==13 && pfpCluPlane==0");
    canvas->SaveAs("Cluster.root");
    canvas->SaveAs("Cluster.pdf");

}


void TrackEnds() {

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

    Int_t j_total=0;

    vector<vector<double>*> TrackStart = {nullptr,nullptr,nullptr};
    vector<vector<double>*> TrackEnd = {nullptr,nullptr,nullptr};


    for (string filename : filelist) {

        TFile file(filename.c_str());
        TTree *Reco = (TTree*) file.Get("LauraPDumper/Reco");
        Int_t n_event=Reco->GetEntries();

        cout << "Opening: " << filename << "==============" << endl;

        // if (i_event<0 || i_event>n_event) {
        //     cout << "event index out of bound" << endl; 
        //     file.Close();
        //     return 1;
        // }

        Reco->SetBranchAddress("pfpTrackStartX", &(TrackStart[0]));
        Reco->SetBranchAddress("pfpTrackStartY", &(TrackStart[1]));
        Reco->SetBranchAddress("pfpTrackStartZ", &(TrackStart[2]));
        Reco->SetBranchAddress("pfpTrackEndX",   &(TrackEnd[0]));
        Reco->SetBranchAddress("pfpTrackEndY",   &(TrackEnd[1]));
        Reco->SetBranchAddress("pfpTrackEndZ",   &(TrackEnd[2]));

        for (Int_t i_event=0; i_event < n_event; i_event++) {

            Reco->GetEntry(i_event);

            for (size_t j=0; j< TrackStart[0]->size(); j++) {        

                graph[0]->SetPoint(
                    j_total++,
                    TrackStart[0]->at(j),
                    TrackStart[1]->at(j),
                    TrackStart[2]->at(j)
                );
                graph[1]->SetPoint(
                    j_total++,
                    TrackEnd[0]->at(j),
                    TrackEnd[1]->at(j),
                    TrackEnd[2]->at(j)
                );
            }
        }

   

        file.Close();
    }

    stringstream title;
    // title << "Tracks Ends of Event #" << i_event;
    title << "Tracks Ends";
    TCanvas* canvas = new TCanvas("c",title.str().c_str());
    canvas->cd();
    graph[0]->Draw("AP");
    graph[1]->Draw("P");  
}
