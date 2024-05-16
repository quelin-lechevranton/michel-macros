// #include "../includes.h"
// #include <vector>
// #include <string>
// #include <sstream>
// #include <iostream>

// vector<string> filelist = {
//     "/afs/cern.ch/work/j/jquelinl/out/protodunevd_10_muon_500MeV_dumped.root",
//     "/afs/cern.ch/work/j/jquelinl/out/pdvd_100_muon_1GeV_dumped.root",
//     "/afs/cern.ch/work/j/jquelinl/out/protodunevd_100_muon_800MeV_dumped.root",
//     "/afs/cern.ch/work/j/jquelinl/out/pdvd_1k_muon_1500MeV_dumped.root",
//     "/afs/cern.ch/work/j/jquelinl/out/pdvd_1k_muon_2GeV_allangles_dumped.root"
// };

vector<string> filelist = {
    "/silver/DUNE/quelin-lechevranton/out/protodunevd_10_muon_500MeV_LPdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/pdvd_100_muon_1GeV_LPdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/protodunevd_100_muon_800MeV_LPdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/pdvd_1k_muon_1500MeV_LPdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/pdvd_1k_muon_2GeV_allangles_LPdumped.root"
};

void TrackEnds();
void Clusters();
void TrackEff(int save=0);
void TrackEff2(int save=0);
void TrackAngle();
void TrackLength(int save=0);

void MacroOnLaura() {
    // TrackEnds();
    // Clusters();
    // TrackEff();
    // TrackEff2();
    // TrackAngle();
    TrackLength();
}


void TrackLength(int save=0) {

    vector<int> *rTrID = nullptr;

    vector<double> *tTrueEnergy = nullptr,
        *tTrueMomentum = nullptr,
        *rTrLen = nullptr;

    vector<vector<double>*> tStartMomentum = {nullptr, nullptr, nullptr};

    // TH2D* h = new TH2D("h",";Momentum (GeV);Total Track Length (cm);counts");
    TGraph* graph = new TGraph();
    graph->GetXaxis()->SetTitle("Momentum (GeV)");
    graph->GetYaxis()->SetTitle("Total Track Length (cm)");
    // graph->GetYaxis()->SetTitle("Max Track Length (cm)");
    graph->SetMarkerStyle(20);
    // graph->SetMarkerSize(1);

    int i_graph=0;

    for (int i_file=0; i_file < filelist.size(); i_file++) {
        string filename = filelist[i_file];

        cout << "\e[3mOpening file #" << i_file+1 << ": " << filename << "\e[0m" << endl;

        TFile file(filename.c_str());
        TTree *Reco = file.Get<TTree>("LauraPDumper/Reco");
        TTree *Truth = file.Get<TTree>("LauraPDumper/Truth");
        int n_ev=Reco->GetEntries();
        
        Reco->SetBranchAddress("pfpTrackID", &rTrID);

        // Truth->SetBranchAddress("genTrueEnergy", &tTrueEnergy);
        Truth->SetBranchAddress("genTrueMomentum", &tTrueMomentum);
        Reco->SetBranchAddress("pfpTrackLength", &rTrLen);

        Truth->SetBranchAddress("genStartMomentumX", &(tStartMomentum[0]));
        Truth->SetBranchAddress("genStartMomentumY", &(tStartMomentum[1]));
        Truth->SetBranchAddress("genStartMomentumZ", &(tStartMomentum[2]));

        for (int i_ev=0; i_ev < n_ev; i_ev++) {

            Reco->GetEntry(i_ev);
            Truth->GetEntry(i_ev);

            int n_gen = tTrueMomentum->size();
            if(!n_gen) continue;
            int i_gen =0;
            
            double MomentumZen=TMath::RadToDeg()*TMath::ACos(- tStartMomentum[0]->at(i_gen) / tTrueMomentum->at(i_gen));

            int n_tr=rTrID->size();
            double TrackLength=0;

            for(int i_tr=0; i_tr < n_tr; i_tr++) {
                double len = rTrLen->at(i_tr);
                // TrackLength = len > TrackLength ? len : TrackLength;
                TrackLength += len;
            } //end of track loop

            graph->SetPoint(i_graph++,tTrueMomentum->at(i_gen),TrackLength);
            // graph->SetPoint(i_graph++,MomentumZen,TrackLength);
        } //end of event loop
        file.Close();
    } //end of file loop

    TCanvas* canvas = new TCanvas("c","");
    canvas->cd();
    graph->Draw("AP");

    if(!save) return;

    canvas->SaveAs("TrackLength_Momentum.pdf");
    canvas->SaveAs("TrackLength_Momentum.root");
}

void TrackAngle() {

    vector<int> *rTrID = nullptr;

    vector<double> *tTrueEnergy = nullptr,
        *tTrueMomentum = nullptr,
        *rTrLen = nullptr,
        *rTrThe = nullptr,
        *rTrPhi = nullptr,
        *rTrZen = nullptr,
        *rTrAzi = nullptr;

    vector<vector<double>*> tStartMomentum = {nullptr, nullptr, nullptr},
        rTrStart = {nullptr, nullptr, nullptr},
        rTrEnd = {nullptr, nullptr, nullptr},
        rTrStartDir = {nullptr, nullptr, nullptr};

    // TGraph* graph = new TGraph();

    for (int i_file=0; i_file < filelist.size(); i_file++) {
        string filename = filelist[i_file];
        
        cout << "\e[3mOpening file #" << i_file+1 << ": " << filename << "\e[0m" << endl;

        TFile file(filename.c_str());
        TTree *Reco = (TTree*) file.Get("LauraPDumper/Reco");
        TTree *Truth = (TTree*) file.Get("LauraPDumper/Truth");
        int n_ev=Reco->GetEntries();
        
        Reco->SetBranchAddress("pfpTrackID", &rTrID);

        // Truth->SetBranchAddress("genTrueEnergy", &tTrueEnergy);
        Truth->SetBranchAddress("genTrueMomentum", &tTrueMomentum);
        Reco->SetBranchAddress("pfpTrackLength", &rTrLen);
        Reco->SetBranchAddress("pfpTrackTheta", &rTrThe);
        Reco->SetBranchAddress("pfpTrackPhi", &rTrPhi);
        Reco->SetBranchAddress("pfpTrackZenithAngle", &rTrZen);
        Reco->SetBranchAddress("pfpTrackAzimuthAngle", &rTrAzi);

        Truth->SetBranchAddress("genStartMomentumX", &(tStartMomentum[0]));
        Truth->SetBranchAddress("genStartMomentumY", &(tStartMomentum[1]));
        Truth->SetBranchAddress("genStartMomentumZ", &(tStartMomentum[2]));
        Reco->SetBranchAddress("pfpTrackStartX", &(rTrStart[0]));
        Reco->SetBranchAddress("pfpTrackStartY", &(rTrStart[1]));
        Reco->SetBranchAddress("pfpTrackStartZ", &(rTrStart[2]));
        Reco->SetBranchAddress("pfpTrackEndX", &(rTrEnd[0]));
        Reco->SetBranchAddress("pfpTrackEndY", &(rTrEnd[1]));
        Reco->SetBranchAddress("pfpTrackEndZ", &(rTrEnd[2]));
        Reco->SetBranchAddress("pfpTrackStartDirectionX", &(rTrStartDir[0]));
        Reco->SetBranchAddress("pfpTrackStartDirectionY", &(rTrStartDir[1]));
        Reco->SetBranchAddress("pfpTrackStartDirectionZ", &(rTrStartDir[2]));

        cout << string(100,'-') << endl;
        for (int i_ev=0; i_ev < n_ev; i_ev++) {

            cout << "Event #" << i_ev+1 << endl;
            Reco->GetEntry(i_ev);
            Truth->GetEntry(i_ev);

            int n_gen = tTrueMomentum->size();

            for(int i_gen=0; i_gen < n_gen; i_gen++) {

                double MomentumZen=TMath::RadToDeg()*TMath::ACos(- tStartMomentum[0]->at(i_gen) / tTrueMomentum->at(i_gen));
                // double MomentumMag=TMath::Sqrt(
                //     TMath::Power(tStartMomentum[0]->at(i_gen),2) + TMath::Power(tStartMomentum[1]->at(i_gen),2) + TMath::Power(tStartMomentum[2]->at(i_gen),2)
                // );
                // double Energy = tTrueEnergy->at(i_gen);
                double Momentum = tTrueMomentum->at(i_gen);
                // double Mass= TMath::Sqrt(TMath::Power(tTrueEnergy->at(i_gen),2)- TMath::Power(tTrueMomentum->at(i_gen),2));
                
                cout << "\tGen #" << i_gen << endl;
                // cout << "\t\tTrueEnergy =\t\t" << Energy << "\tGeV" << endl;
                cout << "\t\tTrueMomentum =\t\t" << Momentum << "\tGeV" << endl;
                // cout << "\t\tTrueMass =\t\t" << Mass << "\tGeV" << endl;
                // cout << "\t\tTrueMomentumMag =\t" << MomentumMag << "\tGeV" << endl;
                cout << "\t\tTrueMomentumZenithX =\t" << MomentumZen << "\tdeg" << endl;
            }

            int n_tr=rTrID->size();
            for(int i_tr=0; i_tr < n_tr; i_tr++) {

                // double DirMag=TMath::Sqrt(
                //     TMath::Power(rTrStartDir[0]->at(i_tr),2) + TMath::Power(rTrStartDir[1]->at(i_tr),2) + TMath::Power(rTrStartDir[2]->at(i_tr),2)
                // );
                double DirZenith=TMath::RadToDeg()*TMath::ACos(- rTrStartDir[0]->at(i_tr) );

                double Length = rTrLen->at(i_tr);
                double Zenith = TMath::RadToDeg()*rTrZen->at(i_tr);

                cout << "\tTrack #" << i_tr+1 << endl;
                cout << "\t\tRecoTrackLength =\t" << Length << "\tcm" << endl;
                // cout << "\t\tRecoTrackDirMag =\t" << DirMag << "\tcm" << endl;
                cout << "\t\tRecoTrackZenithY =\t" << Zenith << "\tdeg" << endl;
                cout << "\t\tRecoTrackZenithX =\t" << DirZenith << "\tdeg" << endl;
            } //end of track loop
        cout << string(100,'-') << endl;
        } //end of event loop
        file.Close();
    } //end of file loop

    // TCanvas* canvas = new TCanvas("c","");
    // canvas->cd();
    // graph->Draw("AP");
}

void TrackEff2(int save=0) {

    int N_ev=0;
    int n_ev=0;
    int n_tr=0;
    int N_trless=0;
    int N_trful=0;
    // int N_gen;
    int n_gen;

    vector<int> *rTrackID = nullptr;
    vector<vector<double>*> tStartMomentum = {nullptr, nullptr, nullptr};
    vector<double> *tTrueMomentum = nullptr;


    int n_xbin=60, x_min=0, x_max=90, n_ybin=30, y_min=0, y_max=4; 
    vector<TH2D*> h(3);
    h[0] = new TH2D("h","incoming muon;Zenith Angle (deg);Energy (GeV);count",n_xbin,x_min,x_max,n_ybin,y_min,y_max);
    h[0]->SetLineColor(kBlack);

    h[1] = new TH2D("hTrless","incoming muon;Zenith Angle (deg);Energy (GeV);count",n_xbin,x_min,x_max,n_ybin,y_min,y_max);
    h[1]->SetLineColor(kRed+1);
    h[1]->SetLineWidth(2);

    h[2] = new TH2D("hTrful","incoming muon;Zenith Angle (deg);Energy (GeV);count",n_xbin,x_min,x_max,n_ybin,y_min,y_max);
    h[2]->SetLineColor(kBlue-4);
    h[2]->SetLineWidth(2);

    for (int i_file=0; i_file < filelist.size(); i_file++) {
        string filename = filelist[i_file];

        cout << "\e[3mOpening file #" << i_file+1 << ": " << filename << "\e[0m" << endl;

        TFile file(filename.c_str());
        TTree* Truth = (TTree*) file.Get("LauraPDumper/Truth");
        TTree* Reco = (TTree*) file.Get("LauraPDumper/Reco");
        n_ev = Truth->GetEntries();
        N_ev += n_ev;

        Truth->SetBranchAddress("genStartMomentumX", &(tStartMomentum[0]));
        Truth->SetBranchAddress("genStartMomentumY", &(tStartMomentum[1]));
        Truth->SetBranchAddress("genStartMomentumZ", &(tStartMomentum[2]));
        Truth->SetBranchAddress("genTrueMomentum", &tTrueMomentum);

        Reco->SetBranchAddress("pfpTrackID", &rTrackID);

        for (int i_ev=0; i_ev<n_ev; i_ev++) {

            Reco->GetEntry(i_ev);
            Truth->GetEntry(i_ev);

            n_tr = rTrackID->size();
            n_gen = tTrueMomentum->size();

            // cout << "Event #" << i_ev << "\tn_gen=" << n_gen << endl;

            if (!n_gen) continue;
            // N_gen++;

            double momentum=tTrueMomentum->at(0);
            double zenith=TMath::RadToDeg()*TMath::ACos(- tStartMomentum[0]->at(0) / tTrueMomentum->at(0));

            h[0]->Fill(zenith,momentum);

            if (!n_tr) {
                N_trless++;
                h[1]->Fill(zenith,momentum);
            }
            if (n_tr>1) {
                N_trful++;
                h[2]->Fill(zenith,momentum);
            }
            

        }

    }   

    float completness=1 - (float) N_trless/N_ev;
    float purity=1 - (float) N_trful/N_ev;

    cout << "total event number:\t" << N_ev << endl; 
    // cout << "total genless event number: " << N_ev-N_gen << endl;
    cout << "total trackless event number:\t" << N_trless << endl;
    cout << "completness (1-Ntrless/Nev):\t" << completness << endl;

    cout << "total trackful event number:\t" << N_trful << endl;
    cout << "purity (1-Ntrful/Nev):\t" << purity << endl;

    auto canvas = new TCanvas("c1","");
    canvas->cd();
    h[0]->Draw("box");
    h[1]->Draw("samebox");
    h[2]->Draw("samebox");

    if (!save) return;
    canvas->SaveAs("TrackEff2.root");
    canvas->SaveAs("TrackEff2.pdf");
}


void TrackEff(int save=0) {

    int N_ev=0;
    int n_ev=0;
    int n_tr=0;
    int N_trless=0;
    int N_trful=0;
    // int N_gen;
    int n_gen;

    vector<int> *rTrackID = nullptr;
    vector<vector<double>*> tStartMomentum = {nullptr, nullptr, nullptr};
    vector<double> *tTrueMomentum = nullptr;


    int n_bin=30, x_min=0, x_max=90;
    vector<TH1D*> hZenith(3);
    hZenith[0] = new TH1D("hZenith","Zenith Angle of incoming muon;Zenith Angle (deg);count",n_bin,x_min,x_max);
    hZenith[0]->SetLineColor(kBlack);

    hZenith[1] = new TH1D("hTrlessZenith","Zenith Angle of Trackless Events;Zenith Angle (deg);count",n_bin,x_min,x_max);
    hZenith[1]->SetLineColor(kRed+1);
    hZenith[1]->SetLineWidth(2);

    hZenith[2] = new TH1D("hTrfulZenith","Zenith Angle of Trackful Events;Zenith Angle (deg);count",n_bin,x_min,x_max);
    hZenith[2]->SetLineColor(kBlue-4);
    hZenith[2]->SetLineWidth(2);

    n_bin=30; x_min=0; x_max=4; 
    vector<TH1D*> hEnergy(3);
    hEnergy[0] = new TH1D("hEnergy","Energy of incoming muon;Momentum (GeV);count",n_bin,x_min,x_max);
    hEnergy[0]->SetLineColor(kBlack);
    // hEnergy[0]->SetLineWidth(2);

    hEnergy[1] = new TH1D("hTrlessEnergy","Energy of Trackless Events;Momentum (GeV);count",n_bin,x_min,x_max);
    hEnergy[1]->SetLineColor(kRed+1);
    hEnergy[1]->SetLineWidth(2);

    hEnergy[2] = new TH1D("hTrfulEnergy","Energy of Trackful Events;Momentum (GeV);count",n_bin,x_min,x_max);
    hEnergy[2]->SetLineColor(kBlue-4);
    hEnergy[2]->SetLineWidth(2);

    for (int i_file=0; i_file < filelist.size(); i_file++) {
        string filename = filelist[i_file];

        cout << "\e[3mOpening file #" << i_file+1 << ": " << filename << "\e[0m" << endl;

        TFile file(filename.c_str());
        TTree* Truth = (TTree*) file.Get("LauraPDumper/Truth");
        TTree* Reco = (TTree*) file.Get("LauraPDumper/Reco");
        n_ev = Truth->GetEntries();
        N_ev += n_ev;

        Truth->SetBranchAddress("genStartMomentumX", &(tStartMomentum[0]));
        Truth->SetBranchAddress("genStartMomentumY", &(tStartMomentum[1]));
        Truth->SetBranchAddress("genStartMomentumZ", &(tStartMomentum[2]));
        Truth->SetBranchAddress("genTrueMomentum", &tTrueMomentum);

        Reco->SetBranchAddress("pfpTrackID", &rTrackID);

        for (int i_ev=0; i_ev<n_ev; i_ev++) {

            Reco->GetEntry(i_ev);
            Truth->GetEntry(i_ev);

            n_tr = rTrackID->size();
            n_gen = tTrueMomentum->size();

            // cout << "Event #" << i_ev << "\tn_gen=" << n_gen << endl;

            if (!n_gen) continue;
            // N_gen++;

            double momentum=tTrueMomentum->at(0);
            double zenith=TMath::RadToDeg()*TMath::ACos(- tStartMomentum[0]->at(0) / tTrueMomentum->at(0));

            hZenith[0]->Fill(zenith);
            hEnergy[0]->Fill(momentum);

            if (!n_tr) {
                N_trless++;
                hZenith[1]->Fill(zenith);
                hEnergy[1]->Fill(momentum);
            }
            if (n_tr>1) {
                N_trful++;
                hZenith[2]->Fill(zenith);
                hEnergy[2]->Fill(momentum);
            }
            

        }

    }   

    float completness=1 - (float) N_trless/N_ev;
    float purity=1 - (float) N_trful/N_ev;

    cout << "total event number:\t" << N_ev << endl; 
    // cout << "total genless event number: " << N_ev-N_gen << endl;
    cout << "total trackless event number:\t" << N_trless << endl;
    cout << "completness (1-Ntrless/Nev):\t" << completness << endl;

    cout << "total trackful event number:\t" << N_trful << endl;
    cout << "purity (1-Ntrful/Nev):\t" << purity << endl;

    auto canvas = new TCanvas("c1","");
    canvas->Divide(2,1);
    canvas->cd(1);
    hZenith[0]->Draw("hist");
    hZenith[1]->Draw("samehist");
    hZenith[2]->Draw("samehist");

    canvas->cd(2);
    hEnergy[0]->Draw("hist");
    hEnergy[1]->Draw("samehist");
    hEnergy[2]->Draw("samehist");

    if (!save) return;
    canvas->SaveAs("TrackEff.root");
    canvas->SaveAs("TrackEff.pdf");
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

            cout << "Event #" << i_event+1 << endl;

            Reco->GetEntry(i_event);
            // cout << "\tn_track=" << TrackID->size();
            // cout << "\tn_particle=" << nParticles;
            // cout << "\tn_particle=" << nClusters->size();

            // for (i_track) {

            // }

            for(Int_t i_part=0; i_part < nParticles; i_part++) {
                
                // cout << "\tPart #" << i_part+1 << endl;

                if (PdgCode->at(i_part)!=13) continue;

                cout << "\tMuon #" << i_part+1 << endl;

                for (Int_t i_clu=0; i_clu < CluPlane->at(i_part).size(); i_clu++) {

                    // cout << "\t\tClu #" << i_clu+1;

                    if(CluPlane->at(i_part)[i_clu]!=0) continue;

                    cout << "\t\tCollec.Clu #" << i_clu+1;

                    Sum=CluSummedADC->at(i_part)[i_clu];
                    Int=CluIntegral->at(i_part)[i_clu];
                    Width=CluWidth->at(i_part)[i_clu];

                    cout << "\t|\tSum=" << Sum << "\tInt=" << Int << "\tWidth=" << Width << endl;

                    histo[0]->Fill(Sum/Width);
                    histo[1]->Fill(Int/Width);
                    // cout << histo[0]->Fill(Sum) << " " << histo[1]->Fill(Int) << endl;
                } //end of cluster loop
            } //end of particule loop

            // for(int i_track : *TrackID) {


            // } //end of track loop
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

    // cout << histo[0]->GetEntries() << endl;
    // cout << histo[1]->GetEntries() << endl;

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

        cout << "\e[3mOpening file: " << filename << "\e[0m" << endl;

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

        cout << string(100,'-') << endl; 
        for (Int_t i_event=0; i_event < n_event; i_event++) {

            cout << "Event #" << i_event << endl; 
            Reco->GetEntry(i_event);

            for (size_t j=0; j< TrackStart[0]->size(); j++) {        

                graph[0]->SetPoint(
                    j_total,
                    TrackStart[0]->at(j),
                    TrackStart[1]->at(j),
                    TrackStart[2]->at(j)
                );

                cout << "\tTrack #" << j << endl; 
                cout << "\t\t( " << TrackStart[0]->at(j) << " , " <<
                                TrackStart[1]->at(j) << " , " <<
                                TrackStart[2]->at(j) << " ) --> ( " <<
                                TrackEnd[0]->at(j) << " , " <<
                                TrackEnd[1]->at(j) << " , " <<
                                TrackEnd[2]->at(j) << " )" << endl;

                graph[1]->SetPoint(
                    j_total++,
                    TrackEnd[0]->at(j),
                    TrackEnd[1]->at(j),
                    TrackEnd[2]->at(j)
                );
            }
            cout << string(100,'-') << endl;
        }
        file.Close();
    }

    cout << j_total << endl;
    cout << graph[0]->GetN() << endl;
    cout << graph[1]->GetN() << endl;

    TCanvas* canvas = new TCanvas("c","Track Ends");
    canvas->cd();
    graph[0]->Draw("AP");
    graph[1]->Draw("P");  
}

/* TEMPLATE */

// void f() {

//     vector<double> *X=nullptr,
//         *Y=nullptr;

//     // TGraph* graph = new TGraph();

//     for (int i_file=0; i_file < filelist.size(); i_file++) {
//         string filename = filelist[i_file];

//         cout << "\e[3mOpening file #" << i_file+1 << ": " << filename << "\e[0m" << endl;

//         TFile file(filename.c_str());
//         TTree *Reco = (TTree*) file.Get("LauraPDumper/Reco");
//         TTree *Truth = (TTree*) file.Get("LauraPDumper/Truth");
//         int n_ev=Reco->GetEntries();
        
//         Reco->SetBranchAddress("", &X);
//         Truth->SetBranchAddress("", &Y);

//         for (int i_ev=0; i_ev < n_ev; i_ev++) {

//             Reco->GetEntry(i_ev);
//             Truth->GetEntry(i_ev);

//             // int n_part=nParticles;

//             for(int i_part=0; i_part < n_part; i_part++) {

//                 // int n_clu = CluPlane->at(i_part).size();

//                 for (int i_clu=0; i_clu < n_clu; i_clu++) {

//                 } //end of cluster loop
//             } //end of particule loop

//             // int n_tr=TrackID->size();

//             for(int i_tr=0; i_tr < n_tr; i_tr++) {

//             } //end of track loop

//             // int n_sh=ShowerID->size();

//             for(int i_sh=0; i_sh < n_sh; i_sh++) {

//             } //end of shower loop
//         } //end of event loop
//         file.Close();
//     } //end of file loop

//     // TCanvas* canvas = new TCanvas("c","");
//     // canvas->cd();
//     // graph->Draw("AP");
// }
