#include "graph.h"

Binning binZ = {207, -310, 311};


void dz() {

    std::string name;
    std::cout << "name: ";
    std::getline(std::cin, name);
    if (!name.empty()) name = "_" + name;

    TString in = "data/pdvd_Muchecks_out" + name + ".root";
    TString out = "out/dz" + name + ".root";

    std::cout << "opening " << in << "..." << std::endl;

    TFile* fin = new TFile(in);
    TFile* fout = new TFile(out, "RECREATE");

    TTree* event = fin->Get<TTree>("checks/Event");
    TTree* muon = fin->Get<TTree>("checks/Muon");

    unsigned int n_event = event->GetEntries();

    size_t EventNMuon;
    std::vector<size_t> *EventiMuon = nullptr;
    event->SetBranchAddress("NMuon", &EventNMuon);
    event->SetBranchAddress("iMuon", &EventiMuon);

    float MuonEndHitZ;
    float MuonEndPointZ;
    float MuonEndSpacePointZ;
    muon->SetBranchAddress("EndHitZ", &MuonEndHitZ);
    muon->SetBranchAddress("EndPointZ", &MuonEndPointZ);
    muon->SetBranchAddress("EndSpacePointZ", &MuonEndSpacePointZ);

    int MuonDoesDecay;
    int MuonHasMichel;
    muon->SetBranchAddress("DoesDecay", &MuonDoesDecay);
    muon->SetBranchAddress("HasMichel", &MuonHasMichel);

    int MuonEndIsInWindow = 0;
    int MuonEndIsInVolumeYZ = 0;
    muon->SetBranchAddress("EndIsInWindow", &MuonEndIsInWindow);
    muon->SetBranchAddress("EndIsInVolumeYZ", &MuonEndIsInVolumeYZ);


    TH1F* hdz = new TH1F("hdz", "DZ;DZ (cm);count", binZ.n, binZ.min, binZ.max);
    TH1F* hdz2 = new TH1F("hdz2", "", binZ.n, binZ.min, binZ.max);

    for (unsigned int i_event=0; i_event<n_event; i_event++) {

        // std::cout << "evt " << i_event << "\r" << std::flush;
        
        event->GetEntry(i_event);

        for (size_t mu=0; mu<EventNMuon; mu++) {

            // std::cout << "\t" << "mu " << mu << "\r" << std::flush;

            muon->GetEntry(EventiMuon->at(mu));

            // if (!MuonDoesDecay) continue;
            // if (MuonHasMichel != 2) continue;
            if (!MuonEndIsInWindow || !MuonEndIsInVolumeYZ) continue;


            float dz = MuonEndHitZ - MuonEndPointZ;
            float dz2 = MuonEndHitZ - MuonEndSpacePointZ;

            if (abs(dz2) > 50) std::cout << MuonEndSpacePointZ << std::endl;

            // std::cout << "\t\t" << "dz: " << dz << " dz2: " << dz2 << std::endl;

            hdz->Fill(dz);

            if (MuonEndSpacePointZ == 0) continue;

            hdz2->Fill(dz2);
        }
    }

    TCanvas* cdz = new TCanvas("cdz", "cdz no associated space point");

    cdz->cd();
    gPad->SetLogy();

    hdz->Draw();

    hdz2->SetLineColor(kOrange+7);
    hdz2->Draw("same");

    cdz->Write();

    std::cout << "writing " << out << "..." << std::endl;

    fin->Close();
    fout->Close();
}