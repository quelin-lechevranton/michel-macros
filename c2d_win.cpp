#include "pdvd_event_display.h"

void c2d_win(void) {

    // Loading from command line inputs

    TString name = CLInput<TString>("name", "data/pdvd_Muchecks_out");

    TString in = name + ".root";

    std::cout << "opening " << in << "..." << std::endl;

    TFile* fin = new TFile(in);

    TTree* event = fin->Get<TTree>("checks/Event");
    TTree* muon = fin->Get<TTree>("checks/Muon");

    float space_radius = CLInput<float>("space_radius", 20.);
    float tick_radius = space_radius / drift_velocity / sampling_rate; // ticks
    float y_radius = cm ? space_radius : tick_radius;
    std::cout << "tick_radius: " << tick_radius << std::endl;

    unsigned i_event = CLInput<unsigned>("i_event", 0);
    int i_mu = CLInput<int>("i_mu", -1);

    
    // Loading file contents

    // Event Tree
    Hits EventHits;
    event->SetBranchAddress("NHit", &EventHits.N);
    event->SetBranchAddress("HitSlice", &EventHits.slice);
    event->SetBranchAddress("HitZ", &EventHits.Z);
    event->SetBranchAddress("HitChannel", &EventHits.channel);
    event->SetBranchAddress("HitTick", &EventHits.tick);
    event->SetBranchAddress("HitADC", &EventHits.adc);

    size_t EventNMuon;
    std::vector<size_t> *EventiMuon = nullptr;
    event->SetBranchAddress("NMuon", &EventNMuon);
    event->SetBranchAddress("iMuon", &EventiMuon);

    // Muon Tree
    int MuonDoesDecay;
    int MuonHasMichel;
    muon->SetBranchAddress("DoesDecay", &MuonDoesDecay);
    muon->SetBranchAddress("HasMichel", &MuonHasMichel);

    Hits MuonHits;
    muon->SetBranchAddress("NHit", &MuonHits.N);
    muon->SetBranchAddress("HitSlice", &MuonHits.slice);
    muon->SetBranchAddress("HitZ", &MuonHits.Z);
    muon->SetBranchAddress("HitChannel", &MuonHits.channel);
    muon->SetBranchAddress("HitTick", &MuonHits.tick);
    muon->SetBranchAddress("HitADC", &MuonHits.adc);

    Hit MuonEnd;
    muon->SetBranchAddress("EndHitSlice", &MuonEnd.slice);
    muon->SetBranchAddress("EndHitZ", &MuonEnd.Z);
    muon->SetBranchAddress("EndHitChannel", &MuonEnd.channel);
    muon->SetBranchAddress("EndHitTick", &MuonEnd.tick);

    int MuonEndIsInWindow = 0;
    int MuonEndIsInVolumeZ = 0;
    int MuonEndIsInVolumeYZ = 0;
    muon->SetBranchAddress("EndIsInWindow", &MuonEndIsInWindow);
    muon->SetBranchAddress("EndIsInVolumeZ", &MuonEndIsInVolumeZ);
    muon->SetBranchAddress("EndIsInVolumeYZ", &MuonEndIsInVolumeYZ);

    float MichelTrueEnergy = 0;
    float MichelHitEnergy = 0;
    float MichelSphereEnergy = 0;
    muon->SetBranchAddress("MichelTrueEnergy", &MichelTrueEnergy);
    muon->SetBranchAddress("MichelHitEnergy", &MichelHitEnergy);
    muon->SetBranchAddress("SphereEnergy", &MichelSphereEnergy);

    Hits MichelHits;
    muon->SetBranchAddress("MichelNHit", &MichelHits.N);
    muon->SetBranchAddress("MichelHitSlice", &MichelHits.slice);
    muon->SetBranchAddress("MichelHitZ", &MichelHits.Z);
    muon->SetBranchAddress("MichelHitChannel", &MichelHits.channel);
    muon->SetBranchAddress("MichelHitTick", &MichelHits.tick);
    muon->SetBranchAddress("MichelHitADC", &MichelHits.adc);

    Hits SphereHits;
    muon->SetBranchAddress("SphereNHit", &SphereHits.N);
    muon->SetBranchAddress("SphereHitSlice", &SphereHits.slice);
    muon->SetBranchAddress("SphereHitZ", &SphereHits.Z);
    muon->SetBranchAddress("SphereHitChannel", &SphereHits.channel);
    muon->SetBranchAddress("SphereHitTick", &SphereHits.tick);
    muon->SetBranchAddress("SphereHitADC", &SphereHits.adc);


    std::cout << "e" << i_event << "\r" << std::flush;

    event->GetEntry(i_event);

    unsigned sl;
    Binning bt, bz;

    if (i_mu == -1) {
        sl = CLInput<unsigned>("slice", 0);
        bz = CLInput<Binning>("Z", binZ);
        if (cm)
            bt = CLInput<Binning>("tick_cm", binX);
        else
            bt = CLInput<Binning>("tick", binT);
    } else {
        if (i_mu < 0 or i_mu >= EventNMuon) {
            std::cout << "no such muon" << std::endl;
            fin->Close();
            return;
        }
        muon->GetEntry(EventiMuon->at(i_mu));
        if (!MuonDoesDecay) {
            std::cout << "no such decaying muon" << std::endl;
            fin->Close();
            return;
        }

        sl = MuonEnd.slice;
        bz = Binning(MuonEnd.Z, 2*space_radius, binZ);
        bt = Binning(MuonEnd.display_tick(), 2*y_radius, binY);
    }

    TCanvas *c = new TCanvas("c", Form("e_%u", i_event));

    TPad* p = drawFrame(c, sl, bz, bt);

    drawTH2F(p, Form("h2_e%u", i_event), EventHits, sl, bz, bt);

    for (unsigned mu=0; mu<EventNMuon; mu++) {

        muon->GetEntry(EventiMuon->at(mu));
        if (!MuonDoesDecay) continue;

        TGraph* gmu = new TGraph();
        gmu->SetName(Form("gmu_e%u_mu%u", i_event, mu));
        gmu->SetLineColor(mu < 20 ? kPink-9 + mu : kPink-9 + mu - 20);
        gmu->SetLineWidth(2);
        gmu->SetMarkerStyle(20);
        drawTGraph(p, gmu, "", MuonHits, sl, bz, bt);

        if (MuonEnd.slice == sl) {
            TEllipse* el = new TEllipse(
                MuonEnd.Z, 
                MuonEnd.display_tick(),
                space_radius, 
                cm ? space_radius : tick_radius
            );
            if (!MuonEndIsInWindow or !MuonEndIsInVolumeYZ) 
                el->SetLineStyle(kDashed);
            el->SetFillStyle(0);
            el->SetLineColor(kViolet+5);
            el->SetLineWidth(2);
            p->cd();
            el->Draw();
        }

        TGraph* gm = new TGraph();
        gm->SetName(Form("gm_e%u_mu%u", i_event, mu));
        gm->SetMarkerStyle(MuonHasMichel == 2 ? kFullDoubleDiamond : kOpenDoubleDiamond);
        gm->SetMarkerSize(1.5);
        gm->SetMarkerColor(kAzure-5);
        drawTGraph(p, gm, "P", MichelHits, sl, bz, bt);

        TGraph* gs = new TGraph();
        gs->SetName(Form("gs_e%u_mu%u", i_event, mu));
        gs->SetMarkerStyle(kFullFourTrianglesX);
        gs->SetMarkerSize(1.5);
        gs->SetMarkerColor(kViolet+5);
        drawTGraph(p, gs, "P", SphereHits, sl, bz, bt);

    }

    TString out = "out/" + (i_mu == -1 ? "c2dw_" : "c2dmu_") + name + ".root";
    std::cout << "writing " << out << "..." << std::endl;
    TFile* fout = new TFile(out, "RECREATE");
    c->Write();
    fout->Close();
}