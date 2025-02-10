#include "pdvd_event_display.h"

void c2d(void) {

    // Loading from command line inputs

    TString name = CLInput<TString>("name", "data/pdvd_Muchecks_out");

    TString in = name.Contains(".root") ? name : name += ".root";
    TString out = "out/c2d_" + removePath(name);

    std::cout << "opening " << in << "..." << std::endl;

    TFile* fin = new TFile(in);
    TFile* fout = new TFile(out, "RECREATE");

    TTree* event = fin->Get<TTree>("checks/Event");
    TTree* muon = fin->Get<TTree>("checks/Muon");

    float space_radius = CLInput<float>("space_radius", 20.); 

    int ev_input = CLInput<int>("n_event", (int) event->GetEntries());
    unsigned first_event = 0;
    unsigned n_event = (unsigned) ev_input;
    if (ev_input < 0) {
        first_event = (unsigned) -ev_input;
        n_event = first_event + 1;
    }

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


    for (unsigned i_event=first_event; i_event<n_event; i_event++) {

        std::cout << "e" << i_event << "\r" << std::flush;

        event->GetEntry(i_event);

        TCanvas *c = new TCanvas(
            Form("e%u", i_event), 
            Form("e%u w/ %u muons", i_event, (unsigned) EventNMuon),
            20, 20, 1200, 700
        );

        std::vector<TPad*> ps = drawFrame(c);

        drawTH2F(ps, Form("h2_e%u", i_event), EventHits);

        unsigned cmu=0;
        for (size_t mu=0; mu<EventNMuon; mu++) {

            muon->GetEntry(EventiMuon->at(mu));

            // if (!MuonDoesDecay) continue;

            std::cout << "\t" << "mu" << std::setw(3) << mu
                      << " end@ sl" << MuonEnd.slice
                      << " : z" << std::setw(4) << (int) MuonEnd.Z
                      << " : t" << std::setw(4) << (int) MuonEnd.tick
                      << "  Nhit: " << std::setw(3) << MichelHits.N
                      << "  TE: " << std::setw(4) << std::setprecision(3) << MichelTrueEnergy
                      << "  HE: " << std::setw(4) << std::setprecision(3) << MichelHitEnergy
                      << "  SE: " << std::setw(4) << std::setprecision(3) << MichelSphereEnergy
                      << std::endl;


            TGraph* gmu = new TGraph();
            gmu->SetName(Form("gmu_e%u_mu%ld", i_event, mu));
            gmu->SetLineColor(kOrange-3 + cmu++);
            gmu->SetLineWidth(2);
            gmu->SetMarkerStyle(20);
            drawTGraph(ps, gmu, "same L", MuonHits);

            TEllipse* el = new TEllipse();
            if (!MuonEndIsInWindow or !MuonEndIsInVolumeYZ) 
                el->SetLineStyle(kDashed);
            el->SetFillStyle(0);
            el->SetLineColor(kRed-7);
            el->SetLineWidth(2);
            drawTEllipse(ps, el, MuonEnd, space_radius);

            TGraph* gm = new TGraph();
            gm->SetName(Form("gm_e%u_mu%ld", i_event, mu));
            gm->SetMarkerStyle(MuonHasMichel == 2 ? kFullCircle : kOpenCircle);
            gm->SetMarkerSize(1);
            gm->SetMarkerColor(kAzure-4);
            drawTGraph(ps, gm, "same P", MichelHits);

            TGraph* gs = new TGraph();
            gs->SetName(Form("gs_e%u_mu%ld", i_event, mu));
            // gs->SetMarkerStyle(71); // OpenCircle
            gs->SetMarkerStyle(70); // Cross
            gs->SetMarkerSize(2);
            gs->SetMarkerColor(kRed-7);
            drawTGraph(ps, gs, "same P", SphereHits);
        }

        c->Write();
    }

    std::cout << "writing " << out << "..." << std::endl;

    fin->Close();
    fout->Close();
}