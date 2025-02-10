#include "pdvd_event_display.h"

void c2d4_bragg(void) {

    // Loading from command line inputs

    TString name = CLInput<TString>("name", "data/bragg");

    TString in = name.Contains(".root") ? name : name += ".root";
    TString out = "out/c2d4_bragg_" + removePath(name);

    std::cout << "opening " << in << "..." << std::endl;

    TFile* fin = new TFile(in);
    TFile* fout = new TFile(out, "RECREATE");

    // TTree* event = fin->Get<TTree>("checks/Event");
    // TTree* muon = fin->Get<TTree>("checks/Muon");
    TTree* event = fin->Get<TTree>("event");
    TTree* muon = fin->Get<TTree>("muon");

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

    // Hits MuonHits;
    // muon->SetBranchAddress("NHit", &MuonHits.N);
    // muon->SetBranchAddress("HitSlice", &MuonHits.slice);
    // muon->SetBranchAddress("HitZ", &MuonHits.Z);
    // muon->SetBranchAddress("HitChannel", &MuonHits.channel);
    // muon->SetBranchAddress("HitTick", &MuonHits.tick);
    // muon->SetBranchAddress("HitADC", &MuonHits.adc);

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


    // BRAGG

    Hits MuonHits;
    muon->SetBranchAddress("NHit", &MuonHits.N);
    muon->SetBranchAddress("HitSliceSort", &MuonHits.slice);
    muon->SetBranchAddress("HitZSort", &MuonHits.Z);
    muon->SetBranchAddress("HitTickSort", &MuonHits.tick);
    muon->SetBranchAddress("HitChannelSort", &MuonHits.channel);


    double BodySlope, BodyOrigin;
    float BodyMinZ, BodyMaxZ;
    muon->SetBranchAddress("BodySlope", &BodySlope);
    muon->SetBranchAddress("BodyOrigin", &BodyOrigin);
    muon->SetBranchAddress("BodyMinZ", &BodyMinZ);
    muon->SetBranchAddress("BodyMaxZ", &BodyMaxZ);

    BasicHits BodyHits;
    muon->SetBranchAddress("BodyNHit", &BodyHits.N);
    muon->SetBranchAddress("BodyHitTick", &BodyHits.tick);
    muon->SetBranchAddress("BodyHitChannel", &BodyHits.channel);

    BasicHits BraggHits;
    muon->SetBranchAddress("BraggNHit", &BraggHits.N);
    muon->SetBranchAddress("BraggHitTick", &BraggHits.tick);
    muon->SetBranchAddress("BraggHitChannel", &BraggHits.channel);


    for (unsigned i_event=first_event; i_event<n_event; i_event++) {

        std::cout << "e" << i_event << "\r" << std::flush;

        event->GetEntry(i_event);

        TCanvas *c = new TCanvas(
            Form("e%u", i_event), 
            Form("e%u w/ %u muons", i_event, (unsigned) EventNMuon),
            20, 20, 1200, 700
        );

        std::vector<TPad*> ps = draw4Frame(c);

        draw4TH2F(ps, Form("h2_e%u", i_event), EventHits);

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
            draw4TGraph(ps, gmu, "same L", MuonHits);

            TEllipse* el = new TEllipse();
            if (!MuonEndIsInWindow or !MuonEndIsInVolumeYZ) 
                el->SetLineStyle(kDashed);
            el->SetFillStyle(0);
            el->SetLineColor(kRed-7);
            el->SetLineWidth(2);
            draw4TEllipse(ps, el, MuonEnd, space_radius);

            TGraph* gm = new TGraph();
            gm->SetName(Form("gm_e%u_mu%ld", i_event, mu));
            gm->SetMarkerStyle(MuonHasMichel == 2 ? kFullCircle : kOpenCircle);
            gm->SetMarkerSize(1);
            gm->SetMarkerColor(kAzure-4);
            draw4TGraph(ps, gm, "same P", MichelHits);

            TGraph* gs = new TGraph();
            gs->SetName(Form("gs_e%u_mu%ld", i_event, mu));
            // gs->SetMarkerStyle(71); // OpenCircle
            gs->SetMarkerStyle(70); // Bold Cross
            gs->SetMarkerSize(2);
            gs->SetMarkerColor(kRed-7);
            draw4TGraph(ps, gs, "same P", SphereHits);

            // BRAGG

            TF1* lr = new TF1(
                Form("lr_e%u_mu%ld", i_event, mu),
                "[0]*x+[1]",
                BodyMinZ, BodyMaxZ
            );
            lr->SetParameters(BodySlope, BodyOrigin);
            ps[MuonEnd.slice%4]->cd();
            lr->SetLineColor(kGreen-3 + cmu);
            lr->SetLineWidth(2);
            lr->Draw("same");

            TGraph* gb = new TGraph();
            gb->SetName(Form("gb_e%u_mu%ld", i_event, mu));
            gb->SetMarkerStyle(68);
            gb->SetMarkerSize(1);
            gb->SetMarkerColor(kGreen-3 + cmu);
            draw4TGraph(ps, gb, "same P", BodyHits);

            TGraph* gbr = new TGraph();
            gbr->SetName(Form("gbr_e%u_mu%ld", i_event, mu));
            gbr->SetMarkerStyle(68);
            gbr->SetMarkerSize(1);
            gbr->SetMarkerColor(kViolet-3 + cmu);
            draw4TGraph(ps, gbr, "same P", BraggHits);
        }

        c->Write();
    }

    std::cout << "writing " << out << "..." << std::endl;

    fin->Close();
    fout->Close();
}