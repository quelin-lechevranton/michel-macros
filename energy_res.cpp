#include "graph.h"

Binning binEnergy = {30U, 0, 100};


void energy_res(void) {

    TString out = "out/energy_res.root";
    
    // std::cout << "opening " << in << "..." << std::endl;

    TFile* f = new TFile(out, "RECREATE");

    ROOT::RDataFrame muon("checks/Muon", "data/pdvd_Muchecks_out_100_r20_in.root");
    auto muon_in = muon.Filter("EndIsInWindow && EndIsInVolumeYZ");
    auto muon_decay = muon_in.Filter("DoesDecay");

    auto hTrueE = muon_decay.Histo1D({
        "hTrueE", "Michel True Energy;Energy (MeV);count", 
        binEnergy.n, binEnergy.min, binEnergy.max
        }, "MichelTrueEnergy");

    auto hHitE = muon_decay.Histo1D({
        "hHitE", "Michel Hit Energy;ADC (~MeV);Count", 
        binEnergy.n, binEnergy.min, binEnergy.max
        }, "MichelHitEnergy");
    
    auto c = new TCanvas();

    c->Divide(2,1);
    c->cd(1);
    hTrueE->Draw();

    // auto hCumulTrueE = hTrueE->GetCumulative();
    // hCumulTrueE->Scale(1./hCumulTrueE->GetMaximum());
    // hCumulTrueE->Draw();

    c->cd(2);
    TH1F* hConvolved = new TH1F("hConvolved", "Convolved with 20% resolution;Convolved Energy (MeV);count", binEnergy.n, binEnergy.min, binEnergy.max);
    for (int i =0; i<hTrueE->GetEntries(); i++) {
        auto x = hTrueE->GetRandom();
        x = gRandom->Gaus(x, 20./100.* x);
        hConvolved->Fill(x);
    }
    hConvolved->Draw();

    c->Write();

    std::cout << "writing " << out << "..." << std::endl;

    f->Close();

}
