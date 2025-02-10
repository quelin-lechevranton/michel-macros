#include "rvec_utilities.h"

#include <TCanvas.h>

Binning binR = {100U, 0, 150};
Binning binE = {30U, 0, 100};
Binning binADC = {100U, 0, 3000};
Binning binN = {100U, 0, 4};


float bragg_length = 20;
float body_length = 20;

ROOT::RDataFrame event("checks/Event","data/pdvd_Muchecks_out_100_r20_in.root");
ROOT::RDataFrame muon("checks/Muon","data/pdvd_Muchecks_out_100_r20_in.root");

size_t prev_e=0, prev_m=0;

auto df = muon.Define("iMuonInEvent", [](unsigned long e, unsigned long m) {
    if (e) {
        if (e != prev_e) {
            prev_e = e;
            prev_m = m;
        }
        return m - prev_m;
    }
    return m;
}, {"iEvent", "iMuon"})

// auto df = muon

.Define("HitCorrTick", Form(
    "(RVec<float>) Where(HitSlice >= 4, 2 * Max(HitTick[HitSlice >= 4]) + %fF - HitTick, HitTick)", 
    cathode_height / drift_velocity / sampling_rate
))

.Define("HitDist", Form(
    "(RVec<float>) sqrt(pow((HitCorrTick - EndHitTick) * %fF, 2) + pow(HitZ - EndHitZ, 2))",
    drift_velocity * sampling_rate
))

.Define("HitIdxSort", "Argsort(HitDist)")
.Define("HitDistSort", "Take(HitDist, HitIdxSort)")

.Define("HitSliceSort", "Take(HitSlice, HitIdxSort)")
.Define("HitChannelSort", "Take(HitChannel, HitIdxSort)")
.Define("HitZSort", "Take(HitZ, HitIdxSort)")
.Define("HitTickSort", "Take(HitTick, HitIdxSort)")
.Define("HitADCSort", "Take(HitADC, HitIdxSort)")

.Define("BodyHitSortIdx", Form(
    "Range(NHit)[%fF < HitDistSort && HitDistSort < %fF]", 
    bragg_length, body_length + bragg_length
))
.Define("BodyNHit", "(unsigned) BodyHitSortIdx.size()")
.Define("BodyHitSlice", "Take(HitSliceSort, BodyHitSortIdx)")
.Define("BodyHitChannel", "Take(HitChannelSort, BodyHitSortIdx)")
.Define("BodyHitZ", "Take(HitZSort, BodyHitSortIdx)")
.Define("BodyHitTick", "Take(HitTickSort, BodyHitSortIdx)")
.Define("BodyHitADC", "Take(HitADCSort, BodyHitSortIdx)")

.Define("BodyMinZ", "Min(BodyHitZ)")
.Define("BodyMaxZ", "Max(BodyHitZ)")

.Define("BodyCovTZ", "Cov(BodyHitZ, BodyHitTick)")
.Define("BodyVarZ", "Var(BodyHitZ)")
.Define("BodyVarT", "Var(BodyHitTick)")

.Define("BodyR2", "pow(BodyCovTZ, 2) / (BodyVarZ * BodyVarT)")
.Define("BodySlope", "BodyCovTZ / BodyVarZ")
.Define("BodyOrigin", "Mean(BodyHitTick) - BodySlope * Mean(BodyHitZ)")
// .Define("BodyMaxErr", "Max(abs(HitCorrTick - (BodySlope * HitZ + BodyOrigin)))")

.Define("BraggHitSortIdx", Form(
    "Range(NHit)[HitDistSort < %fF]", 
    bragg_length
))
.Define("BraggNHit", "(unsigned) BraggHitSortIdx.size()")
.Define("BraggHitSlice", "Take(HitSliceSort, BraggHitSortIdx)")
.Define("BraggHitChannel", "Take(HitChannelSort, BraggHitSortIdx)")
.Define("BraggHitZ", "Take(HitZSort, BraggHitSortIdx)")
.Define("BraggHitTick", "Take(HitTickSort, BraggHitSortIdx)")
.Define("BraggHitADC", "Take(HitADCSort, BraggHitSortIdx)")

.Define("HitdxSort", "sqrt(pow(DiffNext(HitZSort), 2) + pow(DiffNext(HitTickSort), 2))")
// .Define("HitResRangeSort", "Acc(sqrt(pow(DiffPrev(HitZSort), 2) + pow(DiffPrev(HitTickSort), 2)))")
.Define("HitResRangeSort", "Acc(HitdxSort)")

.Define("HitdEdxSort", "HitADCSort / HitdxSort")
.Define("ResRangeAtMaxADC", "Take(HitResRangeSort, ArgMax(HitADCSort))")

.Define("BodydEdx", "Mean(Take(HitdEdxSort, BodyHitSortIdx))")
.Define("BraggdEdx", "Mean(Take(HitdEdxSort, BraggHitSortIdx))")

.Define("Bragg", "BraggdEdx/BodydEdx")




.Define("LocalBragg", [](unsigned n, RVec<float> z, RVec<float> t, RVec<unsigned long> i) {
    std::vector<float> local_bragg;
    for (unsigned i=1; i<n; i++) {
        if (z[i] < bragg_length || z[i] > body_length + bragg_length) continue;
        local_bragg.push_back(abs(t[i] - t[i-1]));
    }
    return local_bragg;
}, {"NHit", "HitZSort", "HitTickSort", "HitIdxSort"})


void bragg() {
    TString out = "data/bragg.root";
    
    std::cout << "writing " << out << "..." << std::endl;

    ROOT::RDF::RSnapshotOptions opt;
    opt.fMode = "UPDATE";
    event.Snapshot("event", out);
    df.Snapshot("muon", out, ".*", opt); // Regex mathch all columns


    auto h2dEdxRR_d = df.Filter("DoesDecay").Histo2D({
        "h2_dEdx_RR_d",
        "dEdx vs RR (decaying);RR [cm];dEdx [ADC.tick/cm]",
        binR.n, binR.min, binR.max,
        binADC.n, binADC.min, binADC.max
    }, "HitResRangeSort", "HitdEdxSort");

    auto h2dEdxRR_n = df.Filter("!DoesDecay").Histo2D({
        "h2_dEdx_RR_n",
        "dEdx vs RR (non-decaying);RR [cm];dEdx [ADC.tick/cm]",
        binR.n, binR.min, binR.max,
        binADC.n, binADC.min, binADC.max
    }, "HitResRangeSort", "HitdEdxSort");

    auto hRRatMax_d = df.Filter("DoesDecay").Histo1D({
        "hRRatMax_d",
        ";RR [cm];count",
        binR.n, binR.min, binR.max
    }, "ResRangeAtMaxADC");

    auto hRRatMax_n = df.Filter("!DoesDecay").Histo1D({
        "hRRatMax_n",
        ";RR [cm];count",
        binR.n, binR.min, binR.max
    }, "ResRangeAtMaxADC");

    auto hBraggdEdx_d = df.Filter("DoesDecay").Histo1D({
        "hBraggdEdx_d",
        ";;count",
        binN.n, binN.min, binN.max
    }, "Bragg");

    auto hBraggdEdx_n = df.Filter("!DoesDecay").Histo1D({
        "hBraggdEdx_n",
        ";;count",
        binN.n, binN.min, binN.max
    }, "Bragg");

    out = "out/bragg_v3_hist.root";
    TFile* f = new TFile(out, "RECREATE");

    TCanvas* c = new TCanvas("c","c");

    c->Divide(2,2);
    c->cd(1);
    h2dEdxRR_d->Draw("colz");
    c->cd(2);
    h2dEdxRR_n->Draw("colz");
    c->cd(3);
    hRRatMax_n->Draw();
    hRRatMax_d->SetLineColor(kRed);
    hRRatMax_d->Draw("same");
    gPad->SetTitle("RR at max ADC");
    c->cd(4);
    hBraggdEdx_n->Draw();
    hBraggdEdx_d->SetLineColor(kRed);
    hBraggdEdx_d->Draw("same");
    gPad->SetTitle("Bragg Mean dEdx / MIP");

    std::cout << "writing " << out << "..." << std::endl;
    c->Write();
    f->Close();

    df.Range(500).Foreach([](
        unsigned long e, 
        unsigned long m, 
        double b, 
        int d, 
        unsigned sl, 
        float z, 
        float t, 
        double mipdedx, 
        double braggdedx, 
        unsigned nmip, 
        unsigned nbragg
    ) {
        if (d) return;
        if (b < 2) return;
        std::cout << "e" << e << "µ" << m << ": bragg of " << b << " @ sl" << sl << " z" << z << " t" << t << " mip" << mipdedx << "(" << nmip << ") bragg" << braggdedx << "(" << nbragg << ")" << std::endl;

        // if (e==8 and m==5) {
        //     std::cout << "e" << e << "µ" << m << ":" << std::endl;
        //     for (unsigned i=0; i<braggz.size(); i++) {
        //         std::cout << "  " << i << ": z" << braggz[i] << " t" << braggt[i] << " a" << bragga[i] << " dedx" << braggdd[i] << std::endl;
        //     }
        // }
    }, {"iEvent", "iMuonInEvent", "Bragg", "DoesDecay", "EndHitSlice", "EndHitZ", "EndHitTick", "BodydEdx", "BraggdEdx", "BodyNHit", "BraggNHit"});



    // auto m = df.Take<double>("BodySlope").GetValue();
    // auto q = df.Take<double>("BodyOrigin").GetValue();
    // auto corr = df.Take<double>("BodyR2").GetValue();
    // auto events = df.Take<unsigned long>("iEvent").GetValue();
    // auto muons = df.Take<unsigned long>("iMuon").GetValue();
    // auto end_slice = df.Take<unsigned>("EndHitSlice").GetValue();
    // auto end_tick = df.Take<float>("EndHitTick").GetValue();
    // auto end_z = df.Take<float>("EndHitZ").GetValue();

    // std::map<unsigned long, unsigned long> event_muon;
    // for (unsigned long i=0; i<events.size(); i++) {
    //     if (event_muon[events[i]] == 0) {
    //         event_muon[events[i]] = muons[i];
    //     }
    // }
    // event_muon[0] = 0;

    // for (unsigned int i=0; i<corr.size(); i++) {
    //     if (corr[i] > 0.8) continue;
    //     std::cout << "e" << events[i] << " mu" << muons[i] - event_muon[events[i]] << " @ sl" << end_slice[i] << " t" << end_tick[i] << " z" << end_z[i] << " lr: " << m[i] << " * x + " << q[i] << std::endl;
    // }
}
