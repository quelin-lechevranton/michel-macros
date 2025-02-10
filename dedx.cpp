#include "graph.h"

Binning binTick = {200, 0, 100};
Binning binEnergy = {30, 0, 100};
Binning binADC = {100, 0, 3000};
Binning binDist = {100, 0, 600};
Binning binNorm = {100, 0, 4};

unsigned int mean_r = 8;


void dedx() {

    std::string name;
    std::cout << "name: ";
    std::getline(std::cin, name);
    if (!name.empty()) name = "_" + name;

    TString in = "data/pdvd_Muchecks_out" + name + ".root";
    TString out = "out/dedx" + name + ".root";
    
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

    Hit MuonEnd;
    muon->SetBranchAddress("EndHitSlice", &MuonEnd.slice);
    muon->SetBranchAddress("EndHitZ", &MuonEnd.Z);
    muon->SetBranchAddress("EndHitChannel", &MuonEnd.channel);
    muon->SetBranchAddress("EndHitTick", &MuonEnd.tick);

    Hits MuonHits;
    muon->SetBranchAddress("NHit", &MuonHits.N);
    muon->SetBranchAddress("HitSlice", &MuonHits.slice);
    muon->SetBranchAddress("HitZ", &MuonHits.Z);
    muon->SetBranchAddress("HitChannel", &MuonHits.channel);
    muon->SetBranchAddress("HitTick", &MuonHits.tick);
    muon->SetBranchAddress("HitADC", &MuonHits.adc);

    int MuonDoesDecay;
    int MuonHasMichel;
    muon->SetBranchAddress("DoesDecay", &MuonDoesDecay);
    muon->SetBranchAddress("HasMichel", &MuonHasMichel);

    int MuonEndIsInWindow = 0;
    int MuonEndIsInVolumeYZ = 0;
    muon->SetBranchAddress("EndIsInWindow", &MuonEndIsInWindow);
    muon->SetBranchAddress("EndIsInVolumeYZ", &MuonEndIsInVolumeYZ);

    TH2F* hdEdx = new TH2F(
        "hdEdx", "dE/dx;Hit distance to end (cm);Hit deposit (ADC.tick)",
        binDist.n, binDist.min, binDist.max,
        binADC.n, binADC.min, binADC.max
    );


    TH2F* hmean = new TH2F(
        "hmean", Form("truncated mean dE/dx over %u hits;Hit distance to end (cm);Hit normalized truncated mean deposit", 2*mean_r+1),
        binDist.n, binDist.min, binDist.max,
        binNorm.n, binNorm.min, binNorm.max
    );

    for (unsigned int i_event=0; i_event<n_event; i_event++) {

        event->GetEntry(i_event);

        for (size_t µ=0; µ<EventNMuon; µ++) {

            muon->GetEntry(EventiMuon->at(µ));

            // if (!MuonDoesDecay) continue;
            if (MuonHasMichel != 0) continue;
            if (!MuonEndIsInWindow || !MuonEndIsInVolumeYZ) continue;

            float mip_adc = 0;
            unsigned int mip_n = 0;

            std::vector<std::pair<float, float>> pt;

            for (unsigned int i=0; i<MuonHits.N; i++) {
                Hit hit = MuonHits.at(i);

                float delta_tick = drift_velocity * abs(hit.tick - MuonEnd.tick) * sampling_rate;
                float dist = pow(delta_tick, 2) + pow(hit.Z - MuonEnd.Z, 2);

                hdEdx->Fill(dist, hit.adc);

                float mean_adc = 0;
                std::vector<float> adcs;
                unsigned int mean_n = 0;
                for (unsigned int j = i-mean_r>0 ? i-mean_r : 0; j<=i+mean_r && j<MuonHits.N; j++) {
                    Hit hit2 = MuonHits.at(j);
                    adcs.push_back(hit2.adc);
                    // mean_n++;
                }
                if (adcs.size() == 0) continue;

                std::sort(adcs.begin(), adcs.end()); 

                // std::cout << "adcs: " << "\r" << std::flush << "\t\t\t";
                // for (float a : adcs) std::cout << a << " ";
                // std::cout << std::endl;

                std::vector<float> truncated_adcs(adcs.begin()+adcs.size()/4, adcs.end()-adcs.size()/4);

                // std::cout << "truncated adcs: " << "\r" << std::flush << "\t\t\t";
                // for (float a : truncated_adcs) std::cout << std::setprecision(3) << a << " ";
                // std::cout << std::endl;

                mean_adc = std::accumulate(truncated_adcs.begin(), truncated_adcs.end(), 0.0);
                mean_adc /= truncated_adcs.size();
                // std::cout << "mean adc: " << mean_adc << std::endl << std::endl;

                if (dist > 10) {
                    mip_adc += mean_adc;
                    mip_n++;
                }

                pt.push_back({dist, mean_adc});
            }

            mip_adc /= mip_n;

            for (auto [dist, mean_adc] : pt) {
            //     std::cout << "new point at " << dist << " cm: " << mean_adc << " / " << mip_adc << std::endl;
                hmean->Fill(dist, mean_adc/mip_adc);
            }
        }
    }


    TGraph* gdEdx = new TGraph();
    gdEdx->SetTitle("dE/dx;Hit distance to end (cm);Hit mean of projected deposit (ADC.tick)");
    for (int binx=1; binx<=binDist.n; binx++) {
        float x = (float) binx * (binDist.max - binDist.min) / binDist.n;
        TH1D* h = hdEdx->ProjectionY(Form("hdEdx_slice_x%.f", x), binx, binx);
        gdEdx->AddPoint(x, h->GetMean());
    }



    TGraph* gmean = new TGraph();
    gmean->SetTitle("mean of projected normalized truncated mean dE/dx;Hit distance to end (cm);Hit mean projected normalized truncated mean deposit");
    for (int binx=1; binx<=binDist.n; binx++) {
        float x = (float) binx * (binDist.max - binDist.min) / binDist.n;
        TH1D* h = hmean->ProjectionY(Form("hmean_slice_x%.f", x), binx, binx);
        gmean->AddPoint(x, h->GetMean());
    }
        

    TCanvas* cdEdx = new TCanvas("cdEdx", "cdEdx");

    cdEdx->Divide(2, 2);

    cdEdx->cd(1);
    gPad->SetMargin(0.15, 0.15, 0.15, 0.15);
    hdEdx->Draw("colz");

    cdEdx->cd(2);
    gPad->SetMargin(0.15, 0.15, 0.15, 0.15);
    hmean->Draw("colz");

    cdEdx->cd(3);
    gPad->SetMargin(0.15, 0.15, 0.15, 0.15);
    gdEdx->Draw("AL");

    cdEdx->cd(4);
    gPad->SetMargin(0.15, 0.15, 0.15, 0.15);
    gmean->Draw("AL");

    cdEdx->Write();

    std::cout << "writing " << out << "..." << std::endl;

    fin->Close();
    fout->Close();
}