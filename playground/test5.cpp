#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TText.h>

TFile* f = new TFile("pdvd_Michecks.root");
TTree* michel = f->Get<TTree>("checks/Michel");
TTree* event = f->Get<TTree>("checks/Event");

size_t n_event = 40;
size_t n_section = 4;
size_t n_channel = 12288;

struct Binning {
    int n;
    double min, max;
};

std::vector<Binning> binChan = {{1168, 1904, 3071},
                                {1168, 4976, 6143},
                                {1168, 8048, 9215},
                                {1168, 11120, 12287} };

std::vector<Binning> binTPC = { {292, 1904, 2195},
                                {292, 2196, 2487},
                                {292, 2488, 2779},
                                {292, 2780, 3071},
                                {292, 4976, 5267},
                                {292, 5268, 5559},
                                {292, 5560, 5851},
                                {292, 5852, 6143},
                                {292, 8048, 8339},
                                {292, 8340, 8631},
                                {292, 8632, 8923},
                                {292, 8924, 9215},
                                {292, 11120, 11411},
                                {292, 11412, 11703},
                                {292, 11704, 11995},
                                {292, 11996, 12287} };

float tick_radius = 50. / 0.5; // us / (us/tick)
float channel_radius = 20. / 0.5; // cm / (chan/cm)


Binning binTick = {2000, 0, 6000};
Binning binEnergy = {30, 0, 90};

enum Tag {kMCP, kSphere, kMuon, kNTag};
std::vector<short> markerStyle = {kFullDoubleDiamond, kFullFourTrianglesX, kOpenSquare};
std::vector<float> markerSize = {2.5, 2.5, .8, 2};
std::vector<short> markerColor = {kAzure-5, kViolet+5, kPink-5};

size_t GetSection(int ch) {
    return size_t(4.*ch / n_channel);
}

void test5() {

    TFile* f = new TFile("test5_cpp.root", "RECREATE");

    TCanvas* ch = new TCanvas("ch");

    TH1F    *hTrueEnergyE       = new TH1F("hTrueEnergyE", "True Energy;K (MeV);#", binEnergy.n, binEnergy.min, binEnergy.max),
            *hMCPEnergyE        = new TH1F("hMCPEnergyE", "Energy from Hits from MCParticle;ADC.Tick (MeV);#", binEnergy.n, binEnergy.min, binEnergy.max),
            *hSphereEnergyE     = new TH1F("hSphereEnergyE", "Energy from Hits selected in sphere;ADC.Tick (MeV);#", binEnergy.n, binEnergy.min, binEnergy.max),
            *hTrueEnergyP       = new TH1F("hTrueEnergyP", "True Energy;K (MeV);#", binEnergy.n, binEnergy.min, binEnergy.max),
            *hMCPEnergyP        = new TH1F("hMCPEnergyP", "Energy from Hits from MCParticle;ADC.Tick (MeV);#", binEnergy.n, binEnergy.min, binEnergy.max),
            *hSphereEnergyP     = new TH1F("hSphereEnergyP", "Energy from Hits selected in sphere;ADC.Tick (MeV);#", binEnergy.n, binEnergy.min, binEnergy.max);

    std::vector<TH1F*> hEnergyE = {hTrueEnergyE, hMCPEnergyE, hSphereEnergyE};
    std::vector<TH1F*> hEnergyP = {hTrueEnergyP, hMCPEnergyP, hSphereEnergyP};

    for (TH1F* hE : hEnergyE) {
        hE->SetLineColor(kAzure-5);
        hE->SetLineWidth(3);
    }
    for (TH1F* hP : hEnergyP) {
        hP->SetLineColor(kPink-5);
        hP->SetLineWidth(3);
    }

    ch->Divide(2, 2);

    ch->cd(1);
    gPad->SetLogy();
    michel->Draw("TrueEnergy>>hTrueEnergyP","IsPositron==1");
    michel->Draw("TrueEnergy>>hTrueEnergyE","IsPositron==0","same");

    ch->cd(2);
    gPad->SetLogy();
    michel->Draw("MCPEnergy>>hMCPEnergyP","IsPositron==1");
    michel->Draw("MCPEnergy>>hMCPEnergyE","IsPositron==0","same");

    ch->cd(3);
    gPad->SetLogy();
    michel->Draw("SphereEnergy>>hSphereEnergyP","IsPositron==1");
    michel->Draw("SphereEnergy>>hSphereEnergyE","IsPositron==0","same");

    ch->Write();


    TCanvas* cCorr = new TCanvas("cCorr");

    TH2F* hCorrEnergyMCPEnergy = new TH2F("hCorrEnergyMCPEnergy", ";True Energy (MeV);MCP Energy (MeV)", binEnergy.n, binEnergy.min, binEnergy.max, binEnergy.n, binEnergy.min, binEnergy.max);
    TH2F* hCorrADCNHit = new TH2F("hCorrADCNHit", ";MCP ADC.Tick (MeV);MCP hit number (#)", binEnergy.n, binEnergy.min, binEnergy.max, 20, 0, 60);
    TH2F* hCorrEnergyNHit = new TH2F("hCorrEnergyNHit", ";True Energy (MeV);MCP hit number (#)", binEnergy.n, binEnergy.min, binEnergy.max, 20, 0, 60);

    cCorr->Divide(2, 2);

    cCorr->cd(1);
    // gPad->SetLogz();
    gStyle->SetOptStat(kFALSE);
    michel->Draw("MCPEnergy:TrueEnergy>>hCorrEnergyMCPEnergy", "", "colz");

    cCorr->cd(2);
    // gPad->SetLogz();
    gStyle->SetOptStat(kFALSE);
    michel->Draw("MCPNHit:MCPEnergy>>hCorrADCNHit", "", "colz");

    cCorr->cd(3);
    // gPad->SetLogz();
    gStyle->SetOptStat(kFALSE);
    michel->Draw("MCPNHit:TrueEnergy>>hCorrEnergyNHit", "", "colz");

    cCorr->Write();

    std::vector<TCanvas*> c2d(n_event);

    size_t nhit;
    std::vector<int> *channels = nullptr;
    std::vector<float> *ticks = nullptr,
                  *adcs = nullptr;

    size_t nmichel;
    std::vector<size_t> *event_imichel = nullptr;

    event->SetBranchAddress("NHit", &nhit);
    event->SetBranchAddress("HitChannel", &channels);
    event->SetBranchAddress("HitTick", &ticks);
    event->SetBranchAddress("HitADC", &adcs);

    event->SetBranchAddress("NMichel", &nmichel);
    event->SetBranchAddress("iMichel", &event_imichel);

    std::vector<size_t>                 michel_nhit(kNTag);
    std::vector<std::vector<int>*>      michel_channels(kNTag);
    std::vector<std::vector<float>*>    michel_ticks(kNTag),
                                        michel_adcs(kNTag);


    michel->SetBranchAddress("MCPNHit", &michel_nhit[kMCP]);
    michel->SetBranchAddress("MCPHitChannel", &michel_channels[kMCP]);
    michel->SetBranchAddress("MCPHitTick", &michel_ticks[kMCP]);
    michel->SetBranchAddress("MCPHitADC", &michel_adcs[kMCP]);

    michel->SetBranchAddress("SphereNHit", &michel_nhit[kSphere]);
    michel->SetBranchAddress("SphereHitChannel", &michel_channels[kSphere]);
    michel->SetBranchAddress("SphereHitTick", &michel_ticks[kSphere]);
    michel->SetBranchAddress("SphereHitADC", &michel_adcs[kSphere]);

    michel->SetBranchAddress("MuonNHit", &michel_nhit[kMuon]);
    michel->SetBranchAddress("MuonHitChannel", &michel_channels[kMuon]);
    michel->SetBranchAddress("MuonHitTick", &michel_ticks[kMuon]);
    michel->SetBranchAddress("MuonHitADC", &michel_adcs[kMuon]);

    int michel_muon_channel;
    float michel_muon_tick;

    michel->SetBranchAddress("MuonEndChannel", &michel_muon_channel);
    michel->SetBranchAddress("MuonEndTick", &michel_muon_tick);

    for (size_t i_event=0; i_event<n_event; i_event++) {

        event->GetEntry(i_event);

        c2d[i_event] = new TCanvas(
            Form("c2d_e%ld", i_event),
            Form("Event#%ld w/ %ld michels", i_event, nmichel)
        );

        std::vector<TH2F*> h2(n_section);
        for (size_t s=0; s<n_section; s++) {
            h2[s] = new TH2F(
                Form("h2_e%ld_s%ld", i_event, s),
                "",
                binTick.n, binTick.min, binTick.max,
                binChan[s].n, binChan[s].min, binChan[s].max
            );
        }

        for (size_t i=0; i<nhit; i++) {
            if (GetSection(channels->at(i)) >= n_section) continue;
            h2[GetSection(channels->at(i))]->Fill(ticks->at(i), channels->at(i), adcs->at(i));
        }

        std::vector<std::vector<TGraph*>> g(kNTag, std::vector<TGraph*>(nmichel));
        for (size_t t=0; t<kNTag; t++) {

            size_t m = 0;
            for (size_t imichel : *event_imichel) {
                michel->GetEntry(imichel);
                
                // g[t][m] = new TGraph(michel_nhit[t], &michel_ticks[t]->at(0), &michel_channels[t]->at(0));
                g[t][m] = new TGraph();
                for (size_t i=0; i<michel_nhit[t]; i++) {
                    g[t][m]->AddPoint(michel_ticks[t]->at(i), michel_channels[t]->at(i));
                }
                m++;
            }
        }

        c2d[i_event]->Divide(2, 2);
        for (size_t s=0; s<n_section; s++) {
            c2d[i_event]->cd(s+1);

            gPad->DrawFrame(
                binTick.min, binChan[s].min,
                binTick.max, binChan[s].max,
                Form("S%ld, W;Tick;Channel",s)
            );
            gPad->SetLogz();

            for (size_t i=0; i<4; i++) {
                size_t tpc = s*n_section+i;
                TLine* line = new TLine(binTick.min, binTPC[tpc].max, binTick.max, binTPC[tpc].max);
                line->Draw();

                TText* text = new TText(binTick.min + 200, binTPC[tpc].max - 20,Form("TPC %ld", tpc));
                text->SetTextAlign(13); // left - top adjusted
                text->SetTextFont(42); // Helvetica thin
                text->SetTextSize(.04); // Size is scaled by the pad size
                text->Draw();
            }

            h2[s]->SetMinimum(.1);
            h2[s]->Draw("colz same");

            for (size_t t=0; t<kNTag; t++) {
                for (size_t m=0; m<(nmichel); m++) {
                    if (g[t][m]->GetN() == 0) continue; 
                    g[t][m]->SetName(Form("g_e%ld_t%ld_m%ld", i_event, t, m));
                    g[t][m]->SetEditable(kFALSE);
                    g[t][m]->SetMarkerStyle(markerStyle[t]);
                    g[t][m]->SetMarkerSize(markerSize[t]);
                    g[t][m]->SetMarkerColor(markerColor[t]+m);
                    g[t][m]->Draw("P same");


                    TEllipse* ellipse = new TEllipse(michel_muon_tick, michel_muon_channel, tick_radius, channel_radius);
                    ellipse->SetFillStyle(0);
                    ellipse->SetLineWidth(2);
                    ellipse->SetLineColor(markerColor[kSphere]);
                    ellipse->Draw();

                }
            }
        }

        c2d[i_event]->Write();
    }


    f->Close();
}