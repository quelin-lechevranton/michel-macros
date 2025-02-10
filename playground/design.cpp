#include "graph.h"

Binning binTick(6000U/20, 0, 6000);
Binning binZ(1168U/8, 0, 300);
Binning bindE(0U, 0, 30);

int font = 82;
float font_size = 0.05;

float canvas_margin = 0.F;
float pad_margin = 0.3F;

std::vector<TH2F*> plotTH2F(std::vector<TPad*> ps, TH2F *h2, const char *draw_opt, Hits h) {
    std::vector<TH2F*> vh2(n_slice);
    for (unsigned int sl=0; sl<n_slice; sl++) {
        vh2[sl] = new TH2F(*h2);
        vh2[sl]->SetName(Form("%s_sl%u", h2->GetName(), sl));
    }
    for (unsigned int i=0; i<h.N; i++) {
        Hit hit = h.at(i);
        vh2[hit.slice]->Fill(hit.Z, hit.display_tick(), ADC_to_E * hit.adc);
    }
    for (unsigned int sl=0; sl<n_slice; sl++) {
        ps[sl]->cd();
        // vh2[sl]->SetMaximum(bindE.max);
        // vh2[sl]->SetMinimum(bindE.min);
        vh2[sl]->Draw(draw_opt);
    }
    return vh2;
}

bool up(unsigned int sl) {return sl >= 4;}
bool left(unsigned int sl) {return sl % 4 == 0;}
bool right(unsigned int sl) {return sl % 4 == 3;}
bool side(unsigned int sl) {return left(sl) || right(sl);}



void design(void) {

    TFile* fin = new TFile("pdvd_Muchecks_out_100_r20_in.root");
    TTree* event = fin->Get<TTree>("checks/Event");

    Hits EventHits;
    event->SetBranchAddress("NHit", &EventHits.N);
    event->SetBranchAddress("HitSlice", &EventHits.slice);
    event->SetBranchAddress("HitZ", &EventHits.Z);
    event->SetBranchAddress("HitChannel", &EventHits.channel);
    event->SetBranchAddress("HitTick", &EventHits.tick);
    event->SetBranchAddress("HitADC", &EventHits.adc);

    event->GetEntry(0);

    TCanvas *c = new TCanvas("cn", "ct", 20, 20, 1200, 800);

    // c->Divide(4, 2, 0, 0);



    float canvas_size = (1.F - 2*canvas_margin);
    float pad_factor = 1 - pad_margin;
    float side_pad_xsize = canvas_size / (2 + 2*pad_factor);
    float center_pad_xsize = pad_factor * side_pad_xsize;

    std::vector<float> xborder = {
        canvas_margin,
        canvas_margin + side_pad_xsize,
        canvas_margin + side_pad_xsize + center_pad_xsize,
        canvas_margin + side_pad_xsize + 2*center_pad_xsize,
        1 - canvas_margin
    };

    std::vector<TPad*> ps(n_slice);
    for (unsigned int sl=0; sl<n_slice; sl++) {

        float xlow = xborder[sl % 4];
        float xup = xborder[sl % 4 + 1];
        float ylow = sl < 4 ? canvas_margin : 0.5;
        float yup = sl < 4 ? 0.5 : 1-canvas_margin;

        c->cd();
        ps[sl] = new TPad(Form("pad_%d", sl), Form("pad_%d", sl), xlow, ylow, xup, yup);
        ps[sl]->Draw();
        ps[sl]->cd();
        
        TH1F* h = gPad->DrawFrame(
            binZ.min, binTick.min,
            binZ.max, binTick.max,
            // sl >= 4 ? ";Z (cm);6000 - ticks"
            //         : ";Z (cm);ticks"
            ";Z (cm); ticks"
        );

        // gPad->SetFillColorAlpha(sl+1, .3);
        gPad->SetLogz();
        TAxis* xa = h->GetXaxis();
        TAxis* ya = h->GetYaxis();
        if (up(sl)) {
            xa->SetTitleSize(0);
            xa->SetLabelSize(0);
        } else {
            xa->SetTitleFont(font);
            xa->SetLabelFont(font);
            if (side(sl)) {
                xa->SetLabelSize(font_size);
                xa->SetTitleSize(font_size);
                xa->SetTitleOffset(1);
                xa->SetLabelOffset(0.0);
            } else {
                xa->SetLabelSize(font_size / pad_factor);
                xa->SetTitleSize(font_size / pad_factor);
                xa->SetTitleOffset(1 * pad_factor);
                xa->SetLabelOffset( - 0.02); //by hand
            }
        }
        if (sl) {
            ya->SetTitleSize(0);
            ya->SetLabelSize(0);
        } else {
            ya->SetTitleFont(font);
            ya->SetTitleSize(font_size);
            ya->SetLabelFont(font);
            ya->SetLabelSize(font_size);
            std::cout << ya->GetLabelOffset() << std::endl;
        }
        if (sl == 4) {
            ya->SetTickSize(0);
        }

        gPad->Update();

        switch (sl) {
            case 0:
                gPad->SetMargin(pad_margin, 0, pad_margin, 0); break;
            case 1: case 2:
                gPad->SetMargin(0, 0, pad_margin, 0); break;
            case 3:
                gPad->SetMargin(0, pad_margin, pad_margin, 0); break;
            case 4:
                gPad->SetMargin(pad_margin, 0, 0, pad_margin); break;
            case 5: case 6:
                gPad->SetMargin(0, 0, 0, pad_margin); break;
            case 7:
                gPad->SetMargin(0, pad_margin, 0, pad_margin); break;
        }
    }
    ps[4]->cd();
    TGaxis *ya = new TGaxis(
        gPad->GetUxmin(),
        gPad->GetUymax(),
        gPad->GetUxmin(),
        gPad->GetUymin(),
        binTick.min,
        binTick.max);
    ya->SetLabelOffset(-0.025);
    ya->SetLabelFont(font);
    ya->SetLabelSize(font_size);
    ya->Draw();



    TH2F* h2 = new TH2F(
        "h2",
        "",
        binZ.n, binZ.min, binZ.max,
        binTick.n, binTick.min, binTick.max
    );
    // h2->SetFillColor(kGray);
    std::vector<TH2F*> vh2 = plotTH2F(ps, h2, "same col", EventHits);
}