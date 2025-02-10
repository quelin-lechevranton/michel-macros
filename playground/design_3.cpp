#include "graph.h"
#include <TExec.h>

Binning binTick(6000U/20, 0, 6000);
Binning binZ(1168U/8, 0, 300);
Binning bindE(0U, 0, 30);
Binning binTickCM(6000U/20, 0, 475);

int font = 82;
float font_size = 0.05;
int font_color = kWhite;
int background_color = kBlack;
TColor c;
int low_palette[] = {c.GetColor(95,95,255), c.GetColor(127,127,255), c.GetColor(159,159,255), c.GetColor(191,191,255), c.GetColor(223,223,255), c.GetColor(255,255,255)};
int upp_palette[] = {c.GetColor(255,95,95), c.GetColor(255,127,127), c.GetColor(255,159,159), c.GetColor(255,191,191), c.GetColor(255,223,223), c.GetColor(255,255,255)};

unsigned int n_pad = n_slice /2;

float canvas_margin = 0.F;
float pad_margin_outer = 0.2F;
float pad_margin_inner = 0.05F;

bool up(unsigned int sl) {return sl >= 4;}
bool left(unsigned int p) {return p == 0;}
bool right(unsigned int p) {return p == 3;}
bool side(unsigned int p) {return left(p) || right(p);}

std::vector<TH2F*> plotTH2F(std::vector<TPad*> ps, TH2F *h2, const char *draw_opt, Hits h) {
    std::vector<TH2F*> vh2(n_slice);

    TExec *low_style = new TExec("low_style", "gStyle->SetPalette(7, low_palette);");
    TExec *upp_style = new TExec("upp_style", "gStyle->SetPalette(7, upp_palette);");

    for (unsigned int sl=0; sl<n_slice; sl++) {
        vh2[sl] = new TH2F(*h2);
        vh2[sl]->SetName(Form("%s_sl%u", h2->GetName(), sl));
    }
    for (unsigned int i=0; i<h.N; i++) {
        Hit hit = h.at(i);
        vh2[hit.slice]->Fill(hit.Z, hit.tick_cm(), ADC_to_E * hit.adc);
    }
    for (unsigned int sl=0; sl<n_slice; sl++) {
        ps[sl % 4]->cd();
        if (up(sl)) upp_style->Draw();
        else low_style->Draw();
        vh2[sl]->Draw(draw_opt);
    }
    return vh2;
}


// TString in = "pdvd_Muchecks_out_100_r20_in.root";

// ROOT::RDataFrame event("checks/Event", in);
// ROOT::RDataFrame muon("checks/Muon", in);

void design_3(void) {

    TFile* fin = new TFile("pdvd_Muchecks_out_100_r20_in.root");
    TTree* event = fin->Get<TTree>("checks/Event");

    Hits EventHits;
    event->SetBranchAddress("NHit", &EventHits.N);
    event->SetBranchAddress("HitSlice", &EventHits.slice);
    event->SetBranchAddress("HitZ", &EventHits.Z);
    event->SetBranchAddress("HitChannel", &EventHits.channel);
    event->SetBranchAddress("HitTick", &EventHits.tick);
    event->SetBranchAddress("HitADC", &EventHits.adc);

    event->GetEntry(1);

    TCanvas *c = new TCanvas("cn", "ct", 20, 20, 1200, 700);

    float canvas_size = (1.F - 2*canvas_margin);
    float pad_factor = 1 - pad_margin_outer + pad_margin_inner;
    float side_pad_xsize = canvas_size / (2 + 2*pad_factor);
    float center_pad_xsize = pad_factor * side_pad_xsize;

    std::vector<float> xborder = {
        canvas_margin,
        canvas_margin + side_pad_xsize,
        canvas_margin + side_pad_xsize + center_pad_xsize,
        canvas_margin + side_pad_xsize + 2*center_pad_xsize,
        1 - canvas_margin
    };

    std::vector<TPad*> ps(n_pad);
    for (unsigned int p=0; p<n_pad; p++) {

        float xlow = xborder[p];
        float xup = xborder[p+1];
        float ylow = canvas_margin;
        float yup = 1-canvas_margin;

        c->cd();
        ps[p] = new TPad(Form("pad_%d", p), Form("pad_%d", p), xlow, ylow, xup, yup);
        ps[p]->Draw();
        ps[p]->cd();

        TH1F* f = gPad->DrawFrame(
            binZ.min, binTickCM.min,
            binZ.max, binTickCM.max,
            ";Z (cm); T (cm)"
        );

        TAxis* x = f->GetXaxis();
        TAxis* y = f->GetYaxis();

        x->SetTitleFont(font);
        x->SetTitleSize(font_size);
        x->SetTitleOffset(1);

        x->SetLabelFont(font);
        x->SetLabelSize(font_size);
        x->SetTickLength(0.03);

        y->SetTitleFont(font);
        y->SetTitleSize(font_size);
        y->SetTitleOffset(2.2);

        y->SetLabelFont(font);
        y->SetLabelSize(font_size);
        y->SetLabelOffset(0.05);
        y->SetTickLength(0.06);

        gPad->SetMargin(pad_margin_inner, pad_margin_inner, pad_margin_outer, pad_margin_outer);
        if (left(p)) {
            gPad->SetLeftMargin(pad_margin_outer);
        } else {
            if (right(p)) {
                gPad->SetRightMargin(pad_margin_outer);
            } else {
                x->SetTitleSize(font_size / pad_factor);
                x->SetLabelSize(font_size / pad_factor);
                x->SetTitleOffset(1 * pad_factor);
                x->SetLabelOffset( - 0.002); //by hand
                x->SetTickLength(0.03 * pad_factor);
            }
            y->SetTitleSize(0);
            y->SetLabelSize(0);
        }
    }


    TH2F* h2 = new TH2F(
        "h2",
        "",
        binZ.n, binZ.min, binZ.max,
        binTickCM.n, binTickCM.min, binTickCM.max
    );
    std::vector<TH2F*> vh2 = plotTH2F(ps, h2, "same pfc col", EventHits);

    c->SaveAs("design_3.pdf");
}