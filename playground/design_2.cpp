#include "graph.h"

Binning binTick(6000U/20, 0, 6000);
Binning binZ(1168U/8, 0, 300);
Binning bindE(0U, 0, 30);

int font = 82;
float font_size = 0.05;
int font_color = kBlack;
int background_color = kWhite;
TColor c;
int palette[] = {background_color, c.GetColor(95,95,255), c.GetColor(127,127,255), c.GetColor(159,159,255), c.GetColor(191,191,255), c.GetColor(223,223,255), c.GetColor(255,255,255)};

    
float canvas_margin = 0.F;
float pad_margin_outer = 0.2F;
float pad_margin_inner = 0.05F;

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
        gStyle->SetPalette(7, palette);
        vh2[sl]->Draw(draw_opt);
    }
    return vh2;
}

bool up(unsigned int sl) {return sl >= 4;}
bool left(unsigned int sl) {return sl % 4 == 0;}
bool right(unsigned int sl) {return sl % 4 == 3;}
bool side(unsigned int sl) {return left(sl) || right(sl);}



void design_2(void) {

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
        
        TH1F* f = gPad->DrawFrame(
            binZ.min, binTick.min,
            binZ.max, binTick.max,
            ";Z (cm); ticks"
        );

        // gPad->SetFillColorAlpha(sl+1, .3);
        gPad->SetFillColor(background_color);
        // f->SetFillColor(background_color); // now working
        // f->SetLineColor(font_color);


        gPad->SetLogz();
        TAxis* xa = f->GetXaxis();
        TAxis* ya = f->GetYaxis();
        xa->SetAxisColor(font_color);
        ya->SetAxisColor(font_color);

        // xa->SetLimits(binZ.min-1, binZ.max+2); // show limits labels

        gPad->SetMargin(pad_margin_inner, pad_margin_inner, pad_margin_inner, pad_margin_inner);
        if (up(sl)) gPad->SetTopMargin(pad_margin_outer);
        else gPad->SetBottomMargin(pad_margin_outer);
        if (left(sl)) gPad->SetLeftMargin(pad_margin_outer);
        else if (right(sl)) gPad->SetRightMargin(pad_margin_outer);


        // WIP

        // xa->SetTitleSize(0);
        // xa->SetLabelSize(0);
        // xa->SetTickLength(0);

        // ya->SetTitleSize(0);
        // ya->SetLabelSize(0);
        // ya->SetTickLength(0);

        // if (left(sl)) {
        //     TGaxis *gya = new TGaxis();
        //     gya->SetTitleFont(font);
        //     gya->SetTitleSize(font_size);
        //     gya->SetTitleColor(font_color);

        //     gya->SetLabelFont(font);
        //     gya->SetLabelSize(font_size);
        //     gya->SetLabelColor(font_color);

        //     if (up(sl)) {
        //         gya->SetLabelOffset(-0.029); //by hand
        //         gya->DrawAxis(
        //             gPad->GetUxmin(),
        //             gPad->GetUymax(),
        //             gPad->GetUxmin(),
        //             gPad->GetUymin(),
        //             binTick.min,
        //             binTick.max);
        //     } else {
        //         gya->DrawAxis(
        //             gPad->GetUxmin(),
        //             gPad->GetUymin(),
        //             gPad->GetUxmin(),
        //             gPad->GetUymax(),
        //             binTick.min,
        //             binTick.max);
        //     }
        // }
        // if (!up(sl)) {
        //     TGaxis *gxa = new TGaxis();
        //     gxa->SetTitleFont(font);
        //     gxa->SetTitleColor(font_color);

        //     gxa->SetLabelFont(font);
        //     gxa->SetLabelColor(font_color);

        //     if (side(sl)) {
        //         gxa->SetTitleSize(font_size);
        //         gxa->SetLabelSize(font_size);
        //     } else {
        //         gxa->SetTitleSize(font_size / pad_factor);
        //         gxa->SetLabelSize(font_size / pad_factor);
        //     }

        //     gxa->DrawAxis(
        //         gPad->GetUxmin(),
        //         gPad->GetUymin(),
        //         gPad->GetUxmax(),
        //         gPad->GetUymin(),
        //         binZ.min,
        //         binZ.max);
        // }


        // END WIP


        if (up(sl)) {
            xa->SetTitleSize(0);
            xa->SetLabelSize(0);
        } else {
            xa->SetTitleFont(font);
            xa->SetLabelFont(font);
            xa->SetTitleColor(font_color);
            xa->SetLabelColor(font_color);
            if (side(sl)) {
                xa->SetLabelSize(font_size);
                xa->SetTitleSize(font_size);
                xa->SetTitleOffset(1);
            } else {
                xa->SetLabelSize(font_size / pad_factor);
                xa->SetTitleSize(font_size / pad_factor);
                xa->SetTitleOffset(1 * pad_factor);
                xa->SetLabelOffset( - 0.002); //by hand
            }
        }
        if (sl) {
            ya->SetTitleSize(0);
            ya->SetLabelSize(0);
        } else {
            ya->SetTitleFont(font);
            ya->SetTitleSize(font_size);
            ya->SetTitleColor(font_color);
            ya->SetLabelFont(font);
            ya->SetLabelSize(font_size);
            ya->SetLabelColor(font_color);
            ya->SetTitleOffset(2);
            ya->SetLabelOffset(0.02);
        }
        if (sl == 4) {
            ya->SetTickSize(0);
            ya->SetTitleFont(font);
            ya->SetTitleSize(font_size);
            ya->SetTitleColor(font_color);
            ya->SetTitleOffset(2);
        }

        if (sl == 4) {
            TGaxis *ya = new TGaxis(
                gPad->GetUxmin(),
                gPad->GetUymax(),
                gPad->GetUxmin(),
                gPad->GetUymin(),
                binTick.min,
                binTick.max);
            ya->SetLabelOffset(-0.040);
            ya->SetLabelFont(font);
            ya->SetLabelSize(font_size);
            ya->SetLabelColor(font_color);
            ya->SetLineColor(font_color);
            ya->Draw();
        }
    }



    TH2F* h2 = new TH2F(
        "h2",
        "",
        binZ.n, binZ.min, binZ.max,
        binTick.n, binTick.min, binTick.max
    );
    // h2->SetFillColor(kGray);
    std::vector<TH2F*> vh2 = plotTH2F(ps, h2, "same col", EventHits);

    c->SaveAs("design_2.pdf");
}