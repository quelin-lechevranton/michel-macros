#include "utilities.h"

#include <TColor.h>

#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TText.h>

#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <TEllipse.h>

#include <TExec.h>

enum Slices { kLowLeft, kLowMidLeft, kLowMidRight, kLowRight, kUpLeft, kUpMidLeft, kUpMidRight, kUpRight };

// TColor col;
// Int_t low_palette[6] = {col.GetColor(95,95,255), col.GetColor(127,127,255), col.GetColor(159,159,255), col.GetColor(191,191,255), col.GetColor(223,223,255), col.GetColor(255,255,255)};
// Int_t upp_palette[6] = {col.GetColor(255,95,95), col.GetColor(255,127,127), col.GetColor(255,159,159), col.GetColor(255,191,191), col.GetColor(255,223,223), col.GetColor(255,255,255)};

bool rev = true;
bool cm = false;

float step = channel_pitch; // cm
Binning binT(0.F, 6000.F, 2*step / sampling_rate / drift_velocity);
Binning binZ(0.F, 300.F, step);
Binning binX(0.F, 480.F, 2*step); //ticks in cm
Binning binY = cm ? binX : binT; // Y axis

int font = 82;
float font_size = 0.05F;
// int font_color = kBlack;
// int background_color = kWhite;

float canvas_margin = 0.F;
float pad_margin_outer = 0.2F;
float pad_margin_inner = 0.05F;

float canvas_width = (1.F - 2*canvas_margin);
float pad_factor = 1 - pad_margin_outer + pad_margin_inner;
float side_pad_width = canvas_width / (2 + 2*pad_factor);
float center_pad_width = pad_factor * side_pad_width;

bool up(unsigned int sl) {return sl >= kUpLeft;}
bool left(unsigned int sl) {return sl % 4 == kLowLeft;}
bool right(unsigned int sl) {return sl % 4 == kLowRight;}
bool side(unsigned int sl) {return left(sl) || right(sl);}

float Hit::display_tick() {
    return rev && up(slice) ? (
        cm ? rev_tick_cm() : rev_tick()
    ) : (
        cm ? tick_cm() : tick
    );
}
float Hit::display4_tick() {
    return cm ? tick_cm() : tick;
}

// enum { kADC, kCVN};

// void SetPalette(int opt) {
//     int n_cont = 125;
//     int n_stop = opt == kADC ? 6
//                 : (opt == kCVN ? 3
//                 : 0);

//     double stop[n_stop], R[n_stop], G[n_stop], B[n_stop];

//     if (opt == kADC) {
//         double stop[]   = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
//         double R[]     = { 0.2, 0.35, 0.5, 0.65, 0.8, 0.95 };
//         double G[]     = { 0.2, 0.45, 0.7, 0.9, 0.8, 0.7 };
//         double B[]     = { 1.0, 0.95, 0.85, 0.8, 0.75, 0.7 };
//     } else if (opt == kCVN) {
//         double stop[] = { 0.0, 0.5, 1.0 };
//         double R[] = { 1.0, 1.0, 0.0 };
//         double G[] = { 0.0, 1.0, 0.0 };
//         double B[] = { 0.0, 1.0, 1.0 };
//     }

//     TColor::CreateGradientColorTable(n_stop, stop, R, G, B, n_cont);
//     gStyle->SetNumberContours(n_cont);
// }

std::vector<TPad*> drawFrame(TCanvas* c) {
    std::vector<TPad*> ps(n_slice);

    std::vector<float> xborder = {
        canvas_margin,
        canvas_margin + side_pad_width,
        canvas_margin + side_pad_width + center_pad_width,
        canvas_margin + side_pad_width + 2*center_pad_width,
        1 - canvas_margin
    };
    std::vector<float> yborder = {
        canvas_margin,
        0.5,
        1 - canvas_margin
    };

    for (unsigned int sl=0; sl<n_slice; sl++) {

        TPad* p = ps[sl] = new TPad(
            Form("pad_%d", sl),
            Form("pad_%d", sl),
            xborder[sl % 4], yborder[sl / 4],
            xborder[sl % 4 + 1], yborder[sl / 4 + 1]
        );

        c->cd();
        p->Draw();
        p->cd();

        const char* yt = Form("%sT (%s)", rev && up(sl) ? "reverse " : "", cm ? "cm" : "tick");

        TH1F* f = gPad->DrawFrame(
            binZ.min, binY.min,
            binZ.max, binY.max,
            Form(";Z [cm];%s",yt)
        );
        gPad->SetLogz();
        
        gPad->SetMargin(
            pad_margin_inner, 
            pad_margin_inner, 
            pad_margin_inner, 
            pad_margin_inner
        );
        if (left(sl))
            gPad->SetLeftMargin(pad_margin_outer);
        if (right(sl))
            gPad->SetRightMargin(pad_margin_outer);
        if (up(sl))
            gPad->SetTopMargin(pad_margin_outer);
        else
            gPad->SetBottomMargin(pad_margin_outer);

        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);

        float fs = side(sl) ? font_size : font_size / pad_factor;

        f->SetTitleFont(font, "xyz");
        f->SetLabelFont(font, "xyz");
        f->SetTitleSize(fs, "xyz");
        f->SetLabelSize(fs, "xyz");

        TAxis* xa = f->GetXaxis();
        TAxis* ya = f->GetYaxis();

        if (side(sl))
            xa->SetTitleOffset(1);
        else {
            xa->SetTitleOffset(1 * pad_factor);
            xa->SetLabelOffset( - 0.002); //by hand
        }
        ya->SetTitleOffset(2);
        ya->SetLabelOffset(0.01);

        if (!left(sl)) {
            ya->SetTitleSize(0);
            ya->SetLabelSize(0);
        }
        if (up(sl)) {
            xa->SetTitleSize(0);
            xa->SetLabelSize(0);
        }

        float title_offset = 0.015;
        TText* t = new TText(
            left(sl) ? pad_margin_outer : pad_margin_inner,
            (up(sl) ? 1 - pad_margin_outer : 1 - pad_margin_inner) + title_offset,
            Form("slice %u", sl)
        );
        t->SetNDC(); // coordinate in pad size
        t->SetTextFont(font);
        t->SetTextSize(fs);
        t->SetTextAlign(kHAlignLeft + kVAlignBottom);
        t->Draw();
    }
    return ps;
}

template <typename H>
void drawTH2F(std::vector<TPad*> ps, const char *name, H h, int pal = kDeepSea, unsigned pal_sl=kUpRight) {

    TExec* ex_pal = new TExec("ex_pal", Form("gStyle->SetPalette(%d)", pal));

    std::vector<TH2F*> vh2(n_slice);
    for (unsigned sl=0; sl<n_slice; sl++) {
        vh2[sl] = new TH2F(
            Form("%s_sl%u", name, sl),
            "",
            binZ.n, binZ.min, binZ.max,
            binY.n, binY.min, binY.max
        );
    }
    for (unsigned i=0; i<h.N; i++) {
        Hit hit = h.at(i);
        vh2[hit.slice]->Fill(hit.Z, hit.display_tick(), hit.adc);
    }
    double max=0;
    for (TH2F* h2 : vh2)
        max = h2->GetMaximum() > max ? h2->GetMaximum() : max;
    for (TH2F* h2 : vh2) {
        h2->SetMinimum(0);
        h2->SetMaximum(max);
    }
    for (unsigned sl=0; sl<n_slice; sl++) {
        ps[sl]->cd();

        ex_pal->Draw();

        if (sl == pal_sl) {
            TAxis* za = vh2[sl]->GetZaxis();
            za->SetTitle("dE [ADC.tick]");
            za->SetTitleFont(font);
            za->SetTitleSize(font_size);
            za->SetTitleOffset(1.5);

            za->SetLabelFont(font);
            za->SetLabelSize(font_size);
            za->SetLabelOffset(0);

            vh2[sl]->Draw("same pfc colz");
        } else {
            vh2[sl]->Draw("same pfc col");
        }
    }
}

template <typename H>
void drawTGraph(std::vector<TPad*> ps, TGraph *g, const char *draw_opt, H h) {
    std::vector<TGraph*> vg(n_slice);
    for (unsigned int sl=0; sl<n_slice; sl++) {
        vg[sl] = new TGraph(*g);
        vg[sl]->SetName(Form("%s_sl%u", g->GetName(), sl));
        vg[sl]->SetEditable(kFALSE);
    }
    // for (unsigned int i=0; i<h.N; i++) {
    for (unsigned int i=0; i<h.N; i++) {
        Hit hit = h.at(i);
        vg[hit.slice]->AddPoint(hit.Z, hit.display_tick());
    }
    for (unsigned int sl=0; sl<n_slice; sl++) {
        if (vg[sl]->GetN() == 0) continue;
        ps[sl]->cd();
        vg[sl]->Draw(draw_opt);
    }
}

void drawTEllipse(std::vector<TPad*> ps, TEllipse *el, Hit h, float rz) {
    TEllipse *e = new TEllipse(*el);
    e->SetX1(h.Z);
    e->SetY1(h.display_tick());
    e->SetR1(rz);
    e->SetR2(cm ? rz : rz / sampling_rate / drift_velocity);
    ps[h.slice]->cd();
    e->Draw();
}















// Window

TPad* drawFrame(TCanvas *c, unsigned int sl, Binning bz, Binning bt) {

    TPad* p = new TPad(
        Form("pad_%d", sl),
        Form("pad_%d", sl),
        canvas_margin, canvas_margin,
        1 - canvas_margin, 1 - canvas_margin
    );

    c->cd();
    p->Draw();
    p->cd();

    const char* yt = Form("%sT [%s]", rev && up(sl) ? "reverse " : "", cm ? "cm" : "tick");

    TH1F* f = gPad->DrawFrame(
        bz.min, bt.min,
        bz.max, bt.max,
        Form(";Z [cm];%s",yt)
    );
    gPad->SetLogz();

    gPad->SetMargin(
        pad_margin_outer,
        pad_margin_outer,
        pad_margin_outer,
        pad_margin_outer
    );

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    f->SetTitleFont(font, "xyz");
    f->SetTitleSize(font_size, "xyz");
    f->SetLabelFont(font, "xyz");
    f->SetLabelSize(font_size, "xyz");

    TAxis* xa = f->GetXaxis();
    TAxis* ya = f->GetYaxis();

    xa->SetTitleOffset(1);
    ya->SetTitleOffset(1.5);
    ya->SetLabelOffset(0.01);

    return p;
}

template <typename H>
void drawTH2F(TPad *p, const char *name, H h, unsigned int sl, Binning bz, Binning bt, int pal = kDeepSea) {

    TExec* ex_pal = new TExec("ex_pal", Form("gStyle->SetPalette(%d)", pal));

    TH2F* h2 = new TH2F(
        Form("%s_sl%u", name, sl),
        "",
        bz.n, bz.min, bz.max,
        bt.n, bt.min, bt.max
    );
    for (unsigned int i=0; i<h.N; i++) {
        Hit hit = h.at(i);
        if (hit.slice != sl) continue;
        if (hit.Z < bz.min or hit.Z > bz.max) continue;
        h2->Fill(hit.Z, hit.display_tick(), hit.adc);
    }
    p->cd();
    ex_pal->Draw();
    h2->Draw("same pfc colz");
}

template <typename H>
void drawTGraph(TPad *p, TGraph *gr, const char *draw_opt, H h, unsigned int sl, Binning bz, Binning bt) {
    TGraph* g = new TGraph(*gr);
    g->SetEditable(kFALSE);
    for (unsigned int i=0; i<h.N; i++) {
        Hit hit = h.at(i);
        if (hit.slice != sl) continue;
        if (hit.Z < bz.min or hit.Z > bz.max) continue;
        g->AddPoint(hit.Z, hit.display_tick());
    }
    if (g->GetN() == 0) return;
    p->cd();
    g->Draw(draw_opt);
}




// Four Frames

std::vector<TPad*> draw4Frame(TCanvas *c) {
    std::vector<TPad*> ps(n_slice/2);

    std::vector<float> xborder = {
        canvas_margin,
        canvas_margin + side_pad_width,
        canvas_margin + side_pad_width + center_pad_width,
        canvas_margin + side_pad_width + 2*center_pad_width,
        1 - canvas_margin
    };

    for (unsigned int sl=0; sl<n_slice/2; sl++) {

        TPad* p = ps[sl] = new TPad(
            Form("pad_%d", sl),
            Form("pad_%d", sl),
            xborder[sl % 4], canvas_margin,
            xborder[sl % 4 + 1], 1 - canvas_margin
        );

        c->cd();
        p->Draw();
        p->cd();

        TH1F* f = gPad->DrawFrame(
            binZ.min, binY.min,
            binZ.max, binY.max,
            Form(";Z [cm];T [%s]", cm ? "cm" : "tick")
        );
        gPad->SetLogz();

        gPad->SetMargin(
            pad_margin_inner,
            pad_margin_inner,
            pad_margin_outer,
            pad_margin_outer
        );
        if (left(sl))
            gPad->SetLeftMargin(pad_margin_outer);
        if (right(sl))
            gPad->SetRightMargin(pad_margin_outer);
        
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);

        float fs = side(sl) ? font_size : font_size / pad_factor;

        f->SetTitleFont(font, "xyz");
        f->SetLabelFont(font, "xyz");
        f->SetTitleSize(fs, "xyz");
        f->SetLabelSize(fs, "xyz");

        TAxis* xa = f->GetXaxis();
        TAxis* ya = f->GetYaxis();

        if (side(sl))
            xa->SetTitleOffset(1);
        else {
            xa->SetTitleOffset(1 * pad_factor);
            xa->SetLabelOffset( - 0.002); //by hand
        }
        ya->SetTitleOffset(2);
        ya->SetLabelOffset(0.01);

        if (!left(sl)) {
            ya->SetTitleSize(0);
            ya->SetLabelSize(0);
        }

        float title_offset = 0.015;
        TText* t = new TText(
            left(sl) ? pad_margin_outer : pad_margin_inner,
            1 - pad_margin_outer + title_offset,
            Form("slices %u & %u", sl, sl+4)
        );
        t->SetNDC(); // coordinate in pad size
        t->SetTextFont(font);
        t->SetTextSize(fs);
        t->SetTextAlign(kHAlignLeft + kVAlignBottom);
        t->Draw();
    }
    return ps;
}

template <typename H>
void draw4TH2F(std::vector<TPad*> ps, const char *name, H h) {

    TExec* ex_up_pal = new TExec("ex_up_pal", "gStyle->SetPalette(kCopper);");
    TExec* ex_low_pal = new TExec("ex_low_pal", "gStyle->SetPalette(kDeepSea);");

    std::vector<TH2F*> vh2(n_slice);
    for (unsigned int sl=0; sl<n_slice; sl++) {
        vh2[sl] = new TH2F(
            Form("%s_sl%u", name, sl),
            "",
            binZ.n, binZ.min, binZ.max,
            binY.n, binY.min, binY.max
        );
    }
    for (unsigned int i=0; i<h.N; i++) {
        Hit hit = h.at(i);
        vh2[hit.slice]->Fill(hit.Z, hit.display4_tick(), hit.adc);
    }
    double max=0;
    for (TH2F* h2 : vh2)
        max = h2->GetMaximum() > max ? h2->GetMaximum() : max;
    for (TH2F* h2 : vh2) {
        h2->SetMinimum(0);
        h2->SetMaximum(max);
    }
    for (unsigned int sl=0; sl<n_slice; sl++) {
        ps[sl % 4]->cd();

        if (up(sl))
            ex_up_pal->Draw();
        else 
            ex_low_pal->Draw();
        
        if (sl == 3) {
            TAxis* za = vh2[sl]->GetZaxis();
            za->SetTitle("ADC");
            za->SetTitleFont(font);
            za->SetTitleSize(font_size);
            za->SetTitleOffset(1.5);

            za->SetLabelFont(font);
            za->SetLabelSize(font_size);
            za->SetLabelOffset(0);

            vh2[sl]->Draw("same pfc colz");
        } else {
            vh2[sl]->Draw("same pfc col");
        }
    }
}

template <typename H>
void draw4TGraph(std::vector<TPad*> ps, TGraph *g, const char *draw_opt, H h) {
    std::vector<TGraph*> vg(n_slice/2);
    for (unsigned int sl=0; sl<n_slice/2; sl++) {
        vg[sl] = new TGraph(*g);
        vg[sl]->SetName(Form("%s_sl%u", g->GetName(), sl));
        vg[sl]->SetEditable(kFALSE);
    }
    // for (unsigned int i=0; i<h.N; i++) {
    for (unsigned int i=0; i<h.N; i++) {
        Hit hit = h.at(i);
        vg[hit.slice%4]->AddPoint(hit.Z, hit.display4_tick());
    }
    for (unsigned int sl=0; sl<n_slice/2; sl++) {
        if (vg[sl]->GetN() == 0) continue;
        ps[sl]->cd();
        vg[sl]->Draw(draw_opt);
    }
}

void draw4TEllipse(std::vector<TPad*> ps, TEllipse *el, Hit h, float rz) {
    TEllipse *e = new TEllipse(*el);
    e->SetX1(h.Z);
    e->SetY1(h.display4_tick());
    e->SetR1(rz);
    e->SetR2(cm ? rz : rz / sampling_rate / drift_velocity);
    ps[h.slice%4]->cd();
    e->Draw();
}