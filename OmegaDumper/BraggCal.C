#include "OmegaLight_tools.h"

const size_t n_file=3;
const vector<string> filelist = omega::ReadFileList(n_file,"list/light.list");

const bool v = true; //verbose
const struct {
    double  Xmin=-320, //cm
            Xmax=350,
            Ymin=-320,
            Ymax=320,
            Zmin=20,
            Zmax=280;
} det;

typedef struct {
    size_t  n,
            min,
            max;
} Binning;
const Binning bRR = {100,0,300}; //cm
const Binning bdEdx = {50,0,5}; //MeV/cm
const Binning bdQdx = {50,0,800};
const Binning bBragg = {200,0,100}; //MeV/cm
const Binning bBraggIntRatio = {20,0,1}; 

const double length_mu_min = 20; //cm

const double dEdx_MIP = 2; //MeV/cm
const double dEdx_min_ratio = 1;
// const double dEdx_min = dEdx_MIP*dEdx_min_ratio;

const size_t n_cal_body_min = 10;

const double bragg_length = 10; //cm
const double bragg_min_ratio_per_int = 3; //MeV/cm
const double bragg_int_ratio_min = 0.5;


void BraggCal(size_t i=0) {
    clock_t start_time=clock();

    // TH2D* hdQdx = new TH2D(
    //     "hdQdx",
    //     ";residual range (cm);dQ/dx (/cm)",
    //     bRR.n,bRR.min,bRR.max,
    //     bdQdx.n,bdQdx.min,bdQdx.max
    // );
    TH2D* hdEdx = new TH2D(
        "hdEdx",
        "before bragg selection;residual range (cm);dE/dx (MeV/cm)",
        bRR.n,bRR.min,bRR.max,
        bdEdx.n,bdEdx.min,bdEdx.max
    );
    TH2D* hdEdx2 = new TH2D(
        "hdEdx",
        "after bragg selection;residual range (cm);dE/dx (MeV/cm)",
        bRR.n,bRR.min,bRR.max,
        bdEdx.n,bdEdx.min,bdEdx.max
    );
    TH1D* hBragg = new TH1D(
        "hBragg",
        ";Bragg dEdx (MeV/cm);#",
        bBragg.n, bBragg.min, bBragg.max
    );
    TH1D* hReverseBragg = new TH1D(
        "hReverseBragg",
        ";Bragg dEdx (MeV/cm);#",
        bBragg.n, bBragg.min, bBragg.max
    );
    TH1D* hBraggIntRatio = new TH1D(
        "hBraggIntRatio",
        ";;#",
        bBraggIntRatio.n,bBraggIntRatio.min,bBraggIntRatio.max
    );

    size_t n_bragg_candidate=0;
    size_t n_bragg=0;

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        omega::Reco R(filename.c_str());

        for (size_t i_evt=0; i_evt < R.N; i_evt++) {
        // for (size_t i_evt=0; i_evt < 20; i_evt++) {

            if (!v) cout << "Event#" << i_evt+1 << "/" << R.N << "\r" << flush;
            if (v) cout << "evt#" << i_evt+1 << "\r" << flush;

            R.GetEvtPfp(i_evt);

            if (R.Trk.N < 1) {
                if (v) cout << "\t\e[91mno track\e[0m " << endl;
                continue;
            }
            if (v) cout << "\tn trk: " << R.Trk.N << endl;

            for (size_t i_pfp=0; i_pfp < R.Pfp.N; i_pfp++) {

                if (!R.Pfp.isTrk[i_pfp]) continue;

                R.GetPfpTrk(i_pfp);

                if (R.Trk.Length < length_mu_min) {
                    if (v) cout << "\t\e[91mtoo short\e[0m (" << R.Trk.Length << ")" << endl;
                    continue;
                }

                if (!omega::IsInside(
                    R.Trk.X,
                    R.Trk.Y,
                    R.Trk.Z,
                    det.Xmin,det.Xmax,
                    det.Ymin,det.Ymax,
                    det.Zmin,det.Zmax
                )) {
                    if (v) cout << " \e[91moutside\e[0m " << endl;
                    continue;
                }

                bool upright = *R.Trk.X[0] > *R.Trk.X.back();
                if (upright) {if (v) cout << "\tuprigth" << endl;}
                else {if (v) cout << "\tupside down" << endl;}

                if (R.Cal.NPt < 1) continue;

                size_t n_cal_head=0;
                while (
                    *R.Cal.ResRange[n_cal_head++] < bragg_length
                    && n_cal_head < R.Cal.NPt
                );
                size_t n_cal_tail=0;
                while (
                    *R.Cal.ResRange[R.Cal.NPt-1-n_cal_tail++] > R.Cal.Range - bragg_length 
                    && n_cal_tail < R.Cal.NPt-1
                );

                size_t n_cal_body = R.Cal.NPt - n_cal_head - n_cal_tail;
                if (v) cout << "\tn cal: " << R.Cal.NPt;
                if (n_cal_body < n_cal_body_min) {
                    if (v) cout << " \e[91mnot enough cal pt\e[0m " << endl;
                    continue;
                }
                if (v) cout << " \u21b4" << endl;

                n_bragg_candidate++;

                double avg_body_dEdx=0;
                for (size_t i_cal=0; i_cal < R.Cal.NPt; i_cal++) {
                    double dEdx = *R.Cal.dEdx[i_cal]; 
                    double resrange;
                    if (upright) {
                        resrange = *R.Cal.ResRange[i_cal];
                    } else {
                        resrange = R.Cal.Range - *R.Cal.ResRange[i_cal];
                    }
                    hdEdx->Fill(resrange,dEdx);

                    if (i_cal >= n_cal_head && i_cal < R.Cal.NPt-n_cal_tail) {
                        avg_body_dEdx += dEdx;
                    }
                } //end calorimetry loop
                avg_body_dEdx /= n_cal_body;
                const double dEdx_min = avg_body_dEdx*dEdx_min_ratio;

                if (v) cout << "\tavg dEdx: " << avg_body_dEdx << endl;

                double bragg_head_int=0;
                size_t n_bragg_head_int=0;
                for (size_t i_cal=0; i_cal < n_cal_head; i_cal++) {
                    double dEdx = *R.Cal.dEdx[i_cal]; 
                    
                    if (dEdx < dEdx_min) continue;
                    bragg_head_int += dEdx;
                    n_bragg_head_int++;
                }
                if (v) cout << "\tBraggHead: " << bragg_head_int << " (" << n_bragg_head_int << "/" << n_cal_head << ")" << endl;

                double normalized_bragg_head = bragg_head_int/avg_body_dEdx;
                if (upright) hBragg->Fill(normalized_bragg_head);
                else hReverseBragg->Fill(normalized_bragg_head);

                // double bragg_min = bragg_min_ratio_per_int*avg_body_dEdx*n_bragg_int;
                // bool is_bragg1 = n_bragg_int > 0 && bragg_int >= bragg_min;
                // double bragg_int_ratio = (double) n_bragg_int/n_cal_head;
                // hBraggIntRatio->Fill(bragg_int_ratio);
                // bool is_bragg1 = bragg_int_ratio > bragg_int_ratio_min;
                // if (!is_bragg1) {if (v) cout << " \e[91mno\e[0m" << endl;}
                // else {if (v) cout << " \e[94myes\e[0m" << endl;}


                double bragg_tail_int=0;
                size_t n_bragg_tail_int=0;
                for (size_t i_cal=R.Cal.NPt-n_cal_tail; i_cal < R.Cal.NPt; i_cal++) {
                    double dEdx = *R.Cal.dEdx[i_cal]; 
                    
                    if (dEdx < dEdx_min) continue;
                    bragg_tail_int += dEdx;
                    n_bragg_tail_int++;
                }
                if (v) cout << "\tBraggTail: " << bragg_tail_int << " (" << n_bragg_tail_int << "/" << n_cal_tail << ")" << endl;

                double normalized_bragg_tail = bragg_tail_int/avg_body_dEdx;
                if (!upright) hBragg->Fill(normalized_bragg_tail);
                else hReverseBragg->Fill(normalized_bragg_tail);

                // bragg_min = bragg_min_ratio_per_int*avg_body_dEdx*n_bragg_int;
                // bool is_bragg2 = n_bragg_int > 0 && bragg_int >= bragg_min;
                // bragg_int_ratio = (double) n_bragg_int/n_cal_tail;
                // hBraggIntRatio->Fill(bragg_int_ratio);
                // bool is_bragg2 = bragg_int_ratio > bragg_int_ratio_min;
                // if (!is_bragg2) {if (v) cout << " \e[91mno\e[0m" << endl;}
                // else {if (v) cout << " \e[94myes\e[0m" << endl;}

                /*
                if (!is_bragg1 && !is_bragg2) continue;
                if (is_bragg1 && is_bragg2) {
                    if (v) cout << "\t\t\e[92mdouble bragg\e[0m" << endl;
                } else if (is_bragg1 == !upright) {
                    if (v) cout << "\t\t\e[92mreverse bragg\e[0m" << endl;
                }

                for (size_t i_cal=0; i_cal < R.Cal.NPt; i_cal++) {
                    double dEdx = *R.Cal.dEdx[i_cal]; 
                    double resrange;
                    if (is_bragg1) {
                        resrange = *R.Cal.ResRange[i_cal];
                    } else {
                        resrange = R.Cal.Range - *R.Cal.ResRange[i_cal];
                    }

                    hdEdx2->Fill(resrange,dEdx);
                } //end calorimetry look

                n_bragg++;
                */
            } //end pfparticle loop
        } //end event loop
    } //end file loop
    cout << endl;

    TCanvas* c1 = new TCanvas("c1","BraggCalo");
    c1->Divide(2,2);
    // c1->Divide(2,1);
    // c1->cd(1);
    // hdQdx->SetMinimum(1);
    // gPad->SetLogz();
    // hdQdx->Draw("colz");
    c1->cd(1);
    hdEdx->SetMinimum(1);
    gPad->SetLogz();
    hdEdx->Draw("colz");
    c1->cd(2);
    hdEdx2->SetMinimum(1);
    gPad->SetLogz();
    hdEdx2->Draw("colz");
    c1->cd(4);
    hBraggIntRatio->Draw("hist");

    TCanvas* c2 = new TCanvas("c2","BraggCal");
    c2->Divide(2,2);
    c2->cd(1);
    hBragg->Draw("hist");
    c2->cd(3);
    hReverseBragg->Draw("hist");

    cout << "nbr of bragg/candidate: " << n_bragg << "/" << n_bragg_candidate << " (" << 100.*n_bragg/n_bragg_candidate << "%)" << endl;

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}

void PlotBragg(size_t i_file=1, size_t i_evt=1, size_t i_pfp=1) {
    i_file--; i_evt--; i_pfp--;

    clock_t start_time=clock();

    TGraph* gdEdx = new TGraph();    gdEdx->SetTitle("muon Energy Loss;residual range (cm);dEdx (MeV/cm)");    

    string filename = filelist[i_file];
    cout << "\e[3mOpening file #" << i_file+1 << "/" << n_file << ": " << filename << "\e[0m" << endl;
    omega::Reco R(filename.c_str());
    // cout << "Event#" << i_evt+1 << "\r" << flush;

    R.GetEvtPfp(i_evt);

    if (!R.Pfp.isTrk[i_pfp]) return;
    R.GetPfpTrk(i_pfp);

    for (size_t i_cal=0; i_cal < R.Cal.NPt; i_cal++) {
        double dEdx = *R.Cal.dEdx[i_cal]; 
        double resrange = *R.Cal.ResRange[i_cal];

        gdEdx->AddPoint(resrange,dEdx);
    } //end calorimetry loop

    TCanvas* c2 = new TCanvas("c2","BraggCal");
    c2->cd(1);
    gdEdx->Draw();

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}