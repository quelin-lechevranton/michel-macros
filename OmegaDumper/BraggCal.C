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
const double n_cal_min = 1; //????

const double dEdx_MIP = 2; //MeV/cm
const double dEdx_min_ratio = 1;
// const double dEdx_min = dEdx_MIP*dEdx_min_ratio;

const double bragg_length = 15; //cm
const double bragg_min_ratio_per_int = 3; //MeV/cm
const double bragg_int_ratio_min = 0.2;


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
            if (v) cout << "\tn trk: " << R.Trk.N << "\u21b4" << endl;
            // if (R.NTrk != 1) {
            //     if (v) cout << " \e[91mmore than one track\e[0m " << endl;
            //     continue;
            // }

            for (size_t i_pfp=0; i_pfp < R.Pfp.N; i_pfp++) {

                if (!R.Pfp.isTrk[i_pfp]) continue;

                R.GetPfpTrk(i_pfp);

                if (R.Trk.Length < length_mu_min) {
                    if (v) cout << "\t\e[91mtoo short\e[0m (" << R.Trk.Length << ")" << endl;
                    continue;
                }

                if (v) cout << "\tn cal: " << R.Cal.NPt;
                if (R.Cal.NPt < n_cal_min) {
                    if (v) cout << " \e[91mnot enough cal pt\e[0m " << endl;
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
                if (v) cout << " \u21b4" << endl;

                n_bragg_candidate++;

                bool upside_down = *R.Trk.X.back() > *R.Trk.X[0];
                if (upside_down) {if (v) cout << "\tupside_down" << endl;}
                else {if (v) cout << "\tuprigth" << endl;}

                // double avg_dQdx=0;
                double avg_dEdx=0;
                for (size_t i_cal=0; i_cal < R.Cal.NPt; i_cal++) {
                    // double dQdx = *R.Cal.dQdx[i_cal]; 
                    double dEdx = *R.Cal.dEdx[i_cal]; 
                    double resrange;
                    if (!upside_down) {
                        resrange = *R.Cal.ResRange[i_cal];
                    } else {
                        resrange = R.Cal.Range - *R.Cal.ResRange[i_cal];
                    }
                    // avg_dQdx += dQdx;
                    avg_dEdx += dEdx;

                    // hdQdx->Fill(resrange,dQdx);
                    hdEdx->Fill(resrange,dEdx);
                } //end calorimetry loop
                // avg_dQdx /= R.Cal.NPt;
                avg_dEdx /= R.Cal.NPt;
                const double dEdx_min = avg_dEdx*dEdx_min_ratio;

                // if (v) cout << "\tavg dQdx: " << avg_dQdx << endl;
                // if (v) cout << "\tavg dEdx: " << avg_dEdx << endl;


                double bragg_int=0;
                size_t n_bragg_int=0;
                size_t n_cal_bragg=0;
                while (*R.Cal.ResRange[n_cal_bragg++] < bragg_length && n_cal_bragg < R.Cal.NPt);

                for (size_t i_cal=0; i_cal < n_cal_bragg; i_cal++) {
                    double dEdx = *R.Cal.dEdx[i_cal]; 
                    
                    if (dEdx < dEdx_min) continue;
                    bragg_int += dEdx;
                    n_bragg_int++;
                }
                hBragg->Fill(bragg_int*dEdx_MIP/avg_dEdx);
                double rrmin = *R.Cal.ResRange[0];
                double rrmax = *R.Cal.ResRange[n_cal_bragg-1];
                if (v) cout << "\tBragg?: " << bragg_int << " (" << n_bragg_int << "/" << n_cal_bragg << " rr:" << rrmin << "-" << rrmax << ")";

                // double bragg_min = bragg_min_ratio_per_int*avg_dEdx*n_bragg_int;
                // bool is_bragg1 = n_bragg_int > 0 && bragg_int >= bragg_min;
                double bragg_int_ratio = (double) n_bragg_int/n_bragg;
                hBraggIntRatio->Fill(bragg_int_ratio);
                bool is_bragg1 = bragg_int_ratio > bragg_int_ratio_min;
                if (!is_bragg1) {if (v) cout << " \e[91mno\e[0m" << endl;}
                else {if (v) cout << " \e[94myes\e[0m" << endl;}

                bragg_int=0;
                n_bragg_int=0;
                n_cal_bragg=0;
                while (*R.Cal.ResRange[R.Cal.NPt-1-n_cal_bragg++] > R.Cal.Range - bragg_length && n_cal_bragg < R.Cal.NPt-1);

                for (size_t i_cal=R.Cal.NPt-n_cal_bragg; i_cal < R.Cal.NPt; i_cal++) {
                    double dEdx = *R.Cal.dEdx[i_cal]; 
                    
                    if (dEdx < dEdx_min) continue;
                    bragg_int += dEdx;
                    n_bragg_int++;
                }
                hBragg->Fill(bragg_int*dEdx_MIP/avg_dEdx);
                rrmin = *R.Cal.ResRange[R.Cal.NPt-n_cal_bragg];
                rrmax = *R.Cal.ResRange[R.Cal.NPt-1];
                if (v) cout << "\tBragg?: " << bragg_int << " (" << n_bragg_int << "/" << n_cal_bragg << " rr:" << rrmin << "-" << rrmax << ")";

                // bragg_min = bragg_min_ratio_per_int*avg_dEdx*n_bragg_int;
                // bool is_bragg2 = n_bragg_int > 0 && bragg_int >= bragg_min;
                bragg_int_ratio = (double) n_bragg_int/n_bragg;
                hBraggIntRatio->Fill(bragg_int_ratio);
                bool is_bragg2 = bragg_int_ratio > bragg_int_ratio_min;
                if (!is_bragg2) {if (v) cout << " \e[91mno\e[0m" << endl;}
                else {if (v) cout << " \e[94myes\e[0m" << endl;}

                if (!is_bragg1 && !is_bragg2) continue;
                if (is_bragg1 && is_bragg2) {
                    if (v) cout << "\t\t\e[92mdouble bragg\e[0m" << endl;
                } else if (is_bragg1 == upside_down) {
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
    c1->cd(3);
    hBragg->Draw("hist");
    c1->cd(4);
    hBraggIntRatio->Draw("hist");

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