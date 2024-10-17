#include "OmegaLight_tools.h"

const size_t n_file=1;
const vector<string> filelist = omega::ReadFileList(n_file,"list/pdvd_cosmics.list");

const bool v = true; //verbose

const omega::Limits det = omega::fiducial;

const omega::Binning bRR = {100,0,300}; //cm
const omega::Binning bdEdx = {50,0,5}; //MeV/cm
const omega::Binning bdQdx = {50,0,800};
const omega::Binning bBragg = {62,-1,30}; //#avg dEdx
const omega::Binning bBraggIntRatio = {20,0,1}; 

const double length_mu_min = 20; //cm

// const double dEdx_MIP = 2; //MeV/cm
const double dEdx_min_ratio = 1.5;
// const double dEdx_min = dEdx_MIP*dEdx_min_ratio;

const size_t n_cal_body_min = 20;

const double bragg_length = 2; //cm
// const double bragg_min_ratio = 10; //cm/MeV?

// const double bragg_int_ratio_min = 0.5;

// const vector<double> bragg_razor = {7,5.5,7,5.5};
const vector<double> bragg_razor = {10,10,10,10};


void BraggCal_cosmics(size_t i=0) {
    clock_t start_time=clock();


    TH2D* hdEdx = new TH2D(
        "hdEdx",
        "before bragg selection;residual range (cm);dE/dx (MeV/cm)",
        bRR.n,bRR.min,bRR.max,
        bdEdx.n,bdEdx.min,bdEdx.max
    );
    TH2D* hdEdx2 = new TH2D(
        "hdEdx2",
        "after bragg selection;residual range (cm);dE/dx (MeV/cm)",
        bRR.n,bRR.min,bRR.max,
        bdEdx.n,bdEdx.min,bdEdx.max
    );
    
    vector<size_t> nBragg={3,4,2};
    vector<vector<vector<TH1D*>>> hBragg;
    vector<string> line = {"Bragg ","Reverse ","Outside "};
    vector<string> col = {"Normalized","Treshold Normalized","Non-normalized","Treshold Non-normalized"};
    vector<int> color={kRed,kBlue};
    for (int i=0; i<nBragg[0]; i++) {
        vector<vector<TH1D*>> tpBragg;
        for (int j=0; j<nBragg[1]; j++) {
            vector<TH1D*> tptpBragg;
            for (int k=0; k<nBragg[2]; k++) {
                stringstream t; t << "h" << i << j <<k;
                stringstream n; 
                n << line[i] << col[j] << ";#AvgdEdx;#";
                TH1D* h = new TH1D(
                    t.str().c_str(), n.str().c_str(),
                    bBragg.n,bBragg.min,bBragg.max
                );
                h->SetLineColor(color[k]);
                tptpBragg.push_back(h);
            }
            tpBragg.push_back(tptpBragg);
        }
        hBragg.push_back(tpBragg);
    }

    vector<size_t> true_positive;
    vector<size_t> false_positive;
    vector<size_t> true_negative;
    vector<size_t> false_negative;

    for (int j=0; j<nBragg[1]; j++) {
        true_positive.push_back(0);
        false_positive.push_back(0);
        true_negative.push_back(0);
        false_negative.push_back(0);
    }


    size_t n_bragg_candidate=0;
    size_t n_bragg=0;

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << i_file+1 << "/" << n_file << ": " << filename << "\e[0m" << endl;

        omega::Reco R(filename.c_str());


        for (size_t i_evt=0; i_evt < R.N; i_evt++) {
        // for (size_t i_evt=0; i_evt < 20; i_evt++) {

            if (!v) cout << "Event#" << i_evt+1 << "/" << R.N << "\r" << flush;
            if (v) cout << "evt#" << i_evt+1 << "\r" << flush;

            R.GetEvtPfp(i_evt);

            if (R.Trk.N != 1) {
                if (v) cout << "\t\e[91mno track\e[0m " << endl;
                continue;
            }
            if (v) cout << "\tn trk: " << R.Trk.N << endl;

            for (size_t i_pfp=0; i_pfp < R.Pfp.N; i_pfp++) {

                if (!R.Pfp.isTrk[i_pfp]) continue;

                R.GetPfpTrk(i_pfp);

                // if (R.Trk.Length < length_mu_min) {
                //     if (v) cout << "\t\e[91mtoo short\e[0m (" << R.Trk.Length << ")" << endl;
                //     continue;
                // }

                if (R.Cal.NPt < 1) continue;

                bool inside = omega::IsInside(
                    R.Trk.X,
                    R.Trk.Y,
                    R.Trk.Z,
                    det
                );
                // if (!inside) {
                //     if (v) cout << " \e[91moutside\e[0m " << endl;
                //     continue;
                // }
                // if (inside) n_bragg_candidate++;

                bool upright = *R.Trk.X[0] > *R.Trk.X.back();
                if (upright) {if (v) cout << "\tuprigth" << endl;}
                else {if (v) cout << "\tupside down" << endl;}

                size_t n_cal_tail=0;
                while (
                    *R.Cal.ResRange[n_cal_tail++] < bragg_length
                    && n_cal_tail < R.Cal.NPt
                );
                size_t n_cal_head=0;
                while (
                    *R.Cal.ResRange[R.Cal.NPt-1-n_cal_head++] > R.Cal.Range - bragg_length 
                    && n_cal_head < R.Cal.NPt-1
                );

                size_t n_cal_body = R.Cal.NPt - n_cal_tail - n_cal_head;
                if (v) cout << "\tn cal: " << R.Cal.NPt;
                if (n_cal_body < n_cal_body_min) {
                    if (v) cout << " \e[91mnot enough cal pt\e[0m " << endl;
                    continue;
                }
                if (v) cout << " \u21b4" << endl;


                double avg_body_dEdx=0;
                for (size_t i_cal=0; i_cal < R.Cal.NPt; i_cal++) {
                    double dEdx = *R.Cal.dEdx[i_cal]; 
                    double resrange;
                    if (upright) {
                        resrange = *R.Cal.ResRange[i_cal];
                    } else {
                        resrange = R.Cal.Range - *R.Cal.ResRange[i_cal];
                    }
                    // if (inside) hdEdx->Fill(resrange,dEdx);
                    hdEdx->Fill(resrange,dEdx);

                    if (i_cal >= n_cal_tail && i_cal < R.Cal.NPt-n_cal_head) {
                        avg_body_dEdx += dEdx;
                    }
                } //end calorimetry loop
                avg_body_dEdx /= n_cal_body;
                const double dEdx_min = avg_body_dEdx*dEdx_min_ratio;
                // const double dEdx_min = 2*dEdx_min_ratio;

                if (v) cout << "\tavg body dEdx: " << avg_body_dEdx << " (" << n_cal_body << ")" << endl;

                if(avg_body_dEdx>4 || avg_body_dEdx<1) continue;

                double bragg_tail_int=0;
                double bragg_tail_treshold_int=0;
                size_t n_bragg_tail_int=0;
                for (size_t i_cal=0; i_cal < n_cal_tail; i_cal++) {
                    double dEdx = *R.Cal.dEdx[i_cal]; 
                    
                    bragg_tail_int += dEdx;
                    if (dEdx < dEdx_min) continue;
                    bragg_tail_treshold_int += dEdx;
                    n_bragg_tail_int++;
                }
                if (v) cout << "\tBraggHead: " << bragg_tail_int << " (" << n_bragg_tail_int << "/" << n_cal_tail << ")" << endl;

                vector<double> bragg_tail = {
                    bragg_tail_int/avg_body_dEdx,
                    bragg_tail_treshold_int/avg_body_dEdx,
                    bragg_tail_int/2,
                    bragg_tail_treshold_int/2
                };
                // vector<double> bragg_tail = {
                //     bragg_tail_int/avg_body_dEdx/n_cal_tail,
                //     bragg_tail_treshold_int/avg_body_dEdx/n_bragg_tail_int,
                //     bragg_tail_int/2/n_cal_tail,
                //     bragg_tail_treshold_int/2/n_bragg_tail_int
                // };

                // double bragg_min = bragg_min_ratio*avg_body_dEdx;
                // bool is_bragg1 = n_bragg_int > 0 && bragg_int >= bragg_min;
                // double bragg_int_ratio = (double) n_bragg_int/n_cal_tail;
                // hBraggIntRatio->Fill(bragg_int_ratio);
                // bool is_bragg1 = bragg_int_ratio > bragg_int_ratio_min;
                // if (!is_bragg1) {if (v) cout << " \e[91mno\e[0m" << endl;}
                // else {if (v) cout << " \e[94myes\e[0m" << endl;}


                double bragg_head_int=0;
                double bragg_head_treshold_int=0;
                size_t n_bragg_head_int=0;
                for (size_t i_cal=R.Cal.NPt-n_cal_head; i_cal < R.Cal.NPt; i_cal++) {
                    double dEdx = *R.Cal.dEdx[i_cal]; 
                    
                    bragg_head_int += dEdx;
                    if (dEdx < dEdx_min) continue;
                    bragg_head_treshold_int += dEdx;
                    n_bragg_head_int++;
                }
                // if (v) cout << "\tBraggTail: " << bragg_head_int << " (" << n_bragg_head_int << "/" << n_cal_head << ")" << endl;

                vector<double> bragg_head = {
                    bragg_head_int/avg_body_dEdx,
                    bragg_head_treshold_int/avg_body_dEdx,
                    bragg_head_int/2,
                    bragg_head_treshold_int/2
                };
                // vector<double> bragg_head = {
                //     bragg_head_int/avg_body_dEdx/n_cal_head,
                //     bragg_head_treshold_int/avg_body_dEdx/n_bragg_head_int,
                //     bragg_head_int/2/n_cal_head,
                //     bragg_head_treshold_int/2/n_bragg_head_int
                // };

                vector<double> bragg, bragg_reverse;
                if (upright) {
                    bragg = bragg_tail;
                    bragg_reverse = bragg_head;
                } else {
                    bragg = bragg_head;
                    bragg_reverse = bragg_tail;
                }

                size_t tru = 0;
                // if (inside) {
                    for (int j=0; j<nBragg[1]; j++) {
                        hBragg[0][j][tru]->Fill(bragg[j]);
                        hBragg[1][j][tru]->Fill(bragg_reverse[j]);
                    }
                // }
                // } else {
                    // for (int j=0; j<nBragg[1]; j++) {
                        // hBragg[2][j][tru]->Fill(bragg[j]);
                    // }
                // }

                // if (!inside) continue;
                n_bragg_candidate++;

                vector<bool> is_bragg;
                for (int j=0; j<nBragg[1]; j++) {
                    is_bragg.push_back (bragg[j] > bragg_razor[j]);
                    if (bragg[j] > bragg_razor[j]) {
                        if (tru) true_positive[j]++;
                        else     false_positive[j]++;
                    } else {
                        if (tru) false_negative[j]++;
                        else     true_negative[j]++;
                    }
                }
                
                if (!is_bragg[1]) continue;

                for (size_t i_cal=0; i_cal < R.Cal.NPt; i_cal++) {
                    double dEdx = *R.Cal.dEdx[i_cal]; 
                    double resrange;
                    if (upright) {
                        resrange = *R.Cal.ResRange[i_cal];
                    } else {
                        resrange = R.Cal.Range - *R.Cal.ResRange[i_cal];
                    }

                    hdEdx2->Fill(resrange,dEdx);
                } //end calorimetry look

                n_bragg++;
            } //end pfparticle loop
        } //end event loop
        i_file++;
    } //end file loop
    cout << endl;

    TCanvas* c1 = new TCanvas("c1","BraggCal_cosmics");
    // c1->Divide(2,2);
    c1->Divide(2,1);
    // // c1->cd(1);
    // // hdQdx->SetMinimum(1);
    // // gPad->SetLogz();
    // // hdQdx->Draw("colz");
    c1->cd(1);
    hdEdx->SetMinimum(1);
    gPad->SetLogz();
    hdEdx->Draw("colz");
    c1->cd(2);
    hdEdx2->SetMinimum(1);
    gPad->SetLogz();
    hdEdx2->Draw("colz");
    // c1->cd(4);
    // hBraggIntRatio->Draw("hist");

    vector<TLine*> razor(nBragg[1]);
    for (int j=0; j<nBragg[1]; j++) {
        razor[j]= new TLine(bragg_razor[j],0,bragg_razor[j],hBragg[0][j][1]->GetMaximum());
        razor[j]->SetLineColor(kViolet);
        razor[j]->SetLineWidth(2);
    }

    TCanvas* c2 = new TCanvas("c2","BraggCal_cosmics");
    c2->Divide(2,2);
    int k=0;
    for (int j=0; j<nBragg[1]; j++) {
        c2->cd(++k);
        gPad->SetLogy();
        hBragg[0][j][0]->Draw("hist");
        hBragg[0][j][1]->Draw("samehist");
        razor[j]->Draw();
    }
    TCanvas* c3 = new TCanvas("c3","BraggCal_cosmics");
    c3->Divide(2,2);
    k=0;
    for (int j=0; j<nBragg[1]; j++) {
        c3->cd(++k);
        gPad->SetLogy();
        hBragg[1][j][0]->Draw("hist");
        hBragg[1][j][1]->Draw("samehist");
        razor[j]->Draw();
    }


    cout << "criterium results: " << endl;
    for (int j=0; j<nBragg[1]; j++) {
        cout << "\t" << col[j] << ": " << endl;
        cout << "\t\ttrue_positive: " << true_positive[j] << " (" << 100.*true_positive[j]/n_bragg_candidate << "%)" << endl;
        cout << "\t\tfalse_positive: " << false_positive[j] << " (" << 100.*false_positive[j]/n_bragg_candidate << "%)" << endl;
        cout << "\t\ttrue_negative: " << true_negative[j] << " (" << 100.*true_negative[j]/n_bragg_candidate << "%)" << endl;
        cout << "\t\tfalse_negative: " << false_negative[j] << " (" << 100.*false_negative[j]/n_bragg_candidate << "%)" << endl;
        cout << "\t\t----------------------------------------" << endl;
        cout << "\t\tefficiency (tp/tp+fn): " << 100.*true_positive[j]/(true_positive[j]+false_negative[j]) << "%" << endl;
        cout << "\t\tpurity (tp/tp+fp): " << 100.*true_positive[j]/(true_positive[j]+false_positive[j]) << "%" << endl;
    }

    c1->SaveAs("out/BraggCal_cosmics_out1.pdf");
    c2->SaveAs("out/BraggCal_cosmics_out2.pdf");
    c3->SaveAs("out/BraggCal_cosmics_out3.pdf");


    // cout << "nbr of bragg/candidate: " << n_bragg << "/" << n_bragg_candidate << " (" << 100.*n_bragg/n_bragg_candidate << "%)" << endl;

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

    TCanvas* c0 = new TCanvas("c0","BraggCal_cosmics");
    c0->cd(1);
    gdEdx->Draw();

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}
