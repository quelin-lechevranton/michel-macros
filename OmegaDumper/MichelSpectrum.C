//#include "OmegaDumper_tools.h"
#include "OmegaLight_tools.h"

const size_t n_file=37;
const vector<string> filelist = omega::ReadFileList(n_file,"list/muminus.list");

const bool v = true;

const omega::Limits det = omega::fiducial;
const omega::Binning bRR = {100,0,300}; //cm
const omega::Binning bdEdx = {50,0,5}; //MeV/cm
const omega::Binning bE = {100,0,1000};

const double length_mu_min = 20; //cm
const double n_cal_min = 1; //????

const double dEdx_MIP = 2; //MeV/cm
const double dEdx_min_ratio = 1;
// const double dEdx_min = dEdx_MIP*dEdx_min_ratio;

const double bragg_length = 15; //cm
const double bragg_min = 20; //MeV/cm

const double length_el_max = 20; //cm
const double end_radius = 30; //cm

void MichelSpectrum(size_t i=0) {
    clock_t start_time=clock();

    TH2D* hdEdx = new TH2D(
        "hdEdx",
        "after bragg selection;residual range (cm);dE/dx (MeV/cm)",
        bRR.n,bRR.min,bRR.max,
        bdEdx.n,bdEdx.min,bdEdx.max
    );
    TGraph2D* gEnd = new TGraph2D();
    gEnd->SetMarkerStyle(20);
    gEnd->SetMarkerColor(kOrange);

    TGraph2D* gEl = new TGraph2D();
    gEl->SetMarkerStyle(20);
    gEl->SetMarkerSize(0.5);
    gEl->SetMarkerColor(kAzure);

    TH1D* hE = new TH1D("hE",";SumADC;#",bE.n,bE.min,bE.max);

    size_t nc=4;
    vector<size_t> c;
    while (nc--) c.push_back(0);

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        omega::Reco R(filename.c_str());

        for (size_t i_evt=0; i_evt < R.N; i_evt++) {

            if (v) cout << "Event#" << i_evt+1 << "/" << R.N << "\r" << flush;

            R.GetEvtPfp(i_evt);

            if (R.Trk.N < 1) continue;

            struct { double *X,*Y,*Z; } End;
            for (size_t i_pfp=0; i_pfp < R.Pfp.N; i_pfp++) {

                if (!R.Pfp.isTrk[i_pfp]) continue;
                R.GetPfpTrk(i_pfp);
                if (R.Trk.Length < length_mu_min) continue;
                if (R.Cal.NPt < n_cal_min) continue;
                
                if (!omega::IsInside(
                    R.Trk.X,
                    R.Trk.Y,
                    R.Trk.Z,
                    det
                )) continue;

                double avg_dEdx=0;
                for (size_t i_cal=0; i_cal < R.Cal.NPt; i_cal++) {
                    avg_dEdx += *R.Cal.dEdx[i_cal];
                } //end calorimetry loop
                avg_dEdx /= R.Cal.NPt;
                const double dEdx_min = avg_dEdx*dEdx_min_ratio;

                bool upside_down = *R.Trk.X.back() > *R.Trk.X[0];

                double bragg_int=0;
                size_t n_cal_bragg=0;
                if (!upside_down) {
                    while (
                        *R.Cal.ResRange[n_cal_bragg++] < bragg_length 
                        && n_cal_bragg < R.Cal.NPt
                    );

                    for (size_t i_cal=0; i_cal < n_cal_bragg; i_cal++) {
                        double dEdx = *R.Cal.dEdx[i_cal]; 
                        if (dEdx < dEdx_min) continue;
                        bragg_int += dEdx;
                    }
                } else {
                    while (
                        *R.Cal.ResRange[R.Cal.NPt-1-n_cal_bragg++] > R.Cal.Range - bragg_length
                        && n_cal_bragg < R.Cal.NPt-1
                    );

                    for (size_t i_cal=R.Cal.NPt-n_cal_bragg; i_cal < R.Cal.NPt; i_cal++) {
                        double dEdx = *R.Cal.dEdx[i_cal]; 
                        if (dEdx < dEdx_min) continue;
                        bragg_int += dEdx;
                    }
                }
                if (bragg_int < bragg_min) continue;
                c[0]++;

                // for (size_t i_cal=0; i_cal < R.Cal.NPt; i_cal++) {
                //     double dEdx = *R.Cal.dEdx[i_cal]; 
                //     double resrange;
                //     if (!upside_down) resrange = *R.Cal.ResRange[i_cal];
                //     else resrange = R.Cal.Range - *R.Cal.ResRange[i_cal];
                //     else continue;

                //     hdEdx->Fill(resrange,dEdx);
                // } //end calorimetry loop

                if (!upside_down) {
                    // gStart->AddPoint(*R.Trk.Y[0],*R.Trk.Z[0],*R.Trk.X[0]);
                    End.X=R.Trk.X.back();
                    End.Y=R.Trk.Y.back();
                    End.Z=R.Trk.Z.back();
                    c[1]++;
                } else {
                    // gStart->AddPoint(*R.Trk.Y.back(),*R.Trk.Z.back(),*R.Trk.X.back());
                    End.X=R.Trk.X[0];
                    End.Y=R.Trk.Y[0];
                    End.Z=R.Trk.Z[0];
                    c[2]++;
                }
                gEnd->AddPoint(*End.Y,*End.Z,*End.X);
                
                double E=0;
                for (size_t j_pfp=0; j_pfp < R.Pfp.N; j_pfp++) {

                    if (R.Pfp.isTrk[j_pfp]) {
                        R.GetPfpTrk(j_pfp);
                        if (R.Trk.Length > length_el_max) continue;
                    }
                    R.GetPfpSpt(j_pfp);

                    for (size_t j_spt=0; j_spt < R.Spt.N; j_spt++) {
                        if (omega::Distance(
                            End.X,End.Y,End.Z,
                            R.Spt.X[j_spt],R.Spt.Y[j_spt],R.Spt.Z[j_spt]
                        ) > end_radius) continue;

                        gEl->AddPoint(*R.Spt.Y[j_spt],*R.Spt.Z[j_spt],*R.Spt.X[j_spt]);
                        R.GetSptHit(j_spt);

                        for (size_t j_hit=0; j_hit < R.Hit.N; j_hit++) {
                            if (*R.Hit.Plane[j_hit] != 2) continue;
                            E += *R.Hit.SumADC[j_hit];
                        }
                    }
                } //end michel spectrum pfparticle loop
                if (E>1) hE->Fill(E);
		        else c[3]++;
            } //end muon selection pfparticle loop
        } //end event loop
    } //end file loop
    cout << endl;

    for (size_t a : c) cout << a << endl;

    TCanvas* c1 = new TCanvas("c1","MichelSpectrum");
    c1->Divide(2,1);
    // c1->cd(1);
    // hdEdx->SetMinimum(1);
    // gPad->SetLogz();
    // hdEdx->Draw("colz");
    c1->cd(1);
    gEnd->Draw("p");
    gEl->Draw("samep");
    c1->cd(2);
    hE->Draw("hist");

    c1->SaveAs("MichelSpectrum.root");

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}
