#include "OmegaDumper_tools.h"

const bool v = true;
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
const Binning bRR = {100,0,300};
const Binning bdEdx = {50,0,5};
const Binning bdQdx = {50,0,800};

const double dQ_min=200; 
const double dx_min=0.5; //cm

const size_t n_bragg_integration=50;
const size_t n_bragg_tail=1;

const size_t n_file=1;
const vector<string> filelist = omega::ReadFileList(n_file,"list/ijclab.list");

void BraggHit() {
    clock_t start_time=clock();

    TH2D* hdQdx = new TH2D(
        "hdQdx",
        ";residual range (cm);dQ/dx (/cm)",
        bRR.n, bRR.min, bRR.max,
        bdQdx.n, bdQdx.min, bdQdx.max
    );

    size_t N_trk=0;

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        omega::Reco R(filename.c_str());

        for (size_t i_evt=0; i_evt < R.N; i_evt++) {

            if (!v) cout << "Event#" << i_evt+1 << "/" << R.N << "\r" << flush;
            if (v) cout << "Event#" << i_evt+1 << "\r" << flush;

            R.GetEvtPfp(i_evt);

            if (v) cout << "\t#trk" << R.Trk.N << "\r" << flush;

            if (R.Trk.N != 1) continue;
            R.GetPfpTrk(0)
            if (R.Trk.Length < 100) continue;

            for (size_t i_pfp=0; i_pfp < R.NPfp; i_pfp++) {

                if (v) cout << "\t\tpfp#" << i_pfp << "\r" << flush;

                if (R.Pfp.isTrk[i_pfp]<0) continue;

                R.GetPfpSpt(i_pfp);
                if (R.Spt.N<n_bragg_integration+n_bragg_tail) {
                    if (v) cout << "\t\t\t\e[91mnot enough spt\e[0m (" << R.Spt.N << ")" << endl;
                    continue;
                }

                if (!omega::IsInside(
                    R.Spt.X,
                    R.Spt.Y,
                    R.Spt.Z,
                    det.Xmin,det.Xmax,
                    det.Ymin,det.Ymax,
                    det.Zmin,det.Zmax
                )) continue;

                double avg_dQdx=0;
                double range=0;
                vector<double> dQs;
                vector<double> dxs;
                vector<double> ranges;
                double *X0=R.Spt.X[0];
                double *Y0=R.Spt.Y[0];
                double *Z0=R.Spt.Z[0];
                for (size_t i_spt=0; i_spt < R.Spt.N; i_spt++) {
                    double dQ = 0;
                    R.GetSptHit(R.Spt.index[i_spt]);
                    for (size_t i_hit=0; i_hit < R.Hit.N; i_hit++) {

                        if (*R.Hit.Plane[i_hit]!=2) continue;

                        dQ += *R.Hit.SumADC[i_hit];
                    } //end hit loop
                    if (dQ < dQ_min) continue;

                    double *X=R.Spt.X[i_spt];
                    double *Y=R.Spt.Y[i_spt];
                    double *Z=R.Spt.Z[i_spt];

                    double dx = omega::Distance(X0,Y0,Z0,X,Y,Z);
                    if (dx < dx_min) continue;
                    range += dx;
                    avg_dQdx+=dQ/dx;
                    X0=X; Y0=Y; Z0=Z;

                    dQs.push_back(dQ);
                    dxs.push_back(dx);
                    ranges.push_back(range);
                } //end spacepoint loop
                avg_dQdx/=R.Spt.N;

                double length = R.TrkLength->at(R.PfpTrkID->at(i_pfp));
                if (TMath::Abs(length-range)/length > .1) {
                    if (v) cout << "\t\t\t\e[91mstrange spt\e[0m" << endl;
                    continue;
                }

                if (v) cout << "\t\t\tavg dQdx: " << avg_dQdx << endl;

                for (size_t i=0; i<dQs.size(); i++) {
                    hdQdx->Fill(range-ranges[i],dQs[i]/dxs[i]);
                }

                N_trk++;
            } //end pfparticle loop
        } //end event loop
    } //end file loop
    cout << endl;

    TCanvas* c1 = new TCanvas("c1","BraggHit");
    c1->cd(1);
    hdQdx->SetMinimum(1);
    c1->SetLogz();
    hdQdx->Draw("colz");

    cout << "nbr of muon: " << N_trk << endl;

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}

void PlotBragg(size_t i_file,size_t i_evt) {
    i_file--; i_evt--;

    clock_t start_time=clock();

    TGraph* gdQ = new TGraph();     gdQ->   SetTitle(";range (cm);dQ");
    TGraph* gdx = new TGraph();     gdx->   SetTitle(";range (cm);dx (cm)");
    TGraph* gdQdx = new TGraph();   gdQdx-> SetTitle(";range (cm);dQ/dx (/cm)");

    string filename = filelist[i_file];
    cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;
    omega::Reco R(filename.c_str());
    cout << "Event#" << i_evt+1 << "\r" << flush;

    R.GetEvent(i_evt);

    if (v) cout << "\t#trk" << R.NTrk << "\r" << flush;

    for (size_t i_pfp=0; i_pfp < R.NPfp; i_pfp++) {

        if (v) cout << "\t\tpfp#" << i_pfp << "\r" << flush;

        if (R.PfpTrkID->at(i_pfp)<0) continue;

        R.GetPfpSpt(i_pfp);
        if (R.Spt.N<n_bragg_integration+n_bragg_tail) {
            if (v) cout << "\t\t\t\e[91mnot enough spt\e[0m (" << R.Spt.N << ")" << endl;
            continue;
        }

        if (!omega::IsInside(
            R.Spt.X,
            R.Spt.Y,
            R.Spt.Z,
            det.Xmin,det.Xmax,
            det.Ymin,det.Ymax,
            det.Zmin,det.Zmax
        )) continue;

        double avg_dQdx=0;
        double range=0;
        double *X0=R.Spt.X[0];
        double *Y0=R.Spt.Y[0];
        double *Z0=R.Spt.Z[0];
        for (size_t i_spt=0; i_spt < R.Spt.N; i_spt++) {
            double dQ = 0;
            R.GetSptHit(R.Spt.index[i_spt]);
            for (size_t i_hit=0; i_hit < R.Hit.N; i_hit++) {

                if (*R.Hit.Plane[i_hit]!=2) continue;

                dQ += *R.Hit.SumADC[i_hit];
            } //end hit loop
            if (dQ < dQ_min) continue;

            double *X=R.Spt.X[i_spt];
            double *Y=R.Spt.Y[i_spt];
            double *Z=R.Spt.Z[i_spt];

            double dx = omega::Distance(X0,Y0,Z0,X,Y,Z);
            if (dx < dx_min) continue;
            range += dx;
            avg_dQdx += dQ/dx;
            X0=X; Y0=Y; Z0=Z;

            gdQ->AddPoint(range,dQ);
            gdx->AddPoint(range,dx);
            gdQdx->AddPoint(range,dQ/dx);
        } //end spacepoint loop
        avg_dQdx/=R.Spt.N;

        double length = R.TrkLength->at(R.PfpTrkID->at(i_pfp));
        if (TMath::Abs(length-range)/length > .1) {
            if (v) cout << "\t\t\t\e[91mstrange spt\e[0m" << endl;
            continue;
        }
        if (v) cout << "\t\t\tavg dQdx: " << avg_dQdx << endl;
    } //end pfparticle loop
    cout << endl;

    cout << "nbr of point: " << gdQ->GetN() << "/" << R.Spt.N << endl;

    TCanvas* c2 = new TCanvas("c2","BraggHit");
    c2->Divide(2,2);
    c2->cd(1);
    gdQ->Draw();
    c2->cd(2);
    gdx->Draw();
    c2->cd(3);
    gdQdx->Draw();

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}