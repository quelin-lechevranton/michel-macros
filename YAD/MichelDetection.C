#include "YAD_tools.h"

const size_t n_file=3;
const vector<string> filelist = yad::ReadFileList(n_file,"ijclab.list");

bool is_outside (vector<double> Xs,vector<double> Ys,vector<double> Zs) {
    bool is_inside=true;
    for (size_t i=0; i<Xs.size(); i++) {
        is_inside = -350 < Xs[i] && Xs[i] < 350 &&
                    -350 < Ys[i] && Ys[i] < 350 &&
                       0 < Zs[i] && Zs[i] < 300;
        if (!is_inside) break;
    }
    return !is_inside;
}

void MichelDetection() {

    clock_t start_time=clock();

    TGraph2D *gMu = new TGraph2D();
    size_t igMu = 0;
    gMu->SetMarkerColor(kOrange+7);
    gMu->SetMarkerStyle(20);
    gMu->SetMarkerSize(0.3);

    // TGraph2D *gMux = new TGraph2D();
    // size_t igMux = 0;
    // gMux->SetMarkerColor(kRed);
    // gMux->SetMarkerStyle(20);
    // gMux->SetMarkerSize(0.3);
    
    TGraph2D *gPrt = new TGraph2D();
    size_t igPrt = 0;
    gPrt->SetMarkerColor(kOrange+7);
    gPrt->SetMarkerStyle(20);
    gPrt->SetMarkerSize(0.3);

    TGraph2D *gEl = new TGraph2D();
    size_t igEl = 0;
    gEl->SetMarkerColor(kAzure+2);
    gEl->SetMarkerStyle(20);
    gEl->SetMarkerSize(0.3);

    TH1D* hNDep = new TH1D("hNDep","Electron Deposit Number;nDep;count",100,0,1000);
    TH2D* hdEdx = new TH2D("hdEdx","Muon Energy Loss;Residual Range (cm);dEdx (MeV/cm)",100,0,200,50,0,5);
    TH1D* hElE = new TH1D("hElE","Electron Spectrum;Total Deposited Energy (MeV);count",100,0,100);

    size_t N_evt=0;
    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());
        yad::Reco R(filename.c_str());

        size_t n_evt = R.GetEntries();
        N_evt+=n_evt;
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {
        // size_t i_evt=45; {

            cout << "Event#" << i_evt+1 << "/" << n_evt << endl;

            T.GetEntry(i_evt);
            R.GetEntry(i_evt);

            for (size_t i_prt=0; i_prt < T.NPrt; i_prt++) {


                if (T.PrtPdg->at(i_prt)!=11) continue; //electrons only
                if (T.PrtNDep->at(i_prt) < 20) continue; //enough electron deposits
                // if (T.PrtNDep->at(i_prt) > 300) continue;

                if (is_outside(T.PrtX->at(i_prt),T.PrtY->at(i_prt),T.PrtZ->at(i_prt))) continue;

                int i_mom = T.PrtMomID->at(i_prt);
                if (i_mom==-1) continue; //no orphelin electrons
                if (T.PrtPdg->at(i_mom)!=13) continue; //electrons coming from muons only


                double detectionRadius=1.; //cm

                size_t nElDep = T.PrtNDep->at(i_prt);
                size_t nMuDep = T.PrtNDep->at(i_mom);


                // vector<size_t> nMuInElDep(nElDep);
                vector<size_t> nElInMuDep(nMuDep);

                //computing the number of electron deposit near every deposit of the muon mother
                for (size_t i_dep=0; i_dep < nMuDep; i_dep++) {

                    double xMu = (*T.DepX)[i_mom][i_dep];
                    double yMu = (*T.DepY)[i_mom][i_dep];
                    double zMu = (*T.DepZ)[i_mom][i_dep];

                    for (size_t j_dep=0; j_dep < nElDep; j_dep++) {
                        
                        double xEl = (*T.DepX)[i_prt][j_dep];
                        double yEl = (*T.DepY)[i_prt][j_dep];
                        double zEl = (*T.DepZ)[i_prt][j_dep];

                        if (TMath::Sqrt(TMath::Power(xEl-xMu,2)+TMath::Power(yEl-yMu,2)+TMath::Power(zEl-zMu,2)) < detectionRadius) {
                            nElInMuDep[i_dep]++;
                        }
                    }
                } //end deposit loop

                //taking the index of the muon deposit with the maximum of electron deposits near itself, considering the deposit as the last muon deposit before the electron emission 
                int i_max = distance(nElInMuDep.begin(), max_element(nElInMuDep.begin(), nElInMuDep.end()));


                //compute total range and distances between points of muon track until electron emission
                double muRange=0;
                vector<double> distances;
                double X0 = (*T.PrtX)[i_mom][0];
                double Y0 = (*T.PrtY)[i_mom][0];
                double Z0 = (*T.PrtZ)[i_mom][0];
                for (size_t i_dep=0; i_dep <= i_max ; i_dep++) {

                    double X = (*T.DepX)[i_mom][i_dep];
                    double Y = (*T.DepY)[i_mom][i_dep];
                    double Z = (*T.DepZ)[i_mom][i_dep];

                    // gMu->SetPoint(igMu++,Y,Z,X);

                    double distance = TMath::Sqrt(TMath::Power(X-X0,2)+TMath::Power(Y-Y0,2)+TMath::Power(Z-Z0,2));
                    distances.push_back(distance);
                    muRange+=distance;
                    X0=X,Y0=Y,Z0=Z;
                }

                //check if muon is dying before emission
                size_t nLast=50;
                double avgLastdEdx=0;
                for (size_t i_dep=i_max-nLast; i_dep <= i_max ; i_dep++) {
                    double E = (*T.DepE)[i_mom][i_dep];
                    double distance=distances[i_dep];

                    avgLastdEdx += E/distance/nLast;
                }
                // if (nElDep < 300) {
                //     cout << avgLastdEdx << endl;
                // } else {
                //     cout << "\t" << avgLastdEdx << endl;
                // }
                if (avgLastdEdx<4) continue; //dying muon only

                //plot dEdx vs. ResidualRange until electron emission
                double muResRange = muRange;
                for (size_t i_dep=1; i_dep <= i_max ; i_dep++) {

                    double E = (*T.DepE)[i_mom][i_dep];

                    double distance=distances[i_dep];
                    muResRange-=distance;
                    hdEdx->Fill(muResRange, E/distance);
                }

                //plot all muon deposits before electron emission
                for (size_t i_dep=0; i_dep <= i_max ; i_dep++) {

                    double X = (*T.DepX)[i_mom][i_dep];
                    double Y = (*T.DepY)[i_mom][i_dep];
                    double Z = (*T.DepZ)[i_mom][i_dep];

                    gMu->SetPoint(igMu++,Y,Z,X);
                    // if (i_dep <= i_max-nLast) {
                    //     gMu->SetPoint(igMu++,Y,Z,X);
                    // } else {
                    //     gMux->SetPoint(igMux++,Y,Z,X);
                    // }
                }

                //plot all prt points
                for (size_t i_ppt=0; i_ppt < T.PrtNPt->at(i_mom); i_ppt++) {
                    double X = (*T.PrtX)[i_mom][i_ppt];
                    double Y = (*T.PrtY)[i_mom][i_ppt];
                    double Z = (*T.PrtZ)[i_mom][i_ppt];

                    gPrt->SetPoint(igPrt++,Y,Z,X);
                } //end particlepoint loop


                //plot all electron deposits + total deposited energy
                double totalE=0;
                hNDep->Fill(nElDep);
                for (size_t i_dep=0; i_dep <= nElDep ; i_dep++) {

                    double X = (*T.DepX)[i_prt][i_dep];
                    double Y = (*T.DepY)[i_prt][i_dep];
                    double Z = (*T.DepZ)[i_prt][i_dep];
                    double E = (*T.DepE)[i_prt][i_dep];

                    gEl->SetPoint(igEl++,Y,Z,X);
                    totalE+=E;
                }
                hElE->Fill(totalE);

                cout << "Michel found with " << totalE << " Mev Deposited over " << nElDep << " points" << endl;

            } //end particle loop
        } //end event loop
    } //end file loop

    TCanvas* c1 = new TCanvas("c1","MichelDetection");
    c1->cd();
    gMu->Draw("p");
    gEl->Draw("samep");
    // gMux->Draw("samep");

    TCanvas* c2 = new TCanvas("c2","MichelDetection");
    c2->Divide(2,2);
    c2->cd(1);
    hNDep->Draw("hist");
    c2->cd(2);
    hdEdx->Draw("colZ");
    c2->cd(3);
    hElE->Draw("hist");

    c1->SaveAs("MichelDetectionGraph2D.root");
    c2->SaveAs("MichelDetectionHist.root");

    cout << N_evt << " events treated in " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}