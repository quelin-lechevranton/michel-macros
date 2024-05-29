#include "YAD_tools.h"

// analysis parameter
double coincidence_radius_gros=1.; //cm
size_t n_step_gros=10; //distance per step : 0.03 cm
double coincidence_radius_fin=1.; //cm
size_t n_scan_fin=4*n_step_gros;
size_t n_last_deposits=50;
double threshold_dEdx=4.;


const size_t n_file=3;
const vector<string> filelist = yad::readFileList(n_file,"list/ijclab.list");



void TrueMichel_v2() {

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
    
    // TGraph2D *gPrt = new TGraph2D();
    // size_t igPrt = 0;
    // gPrt->SetMarkerColor(kPink+7);
    // gPrt->SetMarkerStyle(20);
    // gPrt->SetMarkerSize(0.3);

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

            cout << "Event#" << i_evt+1 << "/" << n_evt << "\r" << flush;

            T.GetEntry(i_evt);
            R.GetEntry(i_evt);

            for (size_t i_prt=0; i_prt < T.NPrt; i_prt++) {


                if (T.PrtPdg->at(i_prt)!=11) continue; //electrons only
                if (T.PrtNDep->at(i_prt) < 20) continue; //enough electron deposits
                // if (T.PrtNDep->at(i_prt) > 300) continue;

                if (!yad::isInside(T.DepX->at(i_prt),T.DepY->at(i_prt),T.DepZ->at(i_prt))) continue;

                int i_mom = T.PrtMomID->at(i_prt);
                if (i_mom==-1) continue; //no orphelin electrons
                if (T.PrtPdg->at(i_mom)!=13) continue; //electrons coming from muons only

                size_t n_el_dep = T.PrtNDep->at(i_prt);
                size_t n_mu_dep = T.PrtNDep->at(i_mom);

                //now we want to find the closest muon deposit to the electron track
                //first step is rough-and-ready second is more precise our the result of first step
                
                size_t n_max_gros=0, i_max_gros=0;

                //computing the number of electron deposit near every deposit of the muon mother
                for (size_t i_mu_dep=0; i_mu_dep < n_mu_dep; i_mu_dep+=n_step_gros) {

                    double x_mu = (*T.DepX)[i_mom][i_mu_dep];
                    double y_mu = (*T.DepY)[i_mom][i_mu_dep];
                    double z_mu = (*T.DepZ)[i_mom][i_mu_dep];

                    size_t n_coincidence=0;
                    for (size_t i_el_dep=0; i_el_dep < n_el_dep; i_el_dep+=n_step_gros) {
                        
                        double x_el = (*T.DepX)[i_prt][i_el_dep];
                        double y_el = (*T.DepY)[i_prt][i_el_dep];
                        double z_el = (*T.DepZ)[i_prt][i_el_dep];

                        double dist= yad::distance(x_el,y_el,z_el,x_mu,y_mu,z_mu);
                        if (dist < coincidence_radius_gros) {
                            n_coincidence++;
                        }
                    }

                    if ( n_coincidence > n_max_gros ) {
                        n_max_gros = n_coincidence;
                        i_max_gros = i_mu_dep;
                    }
                } //end deposit loop

                // more precise count
                size_t n_max=0, i_max=0;
                for (
                    size_t i_mu_dep = i_max_gros-n_scan_fin > 0 ? i_max_gros-n_scan_fin : 0;
                    i_mu_dep < (i_max_gros+n_scan_fin > n_mu_dep ? n_mu_dep : i_max_gros+n_scan_fin) ; i_mu_dep++) {
                    

                    double x_mu = (*T.DepX)[i_mom][i_mu_dep];
                    double y_mu = (*T.DepY)[i_mom][i_mu_dep];
                    double z_mu = (*T.DepZ)[i_mom][i_mu_dep];

                    size_t n_coincidence=0;
                    for (size_t i_el_dep=0; i_el_dep < n_el_dep; i_el_dep++) {
                        
                        double x_el = (*T.DepX)[i_prt][i_el_dep];
                        double y_el = (*T.DepY)[i_prt][i_el_dep];
                        double z_el = (*T.DepZ)[i_prt][i_el_dep];

                        double dist= yad::distance(x_el,y_el,z_el,x_mu,y_mu,z_mu);
                        if (dist < coincidence_radius_fin) {
                            n_coincidence++;
                        }
                    }
                    if ( n_coincidence > n_max ) {
                        n_max = n_coincidence;
                        i_max = i_mu_dep;
                    }
                } //end deposit loop
                // cout << "evt#" << i_evt << ":prt#" << i_prt << ":i_max=" << i_max << endl;

                //compute total range and distances between points of muon track until electron emission
                double mu_range=0;
                vector<double> distances;
                double X0 = (*T.DepX)[i_mom][0];
                double Y0 = (*T.DepY)[i_mom][0];
                double Z0 = (*T.DepZ)[i_mom][0];
                for (size_t i_dep=0; i_dep <= i_max ; i_dep++) {

                    double X = (*T.DepX)[i_mom][i_dep];
                    double Y = (*T.DepY)[i_mom][i_dep];
                    double Z = (*T.DepZ)[i_mom][i_dep];

                    // gMu->SetPoint(igMu++,Y,Z,X);

                    double dist = yad::distance(X,Y,Z,X0,Y0,Z0);
                    distances.push_back(dist);
                    mu_range+=dist;
                    X0=X,Y0=Y,Z0=Z;
                }

                //check if muon is dying before emission
                double avg_last_dEdx=0;
                for (size_t i_dep=i_max-n_last_deposits; i_dep <= i_max ; i_dep++) {
                    double E = (*T.DepE)[i_mom][i_dep];
                    double dist=distances[i_dep];

                    avg_last_dEdx += E/dist/n_last_deposits;
                }
                // if (n_el_dep < 300) {
                //     cout << avg_last_dEdx << endl;
                // } else {
                //     cout << "\t" << avg_last_dEdx << endl;
                // }
                if (avg_last_dEdx<threshold_dEdx) continue; 

                //plot dEdx vs. ResidualRange until electron emission
                double mu_res_range = mu_range;
                for (size_t i_dep=1; i_dep <= i_max ; i_dep++) {

                    double E = (*T.DepE)[i_mom][i_dep];

                    double dist=distances[i_dep];
                    mu_res_range-=dist;
                    hdEdx->Fill(mu_res_range, E/dist);
                }

                //plot all muon deposits before electron emission
                for (size_t i_dep=0; i_dep <= i_max ; i_dep++) {

                    double X = (*T.DepX)[i_mom][i_dep];
                    double Y = (*T.DepY)[i_mom][i_dep];
                    double Z = (*T.DepZ)[i_mom][i_dep];

                    gMu->SetPoint(igMu++,Y,Z,X);
                    // if (i_dep <= i_max-n_last_deposits) {
                    //     gMu->SetPoint(igMu++,Y,Z,X);
                    // } else {
                    //     gMux->SetPoint(igMux++,Y,Z,X);
                    // }
                }

                //plot all prt points
                // for (size_t i_ppt=0; i_ppt < T.PrtNPt->at(i_mom); i_ppt++) {
                //     double X = (*T.PrtX)[i_mom][i_ppt];
                //     double Y = (*T.PrtY)[i_mom][i_ppt];
                //     double Z = (*T.PrtZ)[i_mom][i_ppt];

                //     gPrt->SetPoint(igPrt++,Y,Z,X);
                // } //end particlepoint loop

                //plot all electron deposits + total deposited energy
                double totalE=0;
                for (size_t i_dep=0; i_dep <= n_el_dep ; i_dep++) {

                    double X = (*T.DepX)[i_prt][i_dep];
                    double Y = (*T.DepY)[i_prt][i_dep];
                    double Z = (*T.DepZ)[i_prt][i_dep];
                    double E = (*T.DepE)[i_prt][i_dep];

                    gEl->SetPoint(igEl++,Y,Z,X);
                    if (E < 0) cerr << "\e[91mnegative energy deposit\e[0m of " << E << " MeV at evt#" << i_evt << ":prt#" << i_prt << ":dep#" << i_dep << endl;
                    totalE+=E;
                }
                hElE->Fill(totalE);
                hNDep->Fill(n_el_dep);
            } //end particle loop
        } //end event loop
    } //end file loop

    TCanvas* c1 = new TCanvas("c1","TrueMichel_v2");
    c1->cd();
    gMu->Draw("p");
    gEl->Draw("samep");
    // gPrt->Draw("samep");
    // gMux->Draw("samep");

    TCanvas* c2 = new TCanvas("c2","TrueMichel_v2");
    c2->Divide(2,2);
    c2->cd(1);
    hNDep->Draw("hist");
    c2->cd(2);
    hdEdx->Draw("colZ");
    c2->cd(3);
    hElE->Draw("hist");

    // c1->SaveAs("out/TrueMichel_v2Graph2D.root");
    // c2->SaveAs("out/TrueMichel_v2Hist.root");

    cout << N_evt << " events treated in " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}