#include "FMD_tools.h"

const size_t n_file=1;
const vector<string> filelist = fmd::readFileList(n_file,"list/ijclab.list");

void HitSpt() {
    clock_t start_time=clock();

    vector<unsigned int> color = { kBlue-3, kGreen-2, kRed-3 };
    vector<TGraph2D*> gSpt(3);
    for (size_t i=0; i<3; i++) {
        gSpt[i] = new TGraph2D();
        gSpt[i]->SetTitle(";Y (cm);Z (cm);X (cm)");
        gSpt[i]->SetMarkerStyle(20);
        gSpt[i]->SetMarkerSize(1-0.2*i);
        gSpt[i]->SetMarkerColor(color[i]);
    }

    vector<TGraph*> gSum(3);
    for (size_t i=0; i<3; i++) {
        gSum[i] = new TGraph();
        gSum[i]->SetTitle(";X (cm);SumADC");
        // gSum[i]->SetMarkerStyle(20);
        // gSum[i]->SetMarkerSize(1-0.2*i);
        gSum[i]->SetFillColorAlpha(color[i],0.4);
    }

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        // fmd::Truth T(filename.c_str());
        fmd::Reco R(filename.c_str());

        size_t nSpt=0;
        size_t n_evt = R.GetEntries();
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {
        // size_t i_evt=19; {

            cout << "Event#" << i_evt+1 << "/" << n_evt << "\r" << flush;

            // T.GetEntry(i_evt);
            R.GetEntry(i_evt);

            size_t nSptMultiHit=0;
            for (size_t i_pfp=0; i_pfp < R.NPfp; i_pfp++) {
                cout << "\t\tpfp#" << i_pfp << "\r" << flush;

                nSpt+=R.PfpNSpt->at(i_pfp);

                for (size_t i_hit=0; i_hit < R.PfpNSpt->at(i_pfp); i_hit++) {
                cout << "\t\t\thit#" << i_hit << "\r" << flush;

                    size_t nHit=0;
                    for (size_t i_spt=0; i_spt < R.PfpNSpt->at(i_pfp); i_spt++) {
                        cout << "\t\t\t\tspt#" << i_spt << "\r" << flush;

                        if ((*R.SptHitID)[i_pfp][i_spt]!=i_hit) continue;
                        nHit++;

                        vector<double> Spt = { (*R.SptX)[i_pfp][i_spt], (*R.SptY)[i_pfp][i_spt], (*R.SptZ)[i_pfp][i_spt] };

                        size_t Plane = (*R.HitPlane)[i_pfp][i_hit];
                        gSpt[Plane]->AddPoint(Spt[1],Spt[2],Spt[0]);
                        gSum[Plane]->AddPoint(Spt[0],(*R.HitSumADC)[i_pfp][i_hit]);
                    } //end spacepoint loop
                    if (nHit>1) nSptMultiHit++;
                } //end hit loop
            } //end pfp loop
            cout << "\t\t\t\tnSptMultiHit: " << 100.*nSptMultiHit/nSpt << endl;
        } //end event loop
    } //end file loop
    cout << endl;

    TCanvas* c1 = new TCanvas("c1","HitSpt");
    c1->Divide(2,1);
    c1->cd(1);
    gSpt[0]->Draw("p");
    gSpt[1]->Draw("samep");
    gSpt[2]->Draw("samep");
    c1->cd(2);
    gSum[0]->Draw("ab");
    gSum[1]->Draw("sameb");
    gSum[2]->Draw("sameb");

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}