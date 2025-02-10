#include "utilities.h"

void which_muon(void) {
    TString name = CLInput<TString>("name", "data/pdvd_Muchecks_out");

    TString in = name.Contains(".root") ? name : name += ".root";

    std::cout << "opening " << in << "..." << std::endl;

    TFile* fin = new TFile(in);

    unsigned i_mu = CLInput<unsigned>("i_mu", 0);

    TTree* muon = fin->Get<TTree>("checks/Muon");

    unsigned long i_event;
    muon->SetBranchAddress("iEvent", &i_event);

    muon->GetEntry(i_mu);
    unsigned long muon_event = i_event;

    i_event = 0;
    unsigned mu=0;
    while (i_event != muon_event) muon->GetEntry(mu++);
    std::cout << "e" << muon_event << "µ" << i_mu-mu << std::endl;
}