#include <TString.h>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

ROOT::VecOps::RVec<TString> files = {
    "pdvd_Muchecks_out_100_r5_in.root",
    "pdvd_Muchecks_out_100_r10_in.root",
    "pdvd_Muchecks_out_100_r15_in.root",
    "pdvd_Muchecks_out_100_r20_in.root",
    "pdvd_Muchecks_out_100_r25_in.root",
    "pdvd_Muchecks_out_100_r30_in.root",
    "pdvd_Muchecks_out_100_r35_in.root",
    "pdvd_Muchecks_out_100_r40_in.root",
    "pdvd_Muchecks_out_100_r45_in.root",
    "pdvd_Muchecks_out_100_r50_in.root"
};

ROOT::VecOps::RVec<ROOT::RDataFrame> muon = {
    ROOT::RDataFrame("checks/Muon", files[0]),
    ROOT::RDataFrame("checks/Muon", files[1]),
    ROOT::RDataFrame("checks/Muon", files[2]),
    ROOT::RDataFrame("checks/Muon", files[3]),
    ROOT::RDataFrame("checks/Muon", files[4]),
    ROOT::RDataFrame("checks/Muon", files[5]),
    ROOT::RDataFrame("checks/Muon", files[6]),
    ROOT::RDataFrame("checks/Muon", files[7]),
    ROOT::RDataFrame("checks/Muon", files[8]),
    ROOT::RDataFrame("checks/Muon", files[9])
};


ROOT::VecOps::RVec<ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>> muon_in = {
    muon[0].Filter("EndIsInWindow && EndIsInVolumeYZ"),
    muon[1].Filter("EndIsInWindow && EndIsInVolumeYZ"),
    muon[2].Filter("EndIsInWindow && EndIsInVolumeYZ"),
    muon[3].Filter("EndIsInWindow && EndIsInVolumeYZ"),
    muon[4].Filter("EndIsInWindow && EndIsInVolumeYZ"),
    muon[5].Filter("EndIsInWindow && EndIsInVolumeYZ"),
    muon[6].Filter("EndIsInWindow && EndIsInVolumeYZ"),
    muon[7].Filter("EndIsInWindow && EndIsInVolumeYZ"),
    muon[8].Filter("EndIsInWindow && EndIsInVolumeYZ"),
    muon[9].Filter("EndIsInWindow && EndIsInVolumeYZ")
};

ROOT::VecOps::RVec<ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>> muon_decay = {
    muon_in[0].Filter("DoesDecay"),
    muon_in[1].Filter("DoesDecay"),
    muon_in[2].Filter("DoesDecay"),
    muon_in[3].Filter("DoesDecay"),
    muon_in[4].Filter("DoesDecay"),
    muon_in[5].Filter("DoesDecay"),
    muon_in[6].Filter("DoesDecay"),
    muon_in[7].Filter("DoesDecay"),
    muon_in[8].Filter("DoesDecay"),
    muon_in[9].Filter("DoesDecay")
};

void loadrdf() {}
