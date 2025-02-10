#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <ROOT/RDataFrame.hxx>


void endprocess() {

    ROOT::RDataFrame muon("checks/Muon","data/pdvd_Muchecks_out_100_r20.root");

    auto in = muon.Filter("EndIsInWindow && EndIsInVolumeYZ");

    auto in_minus = in.Filter("!IsAnti");
    auto in_plus = in.Filter("IsAnti");

    auto n = in.Count().GetValue();
    auto n_minus = in_minus.Count().GetValue();
    auto n_plus = in_plus.Count().GetValue();

    auto vs_ep_minus = in_minus.Take<std::string>("MuonEndProcess");
    auto vs_ep_plus = in_plus.Take<std::string>("MuonEndProcess");

    std::map<std::string, unsigned int> ep_count_minus;
    for (auto ep : vs_ep_minus.GetValue()) {
        ep_count_minus[ep]++;
    }

    std::map<std::string, unsigned int> ep_count_plus;
    for (auto ep : vs_ep_plus.GetValue()) {
        ep_count_plus[ep]++;
    }

    std::cout << "number of muon inside volume: " << n << std::endl;
    std::cout << "µ- (" << n_minus << ", " << std::setprecision(3) << 100.*n_minus/n << "%)" << std::endl; 
    for (auto [ep, count] : ep_count_minus) {
        std::cout << "\t" << ep << ": " << count << ", " << std::setprecision(3) << 100.*count/n_minus << "%" << std::endl;
    }

    std::cout << "µ+ (" << n_plus << ", " << std::setprecision(3) << 100.*n_plus/n << "%)" << std::endl;
    for (auto [ep, count] : ep_count_plus) {
        std::cout << "\t" << ep << ": " << count << ", " << std::setprecision(3) << 100.*count/n_plus << "%" << std::endl;
    }
}
