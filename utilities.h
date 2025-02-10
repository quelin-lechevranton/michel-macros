#include "detector.h"

#include <iostream>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>

struct Binning {
    int n;
    float min, max;
    float step() { return (max - min) / n; }

    Binning() : n(0), min(0), max(0) {}
    Binning(unsigned n, float min, float max) : n(n), min(min), max(max) {}
    Binning(float min, float max, float step) : n((max - min) / step), min(min), max(max) {}
    Binning(float c, float r, Binning b) {
        *this = Binning(
            c-r > b.min ? c-r : b.min,
            c+r < b.max ? c+r : b.max,
            b.step());
    }
};


unsigned GetSlice(int ch) {
    unsigned tpc=0;
    while (tpc<n_tpc && !(map_tpc_chs[tpc].first <= ch && ch <= map_tpc_chs[tpc].second)) tpc++;
    if (tpc == n_tpc) return 8;
    return unsigned(tpc/4 * 2 + tpc % 2);
}

struct Hit {
    unsigned slice;
    float Z;
    int channel;
    float tick;
    float adc;

    float rev_tick() {
        return slice >= 4 ? tick_window - tick : tick;
    }
    float tick_cm() {
        return tick * sampling_rate * drift_velocity;
    }
    float rev_tick_cm() {
        return rev_tick() * sampling_rate * drift_velocity;
    }
    float display_tick();
    float display4_tick();
};
struct Hits {
    unsigned N;
    std::vector<unsigned> *slice; 
    std::vector<float> *Z;
    std::vector<int> *channel;
    std::vector<float> *tick;
    std::vector<float> *adc;

    Hits() : 
        N(0), 
        slice(nullptr), 
        Z(nullptr), 
        channel(nullptr), 
        tick(nullptr), 
        adc(nullptr) {} 
    Hit at(unsigned i) {
        return Hit{
            slice->at(i), 
            Z->at(i), 
            channel->at(i), 
            tick->at(i), 
            adc == nullptr ? 0.F : adc->at(i)
        };
    }
};
struct BasicHits {
    unsigned N;
    std::vector<int> *channel;
    std::vector<float> *tick;
    std::vector<float> *score;

    BasicHits() : 
        N(0),
        channel(nullptr), 
        tick(nullptr), 
        score(nullptr) {}
    Hit at(unsigned i) {
        return Hit{
            GetSlice(channel->at(i)),
            map_ch_z[channel->at(i)],
            channel->at(i),
            tick->at(i),
            score == nullptr ? 0.F : score->at(i)
        };
    }
};

template <typename T>
T CLInput(const char *prompt, T def) {
    std::string input;
    std::cout << prompt << " (" << def << ") ";
    std::getline(std::cin, input);
    if (input.empty()) 
        return def;
    return T(stof(input));
}
template <>
TString CLInput<TString>(const char *prompt, TString def) {
    std::string input;
    std::cout << prompt << " (" << def << ") ";
    std::getline(std::cin, input);
    if (input.empty()) 
        return def;
    if (input.front() == '_')
        return def + TString(input);
    return TString(input);
}
template <>
Binning CLInput<Binning>(const char *prompt, Binning def) {
    float min = CLInput<float>(Form("%s_min", prompt), def.min);
    float max = CLInput<float>(Form("%s_max", prompt), def.max);
    float step = CLInput<float>(Form("%s_step", prompt), def.step());
    return Binning(min, max, step);
}
TString removePath(TString name) {
    return name.Last('/') == kNPOS ? 
        name : 
        name(name.Last('/') + 1, name.Length() - name.Last('/') - 1);
}