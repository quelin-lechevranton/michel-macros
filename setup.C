R__ADD_INCLUDE_PATH("gallery/Event.h")
R__ADD_INCLUDE_PATH("lardataobj/Simulation/SimEnergyDeposit.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_types.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_vectors.h")
R__ADD_INCLUDE_PATH("nusimdata/v1_27_01/include/nusimdata/SimulationBase")
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "tools.h"


art::InputTag monte_tag("generator::SinglesGen");
art::InputTag depo_tag("largeant:LArG4DetectorServicevolTPCActive");
art::InputTag point_tag("pandora");
art::InputTag track_tag("pandoraTrack");
art::InputTag calo_tag("pandoraGnocalo");

gallery::Event ev({"/eos/user/t/thoudy/pdvd/sims/out/protodunevd_10_muon_reco.root"});

auto const monte_list = ev.getValidHandle<vector<simb::MCTruth>>(monte_tag);
auto const depo_list = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(depo_tag);
auto const point_list = ev.getValidHandle<vector<recob::SpacePoint>>(point_tag);
auto const track_list = ev.getValidHandle<vector<recob::Track>>(track_tag);
auto const calo_list = ev.getValidHandle<vector<anab::Calorimetry>>(calo_tag);


// const simb::MCTruth& monte = monte_list->at(0);
// const sim::SimEnergyDeposit& depo = depo_list->at(0);
// const recob::SpacePoint& point = point_list->at(0);
// const recob::Track& track = track_list->at(0);
// const anab::Calorimetry& calo = calo_list->at(0);