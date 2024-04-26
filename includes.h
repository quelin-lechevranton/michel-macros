// Framework includes
// #include "art/Framework/Core/EDAnalyzer.h"
// #include "art/Framework/Core/ModuleMacros.h"
// // #include "art/Framework/Principal/Event.h"
// #include "art/Framework/Principal/Run.h"
// #include "art/Framework/Principal/SubRun.h"
// #include "art/Framework/Principal/Handle.h"
// #include "art/Framework/Principal/View.h"
// #include "art/Framework/Services/Registry/ServiceHandle.h"
// #include "art_root_io/TFileService.h"
// #include "art_root_io/TFileDirectory.h"
// #include "canvas/Utilities/InputTag.h"
// #include "canvas/Persistency/Common/Ptr.h"
// #include "canvas/Persistency/Common/FindMany.h"
// #include "canvas/Persistency/Common/FindManyP.h"
// #include "canvas/Persistency/Common/PtrVector.h"

// // Utility libraries
// #include "fhiclcpp/ParameterSet.h"
// #include "messagefacility/MessageLogger/MessageLogger.h"

// // LArSoft includes
// #include "nusimdata/SimulationBase/MCTruth.h"
// #include "nusimdata/SimulationBase/MCNeutrino.h"
// #include "nusimdata/SimulationBase/MCParticle.h" // simb::MCParticle
// #include "nusimdata/SimulationBase/MCTrajectory.h"
// #include "lardataobj/Simulation/SimEnergyDeposit.h"
// #include "lardataobj/RawData/raw.h"
// #include "lardataobj/RecoBase/Hit.h"
// #include "lardataobj/RecoBase/Track.h"
// #include "lardataobj/RecoBase/Shower.h"
// #include "lardataobj/RecoBase/Cluster.h"
// #include "lardataobj/RecoBase/PFParticle.h" //top level information form Pandora inLArSoft in recob::PFParticle
// #include "lardataobj/RawData/OpDetWaveform.h"
// #include "lardataobj/AnalysisBase/Calorimetry.h"
// #include "larcore/CoreUtils/ServiceUtil.h"
// #include "larcore/Geometry/Geometry.h" // DetSim Analysis
// #include "larcorealg/Geometry/GeometryCore.h" // DetSim Analysis
// #include "larcoreobj/SimpleTypesAndConstants/geo_types.h"// DetSim Analysis
// #include "lardata/ArtDataHelper/TrackUtils.h"
// #include "lardata/ArtDataHelper/HitCreator.h" // RawDigit
// #include "lardata/DetectorInfoServices/DetectorClocksService.h"
// #include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
// #include "larsim/MCCheater/BackTrackerService.h"
// #include "larsim/MCCheater/PhotonBackTrackerService.h"
// #include "larsim/MCCheater/ParticleInventoryService.h"
// #include "larsim/MCCheater/BackTracker.h"
// #include "larsim/Utils/TruthMatchUtils.h"
// #include "larsim/MCCheater/BackTrackerService.h"

// ROOT includes.
#include "TH1.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraph2D.h"

// C++ includes
#include <TTree.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <vector>
#include <cmath>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

// #include "duneana/DAQSimAna/TriggerPrimitiveFinderTool.h"

// #include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
// #include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
// #include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
// #include "dunereco/AnaUtils/DUNEAnaEventUtils.h"

//#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// #include "tools.h"