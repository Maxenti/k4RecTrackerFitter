//======================================================================
// GenFit2DCHFitter.cpp  -- fit GGTF 3D hits with GenFit2 (minimal working)
//======================================================================

#include <memory>
#include <vector>
#include <string>
#include <cmath>

#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/ISvcLocator.h"

#include "k4FWCore/Transformer.h"

// EDM4hep
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackCollection.h"

// ROOT
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"

// GenFit2
#include "RKTrackRep.h"
#include "Track.h"
#include "KalmanFitterRefTrack.h"
#include "SpacepointMeasurement.h"
#include "TrackPoint.h"
#include "MeasuredStateOnPlane.h"
#include "Exception.h"
#include "FieldManager.h"
#include "ConstField.h"

namespace {

// conservative seed covariance (pos [mm]^2, mom [(GeV/c)^2])
inline TMatrixDSym defaultSeedCov() {
  TMatrixDSym C(6);
  C.Zero();
  for (int i = 0; i < 3; ++i) C(i,i) = 10.*10.;   // 1 cm^2
  for (int i = 3; i < 6; ++i) C(i,i) = 1.0*1.0;   // (1 GeV)^2
  return C;
}

} // namespace

struct GenFit2DCHFitter final
  : k4FWCore::Transformer<edm4hep::TrackCollection (const edm4hep::TrackerHit3DCollection&)> {

  using Traits    = Gaudi::Functional::Traits::use_<>;
  using KeyValues = Gaudi::Functional::details::DataHandleMixin<
                      std::tuple<>, std::tuple<>, Traits>::KeyValues;

  GenFit2DCHFitter(const std::string& name, ISvcLocator* svcLoc)
  : Transformer(name, svcLoc,
      /* inputs  */ std::tuple<KeyValues>{ KeyValues{"inputHits",  std::vector<std::string>{"GGTF_3DHits"}} },
      /* outputs */ std::tuple<KeyValues>{ KeyValues{"outputTracks", std::vector<std::string>{"GenFitTracks"}} }) {

    declareProperty("Bz",  m_Bz,  "Magnetic field Bz [Tesla]");
    declareProperty("PDG", m_pdg, "PDG hypothesis (e.g. 13 for mu-)");
  }

  Gaudi::Property<double> m_Bz {this, "Bz", 2.0, "Bz field [T]"};
  Gaudi::Property<int>    m_pdg{this, "PDG", 13,  "PDG hypothesis"};

  StatusCode initialize() override {
    // Set a constant B field (straightforward & fast)
    genfit::FieldManager::getInstance()->init(new genfit::ConstField(0., 0., m_Bz.value()));

    m_fitter = std::make_unique<genfit::KalmanFitterRefTrack>();
    info() << "GenFit2DCHFitter initialized | Bz=" << m_Bz.value()
           << " | PDG=" << m_pdg.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  edm4hep::TrackCollection operator()(const edm4hep::TrackerHit3DCollection& hits) const override {
    edm4hep::TrackCollection out;
    if (hits.size() < 3u) return out;

    // naive seed from first/last hit
    const auto p0 = hits[0].getPosition();
    const auto pl = hits[hits.size()-1].getPosition();
    TVector3 pos0(p0.x, p0.y, p0.z);
    TVector3 posL(pl.x, pl.y, pl.z);
    TVector3 dir = (posL - pos0);
    if (dir.Mag2() < 1e-6) dir = TVector3(0,0,1);
    dir = dir.Unit();
    const double pMag = 1.0; // GeV/c
    TVector3 mom0 = pMag * dir;

    // rep & track (track takes ownership of rep pointer)
    auto rep_up = std::make_unique<genfit::RKTrackRep>(m_pdg.value());
    genfit::Track fitTrack(rep_up.release(), pos0, mom0);
    fitTrack.setCovSeed(defaultSeedCov());

    // add spacepoint measurements
    int detId = 0, hitId = 0;
    for (const auto& h : hits) {
      const auto q = h.getPosition();
      TVectorD pos(3); pos[0]=q.x; pos[1]=q.y; pos[2]=q.z;

      TMatrixDSym C(3); C.Zero();
      // map edm4hep cov if you have it; for now small diagonal:
      C(0,0)=1e-4; C(1,1)=1e-4; C(2,2)=1e-4;

      auto* tp   = new genfit::TrackPoint(&fitTrack);
      auto* meas = new genfit::SpacepointMeasurement(pos, C, detId, hitId, tp, false, false);
      tp->addRawMeasurement(meas);
      fitTrack.insertPoint(tp);
      ++hitId;
    }

    try {
      m_fitter->processTrack(&fitTrack);
    } catch (const genfit::Exception& e) {
      warning() << "GenFit exception: " << e.what() << endmsg;
      return out;
    } catch (const std::exception& e) {
      warning() << "Std exception: " << e.what() << endmsg;
      return out;
    } catch (...) {
      warning() << "Unknown exception during fit" << endmsg;
      return out;
    }

    if (fitTrack.getNumReps() == 0) {
      warning() << "No track reps after fit" << endmsg;
      return out;
    }
    genfit::AbsTrackRep* rep = fitTrack.getTrackRep(0);
    if (!rep) {
      warning() << "Null track rep" << endmsg;
      return out;
    }

    // get fit status
    const auto* fs = fitTrack.getFitStatus(rep);
    const int lastId = static_cast<int>(fitTrack.getNumPoints()) - 1;
    (void) fitTrack.getFittedState(lastId, rep, /*biased=*/true); // available if you want to export params

    // fill edm4hep track (simple)
    auto trk = out.create();
    trk.setType(m_pdg.value());
    if (fs) {
      trk.setChi2(fs->getChi2());
      trk.setNdf (fs->getNdf());
    }
    for (const auto& h : hits) trk.addToTrackerHits(h);

    return out;
  }

  StatusCode finalize() override { return StatusCode::SUCCESS; }

private:
  std::unique_ptr<genfit::KalmanFitterRefTrack> m_fitter;
};

DECLARE_COMPONENT(GenFit2DCHFitter)
