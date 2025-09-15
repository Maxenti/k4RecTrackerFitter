// Tracking/src/SimTo3DHits.cpp
#include "Gaudi/Algorithm.h"
#include "k4FWCore/DataHandle.h"

#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/MutableTrackerHit3D.h"
#include "edm4hep/CovMatrix3f.h"

class SimTo3DHits final : public Gaudi::Algorithm {
public:
  SimTo3DHits(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm{name, svcLoc}
    , m_in {"InputSimHits", Gaudi::DataHandle::Reader, this}
    , m_out{"OutputHits",   Gaudi::DataHandle::Writer, this} {}

  StatusCode initialize() override {
    return Gaudi::Algorithm::initialize();
  }

  // Gaudi v39 execute signature; method is const
  StatusCode execute(const EventContext& /*ctx*/) const override {
    // DataHandle APIs are non-const -> declare handles 'mutable' (below)
    const auto* in  = m_in.get();
    auto*       out = m_out.createAndPut();

    for (const auto& s : *in) {
      edm4hep::MutableTrackerHit3D h;

      const auto& p = s.getPosition();  // dd4hep::Position (x,y,z) doubles
      h.setPosition({static_cast<float>(p.x),
                     static_cast<float>(p.y),
                     static_cast<float>(p.z)});

      // Packed 3x3 cov (xx, yy, zz, xy, xz, yz) – 6 numbers
      h.setCovMatrix(edm4hep::CovMatrix3f{1e-4f, 1e-4f, 1e-4f, 0.f, 0.f, 0.f});

      // Time if you want it (optional)
      h.setTime(0.f);

      out->push_back(h);
    }
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() override {
    return Gaudi::Algorithm::finalize();
  }

private:
  // Must be 'mutable' because execute() is const but k4FWCore handles are not
  mutable k4FWCore::DataHandle<edm4hep::SimTrackerHitCollection> m_in;
  mutable k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>  m_out;
};

DECLARE_COMPONENT(SimTo3DHits)
