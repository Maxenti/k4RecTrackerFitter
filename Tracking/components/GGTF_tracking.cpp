//======================================================================
// GGTF_tracking.cpp  (tracks + optional 3D hits output)  [diagnostic+mem-lite]
//======================================================================

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <memory>
#include <numeric>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <type_traits>  // for label stamping trait checks

// ONNX & Torch
#include <ATen/ATen.h>
#include <ATen/Parallel.h>               // at::set_num_threads / set_num_interop_threads
#include <torch/torch.h>
#include "onnxruntime_cxx_api.h"

// ROOT (helpers)
#include "TVector3.h"

// Gaudi + k4FWCore
#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartIF.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"

// EDM4hep & extensions
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHit3DCollection.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/TrackCollection.h"

// DD4hep (optional)
#include "DD4hep/Detector.h"
#include "DDSegmentation/BitFieldCoder.h"

// Local
#include "utils.hpp"

namespace {
inline edm4hep::CovMatrix3f small_cov_3d() {
  // order: xx, xy, yy, xz, yz, zz
  return {1e-4f, 0.f, 1e-4f, 0.f, 0.f, 1e-4f};
}

inline std::pair<long,long> readRSSkB() {
  // returns {VmRSS_kB, VmHWM_kB}
  std::ifstream f("/proc/self/status");
  std::string key;
  long rss=0, hwm=0, val=0; std::string unit;
  while (f >> key) {
    if (key == "VmRSS:") { f >> val >> unit; rss = val; }
    else if (key == "VmHWM:") { f >> val >> unit; hwm = val; }
    f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return {rss,hwm};
}

struct StepTimer {
  std::chrono::steady_clock::time_point t0;
  StepTimer() : t0(std::chrono::steady_clock::now()) {}
  double ms() const {
    auto dt = std::chrono::steady_clock::now() - t0;
    return std::chrono::duration_cast<std::chrono::milliseconds>(dt).count();
  }
};

// ---- label stamping helpers for 3D hits (schema-safe) ----
template <typename, typename = void> struct has_setQuality : std::false_type {};
template <typename T>
struct has_setQuality<T, std::void_t<decltype(std::declval<T&>().setQuality(int{}))>> : std::true_type {};

template <typename, typename = void> struct has_setType : std::false_type {};
template <typename T>
struct has_setType<T, std::void_t<decltype(std::declval<T&>().setType(int{}))>> : std::true_type {};

template <typename Hit3D>
inline void set_label_on_3d(Hit3D& h, int label) {
  if constexpr (has_setQuality<Hit3D>::value) {
    h.setQuality(label);
  } else if constexpr (has_setType<Hit3D>::value) {
    h.setType(label);
  } // else: schema lacks both; do nothing
}

} // namespace

struct GGTF_tracking final
  : k4FWCore::MultiTransformer<
        std::tuple<extension::TrackCollection, edm4hep::TrackerHit3DCollection>(
            const std::vector<const edm4hep::TrackerHitPlaneCollection*>&,
            const std::vector<const extension::SenseWireHitCollection*>&)> {

  GGTF_tracking(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc,
                         // inputs
                         {
                             KeyValues("inputPlanarHits", {"inputPlanarHits"}),
                             KeyValues("inputWireHits",   {"inputWireHits"}),
                         },
                         // outputs
                         {
                             KeyValues("outputTracks", {"outputTracks"}),
                             KeyValues("output3DHits", {"GGTF_3DHits"}),
                         }) {
    // Acquire GeoSvc (optional)
    m_geoSvc = serviceLocator()->service(m_geoSvcName.value());
  }

  // ---- properties ----
  Gaudi::Property<std::string> modelPath{this, "modelPath", "", "Path to ONNX model"};
  Gaudi::Property<double>      tbeta{this, "tbeta", 0.6, "clustering beta threshold"};
  Gaudi::Property<double>      td{this, "td", 0.3, "clustering distance threshold"};

  Gaudi::Property<bool>        produce3DHits{this, "produce3DHits", true,
                                            "If false, do not create GGTF_3DHits (returns empty collection)"};

  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "GeoSvc instance name"};
  Gaudi::Property<std::string> m_DCH_name{this, "DCH_name", "DCH_v2", "Drift chamber name"};

  // Safety: cap hits per event (0 = unlimited)
  Gaudi::Property<int>         maxHitsPerEvent{this, "maxHitsPerEvent", 0,
                                             "If >0, hard-limit number of input hits per event"};

  // ---- lifecycle ----
  StatusCode initialize() override {
    // Torch threading (be extra conservative)
    at::set_num_threads(1);
    at::set_num_interop_threads(1);

    // ONNX init (lean settings)
    fInfo = std::make_unique<Ort::MemoryInfo>(
        Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault));
    fEnv = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "ONNX_Runtime");
    fSessionOptions.SetIntraOpNumThreads(1);
    fSessionOptions.SetInterOpNumThreads(1);
    fSessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);
    fSessionOptions.DisableMemPattern();
    fSessionOptions.SetExecutionMode(ExecutionMode::ORT_SEQUENTIAL);

    fSession = std::make_unique<Ort::Session>(*fEnv, modelPath.value().c_str(), fSessionOptions);

    {
      Ort::AllocatorWithDefaultOptions alloc;
      const std::size_t i = 0;
      // Store names as std::string (avoid leaking C-strings)
      fInamesStr.emplace_back(fSession->GetInputNameAllocated(i,  alloc).get());
      fOnamesStr.emplace_back(fSession->GetOutputNameAllocated(i, alloc).get());
      fInames   = { fInamesStr.back().c_str() };
      fOnames   = { fOnamesStr.back().c_str() };
    }

    // Optional: geometry decoder (safe if absent)
    if (m_geoSvc && m_geoSvc->getDetector()) {
      try {
        auto* det = m_geoSvc->getDetector();
        auto sd   = det->sensitiveDetector(m_DCH_name.value());
        if (sd.isValid()) {
          auto ro = sd.readout();
          if (ro.isValid()) {
            dc_decoder = ro.idSpec().decoder();
          }
        }
      } catch (...) { /* continue without decoder */ }
    }

    auto [rss,hwm] = readRSSkB();
    info() << "GGTF_tracking init | modelPath=" << modelPath.value()
           << " | DCH=" << m_DCH_name.value()
           << " | decoder=" << (dc_decoder ? "OK" : "null")
           << " | produce3DHits=" << (produce3DHits.value() ? "true" : "false")
           << " | RSS=" << rss/1024.0 << "MB, Peak=" << hwm/1024.0 << "MB"
           << endmsg;
    return StatusCode::SUCCESS;
  }

  std::tuple<extension::TrackCollection, edm4hep::TrackerHit3DCollection>
  operator()(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& inputPlanar,
             const std::vector<const extension::SenseWireHitCollection*>&   inputWire) const override {
    torch::NoGradGuard nograd; // no autograd anywhere in this algorithm
    ++m_evt;

    extension::TrackCollection      outTracks;
    edm4hep::TrackerHit3DCollection out3D;

    auto logMem = [&](const char* tag) {
      auto [rss,hwm] = readRSSkB();
      info() << "[evt " << m_evt << "] " << tag
             << " | RSS=" << rss/1024.0 << "MB, Peak=" << hwm/1024.0 << "MB" << endmsg;
    };

    StepTimer t_all;

    // ------- flatten inputs + remember (type, coll, idx) -------
    int64_t nPlanar=0, nWire=0;
    for (auto c : inputPlanar) nPlanar += c ? c->size() : 0;
    for (auto c : inputWire)   nWire   += c ? c->size() : 0;
    int64_t nHitsEst = nPlanar + nWire;
    if (maxHitsPerEvent > 0 && nHitsEst > maxHitsPerEvent) {
      warning() << "[evt " << m_evt << "] capping hits " << nHitsEst
                << " -> " << int(maxHitsPerEvent) << endmsg;
      nHitsEst = maxHitsPerEvent;
    }

    std::vector<float>   gInputs;   gInputs.reserve(std::max<int64_t>(nHitsEst*7, 128));
    std::vector<int64_t> tagType;   tagType.reserve(std::max<int64_t>(nHitsEst, 128));   // 0=planar, 1=wire
    std::vector<int64_t> tagIndexA; tagIndexA.reserve(tagType.capacity()); // collection id
    std::vector<int64_t> tagIndexB; tagIndexB.reserve(tagType.capacity()); // hit index

    auto push_planar = [&](int icoll, int ihit, const edm4hep::TrackerHitPlane& h) {
      const auto p = h.getPosition();
      gInputs.insert(gInputs.end(), {float(p.x), float(p.y), float(p.z), 1.f, 0.f, 0.f, 0.f});
      tagType.push_back(0); tagIndexA.push_back(icoll); tagIndexB.push_back(ihit);
    };

    auto push_wire = [&](int icoll, int ihit, const extension::SenseWireHit& h) {
      const auto wp = h.getPosition();
      const double d   = h.getDistanceToWire();
      const double phi = h.getWireAzimuthalAngle();
      const double st  = h.getWireStereoAngle();

      TVector3 wpos(wp.x, wp.y, wp.z);
      TVector3 dir(0,0,1); dir.RotateX(st); dir.RotateZ(phi); dir = dir.Unit();

      TVector3 xprime(1.0, 0.0, -dir.X()/std::max(1e-9, dir.Z())); xprime = xprime.Unit();
      TVector3 yprime = dir.Cross(xprime).Unit();

      const TVector3 l(-d,0,0), r(+d,0,0);
      const auto L = xprime*l.X() + yprime*l.Y() + dir*l.Z() + wpos;
      const auto R = xprime*r.X() + yprime*r.Y() + dir*r.Z() + wpos;
      const auto M = 0.5*(L+R);

      gInputs.insert(gInputs.end(),
                     {float(M.X()), float(M.Y()), float(M.Z()),
                      0.f, float(R.X()-L.X()), float(R.Y()-L.Y()), float(R.Z()-L.Z())});
      tagType.push_back(1); tagIndexA.push_back(icoll); tagIndexB.push_back(ihit);
    };

    {
      StepTimer t_flat;
      int ic = 0;
      for (auto coll : inputPlanar) {
        if (!coll) { ++ic; continue; }
        for (int i=0, n=coll->size(); i<n; ++i) {
          if (maxHitsPerEvent>0 && (int)tagType.size()>=maxHitsPerEvent) break;
          push_planar(ic, i, (*coll)[i]);
        }
        ++ic;
        if (maxHitsPerEvent>0 && (int)tagType.size()>=maxHitsPerEvent) break;
      }
      ic = 0;
      for (auto coll : inputWire) {
        if (!coll) { ++ic; continue; }
        for (int i=0, n=coll->size(); i<n; ++i) {
          if (maxHitsPerEvent>0 && (int)tagType.size()>=maxHitsPerEvent) break;
          push_wire(ic, i, (*coll)[i]);
        }
        ++ic;
        if (maxHitsPerEvent>0 && (int)tagType.size()>=maxHitsPerEvent) break;
      }
      info() << "[evt " << m_evt << "] flatten: "
             << "planar=" << nPlanar << " wire=" << nWire
             << " -> used=" << tagType.size()
             << " in " << t_flat.ms() << " ms" << endmsg;
    }

    const int64_t nHits = static_cast<int64_t>(tagType.size());
    if (nHits == 0) {
      logMem("no-hits (early return)");
      return {std::move(outTracks), std::move(out3D)};
    }
    logMem("after-flatten");

    // ------- ONNX run -------
    std::vector<float> embed;
    {
      StepTimer t_onnx;
      const std::vector<int64_t> shape{nHits, 7};
      std::vector<Ort::Value> inputs;
      inputs.emplace_back(Ort::Value::CreateTensor<float>(*fInfo, gInputs.data(),
                                                          gInputs.size(), shape.data(), shape.size()));
      auto outputs = fSession->Run(Ort::RunOptions{nullptr},
                                   fInames.data(), inputs.data(), fInames.size(),
                                   fOnames.data(), fOnames.size());
      float* outptr = outputs.front().GetTensorMutableData<float>();
      embed.assign(outptr, outptr + 4*nHits); // 4-dim embedding per hit
      // free ONNX tensors and the input buffer ASAP
      outputs.clear();
      inputs.clear();
      gInputs.clear();
      gInputs.shrink_to_fit();

      info() << "[evt " << m_evt << "] onnx: nHits=" << nHits
             << " in " << t_onnx.ms() << " ms" << endmsg;
    }
    logMem("after-onnx");

    // ------- clustering -------
    torch::Tensor clustering;
    {
      StepTimer t_cluster;
      clustering = get_clustering(embed, nHits, tbeta.value(), td.value()); // int64 labels, len=nHits
      info() << "[evt " << m_evt << "] clustering: "
             << " in " << t_cluster.ms() << " ms" << endmsg;
    }
    embed.clear(); embed.shrink_to_fit();
    logMem("after-clustering");

    // Unique labels and inverse index
    torch::Tensor uniques, invIdx;
    {
      StepTimer t_uni;
      std::tie(uniques, invIdx) = at::_unique(clustering, /*sorted=*/true, /*return_inverse=*/true);
      info() << "[evt " << m_evt << "] unique: nLabels=" << uniques.size(0)
             << " in " << t_uni.ms() << " ms" << endmsg;
    }

    // Build groups in one pass using invIdx (AVOIDS (label==..) masks)
    std::vector<std::vector<int64_t>> groups;
    int64_t zeroLabelPos = -1;
    torch::Tensor uniques_cpu;  // keep CPU view for label values later
    {
      StepTimer t_bucket;
      uniques_cpu = uniques.to(torch::kCPU);
      auto inv_cpu = invIdx.to(torch::kCPU).contiguous();
      const int64_t nLabels = uniques_cpu.size(0);
      groups.resize(nLabels);

      // find position of label "0" (noise) if present
      for (int64_t i=0; i<nLabels; ++i) {
        if (uniques_cpu[i].item<int64_t>() == 0) { zeroLabelPos = i; break; }
      }

      auto invAcc = inv_cpu.accessor<int64_t,1>();
      for (int64_t idx=0; idx<nHits; ++idx) {
        const int64_t pos = invAcc[idx]; // index into uniques
        if (pos>=0 && pos<nLabels) groups[pos].push_back(idx);
      }
      info() << "[evt " << m_evt << "] bucket: labels=" << nLabels
             << " in " << t_bucket.ms() << " ms" << endmsg;
    }
    // free big tensors early
    clustering = torch::Tensor();
    invIdx     = torch::Tensor();
    logMem("after-bucket");

    // ------- helpers to emit 3D hits (used only if produce3DHits) -------
    auto make3D_from_planar = [&](const edm4hep::TrackerHitPlane& hp, int label) {
      auto h = out3D.create();
      const auto p = hp.getPosition();
      h.position()  = edm4hep::Vector3d{p.x, p.y, p.z};
      h.covMatrix() = small_cov_3d();
      h.time()      = 0.f;
      set_label_on_3d(h, label); // stamp GGTF label
      return h;
    };
    auto make3D_from_wire = [&](const extension::SenseWireHit& hw, int label) {
      const auto wp = hw.getPosition();
      const double d   = hw.getDistanceToWire();
      const double phi = hw.getWireAzimuthalAngle();
      const double st  = hw.getWireStereoAngle();

      TVector3 wpos(wp.x, wp.y, wp.z);
      TVector3 dir(0,0,1); dir.RotateX(st); dir.RotateZ(phi); dir = dir.Unit();

      TVector3 xprime(1.0, 0.0, -dir.X()/std::max(1e-9, dir.Z())); xprime = xprime.Unit();
      TVector3 yprime = dir.Cross(xprime).Unit();

      const TVector3 l(-d,0,0), r(+d,0,0);
      const auto L = xprime*l.X() + yprime*l.Y() + dir*l.Z() + wpos;
      const auto R = xprime*r.X() + yprime*r.Y() + dir*r.Z() + wpos;
      const auto M = 0.5*(L+R);

      auto h = out3D.create();
      h.position()  = edm4hep::Vector3d{M.X(), M.Y(), M.Z()};
      h.covMatrix() = small_cov_3d();
      h.time()      = 0.f;
      set_label_on_3d(h, label); // stamp GGTF label
      return h;
    };

    // quick helper to attach relation & optionally create 3D
    auto add_hit_to_track = [&](extension::MutableTrack& trk, int64_t flatIdx, int label) {
      const int64_t t  = tagType[flatIdx];
      const int64_t ia = tagIndexA[flatIdx];
      const int64_t ib = tagIndexB[flatIdx];

      if (t == 0) {
        const auto& hp = (*inputPlanar[ia])[int(ib)];
        if (produce3DHits.value()) { (void)make3D_from_planar(hp, label); }
        trk.addToTrackerHits(hp);
      } else {
        const auto& hw = (*inputWire[ia])[int(ib)];
        if (produce3DHits.value()) { (void)make3D_from_wire(hw, label); }
        trk.addToTrackerHits(hw);
      }
    };

    // ------- assemble tracks from groups -------
    {
      StepTimer t_build;
      const int64_t nLabels = static_cast<int64_t>(groups.size());
      int nTracks = 0;
      for (int64_t li = 0; li < nLabels; ++li) {
        // skip noise label "0" position (if any)
        if (zeroLabelPos == li) continue;
        if (groups[li].empty()) continue;

        auto trk = outTracks.create();
        const int labelValue = uniques_cpu[li].item<int>(); // use CPU tensor (avoid device sync)
        trk.setType(labelValue);

        for (int64_t k : groups[li]) add_hit_to_track(trk, k, labelValue); // pass label to 3D maker
        ++nTracks;
      }
      info() << "[evt " << m_evt << "] build: tracks=" << nTracks
             << " in " << t_build.ms() << " ms" << endmsg;
    }
    logMem("after-build");

    info() << "[evt " << m_evt << "] TOTAL " << t_all.ms() << " ms" << endmsg;
    return std::make_tuple(std::move(outTracks), std::move(out3D));
  }

  StatusCode finalize() override {
    auto [rss,hwm] = readRSSkB();
    info() << "GGTF_tracking finalize | events=" << m_evt
           << " | RSS=" << rss/1024.0 << "MB, Peak=" << hwm/1024.0 << "MB"
           << endmsg;
    return StatusCode::SUCCESS;
  }

private:
  // ONNX
  std::unique_ptr<Ort::Env>          fEnv;
  std::unique_ptr<Ort::Session>      fSession;
  Ort::SessionOptions                fSessionOptions;
  std::unique_ptr<Ort::MemoryInfo>   fInfo;
  // names as strings to avoid heap leaks from AllocatedStringPtr::release()
  std::vector<std::string>           fInamesStr, fOnamesStr;
  std::vector<const char*>           fInames, fOnames;

  // Geometry (optional)
  SmartIF<IGeoSvc> m_geoSvc;
  dd4hep::DDSegmentation::BitFieldCoder* dc_decoder{nullptr};

  // stats
  mutable int m_evt{0};
};

DECLARE_COMPONENT(GGTF_tracking)
