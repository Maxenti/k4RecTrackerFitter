// ======================================================================
// GenFit2DCHFitter.cpp  -- fit GGTF 3D hits with GenFit2 (robust)
//   - TGeo material effects (cm units internally)
//   - Group hits per GGTF label (quality/type); DBSCAN fallback if all zero
//   - Diagonal anisotropic measurement covariance (XY/Z) + PD guards
//   - XY-circle pT seeding with guards; per-group sort + optional dedup
//   - Clamp on pT AND |p| to avoid βγ→0 failures
//   - Optional cap on measurements per group (for conditioning)
//   - Retry path: inflate covariances if no FitterInfo after base fit
//   - Canonical SpacepointMeasurement attach (7-arg ABI)
// ======================================================================

#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <exception>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <numeric>
#include <queue>

// Gaudi
#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/ISvcLocator.h"

// k4FWCore
#include "k4FWCore/Transformer.h"

// EDM4hep
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackCollection.h"

// ROOT
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TGeoManager.h"
#include "TDecompChol.h"

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
#include "MaterialEffects.h"
#include "TGeoMaterialInterface.h"

namespace {

// --- trait helpers to fetch a label without breaking older EDM4hep builds
template <typename, typename = void> struct has_getQuality : std::false_type {};
template <typename T> struct has_getQuality<T, std::void_t<decltype(std::declval<T>().getQuality())>> : std::true_type {};
template <typename, typename = void> struct has_getType    : std::false_type {};
template <typename T> struct has_getType<T, std::void_t<decltype(std::declval<T>().getType())>>       : std::true_type {};

template <typename HitT>
inline int hitLabel(const HitT& h) {
  if constexpr (has_getQuality<HitT>::value) return static_cast<int>(h.getQuality());
  else if constexpr (has_getType<HitT>::value) return static_cast<int>(h.getType());
  else return 0;
}

inline void makeDiagonalFloor(TMatrixDSym& C, double eps) {
  for (int i = 0; i < C.GetNrows(); ++i) C(i,i) = std::max(C(i,i), eps);
}

inline void ensurePD(TMatrixDSym& C, double floorDiag, double inflateFactor, int maxIters) {
  makeDiagonalFloor(C, floorDiag);
  for (int it = 0; it < maxIters; ++it) {
    TDecompChol chol(C);
    if (chol.Decompose()) return;
    for (int i = 0; i < C.GetNrows(); ++i) C(i,i) *= inflateFactor;
  }
  for (int i = 0; i < C.GetNrows(); ++i) C(i,i) += floorDiag;
}

inline double projS(const TVector3& origin, const TVector3& dirUnit, const TVector3& p) {
  return dirUnit.Dot(p - origin);
}

// 3-point circle in XY; returns R in the same units as inputs (cm inside here)
inline double circleRadiusXY(const TVector3& A, const TVector3& B, const TVector3& C) {
  const double x1=A.X(), y1=A.Y(), x2=B.X(), y2=B.Y(), x3=C.X(), y3=C.Y();
  const double a = x1*(y2 - y3) - y1*(x2 - x3) + x2*y3 - x3*y2;
  const double d = 2.0 * a;
  if (std::abs(d) < 1e-9) return 1e9; // near-line → huge R
  const double x1s=x1*x1, y1s=y1*y1, x2s=x2*x2, y2s=y2*y2, x3s=x3*x3, y3s=y3*y3;
  const double ux = ((x1s+y1s)*(y2-y3) + (x2s+y2s)*(y3-y1) + (x3s+y3s)*(y1-y2)) / d;
  const double uy = ((x1s+y1s)*(x3-x2) + (x2s+y2s)*(x1-x3) + (x3s+y3s)*(x2-x1)) / d;
  const double dx = x1 - ux, dy = y1 - uy;
  const double R = std::sqrt(dx*dx + dy*dy);
  return (std::isfinite(R) && R>1e-9) ? R : 1e9;
}

struct Group {
  std::vector<TVector3> P;     // positions (internal units)
  std::vector<size_t>   idxs;  // indices back to original hits
};

// Minimal DBSCAN (cm units)
static std::vector<int> dbscan(
  const std::vector<TVector3>& P,
  double eps, unsigned minPts)
{
  const int N = static_cast<int>(P.size());
  const double eps2 = eps*eps;
  std::vector<int> labels(N, -99); // -99=UNVISITED, -1=NOISE, >=0 clusterId
  int clusterId = 0;

  auto regionQuery = [&](int i, std::vector<int>& out){
    out.clear();
    const auto& pi = P[i];
    for (int j=0;j<N;++j){
      if ((P[j]-pi).Mag2() <= eps2) out.push_back(j);
    }
  };

  for (int i=0;i<N;++i){
    if (labels[i] != -99) continue;
    std::vector<int> neigh; regionQuery(i, neigh);
    if (neigh.size() < minPts){ labels[i] = -1; continue; }
    labels[i] = clusterId;
    std::queue<int> q;
    for (int n : neigh) if (labels[n] == -99 || labels[n] == -1){ labels[n] = clusterId; q.push(n); }
    while(!q.empty()){
      int k = q.front(); q.pop();
      std::vector<int> neigh2; regionQuery(k, neigh2);
      if (neigh2.size() >= minPts){
        for (int n : neigh2){
          if (labels[n] == -99 || labels[n] == -1){ labels[n] = clusterId; q.push(n); }
        }
      }
    }
    ++clusterId;
  }
  return labels;
}

} // namespace

struct GenFit2DCHFitter final
  : k4FWCore::Transformer<edm4hep::TrackCollection (const edm4hep::TrackerHit3DCollection&)> {

  using Traits    = Gaudi::Functional::Traits::use_<>;
  using KeyValues = Gaudi::Functional::details::DataHandleMixin<
                      std::tuple<>, std::tuple<>, Traits>::KeyValues;

  GenFit2DCHFitter(const std::string& name, ISvcLocator* svcLoc)
  : Transformer(name, svcLoc,
      std::tuple<KeyValues>{ KeyValues{"inputHits",  std::vector<std::string>{"GGTF_3DHits"}} },
      std::tuple<KeyValues>{ KeyValues{"outputTracks", std::vector<std::string>{"GenFitTracks"}} }) {

    declareProperty("Bz",  m_Bz,  "Magnetic field Bz [Tesla]");
    declareProperty("PDG", m_pdg, "PDG hypothesis (e.g. 13 for mu-)");

    declareProperty("UseMaterialEffects", m_useMatEff,
      "Enable GenFit material effects (MS / dE/dx) using TGeo geometry.");

    // Units / scales (cm internal)
    declareProperty("PositionUnitScale",      m_posScale,
      "Multiply EDM positions by this before passing to GenFit (0.1 converts mm -> cm).");
    declareProperty("InternalLengthToMeters", m_internalLenToM,
      "Conversion from internal length units to meters for pT seeding (0.01 converts cm -> m).");

    // Measurement covariance (anisotropic)
    declareProperty("HitSigmaXYMM", m_hitSigmaXYMM, "Spacepoint sigma in X/Y [mm]");
    declareProperty("HitSigmaZMM",  m_hitSigmaZMM,  "Spacepoint sigma in Z [mm]");

    // Seed / robustness
    declareProperty("SeedPosSigmaMM",  m_seedPosSigmaMM,  "Seed position sigma [mm] (diag)");
    declareProperty("SeedMomSigmaGeV", m_seedMomSigmaGeV, "Seed momentum sigma [GeV] (diag)");
    declareProperty("SeedPTMinGeV",    m_seedPTMinGeV,    "Clamp pT seed min [GeV]");
    declareProperty("SeedPTMaxGeV",    m_seedPTMaxGeV,    "Clamp pT seed max [GeV]");
    declareProperty("SeedPMinGeV",     m_seedPMinGeV,     "Clamp total |p| seed min [GeV]");
    declareProperty("SortHits",        m_sortHits,        "Sort hits along seed direction");
    declareProperty("DeduplicateHits", m_dedup,           "Drop nearly identical consecutive hits");
    declareProperty("DedupTolMM",      m_dedupTolMM,      "Dedup distance tolerance [mm]");

    // Grouping control
    declareProperty("MinGroupSize",          m_minGroupSize,        "Minimum hits per group to attempt a fit");
    declareProperty("UseFallbackClustering", m_useFallbackClust,    "Enable DBSCAN fallback when labels absent");
    declareProperty("FallbackEpsCM",         m_fallbackEpsCM,       "DBSCAN epsilon (cm) for fallback clustering");
    declareProperty("FallbackMinPts",        m_fallbackMinPts,      "DBSCAN minPts for fallback clustering");

    // Retry controls (if base fit yields no FitterInfo)
    declareProperty("RetryIfNoFitterInfo", m_retryIfNoFI,  "Retry with inflated covariances if no FitterInfo");
    declareProperty("RetryMeasInfl",       m_retryMeasInfl, "Measurement variance inflation factor");
    declareProperty("RetrySeedPosInfl",    m_retrySeedPosInfl, "Seed position sigma inflation factor");
    declareProperty("RetrySeedMomInfl",    m_retrySeedMomInfl, "Seed momentum sigma inflation factor");

    // Optional cap on measurements per group
    declareProperty("MaxMeasPerGroup",     m_maxMeasPerGroup, "If >0, downsample measurements to this count");
  }

  // Physics knobs
  Gaudi::Property<double>  m_Bz {this, "Bz", 2.0, "Bz field [T]"};
  Gaudi::Property<int>     m_pdg{this, "PDG", 13,  "PDG hypothesis"};

  // Material effects
  Gaudi::Property<bool>    m_useMatEff{this, "UseMaterialEffects", true,
                                      "Use TGeoMaterialInterface for GenFit MaterialEffects"};

  // Units / scales (defaults chosen to match TGeo's cm)
  Gaudi::Property<double>  m_posScale        {this, "PositionUnitScale", 0.1,  "Multiply positions (0.1: mm->cm)"};
  Gaudi::Property<double>  m_internalLenToM  {this, "InternalLengthToMeters", 0.01, "Length to meters (0.01: cm->m)"};

  // Measurement covariance (anisotropic)
  Gaudi::Property<double>  m_hitSigmaXYMM {this, "HitSigmaXYMM", 0.50, "XY sigma [mm]"};
  Gaudi::Property<double>  m_hitSigmaZMM  {this, "HitSigmaZMM",  2.00, "Z  sigma [mm]"};

  // Seed / robustness
  Gaudi::Property<double>  m_seedPosSigmaMM  {this, "SeedPosSigmaMM", 80.0, "Seed pos sigma [mm]"};
  Gaudi::Property<double>  m_seedMomSigmaGeV {this, "SeedMomSigmaGeV", 5.0,  "Seed mom sigma [GeV]"};
  Gaudi::Property<double>  m_seedPTMinGeV    {this, "SeedPTMinGeV",    0.20, "Min pT [GeV]"};
  Gaudi::Property<double>  m_seedPTMaxGeV    {this, "SeedPTMaxGeV",   200.0, "Max pT [GeV]"};
  Gaudi::Property<double>  m_seedPMinGeV     {this, "SeedPMinGeV",     0.80, "Min |p| [GeV]"};
  Gaudi::Property<bool>    m_sortHits        {this, "SortHits", true,  "Sort hits along seed direction"};
  Gaudi::Property<bool>    m_dedup           {this, "DeduplicateHits", true, "Drop nearly-identical hits"};
  Gaudi::Property<double>  m_dedupTolMM      {this, "DedupTolMM", 0.25, "Dedup tol [mm]"};

  // Grouping
  Gaudi::Property<unsigned> m_minGroupSize     {this, "MinGroupSize", 6u,      "Minimum hits per group"};
  Gaudi::Property<bool>     m_useFallbackClust {this, "UseFallbackClustering", true, "Enable DBSCAN fallback"};
  Gaudi::Property<double>   m_fallbackEpsCM    {this, "FallbackEpsCM", 2.0,     "DBSCAN epsilon in cm"};
  Gaudi::Property<unsigned> m_fallbackMinPts   {this, "FallbackMinPts", 6u,     "DBSCAN minPts"};

  // Retry
  Gaudi::Property<bool>     m_retryIfNoFI      {this, "RetryIfNoFitterInfo", true, "Retry if no FitterInfo"};
  Gaudi::Property<double>   m_retryMeasInfl    {this, "RetryMeasInfl", 4.0, "Measurement variance inflation factor"};
  Gaudi::Property<double>   m_retrySeedPosInfl {this, "RetrySeedPosInfl", 3.0, "Seed position sigma inflation factor"};
  Gaudi::Property<double>   m_retrySeedMomInfl {this, "RetrySeedMomInfl", 3.0, "Seed momentum sigma inflation factor"};

  // Optional cap on measurements
  Gaudi::Property<unsigned> m_maxMeasPerGroup  {this, "MaxMeasPerGroup", 0u, "If >0, downsample measurements per group"};

  StatusCode initialize() override {
    genfit::FieldManager::getInstance()->init(new genfit::ConstField(0., 0., m_Bz.value()));

    if (m_useMatEff.value()) {
      if (!gGeoManager) {
        warning() << "UseMaterialEffects=True but gGeoManager is null; proceeding WITHOUT material effects." << endmsg;
      } else {
        try {
          genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
          info() << "Initialized GenFit MaterialEffects with TGeoMaterialInterface (DD4hep/TGeo)." << endmsg;
        } catch (const std::exception& e) {
          warning() << "Failed to initialize MaterialEffects(TGeo): " << e.what()
                    << " — proceeding WITHOUT material effects." << endmsg;
        }
      }
    } else {
      info() << "UseMaterialEffects = false — running WITHOUT material effects." << endmsg;
    }

    m_fitter = std::make_unique<genfit::KalmanFitterRefTrack>();
    m_fitter->setMaxIterations(12);

    info() << "GenFit2DCHFitter init | Bz=" << m_Bz.value()
           << " | PDG=" << m_pdg.value()
           << " | UseMaterialEffects=" << (m_useMatEff.value() ? "true" : "false")
           << " | posScale=" << m_posScale.value()
           << " | len2m=" << m_internalLenToM.value()
           << " | HitSigmaXY=" << m_hitSigmaXYMM.value() << " mm"
           << " | HitSigmaZ="  << m_hitSigmaZMM.value()  << " mm"
           << " | SeedPosSigma=" << m_seedPosSigmaMM.value() << " mm"
           << " | SeedMomSigma=" << m_seedMomSigmaGeV.value() << " GeV"
           << " | SeedPTMin=" << m_seedPTMinGeV.value() << " | SeedPTMax=" << m_seedPTMaxGeV.value()
           << " | SeedPMin="  << m_seedPMinGeV.value()
           << " | MinGroupSize=" << m_minGroupSize.value()
           << " | UseFallbackClustering=" << (m_useFallbackClust.value() ? "true":"false")
           << " | FallbackEpsCM=" << m_fallbackEpsCM.value()
           << " | FallbackMinPts=" << m_fallbackMinPts.value()
           << " | RetryIfNoFitterInfo=" << (m_retryIfNoFI.value() ? "true":"false")
           << " | RetryMeasInfl=" << m_retryMeasInfl.value()
           << " | RetrySeedPosInfl=" << m_retrySeedPosInfl.value()
           << " | RetrySeedMomInfl=" << m_retrySeedMomInfl.value()
           << " | MaxMeasPerGroup=" << m_maxMeasPerGroup.value()
           << endmsg;
    return StatusCode::SUCCESS;
  }

  // Seed covariance (scaled to internal units)
  TMatrixDSym makeSeedCov(double posScale=1.0, double momScale=1.0) const {
    TMatrixDSym C(6); C.Zero();
    const double sP = (m_seedPosSigmaMM.value() * posScale) * m_posScale.value(); // mm * scale * 0.1 -> cm
    const double sM = (m_seedMomSigmaGeV.value() * momScale);
    for (int i = 0; i < 3; ++i) C(i,i) = sP * sP;
    for (int i = 3; i < 6; ++i) C(i,i) = sM * sM;
    ensurePD(C, 1e-6, 4.0, 6);
    return C;
  }

  edm4hep::TrackCollection operator()(const edm4hep::TrackerHit3DCollection& hits) const override {
    edm4hep::TrackCollection out;
    if (hits.size() < 3u) return out;

    // ---- 1) Group by GGTF label (quality/type); else fallback DBSCAN ----
    std::unordered_map<int, Group> groups;
    bool anyNonZero = false;
    groups.reserve(hits.size()/8 + 1);

    for (size_t i = 0; i < hits.size(); ++i) {
      const auto& h = hits[i];
      const int lbl = hitLabel(h);
      anyNonZero = anyNonZero || (lbl != 0);
      const auto p = h.getPosition();
      TVector3 v(m_posScale.value()*p.x, m_posScale.value()*p.y, m_posScale.value()*p.z);
      auto& g = groups[lbl];
      g.P.push_back(v);
      g.idxs.push_back(i);
    }

    if (!anyNonZero && m_useFallbackClust.value()) {
      std::vector<TVector3> Pall; Pall.reserve(hits.size());
      for (size_t i = 0; i < hits.size(); ++i){
        const auto p = hits[i].getPosition();
        Pall.emplace_back(m_posScale.value()*p.x, m_posScale.value()*p.y, m_posScale.value()*p.z);
      }
      const auto labels = dbscan(Pall, m_fallbackEpsCM.value(), m_fallbackMinPts.value());
      std::unordered_map<int, Group> g2;
      int maxLbl = -1; for (int L : labels) maxLbl = std::max(maxLbl, L);
      if (maxLbl >= 0) {
        for (size_t i=0;i<labels.size();++i){
          int L = labels[i]; if (L < 0) continue;
          auto& g = g2[L]; g.P.push_back(Pall[i]); g.idxs.push_back(i);
        }
        groups.swap(g2);
      }
    }

    info() << "GF2: groups=" << groups.size()
           << " (label=quality/type or DBSCAN; 0 means unlabeled)"
           << endmsg;

    // --- 2) Measurement covariance (anisotropic, cm)
    const double sxy = m_hitSigmaXYMM.value() * m_posScale.value();
    const double sz  = m_hitSigmaZMM.value()  * m_posScale.value();
    TMatrixDSym Cmeas(3); Cmeas.Zero();
    Cmeas(0,0)=sxy*sxy; Cmeas(1,1)=sxy*sxy; Cmeas(2,2)=sz*sz;
    ensurePD(Cmeas, 1e-6, 4.0, 6);

    // ---- 3) Fit each group independently
    for (auto& kv : groups) {
      const int label = kv.first;
      auto& P    = kv.second.P;
      auto& idxs = kv.second.idxs;

      // optional dedup (consecutive)
      if (m_dedup.value() && !P.empty()) {
        const double tol2 = std::pow(m_dedupTolMM.value() * m_posScale.value(), 2);
        std::vector<TVector3> P2; P2.reserve(P.size());
        std::vector<size_t>   I2; I2.reserve(idxs.size());
        TVector3 prev(1e99,1e99,1e99);
        for (size_t k = 0; k < P.size(); ++k) {
          if ((P[k] - prev).Mag2() < tol2) continue;
          P2.push_back(P[k]); I2.push_back(idxs[k]); prev = P[k];
        }
        P.swap(P2); idxs.swap(I2);
      }

      if (P.size() < m_minGroupSize.value()) continue;

      // Initial direction from endpoints
      TVector3 pos0 = P.front();
      TVector3 posL = P.back();
      TVector3 dir  = (posL - pos0);
      if (dir.Mag2() < 1e-12) dir = TVector3(0,0,1);
      dir = dir.Unit();

      // Sort along seed direction (carry indices)
      if (m_sortHits.value()) {
        std::vector<size_t> order(P.size()); std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b){
          return projS(pos0, dir, P[a]) < projS(pos0, dir, P[b]);
        });
        std::vector<TVector3> Psorted; Psorted.reserve(P.size());
        std::vector<size_t>   Isorted; Isorted.reserve(idxs.size());
        for (auto k : order) { Psorted.push_back(P[k]); Isorted.push_back(idxs[k]); }
        P.swap(Psorted); idxs.swap(Isorted);
        pos0 = P.front(); posL = P.back();
        dir  = (posL - pos0).Mag2() < 1e-12 ? TVector3(0,0,1) : (posL - pos0).Unit();
      }

      // Optional cap on measurements per group (downsample evenly)
      if (m_maxMeasPerGroup.value() > 0 && P.size() > m_maxMeasPerGroup.value()) {
        const unsigned cap = m_maxMeasPerGroup.value();
        const size_t N = P.size();
        std::vector<TVector3> P2; P2.reserve(cap);
        std::vector<size_t>   I2; I2.reserve(cap);
        for (size_t j = 0; j < cap; ++j) {
          size_t k = static_cast<size_t>( std::round( double(j) * (N-1) / double(cap-1) ) );
          P2.push_back(P[k]); I2.push_back(idxs[k]);
        }
        P.swap(P2); idxs.swap(I2);
      }

      // Seed pT from XY curvature
      const TVector3 A = P.front();
      const TVector3 B = P[P.size()/2];
      const TVector3 C = P.back();
      const double R_int = circleRadiusXY(A,B,C);            // cm
      const double R_m   = R_int * m_internalLenToM.value(); // m
      double pTseed = 0.3 * m_Bz.value() * R_m;              // GeV
      if (!std::isfinite(pTseed) || pTseed <= 0) pTseed = 1.0;
      pTseed = std::clamp(pTseed, m_seedPTMinGeV.value(), m_seedPTMaxGeV.value());

      const double cosTheta = std::clamp(std::abs(dir.Z()), 0.0, 0.999999);
      const double sinTheta = std::sqrt(std::max(1e-4, 1.0 - cosTheta*cosTheta));
      double pMag = pTseed / sinTheta;
      pMag = std::max(pMag, m_seedPMinGeV.value());
      pMag = std::clamp(pMag, m_seedPTMinGeV.value(), m_seedPTMaxGeV.value());
      TVector3 mom0 = pMag * dir;

      info() << "Seed[label=" << label << "] | n=" << P.size()
             << " R_int(cm)=" << R_int
             << " pTseed(GeV)=" << pTseed
             << " sinTheta=" << sinTheta
             << " pMag(GeV)=" << pMag << endmsg;

      // Build track (base)
      auto rep_up = std::make_unique<genfit::RKTrackRep>(m_pdg.value());
      genfit::Track fitTrack(rep_up.release(), pos0, mom0);
      fitTrack.setCovSeed(makeSeedCov(/*posScale*/1.0, /*momScale*/1.0));

      // Add measurements
      int detId = label; int hitId = 0;
      for (const auto& v : P) {
        TVectorD pos(3); pos[0]=v.X(); pos[1]=v.Y(); pos[2]=v.Z();
        auto* tp   = new genfit::TrackPoint(&fitTrack);
        auto* meas = new genfit::SpacepointMeasurement(pos, Cmeas, detId, hitId, tp, false, false);
        tp->addRawMeasurement(meas);
        fitTrack.insertPoint(tp);
        ++hitId;
      }
      if (fitTrack.getNumPointsWithMeasurement() == 0) continue;

      auto run_fit = [&](genfit::Track& trk, const char* tag)->const genfit::FitStatus*{
        try {
          m_fitter->processTrack(&trk);
        } catch (const genfit::Exception& e) {
          warning() << "GenFit exception during processTrack (label=" << label << ", " << tag << "): " << e.what() << endmsg;
          return nullptr;
        } catch (const std::exception& e) {
          warning() << "Std exception during processTrack (label=" << label << ", " << tag << "): " << e.what() << endmsg;
          return nullptr;
        } catch (...) {
          warning() << "Unknown exception during processTrack (label=" << label << ", " << tag << ")" << endmsg;
          return nullptr;
        }
        genfit::AbsTrackRep* rep = trk.getCardinalRep();
        return rep ? trk.getFitStatus(rep) : nullptr;
      };

      const genfit::FitStatus* fs_base = run_fit(fitTrack, "base");
      if (fs_base) {
        info() << "FitStatus[label=" << label << "][base] | fitted=" << fs_base->isFitted()
               << " converged=" << fs_base->isFitConverged()
               << " Ndf=" << fs_base->getNdf()
               << " Chi2=" << fs_base->getChi2()
               << " pVal=" << fs_base->getPVal() << endmsg;
      }

      // Count FitterInfo presence
      auto has_any_FI = [&](genfit::Track& trk)->bool{
        genfit::AbsTrackRep* rep = trk.getCardinalRep();
        if (!rep) return false;
        const size_t nPts = trk.getNumPoints();
        for (size_t i=0;i<nPts;++i) {
          auto* tp = trk.getPoint(i);
          if (tp && tp->hasFitterInfo(rep)) return true;
        }
        return false;
      };

      bool ok = (fs_base && has_any_FI(fitTrack));

      // Retry path if requested
      if (!ok && m_retryIfNoFI.value()) {
        warning() << "No TrackPoint has FitterInfo (label=" << label << ", base) — will consider retry." << endmsg;

        // Rebuild a track with inflated covariances
        auto rep2_up = std::make_unique<genfit::RKTrackRep>(m_pdg.value());
        genfit::Track fitTrack2(rep2_up.release(), pos0, mom0);

        // Seed cov inflation
        fitTrack2.setCovSeed(makeSeedCov(/*posScale*/m_retrySeedPosInfl.value(),
                                         /*momScale*/m_retrySeedMomInfl.value()));

        // Measurement inflation
        TMatrixDSym Cmeas_retry = Cmeas;
        for (int i=0;i<3;++i) Cmeas_retry(i,i) *= m_retryMeasInfl.value();
        ensurePD(Cmeas_retry, 1e-6, 4.0, 6);

        int det2 = label, hid2 = 0;
        for (const auto& v : P) {
          TVectorD pos(3); pos[0]=v.X(); pos[1]=v.Y(); pos[2]=v.Z();
          auto* tp   = new genfit::TrackPoint(&fitTrack2);
          auto* meas = new genfit::SpacepointMeasurement(pos, Cmeas_retry, det2, hid2, tp, false, false);
          tp->addRawMeasurement(meas);
          fitTrack2.insertPoint(tp);
          ++hid2;
        }

        const genfit::FitStatus* fs_retry = run_fit(fitTrack2, "retry");
        if (fs_retry) {
          info() << "FitStatus[label=" << label << "][retry] | fitted=" << fs_retry->isFitted()
                 << " converged=" << fs_retry->isFitConverged()
                 << " Ndf=" << fs_retry->getNdf()
                 << " Chi2=" << fs_retry->getChi2()
                 << " pVal=" << fs_retry->getPVal() << endmsg;
        }

        if (fs_retry && has_any_FI(fitTrack2)) {
          // Export this one
          auto trk = out.create();
          trk.setType(m_pdg.value());
          try { trk.setChi2(fs_retry->getChi2()); trk.setNdf(fs_retry->getNdf()); } catch (...) {}
          for (auto i : idxs) trk.addToTrackerHits(hits[i]);
          continue;
        } else {
          warning() << "No TrackPoint has FitterInfo (label=" << label << ", retry) — will consider retry." << endmsg;
        }
      }

      // Base export (if base had FI)
      if (ok) {
        auto trk = out.create();
        trk.setType(m_pdg.value());
        if (fs_base) { try { trk.setChi2(fs_base->getChi2()); trk.setNdf(fs_base->getNdf()); } catch (...) {} }
        for (auto i : idxs) trk.addToTrackerHits(hits[i]);
      }
    } // end group loop

    return out;
  }

  StatusCode finalize() override { return StatusCode::SUCCESS; }

private:
  std::unique_ptr<genfit::KalmanFitterRefTrack> m_fitter;
};

DECLARE_COMPONENT(GenFit2DCHFitter)
