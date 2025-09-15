// fit_ggtf_tracks.cpp — Podio Frames (ROOTReader) backend

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

// Podio Frames
#include "podio/ROOTReader.h"
#include "podio/Frame.h"

// extension types
#include "extension/TrackCollection.h"
#include "extension/TrackerHit.h"

// ROOT
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

// ---- traits to extract (x,y,z) from a hit
template <typename T, typename = void>
struct HasGetPosition : std::false_type {};
template <typename T>
struct HasGetPosition<T, std::void_t<decltype(std::declval<const T&>().getPosition())>>
  : std::true_type {};

template <typename T, typename = void>
struct HasXYZ : std::false_type {};
template <typename T>
struct HasXYZ<T, std::void_t<
  decltype(std::declval<const T&>().getX()),
  decltype(std::declval<const T&>().getY()),
  decltype(std::declval<const T&>().getZ())>> : std::true_type {};

template <typename HitT>
inline bool extractXYZ(const HitT& h, double& x, double& y, double& z) {
  if constexpr (HasGetPosition<HitT>::value) {
    const auto p = h.getPosition();
    x = static_cast<double>(p.x);
    y = static_cast<double>(p.y);
    z = static_cast<double>(p.z);
    return true;
  } else if constexpr (HasXYZ<HitT>::value) {
    x = static_cast<double>(h.getX());
    y = static_cast<double>(h.getY());
    z = static_cast<double>(h.getZ());
    return true;
  } else {
    return false;
  }
}

// ---- simple Kåsa circle fit in (x,y)
struct CircleResult { double R{NAN}, xc{NAN}, yc{NAN}, chi2{NAN}; bool ok{false}; };

static CircleResult fitCircleKasa(const std::vector<std::array<double,3>>& xyz) {
  const size_t n = xyz.size();
  CircleResult res;
  if (n < 3) return res;

  long double Sx=0,Sy=0,Sxx=0,Syy=0,Sxy=0,Szz=0,Sxz=0,Syz=0;
  for (auto& p : xyz) {
    const long double x=p[0], y=p[1], z=x*x+y*y;
    Sx+=x; Sy+=y; Sxx+=x*x; Syy+=y*y; Sxy+=x*y; Szz+=z; Sxz+=x*z; Syz+=y*z;
  }
  const long double A11=Sxx, A12=Sxy, A13=Sx;
  const long double A21=Sxy, A22=Syy, A23=Sy;
  const long double A31=Sx,  A32=Sy,  A33=n;
  const long double B1=Sxz,  B2=Syz,  B3=Szz;

  const long double detA =
      A11*(A22*A33 - A23*A32) -
      A12*(A21*A33 - A23*A31) +
      A13*(A21*A32 - A22*A31);
  if (std::fabsl(detA) < 1e-12L) return res;

  auto det = [](long double a11,long double a12,long double a13,
                long double a21,long double a22,long double a23,
                long double a31,long double a32,long double a33){
    return a11*(a22*a33 - a23*a32)
         - a12*(a21*a33 - a23*a31)
         + a13*(a21*a32 - a22*a31);
  };
  const long double a = det(B1,A12,A13, B2,A22,A23, B3,A32,A33)/detA;
  const long double b = det(A11,B1,A13, A21,B2,A23, A31,B3,A33)/detA;
  const long double c = det(A11,A12,B1, A21,A22,B2, A31,A32,B3)/detA;

  const long double xc=a/2, yc=b/2, R2=c+xc*xc+yc*yc;
  if (R2 <= 0) return res;
  const double R = std::sqrt(static_cast<double>(R2));

  long double chi2=0;
  for (auto& p : xyz) {
    const long double dx=p[0]-xc, dy=p[1]-yc, r=std::sqrt(dx*dx+dy*dy);
    const long double dr=r-R; chi2+=dr*dr;
  }
  res.R=R; res.xc=static_cast<double>(xc); res.yc=static_cast<double>(yc);
  res.chi2=static_cast<double>(chi2); res.ok=true; return res;
}

// ---- z vs. arc-length s = R * unwrapped phi
static std::pair<double,double> fitZvsS(const std::vector<std::array<double,3>>& xyz,
                                        double xc, double yc, double Rm) {
  const size_t n = xyz.size();
  if (n < 2 || !(Rm > 0)) return {NAN, NAN};

  std::vector<double> phi(n), z(n);
  for (size_t i=0;i<n;++i) {
    const double dx = xyz[i][0]-xc, dy = xyz[i][1]-yc;
    phi[i] = std::atan2(dy, dx);
    z[i]   = xyz[i][2];
  }
  std::vector<size_t> idx(n); std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a,size_t b){ return phi[a]<phi[b]; });

  std::vector<double> phiu(n), zu(n);
  double last=0;
  for (size_t k=0;k<n;++k) {
    const double p = phi[idx[k]];
    if (k==0) phiu[k]=p;
    else {
      double dp = p - phi[idx[k-1]];
      if (dp >  M_PI) dp -= 2*M_PI;
      if (dp < -M_PI) dp += 2*M_PI;
      phiu[k] = last + dp;
    }
    last = phiu[k];
    zu[k] = z[idx[k]];
  }

  std::vector<double> s(n);
  for (size_t k=0;k<n;++k) s[k] = Rm * phiu[k];

  double S=0, Sz=0, Ss=0, Sss=0, Ssz=0;
  for (size_t k=0;k<n;++k) { S+=1; Sz+=zu[k]; Ss+=s[k]; Sss+=s[k]*s[k]; Ssz+=s[k]*zu[k]; }
  const double denom = (S*Sss - Ss*Ss);
  if (std::fabs(denom) < 1e-12) return {NAN, NAN};
  const double m  = (S*Ssz - Ss*Sz)/denom;
  // const double z0 = (Sz - m*Ss)/S; // not currently used

  double chi2=0;
  for (size_t k=0;k<n;++k) { const double dz = zu[k] - (/*z0 +*/ m*s[k]); chi2 += dz*dz; }
  return {m, chi2};
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
              << " <in.root> <out.root> [--coll CDCHTracks] [--B 2.0]\n";
    return 1;
  }
  std::string inFile  = argv[1];
  std::string outFile = argv[2];
  std::string coll    = "CDCHTracks";
  double      B       = 2.0; // Tesla

  for (int i=3;i<argc;++i) {
    std::string a = argv[i];
    if (a=="--coll" && i+1<argc) { coll=argv[++i]; continue; }
    if (a=="--B"    && i+1<argc) { B=std::stod(argv[++i]); continue; }
  }

  // Podio ROOTReader (Frames): open file, then loop entries by index
  podio::ROOTReader rrd;
  rrd.openFile(inFile);
  const auto nEntries = rrd.getEntries("events");

  // ROOT outputs
  TFile fout(outFile.c_str(), "RECREATE");
  TH1F  h_pt("h_pt","p_{T} [GeV];p_{T} [GeV];Entries",100,0,100);
  TH1F  h_tanl("h_tanl","tan#lambda;tan#lambda;Entries",100,-5,5);
  TTree t("fit","fit results");
  int iev=0, itrk=0, nhit=0; double R=0, pt=0, tanl=0, chi2_xy=0, chi2_z=0;
  t.Branch("iev",&iev,"iev/I"); t.Branch("itrk",&itrk,"itrk/I");
  t.Branch("nhit",&nhit,"nhit/I"); t.Branch("R",&R,"R/D");
  t.Branch("pt",&pt,"pt/D"); t.Branch("tanl",&tanl,"tanl/D");
  t.Branch("chi2_xy",&chi2_xy,"chi2_xy/D"); t.Branch("chi2_z",&chi2_z,"chi2_z/D");

  size_t n_used=0;
  for (size_t i=0; i<nEntries; ++i) {
    iev = static_cast<int>(i);
    // readEntry -> unique_ptr<ROOTFrameData>; move-construct podio::Frame
    auto fdata = rrd.readEntry("events", i);
    podio::Frame f{std::move(fdata)};

    // presence check
    const auto names = f.getAvailableCollections();
    if (std::find(names.begin(), names.end(), coll) == names.end()) continue;

    const auto& trks = f.get<extension::TrackCollection>(coll);
    if (trks.size()==0) continue;

    for (size_t it=0; it<trks.size(); ++it) {
      itrk = static_cast<int>(it);
      const auto trk = trks.at(it);

      std::vector<std::array<double,3>> xyz; xyz.reserve(64);
      for (auto hIt = trk.trackerHits_begin(); hIt != trk.trackerHits_end(); ++hIt) {
        const auto& h = *hIt;
        double x,y,z;
        if (extractXYZ(h,x,y,z)) xyz.push_back({x,y,z});
      }
      nhit = static_cast<int>(xyz.size());
      if (nhit < 3) continue;

      const auto c = fitCircleKasa(xyz);
      if (!c.ok || !(c.R>0)) continue;
      R = c.R;

      const double Rm = R * 1e-3;      // mm -> m
      pt      = 0.3 * B * Rm;          // GeV
      chi2_xy = c.chi2;

      const auto zf = fitZvsS(xyz, c.xc, c.yc, Rm);
      tanl   = zf.first;
      chi2_z = zf.second;

      h_pt.Fill(pt);
      h_tanl.Fill(tanl);
      t.Fill();
      ++n_used;
    }
  }

  std::cout << "Finished. Tracks used: " << n_used << "\n";
  fout.Write(); fout.Close();
  return 0;
}
