import os
import subprocess
import shutil

from Gaudi.Configuration import INFO, DEBUG
from Gaudi.Configuration import ApplicationMgr as GaudiApp
from Configurables import EventDataSvc, UniqueIDGenSvc, RndmGenSvc, GeoSvc
from Configurables import AuditorSvc, ChronoAuditor, MemoryAuditor, MessageSvc

# Preload factories (ensures Gaudi registries are in memory)
import ROOT
ROOT.gSystem.Load("libTracking.so")
ROOT.gSystem.Load("libDCHdigi.so")

from k4FWCore import IOSvc
from k4FWCore.parseArgs import parser

# ----------------- CLI -----------------
parser.add_argument("--inputFile",  default="ddsim_output_edm4hep.root")
parser.add_argument("--outputFile", default="output_digi_tracks.root")
parser.add_argument("--modelPath",  default="", help="ONNX model path or .md5/http/root URL")
parser.add_argument("--tbeta",      type=float, default=0.05)
parser.add_argument("--td",         type=float, default=0.05)
parser.add_argument("--dchSimHits", default="DCHCollection")
parser.add_argument("--compactXML", default="")
parser.add_argument("--dchName",    default="DCH_v2")

# Tracking toggles
parser.add_argument("--produce3DHits", action="store_true", default=False)
parser.add_argument("--ggtfLog", choices=["INFO", "DEBUG"], default="INFO")
parser.add_argument("--maxHitsPerEvent", type=int, default=0)

# Stage control
parser.add_argument("--stage", choices=["digi", "ggtf", "fit"], default="ggtf")

# Fitter toggles
parser.add_argument("--fitter", choices=["none","genfit2"], default="none")
parser.add_argument("--fitOut", default="GenFitTracks")
parser.add_argument("--fitterLog", choices=["INFO","DEBUG"], default="INFO")

# ----------------- GenFit runtime tuning -----------------
# Booleans (material effects, sorting, dedup) with on/off pairs for convenience
parser.add_argument("--gf-useMat", dest="gf_useMat", action="store_true", default=True)
parser.add_argument("--no-gf-useMat", dest="gf_useMat", action="store_false")
parser.add_argument("--gf-sortHits", dest="gf_sortHits", action="store_true", default=True)
parser.add_argument("--no-gf-sortHits", dest="gf_sortHits", action="store_false")
parser.add_argument("--gf-dedup", dest="gf_dedup", action="store_true", default=True)
parser.add_argument("--no-gf-dedup", dest="gf_dedup", action="store_false")

# Units / scales  (use cm internally)
parser.add_argument("--gf-posScale", type=float, default=0.1)
parser.add_argument("--gf-len2m", type=float, default=0.01)

# Dedup tolerance
parser.add_argument("--gf-dedupTol", type=float, default=0.10)

# Measurement covariances (anisotropic)
parser.add_argument("--gf-hitSigmaXY", type=float, default=0.60)
parser.add_argument("--gf-hitSigmaZ",  type=float, default=3.00)

# Seed covariances and clamps
parser.add_argument("--gf-seedPosSigma", type=float, default=30.0)
parser.add_argument("--gf-seedMomSigma", type=float, default=2.0)
parser.add_argument("--gf-seedPTMin", type=float, default=0.30)
parser.add_argument("--gf-seedPTMax", type=float, default=50.0)
# NEW: clamp on total |p|
parser.add_argument("--gf-seedPMin", type=float, default=0.80)

# Grouping & fallback clustering (for fitter)
parser.add_argument("--gf-minGroup",       type=int,   default=6)
parser.add_argument("--gf-useFallback",    dest="gf_useFallback", action="store_true",  default=True)
parser.add_argument("--no-gf-useFallback", dest="gf_useFallback", action="store_false")
parser.add_argument("--gf-fallbackEpsCM",  type=float, default=2.0)
parser.add_argument("--gf-fallbackMinPts", type=int,   default=6)

# Retry controls (if base fit has no FitterInfo)
parser.add_argument("--gf-retry",           dest="gf_retry", action="store_true",  default=True)
parser.add_argument("--no-gf-retry",        dest="gf_retry", action="store_false")
parser.add_argument("--gf-retryMeasInfl",   type=float, default=4.0)
parser.add_argument("--gf-retrySeedPos",    type=float, default=3.0)
parser.add_argument("--gf-retrySeedMom",    type=float, default=3.0)

# Optional: cap number of measurements per group
parser.add_argument("--gf-maxMeasPerGroup", type=int, default=0)

args = parser.parse_args()

# ----------------- Message/Auditing -----------------
MessageSvc().Format = "% F%18W%S%7W%R%T %0W%M"

GaudiApp().PluginDebugLevel = 1
GaudiApp().AuditAlgorithms = True
try:
    GaudiApp().AuditTools = True
except Exception:
    pass

AuditorSvc().Auditors = [ ChronoAuditor(), MemoryAuditor() ]

# ----------------- IO -----------------
svc = IOSvc("IOSvc")
svc.Input  = args.inputFile
svc.Output = args.outputFile
print(f"[IO] input={args.inputFile}  output={args.outputFile}")

# ----------------- Geometry -----------------
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [args.compactXML] if args.compactXML else []
geoservice.OutputLevel = INFO
print(f"[Geo] compactXML={args.compactXML}  DD4hep_XMLPATH={os.environ.get('DD4hep_XMLPATH','')}")

# --------- Model staging helper ----------
def stage_model(spec: str) -> str:
    if not spec:
        raise RuntimeError("modelPath is empty")
    if spec.endswith(".onnx") and os.path.exists(spec):
        return os.path.abspath(spec)
    out = os.path.abspath("model.onnx")
    if spec.endswith(".onnx.md5"):
        with open(spec) as f: md5 = f.read().split()[0]
        url = f"https://key4hep.web.cern.ch/testFiles/k4RecTracker/{md5}"
        print(f"[model] from md5: {md5} -> {url}")
        subprocess.run(["wget", "--no-verbose", "--timeout=180", "--tries=2", "-O", out, url], check=True)
        if not os.path.exists(out) or os.path.getsize(out) == 0:
            raise RuntimeError(f"Downloaded model is missing/empty: {out}")
        return out
    if spec.startswith(("http://","https://")):
        print(f"[model] download {spec} -> {out}")
        subprocess.run(["wget", "--no-verbose", "--timeout=180", "--tries=2", "-O", out, spec], check=True)
        if not os.path.exists(out) or os.path.getsize(out) == 0:
            raise RuntimeError(f"Downloaded model is missing/empty: {out}")
        return out
    if spec.startswith("root://"):
        print(f"[model] xrdcp {spec} -> {out}")
        subprocess.run(["xrdcp", "-f", spec, out], check=True)
        if not os.path.exists(out) or os.path.getsize(out) == 0:
            raise RuntimeError(f"xrdcp model is missing/empty: {out}")
        return out
    if spec.endswith(".onnx"):
        shutil.copy2(spec, out)
        return out
    raise RuntimeError(f"Unrecognized model spec: {spec}")

# ----------------- DCH Digitizer -----------------
from Configurables import DCHdigi_v01
dch_digitizer = DCHdigi_v01(
    "DCHdigi",
    DCH_simhits=[args.dchSimHits],
    DCH_name=args.dchName,
    fileDataAlg="DataAlgFORGEANT.root",
    calculate_dndx=False,
    create_debug_histograms=False,
    zResolution_mm=30.0,
    xyResolution_mm=0.1
)

# ----------------- Track Finder (GGTF) -----------------
try:
    from TrackingConf import GGTF_tracking
except Exception:
    from Configurables import GGTF_tracking  # fallback

GGTF = GGTF_tracking(
    "GGTF_tracking",
    inputWireHits=["DCH_DigiCollection"],
    inputPlanarHits=[],
    outputTracks=["CDCHTracks"],
    output3DHits=["GGTF_3DHits"],
    OutputLevel=INFO,
)

GGTF.modelPath = stage_model(args.modelPath)
GGTF.tbeta     = args.tbeta
GGTF.td        = args.td
try:
    GGTF.produce3DHits = bool(args.produce3DHits)
except Exception:
    pass
try:
    if int(args.maxHitsPerEvent) > 0:
        GGTF.maxHitsPerEvent = int(args.maxHitsPerEvent)
except Exception:
    print("[warn] GGTF.maxHitsPerEvent property not present; ignoring.")
GGTF.OutputLevel = DEBUG if args.ggtfLog == "DEBUG" else INFO

print(f"[GGTF] stage={args.stage} modelPath={GGTF.modelPath} tbeta={GGTF.tbeta} td={GGTF.td} "
      f"produce3DHits={getattr(GGTF, 'produce3DHits', 'n/a')} "
      f"maxHitsPerEvent={getattr(GGTF, 'maxHitsPerEvent', 0)} "
      f"log={args.ggtfLog}")

# ----------------- Ensure cluster-size file for DCHdigi -----------------
cluster_file = "DataAlgFORGEANT.root"
if not os.path.exists(cluster_file):
    url = "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/IDEA/DataAlgFORGEANT.root"
    print(f"[setup] Fetching {cluster_file} from {url}")
    subprocess.run(["wget", "--no-verbose", "--timeout=180", "--tries=2", "--no-clobber", url], check=True)

# ----- GenFit services (field + material) -----
field_svc_name = None
material_svc_name = None
field_svc_obj = None
material_svc_obj = None

try:
    from Configurables import DD4hepFieldSvc
    field_svc_obj = DD4hepFieldSvc("GenFitFieldSvc", GeoSvcName="GeoSvc")
    field_svc_name = field_svc_obj.getName()
    print(f"[genfit] Using DD4hepFieldSvc -> {field_svc_name}")
except Exception as e:
    print(f"[genfit] DD4hepFieldSvc not available: {e}")
    try:
        from Configurables import ConstBFieldSvc
        field_svc_obj = ConstBFieldSvc("GenFitFieldSvc", Bz=args.gf_bz)
        field_svc_name = field_svc_obj.getName()
        print(f"[genfit] Using ConstBFieldSvc Bz={args.gf_bz}T -> {field_svc_name}")
    except Exception as e2:
        print(f"[genfit] No field service configured: {e2}")

try:
    from Configurables import DD4hepMaterialSvc
    material_svc_obj = DD4hepMaterialSvc("GenFitMaterialSvc", GeoSvcName="GeoSvc")
    material_svc_name = material_svc_obj.getName()
    print(f"[genfit] Using DD4hepMaterialSvc -> {material_svc_name}")
except Exception as e:
    print(f"[genfit] DD4hepMaterialSvc not available: {e}")
    try:
        from Configurables import GenFitMaterialSvc
        material_svc_obj = GenFitMaterialSvc("GenFitMaterialSvc", GeoSvcName="GeoSvc")
        material_svc_name = material_svc_obj.getName()
        print(f"[genfit] Using GenFitMaterialSvc -> {material_svc_name}")
    except Exception as e2:
        print(f"[genfit] No material service configured: {e2}")

# ----------------- Optional fitter stage -----------------
def _set_if_has(obj, name, value):
    try:
        if hasattr(obj, name):
            setattr(obj, name, value)
            print(f"[fitter] set {name} = {value}")
            return True
    except Exception as e:
        print(f"[fitter] could not set {name}: {e}")
    return False

fitter_alg = None
if args.stage == "fit" and args.fitter == "genfit2":
    try:
        if not bool(getattr(GGTF, "produce3DHits", False)):
            GGTF.produce3DHits = True
            print("[fitter] Enabling GGTF.produce3DHits = True (required by GenFit2).")
    except Exception:
        pass

    try:
        from Configurables import GenFit2DCHFitter
        fitter_alg = GenFit2DCHFitter("GenFit2DCHFitter")
        fitter_alg.OutputLevel = DEBUG if args.fitterLog == "DEBUG" else INFO

        # Accept either property style (string vs list) across builds
        for prop, val in (("input3DHits", "GGTF_3DHits"), ("inputHits", ["GGTF_3DHits"])):
            if _set_if_has(fitter_alg, prop, val): break
        for prop, val in (("outputTracks", args.fitOut), ("outputTracks", [args.fitOut])):
            if _set_if_has(fitter_alg, prop, val): break

        # Wire optional services if supported
        for prop, val in [
            ("FieldSvc", field_svc_name),
            ("BFieldSvc", field_svc_name),
            ("MaterialSvc", material_svc_name),
            ("MaterialEffectsSvc", material_svc_name),
            ("GeoSvcName", "GeoSvc"),
        ]:
            if val is not None:
                _set_if_has(fitter_alg, prop, val)

        # Physics hypothesis / field
        _set_if_has(fitter_alg, "Bz", args.gf_bz)
        _set_if_has(fitter_alg, "PDG", args.gf_pdg)

        # Material effects & units
        _set_if_has(fitter_alg, "UseMaterialEffects", args.gf_useMat)
        _set_if_has(fitter_alg, "PositionUnitScale", args.gf_posScale)
        _set_if_has(fitter_alg, "InternalLengthToMeters", args.gf_len2m)

        # Measurement & seed
        _set_if_has(fitter_alg, "HitSigmaXYMM", args.gf_hitSigmaXY)
        _set_if_has(fitter_alg, "HitSigmaZMM",  args.gf_hitSigmaZ)
        _set_if_has(fitter_alg, "SeedPosSigmaMM",  args.gf_seedPosSigma)
        _set_if_has(fitter_alg, "SeedMomSigmaGeV", args.gf_seedMomSigma)
        _set_if_has(fitter_alg, "SeedPTMinGeV",    args.gf_seedPTMin)
        _set_if_has(fitter_alg, "SeedPTMaxGeV",    args.gf_seedPTMax)
        _set_if_has(fitter_alg, "SeedPMinGeV",     args.gf_seedPMin)

        # Sorting / dedup
        _set_if_has(fitter_alg, "SortHits",        args.gf_sortHits)
        _set_if_has(fitter_alg, "DeduplicateHits", args.gf_dedup)
        _set_if_has(fitter_alg, "DedupTolMM",      args.gf_dedupTol)

        # Grouping & fallback clustering
        _set_if_has(fitter_alg, "MinGroupSize",          int(args.gf_minGroup))
        _set_if_has(fitter_alg, "UseFallbackClustering", bool(args.gf_useFallback))
        _set_if_has(fitter_alg, "FallbackEpsCM",         float(args.gf_fallbackEpsCM))
        _set_if_has(fitter_alg, "FallbackMinPts",        int(args.gf_fallbackMinPts))

        # Retry controls
        _set_if_has(fitter_alg, "RetryIfNoFitterInfo", bool(args.gf_retry))
        _set_if_has(fitter_alg, "RetryMeasInfl",       float(args.gf_retryMeasInfl))
        _set_if_has(fitter_alg, "RetrySeedPosInfl",    float(args.gf_retrySeedPos))
        _set_if_has(fitter_alg, "RetrySeedMomInfl",    float(args.gf_retrySeedMom))

        # Optional cap on measurements per group
        _set_if_has(fitter_alg, "MaxMeasPerGroup",     int(args.gf_maxMeasPerGroup))

        # Friendly catch-alls (exist in some builds)
        for name, val in [("pdgHypothesis", args.gf_pdg), ("minHitsOnTrack", 4), ("maxChi2", 1e6)]:
            _set_if_has(fitter_alg, name, val)

        print(f"[fitter] GenFit2DCHFitter configured; output -> '{args.fitOut}' | log={args.fitterLog}")
    except Exception as e:
        print(f"[fitter] GenFit2DCHFitter not found/failed to configure: {e}")
        print("[fitter] Continuing WITHOUT fitter.")
        fitter_alg = None

# --- Optional: collection size probe
try:
    from Configurables import EDM4hepCollectionSizePrinter as SizePrinter
    size_printer = SizePrinter("SizePrinter", CollectionsToPrint=["GGTF_3DHits", args.fitOut], OutputLevel=INFO)
except Exception:
    size_printer = None

# ----------------- AppMgr / Stage selection -----------------
top_algs = []
if args.stage in ("digi", "ggtf", "fit"):
    top_algs.append(dch_digitizer)
if args.stage in ("ggtf", "fit"):
    top_algs.append(GGTF)
if args.stage == "fit" and fitter_alg is not None:
    top_algs.append(fitter_alg)
if size_printer is not None:
    top_algs.append(size_printer)

# Try reader/writer if available
try:
    from k4FWCore import getReader, getWriter
    reader = getReader()
    writer = getWriter()
    top_algs = [reader] + top_algs + [writer]
except Exception:
    pass

print(f"[pipeline] TopAlg order: {[alg.getFullName() for alg in top_algs]}")

# ---- Ext services: include field/material if created ----
ext_svcs = [geoservice, EventDataSvc("EventDataSvc"), UniqueIDGenSvc("uidSvc"), RndmGenSvc()]
if field_svc_obj is not None:
    ext_svcs.append(field_svc_obj)
if material_svc_obj is not None:
    ext_svcs.append(material_svc_obj)

mgr = GaudiApp(TopAlg=top_algs, EvtSel="NONE", EvtMax=-1, ExtSvc=ext_svcs, OutputLevel=DEBUG)
