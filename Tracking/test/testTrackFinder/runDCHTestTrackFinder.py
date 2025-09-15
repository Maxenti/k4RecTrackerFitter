import os
import math
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
parser.add_argument("--inputFile",  default="ddsim_output_edm4hep.root",
                    help="Input EDM4hep file (from ddsim or digi)")
parser.add_argument("--outputFile", default="output_digi_tracks.root",
                    help="Output EDM4hep file")
parser.add_argument("--modelPath",  default="",
                    help="ONNX model path for GGTF_tracking; accepts .onnx, .onnx.md5, http(s)://, root://")
parser.add_argument("--tbeta",      type=float, default=0.05,
                    help="GGTF beta threshold")
parser.add_argument("--td",         type=float, default=0.05,
                    help="GGTF distance threshold")
parser.add_argument("--dchSimHits", default="DCHCollection",
                    help="Name of DCH sim-hit collection in the input file")
parser.add_argument("--compactXML", default="",
                    help="Path or URL to compact XML to load in GeoSvc (use the same one passed to ddsim)")
parser.add_argument("--dchName",    default="DCH_v2",
                    help="DD4hep detector name for the DCH (e.g. DCH_v2, CDCH, DCH)")

# Tracking toggles
parser.add_argument("--produce3DHits", action="store_true", default=False,
                    help="If set, also write GGTF3DHits (disabled by default).")
parser.add_argument("--ggtfLog", choices=["INFO", "DEBUG"], default="INFO",
                    help="GGTF_tracking OutputLevel.")
parser.add_argument("--maxHitsPerEvent", type=int, default=0,
                    help="If >0, hard-cap number of input hits per event in GGTF (for diagnostics).")

# Stage control: run only part of the chain for isolation
parser.add_argument("--stage", choices=["digi", "ggtf", "fit"], default="ggtf",
                    help="Which pipeline stage(s) to run: digi | ggtf | fit")

# Fitter toggles
parser.add_argument("--fitter", choices=["none","genfit2"], default="none",
                    help="Optionally run a fitter after GGTF. 'genfit2' consumes GGTF3DHits.")
parser.add_argument("--fitOut", default="GenFitTracks",
                    help="Output collection name for fitted tracks (when fitter is enabled).")

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

# --------- Model staging helper (handles .md5/http/root) ----------
def stage_model(spec: str) -> str:
    """
    Return a local .onnx path for ONNXRuntime.
    Accepts:
      - /path/to/model.onnx
      - /path/to/model.onnx.md5 (downloads by hash)
      - http(s)://... (downloads)
      - root://... (xrdcp)
    """
    if not spec:
        raise RuntimeError("modelPath is empty")

    # Already a real ONNX on disk
    if spec.endswith(".onnx") and os.path.exists(spec):
        return os.path.abspath(spec)

    out = os.path.abspath("model.onnx")

    # .md5 checksum file → resolve to hosted ONNX
    if spec.endswith(".onnx.md5"):
        with open(spec) as f:
            md5 = f.read().split()[0]
        url = f"https://key4hep.web.cern.ch/testFiles/k4RecTracker/{md5}"
        print(f"[model] from md5: {md5} -> {url}")
        subprocess.run(["wget", "--no-verbose", "--timeout=180", "--tries=2", "-O", out, url], check=True)
        if not os.path.exists(out) or os.path.getsize(out) == 0:
            raise RuntimeError(f"Downloaded model is missing/empty: {out}")
        return out

    # http(s) URL
    if spec.startswith("http://") or spec.startswith("https://"):
        print(f"[model] download {spec} -> {out}")
        subprocess.run(["wget", "--no-verbose", "--timeout=180", "--tries=2", "-O", out, spec], check=True)
        if not os.path.exists(out) or os.path.getsize(out) == 0:
            raise RuntimeError(f"Downloaded model is missing/empty: {out}")
        return out

    # xrootd URL
    if spec.startswith("root://"):
        print(f"[model] xrdcp {spec} -> {out}")
        subprocess.run(["xrdcp", "-f", spec, out], check=True)
        if not os.path.exists(out) or os.path.getsize(out) == 0:
            raise RuntimeError(f"xrdcp model is missing/empty: {out}")
        return out

    # Plain path to .onnx we couldn’t find above → try copying into CWD
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
    from Configurables import GGTF_tracking  # fallback to registry-based

GGTF = GGTF_tracking(
    "GGTF_tracking",
    inputWireHits=["DCH_DigiCollection"],  # produced by DCHdigi_v01
    inputPlanarHits=[],
    outputTracks=["CDCHTracks"],
    output3DHits=["GGTF_3DHits"],
    OutputLevel=INFO,
)

# Stage the model now (robust to .onnx.md5/http/root)
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

# Map --ggtfLog
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
        field_svc_obj = ConstBFieldSvc("GenFitFieldSvc", Bz=2.0)  # adjust as needed
        field_svc_name = field_svc_obj.getName()
        print(f"[genfit] Using ConstBFieldSvc Bz=2.0T -> {field_svc_name}")
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
fitter_alg = None
if args.stage == "fit" and args.fitter == "genfit2":
    # GenFit2 consumes 3D spacepoints; force-enable if user forgot
    try:
        if not bool(getattr(GGTF, "produce3DHits", False)):
            GGTF.produce3DHits = True
            print("[fitter] Enabling GGTF.produce3DHits = True (required by GenFit2).")
    except Exception:
        pass

    try:
        from Configurables import GenFit2DCHFitter
        fitter_alg = GenFit2DCHFitter("GenFit2DCHFitter", OutputLevel=INFO)

        # Accept either property style (string vs list) across builds
        # input3DHits / inputHits
        for prop, val in (("input3DHits", "GGTF3DHits"), ("inputHits", ["GGTF3DHits"])):
            try:
                if hasattr(fitter_alg, prop):
                    setattr(fitter_alg, prop, val)
                    print(f"[fitter] set {prop} = {val}")
                    break
            except Exception:
                pass

        # outputTracks (string or list)
        for prop, val in (("outputTracks", args.fitOut), ("outputTracks", [args.fitOut])):
            try:
                if hasattr(fitter_alg, prop):
                    setattr(fitter_alg, prop, val)
                    print(f"[fitter] set {prop} = {val}")
                    break
            except Exception:
                pass

        # Wire services if supported
        for prop, val in [
            ("FieldSvc", field_svc_name),
            ("BFieldSvc", field_svc_name),          # some builds use BFieldSvc name
            ("MaterialSvc", material_svc_name),
            ("MaterialEffectsSvc", material_svc_name),
            ("GeoSvcName", "GeoSvc"),
        ]:
            try:
                if val is not None and hasattr(fitter_alg, prop):
                    setattr(fitter_alg, prop, val)
                    print(f"[genfit] Set {prop} = {val}")
            except Exception as e:
                print(f"[genfit] Could not set {prop}: {e}")

        # Friendly defaults (guarded)
        for name, val in [
            ("pdgHypothesis", 13),     # muon hypothesis
            ("minHitsOnTrack", 4),
            ("maxChi2", 1e6),
        ]:
            try:
                if hasattr(fitter_alg, name):
                    setattr(fitter_alg, name, val)
                    print(f"[fitter] set {name} = {val}")
            except Exception:
                pass

        print(f"[fitter] GenFit2DCHFitter configured; output -> '{args.fitOut}'")
    except Exception as e:
        print(f"[fitter] GenFit2DCHFitter not found/failed to configure: {e}")
        print("[fitter] Continuing WITHOUT fitter.")
        fitter_alg = None

# --- Optional: collection size probe (if available)
try:
    from Configurables import EDM4hepCollectionSizePrinter as SizePrinter
    size_printer = SizePrinter(
        "SizePrinter",
        CollectionsToPrint=["GGTF3DHits", args.fitOut],
        OutputLevel=INFO,
    )
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

mgr = GaudiApp(
    TopAlg=top_algs,
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=ext_svcs,
    OutputLevel=INFO,
)
