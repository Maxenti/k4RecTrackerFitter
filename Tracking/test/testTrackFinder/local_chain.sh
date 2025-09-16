#!/bin/bash
# local_chain.sh — run ONLY the k4run stage locally
set -euo pipefail

########## defaults you can edit ##########
DEFAULT_INPUT="sim_local.root"
DEFAULT_OUTPUT="reco_local.root"
DEFAULT_MODEL_SPEC="/afs/cern.ch/user/c/cglenn/FCCWork/k4RecTracker/Tracking/test/testTrackFinder/model.onnx"
DEFAULT_COMPACT_XML="/eos/user/c/cglenn/FCCWork/GithubRepos/k4geoMax/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03CF_2umAu.xml"
DEFAULT_DCH_SIMHITS="DCHCollection"
DEFAULT_DCH_NAME="DCH_v2"
###########################################

# CLI overrides (optional)
INPUT="${1:-$DEFAULT_INPUT}"
OUTPUT="${2:-$DEFAULT_OUTPUT}"
MODEL_SPEC="${3:-$DEFAULT_MODEL_SPEC}"
COMPACT_XML="${4:-$DEFAULT_COMPACT_XML}"
DCH_SIMHITS="${5:-$DEFAULT_DCH_SIMHITS}"
DCH_NAME="${6:-$DEFAULT_DCH_NAME}"

# ----------------- env knobs (override as needed) -----------------
: "${GGTF_LOG:=INFO}"          # INFO|DEBUG
: "${PRODUCE_3DHITS:=1}"       # 0|1
: "${MAX_HITS:=0}"             # 0=no cap; try 120 for stability tests
: "${FITTER:=genfit2}"         # none|genfit2
: "${FIT_OUT:=GenFitTracks}"
: "${TIMEOUT_K4RUN:=0}"

# GGTF clustering thresholds (tighter helps avoid mixed clusters)
: "${TBETA:=0.10}"             # slightly larger than 0.05
: "${TD:=0.08}"

# GenFit “safe stability profile” (very forgiving first updates)
: "${GF_POS_SCALE:=0.1}"       # mm -> cm internally
: "${GF_LEN2M:=0.01}"          # cm -> m for pT seeding
: "${GF_HIT_SIGMA_XY:=5.0}"    # mm
: "${GF_HIT_SIGMA_Z:=30.0}"    # mm
: "${GF_SEED_POS_SIGMA:=500}"  # mm
: "${GF_SEED_MOM_SIGMA:=30.0}" # GeV
: "${GF_DEDUP_TOL:=1.00}"      # mm
: "${GF_USE_MAT:=0}"           # 0 = disable material effects while stabilising

# Seed clamps
: "${GF_SEED_PT_MIN:=1.0}"     # GeV  (keeps βγ away from 0)
: "${GF_SEED_PT_MAX:=20.0}"    # GeV
: "${GF_SEED_P_MIN:=1.0}"      # GeV  (new: clamp on total |p|)

# Grouping & fallback clustering for fitter
: "${GF_MIN_GROUP:=8}"
: "${GF_USE_FALLBACK:=1}"      # 1=on, 0=off
: "${GF_FALLBACK_EPS_CM:=2.5}"
: "${GF_FALLBACK_MINPTS:=8}"

# Retry controls (if base fit has no FitterInfo)
: "${GF_RETRY:=1}"             # 1=on, 0=off
: "${GF_RETRY_MEAS_INFL:=9.0}" # variance k (C' = k*C)
: "${GF_RETRY_SEED_POS:=5.0}"  # seed pos sigma ×
: "${GF_RETRY_SEED_MOM:=5.0}"  # seed mom sigma ×

# Optional: cap number of measurements per group (helps conditioning)
: "${GF_MAX_MEAS_PER_GROUP:=60}"

echo "[cfg] INPUT=$INPUT"
echo "[cfg] OUTPUT=$OUTPUT"
echo "[cfg] MODEL_SPEC=$MODEL_SPEC"
echo "[cfg] COMPACT_XML=$COMPACT_XML  DCH_SIMHITS=$DCH_SIMHITS  DCH_NAME=$DCH_NAME"
echo "[cfg] GGTF_LOG=$GGTF_LOG PRODUCE_3DHITS=$PRODUCE_3DHITS MAX_HITS=$MAX_HITS FITTER=$FITTER FIT_OUT=$FIT_OUT TIMEOUT_K4RUN=$TIMEOUT_K4RUN"
echo "[cfg] TBETA=$TBETA TD=$TD"
echo "[cfg] GF_POS_SCALE=$GF_POS_SCALE GF_LEN2M=$GF_LEN2M"
echo "[cfg] GF_HIT_SIGMA_XY=$GF_HIT_SIGMA_XY GF_HIT_SIGMA_Z=$GF_HIT_SIGMA_Z"
echo "[cfg] GF_SEED_POS_SIGMA=$GF_SEED_POS_SIGMA GF_SEED_MOM_SIGMA=$GF_SEED_MOM_SIGMA GF_DEDUP_TOL=$GF_DEDUP_TOL"
echo "[cfg] GF_USE_MAT=$GF_USE_MAT GF_SEED_PT_MIN=$GF_SEED_PT_MIN GF_SEED_PT_MAX=$GF_SEED_PT_MAX GF_SEED_P_MIN=$GF_SEED_P_MIN"
echo "[cfg] GF_MIN_GROUP=$GF_MIN_GROUP GF_USE_FALLBACK=$GF_USE_FALLBACK GF_FALLBACK_EPS_CM=$GF_FALLBACK_EPS_CM GF_FALLBACK_MINPTS=$GF_FALLBACK_MINPTS"
echo "[cfg] GF_RETRY=$GF_RETRY GF_RETRY_MEAS_INFL=$GF_RETRY_MEAS_INFL GF_RETRY_SEED_POS=$GF_RETRY_SEED_POS GF_RETRY_SEED_MOM=$GF_RETRY_SEED_MOM"
echo "[cfg] GF_MAX_MEAS_PER_GROUP=$GF_MAX_MEAS_PER_GROUP"

# --- keep memory tame ---
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export TORCH_NUM_THREADS=1
export MALLOC_ARENA_MAX=2
export ORT_DISABLE_MEMORY_ARENA=1
export ORT_ENABLE_MEM_PATTERN=0

export GAUDI_PLUGIN_PATH="${GAUDI_PLUGIN_PATH:-.}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-.}"

# logs
: > k4run.log
: > progress.log

# Build k4run args (Python script will stage the model)
K4_ARGS=(
  ./runDCHTestTrackFinder.py
  --inputFile  "$INPUT"
  --outputFile "$OUTPUT"
  --modelPath  "$MODEL_SPEC"
  --compactXML "${COMPACT_XML}"
  --dchName    "${DCH_NAME}"
  --dchSimHits "${DCH_SIMHITS}"
  --ggtfLog    "${GGTF_LOG}"
  --tbeta      "${TBETA}"
  --td         "${TD}"
  --fitter     "${FITTER}"
  --fitOut     "${FIT_OUT}"
  --stage      "fit"
)

# Ensure GGTF produces 3D hits for the fitter
[[ "${PRODUCE_3DHITS}" == "1" ]] && K4_ARGS+=( --produce3DHits )
[[ "${MAX_HITS}" -gt 0 ]]        && K4_ARGS+=( --maxHitsPerEvent "${MAX_HITS}" )

# --- GenFit stability profile ---
K4_ARGS+=(
  --gf-posScale       "${GF_POS_SCALE}"
  --gf-len2m          "${GF_LEN2M}"
  --gf-hitSigmaXY     "${GF_HIT_SIGMA_XY}"
  --gf-hitSigmaZ      "${GF_HIT_SIGMA_Z}"
  --gf-seedPosSigma   "${GF_SEED_POS_SIGMA}"
  --gf-seedMomSigma   "${GF_SEED_MOM_SIGMA}"
  --gf-dedupTol       "${GF_DEDUP_TOL}"
  --gf-seedPTMin      "${GF_SEED_PT_MIN}"
  --gf-seedPTMax      "${GF_SEED_PT_MAX}"
  --gf-seedPMin       "${GF_SEED_P_MIN}"
  --gf-minGroup       "${GF_MIN_GROUP}"
  --gf-fallbackEpsCM  "${GF_FALLBACK_EPS_CM}"
  --gf-fallbackMinPts "${GF_FALLBACK_MINPTS}"
  --gf-retryMeasInfl  "${GF_RETRY_MEAS_INFL}"
  --gf-retrySeedPos   "${GF_RETRY_SEED_POS}"
  --gf-retrySeedMom   "${GF_RETRY_SEED_MOM}"
  --gf-maxMeasPerGroup "${GF_MAX_MEAS_PER_GROUP}"
)
if [[ "${GF_USE_MAT}" == "1" ]]; then
  K4_ARGS+=( --gf-useMat )
else
  K4_ARGS+=( --no-gf-useMat )
fi
if [[ "${GF_USE_FALLBACK}" == "1" ]]; then
  K4_ARGS+=( --gf-useFallback )
else
  K4_ARGS+=( --no-gf-useFallback )
fi
if [[ "${GF_RETRY}" == "1" ]]; then
  K4_ARGS+=( --gf-retry )
else
  K4_ARGS+=( --no-gf-retry )
fi

echo "[k4run] args: ${K4_ARGS[*]}" | tee -a k4run.log

run_cmd() {
  if [[ "${TIMEOUT_K4RUN}" -gt 0 ]]; then
    timeout --signal=TERM --kill-after=30 "${TIMEOUT_K4RUN}" "$@"
  else
    "$@"
  fi
}

# Full log -> k4run.log ; filtered checkpoints -> progress.log
( run_cmd /usr/bin/time -v stdbuf -oL -eL k4run "${K4_ARGS[@]}" ) 2>&1 \
  | tee -a k4run.log \
  | awk '/GGTF_tracking|RSS=|Peak=|TOTAL|flatten:|onnx:|clustering|unique|bucket|build|Application Manager|tracks=|MemoryAuditor/{print; fflush()}' > progress.log

K4_RC=${PIPESTATUS[0]}

check_edm_root() {
  python3 - "$1" <<'PY'
import sys, ROOT
f = ROOT.TFile.Open(sys.argv[1])
ok = bool(f and not f.IsZombie())
nev = 0
if ok:
    t = f.Get("events")
    nev = int(t.GetEntries()) if t else 0
f.Close() if f else None
print(f"[verify] {sys.argv[1]} OK=${ok} NEV=${nev}")
sys.exit(0 if ok and nev>0 else 1)
PY
}

if [[ $K4_RC -ne 0 ]]; then
  if check_edm_root "${OUTPUT:-reco_local.root}"; then
    echo "[note] k4run rc=$K4_RC but output looks fine; overriding to 0."
    K4_RC=0
  fi
fi

exit $K4_RC
