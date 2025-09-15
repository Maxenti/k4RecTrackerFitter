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

# env knobs (optional)
: "${GGTF_LOG:=INFO}"          # INFO|DEBUG
: "${PRODUCE_3DHITS:=1}"       # 0|1
: "${MAX_HITS:=0}"             # 0=no cap
: "${FITTER:=genfit2}"            # none|genfit2
: "${FIT_OUT:=GenFitTracks}"
: "${TIMEOUT_K4RUN:=0}"

echo "[cfg] INPUT=$INPUT"
echo "[cfg] OUTPUT=$OUTPUT"
echo "[cfg] MODEL_SPEC=$MODEL_SPEC"
echo "[cfg] COMPACT_XML=$COMPACT_XML  DCH_SIMHITS=$DCH_SIMHITS  DCH_NAME=$DCH_NAME"
echo "[cfg] GGTF_LOG=$GGTF_LOG PRODUCE_3DHITS=$PRODUCE_3DHITS MAX_HITS=$MAX_HITS FITTER=$FITTER FIT_OUT=$FIT_OUT TIMEOUT_K4RUN=$TIMEOUT_K4RUN"



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
  --fitter     "${FITTER}"
  --fitOut     "${FIT_OUT}"
  --stage      "fit"
)

[[ "${PRODUCE_3DHITS}" == "1" ]] && K4_ARGS+=( --produce3DHits )
[[ "${MAX_HITS}" -gt 0 ]]        && K4_ARGS+=( --maxHitsPerEvent "${MAX_HITS}" )

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

echo "[DONE] See k4run.log (full) and progress.log (key checkpoints)"
