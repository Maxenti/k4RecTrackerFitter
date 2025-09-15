#!/bin/bash
# run_dch_chain.sh  (DCH-only chain; uses TrackingConf.py from the bundle)
set -euo pipefail
set -o pipefail

# ---- make stdout/stderr line-buffered & timestamped for condor_tail ----
exec > >(stdbuf -oL awk '{ print strftime("[%Y-%m-%dT%H:%M:%S]"), $0; fflush() }') \
     2> >(stdbuf -oL awk '{ print strftime("[%Y-%m-%dT%H:%M:%S]"), $0 > "/dev/stderr"; fflush() }')

echo "[$(date)] START on $(hostname) ; PWD=$PWD"

# --------- args ---------
JOB_ID="${1:?job id}"
EVENTS="${2:?events per job}"
SKIP="${3:?skip events}"
SEED="${4:?random seed}"
OUT_EOS_URL="${5:?xrootd eos dest OR /eos/... path}"
MODEL_SPEC="${6:?can be .onnx, .onnx.md5, https://..., or root://...}"
COMPACT_XML="${7:?absolute path or EOS/xrootd URL to compact xml}"
BUNDLE_URL="${8:-}"   # optional: xrootd/https path to k4rt_runtime.tar.gz

# --------- runtime knobs from environment (JDL 'environment = ...') ---------
: "${GGTF_LOG:=INFO}"                 # INFO|DEBUG
: "${PRODUCE_3DHITS:=0}"              # 0|1
: "${MAX_HITS:=0}"                    # 0 => unlimited; >0 => cap per event
: "${FITTER:=none}"                   # none|genfit2
: "${FIT_OUT:=GenFitTracks}"          # fitter output coll name
: "${STAGE:=ggtf}"                    # digi|ggtf|fit  (matches runDCHTestTrackFinder.py)
: "${SKIP_DDSIM:=0}"                  # 1 => reuse existing sim_*.root if present
: "${TIMEOUT_DDSIM:=0}"               # seconds; 0 => no timeout
: "${TIMEOUT_K4RUN:=0}"               # seconds; 0 => no timeout
: "${MAX_VMEM_MB:=0}"                 # optional soft vmem cap (ulimit -v)

echo "[knobs] GGTF_LOG=$GGTF_LOG PRODUCE_3DHITS=$PRODUCE_3DHITS MAX_HITS=$MAX_HITS FITTER=$FITTER FIT_OUT=$FIT_OUT STAGE=$STAGE SKIP_DDSIM=$SKIP_DDSIM"
echo "[knobs] TIMEOUT_DDSIM=$TIMEOUT_DDSIM TIMEOUT_K4RUN=$TIMEOUT_K4RUN MAX_VMEM_MB=$MAX_VMEM_MB"

# Small helpers
WGET="wget --no-verbose --timeout=180 --tries=2"
XRDX="timeout 300 xrdcp -f"
TMO_NET="timeout 300"

# --------- resolve EOS dest early (so traps can use it) ---------
DEST_URL="${OUT_EOS_URL}"
EOS_DIR=""
if [[ "${OUT_EOS_URL}" == /eos/* ]]; then
  EOS_DIR="${OUT_EOS_URL}"
  DEST_URL="root://eosuser.cern.ch//${OUT_EOS_URL#/}"
elif [[ "${OUT_EOS_URL}" == root://* ]]; then
  TMP="${OUT_EOS_URL#*://}"         # eosuser.cern.ch//eos/...
  EOS_DIR="/${TMP#*/}"              # /eos/...
fi

copy_logs() {
  # run from inside WORKDIR
  if [[ -n "${EOS_DIR}" ]]; then
    echo "[post] shipping logs to ${DEST_URL}"
    ${XRDX} ddsim.log     "${DEST_URL}/ddsim_${JOB_ID}.log"     || true
    ${XRDX} k4run.log     "${DEST_URL}/k4run_${JOB_ID}.log"     || true
    ${XRDX} progress.log  "${DEST_URL}/progress_${JOB_ID}.log"  || true
  fi
}

# Ensure we still ship logs even on failure
trap 'RC=$?; echo "[trap] EXIT with rc=$RC"; (cd "${WORKDIR:-.}" 2>/dev/null && copy_logs) || true; exit $RC' EXIT

# --- load Key4HEP stack ---
load_k4hep () { set +u; source "$1" -r 2025-09-11; set -u; }
if [ -f /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh ]; then
  load_k4hep /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
else
  load_k4hep /cvmfs/sw.hsf.org/key4hep/setup.sh
fi
echo "[env] PYTHON=$(command -v python3)"

# --------- workdir ---------
WORKDIR="job_${JOB_ID}"
mkdir -p "${WORKDIR}" logs

# Create top-level symlinks so Condor can transfer logs from IWD root
: > "${WORKDIR}/ddsim.log";    ln -sf "${WORKDIR}/ddsim.log"    ddsim.log
: > "${WORKDIR}/k4run.log";    ln -sf "${WORKDIR}/k4run.log"    k4run.log
: > "${WORKDIR}/progress.log"; ln -sf "${WORKDIR}/progress.log" progress.log

cd "${WORKDIR}"

# Re-entry guard
if [[ -f .already_started ]]; then
  echo "[guard] Re-entry detected; first start was at: $(cat .already_started)"
else
  date -Ins > .already_started
fi

# Timestamp every traced command
export PS4='+ [${BASHPID}] $(date -Ins) : '
set -x

# --------- fetch & unpack runtime bundle ---------
if [[ -n "${BUNDLE_URL}" ]]; then
  echo "[bundle] fetching ${BUNDLE_URL}"
  if [[ "${BUNDLE_URL}" == http* ]]; then
    ${TMO_NET} ${WGET} -O k4rt_runtime.tar.gz "${BUNDLE_URL}"
  else
    ${XRDX} "${BUNDLE_URL}" k4rt_runtime.tar.gz
  fi
else
  if [[ -f ../k4rt_runtime.tar.gz ]]; then
    cp -v ../k4rt_runtime.tar.gz .
  else
    echo "[FATAL] No bundle provided (BUNDLE_URL empty) and ../k4rt_runtime.tar.gz not found."
    exit 1
  fi
fi

tar xzf k4rt_runtime.tar.gz
echo "[bundle] contents:"


ls -1

# ---- prefer the bundle's setup helper if present
if [[ -f ./setup_runtime_env.sh ]]; then
  # shellcheck disable=SC1091
  source ./setup_runtime_env.sh
else
  export PYTHONPATH="$PWD/python:${PYTHONPATH:-}"
  export GAUDI_PLUGIN_PATH="$PWD:$PWD/Tracking:$PWD/DCHdigi:$PWD/ARCdigi:$PWD/VTXdigi:${GAUDI_PLUGIN_PATH:-}"
  export LD_LIBRARY_PATH="$PWD/Tracking:$PWD/DCHdigi:$PWD/ARCdigi:$PWD/VTXdigi:${LD_LIBRARY_PATH:-}"
  export GAUDI_CONFIGURABLES_PATH="$PWD:${GAUDI_CONFIGURABLES_PATH:-}"
  export GAUDI_CONFDB_PATH="$PWD:${GAUDI_CONFDB_PATH:-}"
fi

# If the bundle contains GenFit2, expose it (harmless if absent)
if [[ -d genfit/lib64 || -d genfit/lib ]]; then
  export LD_LIBRARY_PATH="$PWD/genfit/lib64:$PWD/genfit/lib:${LD_LIBRARY_PATH:-}"
  echo "[env] Added GenFit2 libs to LD_LIBRARY_PATH"
fi

echo "[debug] GAUDI_PLUGIN_PATH=$GAUDI_PLUGIN_PATH"
echo "[debug] PYTHONPATH=$PYTHONPATH"

# ---- preflight import using the fallback (works even w/o confdb)
python3 - <<'PY'
import sys
print("py:", sys.executable)
try:
    import TrackingConf as T
    GG = T.GGTF_tracking
    print("TrackingConf OK; GGTF_tracking =", GG)
except Exception as e:
    print("[FATAL] Cannot import GGTF_tracking via TrackingConf fallback:", e)
    raise
PY

# --------- stage model (.onnx or via .md5) ---------
MODEL_PATH="model.onnx"
if [[ "${MODEL_SPEC}" == *.onnx ]]; then
  if [[ "${MODEL_SPEC}" == http* ]]; then
    ${TMO_NET} ${WGET} -O "${MODEL_PATH}" "${MODEL_SPEC}"
  elif [[ "${MODEL_SPEC}" == root://* ]]; then
    ${XRDX} "${MODEL_SPEC}" "${MODEL_PATH}"
  else
    cp -v "${MODEL_SPEC}" "${MODEL_PATH}"
  fi
elif [[ "${MODEL_SPEC}" == *.onnx.md5 ]]; then
  HASH=$(awk '{print $1}' "${MODEL_SPEC}")
  URL="https://key4hep.web.cern.ch/testFiles/k4RecTracker/${HASH}"
  echo "[model] from md5: ${HASH} -> ${URL}"
  ${TMO_NET} ${WGET} -O "${MODEL_PATH}" "${URL}"
elif [[ "${MODEL_SPEC}" == http* ]]; then
  ${TMO_NET} ${WGET} -O "${MODEL_PATH}" "${MODEL_SPEC}"
elif [[ "${MODEL_SPEC}" == root://* ]]; then
  ${XRDX} "${MODEL_SPEC}" "${MODEL_PATH}"
else
  echo "[FATAL] Unrecognized model spec: ${MODEL_SPEC}"
  exit 1
fi
[[ -s "${MODEL_PATH}" ]] || { echo "[FATAL] model not found/empty: ${MODEL_PATH}"; exit 2; }

# --------- geometry: PREFER bundled tree (geom/src + compact.path) ---------
LOCAL_XML=""
if [[ -f geom/compact.path ]]; then
  REL_COMPACT="$(cat geom/compact.path)"
  CANDIDATE="geom/src/${REL_COMPACT}"
  if [[ -f "${CANDIDATE}" ]]; then
    LOCAL_XML="${CANDIDATE}"
    export DD4hep_XMLPATH="$PWD/geom/src:${DD4hep_XMLPATH:-}"
    echo "[geom] Using bundled geometry: ${LOCAL_XML}"
    echo "[geom] DD4hep_XMLPATH=$DD4hep_XMLPATH"
  else
    echo "[warn] geom/compact.path present but file missing: ${CANDIDATE}"
  fi
fi

# Fallback: use provided COMPACT_XML (stage locally if EOS/xrootd)
if [[ -z "${LOCAL_XML}" ]]; then
  LOCAL_XML="compact_${JOB_ID}.xml"
  if [[ "${COMPACT_XML}" == root://* ]]; then
    ${XRDX} "${COMPACT_XML}" "${LOCAL_XML}"
  elif [[ "${COMPACT_XML}" == /eos/* ]]; then
    ${XRDX} "root://eosuser.cern.ch//${COMPACT_XML#/}" "${LOCAL_XML}"
  else
    cp -v "${COMPACT_XML}" "${LOCAL_XML}"
  fi
  [[ -s "${LOCAL_XML}" ]] || { echo "[FATAL] compact xml not staged: ${LOCAL_XML}"; exit 3; }
  export DD4hep_XMLPATH="$(dirname "$(readlink -f "${LOCAL_XML}")"):${DD4hep_XMLPATH:-}"
  echo "[geom] Using staged compact: ${LOCAL_XML}"
  echo "[geom] DD4hep_XMLPATH=$DD4hep_XMLPATH"
fi

# ---- CPU + allocator knobs ----
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export TORCH_NUM_THREADS=1
export MALLOC_ARENA_MAX=2
export ORT_DISABLE_MEMORY_ARENA=1
export ORT_ENABLE_MEM_PATTERN=0

# Optional soft cap on virtual memory (MB)
if [[ -n "${MAX_VMEM_MB:-}" && "${MAX_VMEM_MB}" -gt 0 ]]; then
  ulimit -v "$(( MAX_VMEM_MB * 1024 ))" || true
  echo "[env] ulimit -v set to ${MAX_VMEM_MB} MB"
fi

# Helper: run with optional timeout (0 => no timeout)
run_with_to() {
  local to="$1"; shift
  if [[ "${to}" -gt 0 ]]; then
    timeout --signal=TERM --kill-after=30 "${to}" "$@"
  else
    "$@"
  fi
}

# --------- produce simulation (DCH-only) ---------
SIM_OUT="sim_${JOB_ID}.root"
if [[ "${SKIP_DDSIM}" == "1" && -s "${SIM_OUT}" ]]; then
  echo "[ddsim] SKIP_DDSIM=1 and ${SIM_OUT} exists; reusing."
else
  echo "[ddsim] events=${EVENTS} skip=${SKIP} seed=${SEED}"
  ( run_with_to "${TIMEOUT_DDSIM}" /usr/bin/time -v stdbuf -oL -eL ddsim --steeringFile ../SteeringFile_IDEA_o1_v03DCH.py \
        --compactFile "${LOCAL_XML}" \
        --numberOfEvents "${EVENTS}" \
        --skipNEvents "${SKIP}" \
        --random.seed "${SEED}" \
        --outputFile "${SIM_OUT}" \
        --runType batch ) 2>&1 | tee -a ddsim.log
  DD_RC=${PIPESTATUS[0]}
  [[ $DD_RC -eq 0 ]] || { echo "[FATAL] ddsim failed (rc=$DD_RC)"; exit $DD_RC; }
fi

# --------- run DCH-only tracking / stages ---------
TRK_OUT="tracks_${JOB_ID}.root"

# Build k4run args with our knobs (including stage)
K4_ARGS=(
  ../runDCHTestTrackFinder.py
  --inputFile "${SIM_OUT}"
  --outputFile "${TRK_OUT}"
  --modelPath "${MODEL_PATH}"
  --tbeta 0.15 --td 0.10
  --dchSimHits DCHCollection
  --compactXML "${LOCAL_XML}"
  --dchName DCH_v2
  --ggtfLog "${GGTF_LOG}"
  --fitter "${FITTER}"
  --fitOut "${FIT_OUT}"
  --stage "${STAGE}"
)
if [[ "${PRODUCE_3DHITS}" == "1" ]]; then
  K4_ARGS+=( --produce3DHits )
fi
if [[ "${MAX_HITS}" -gt 0 ]]; then
  K4_ARGS+=( --maxHitsPerEvent "${MAX_HITS}" )
fi

echo "[k4run] args: ${K4_ARGS[*]}"
# Full log -> k4run.log ; filtered progress -> progress.log
( run_with_to "${TIMEOUT_K4RUN}" /usr/bin/time -v stdbuf -oL -eL k4run "${K4_ARGS[@]}" ) 2>&1 \
  | tee -a k4run.log \
  | awk '/GGTF_tracking|RSS=|Peak=|TOTAL|flatten:|onnx:|clustering|unique|bucket|build|Application Manager|tracks=/{print; fflush()}' > progress.log
K4_RC=${PIPESTATUS[0]}
[[ $K4_RC -eq 0 ]] || { echo "[FATAL] k4run failed (rc=$K4_RC)"; exit $K4_RC; }

# --------- copy products to EOS ---------
if [[ -n "${EOS_DIR}" ]]; then
  echo "[xrdfs] mkdir -p ${EOS_DIR}"
  timeout 60 xrdfs eosuser.cern.ch mkdir -p "${EOS_DIR}" || true
fi

echo "[xrdcp] -> ${DEST_URL}"
${XRDX} "${SIM_OUT}"  "${DEST_URL}/${SIM_OUT}"
${XRDX} "${TRK_OUT}"  "${DEST_URL}/${TRK_OUT}"
ls -lh "${SIM_OUT}" "${TRK_OUT}" || true

# ship logs (also handled by trap, but doing it here too is harmless)
copy_logs

echo "[$(date)] DONE job ${JOB_ID}"
