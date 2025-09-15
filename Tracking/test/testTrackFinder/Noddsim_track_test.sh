#!/usr/bin/env bash
set -euo pipefail

# Figure out repo root from this script's location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"   # .../k4RecTracker
BUILD_DIR="${REPO_ROOT}/build"


# Resolve model path
MODEL_PATH="${1-}"   # first arg if given
if [[ -z "${MODEL_PATH}" ]]; then
  # Try to fetch via md5 like ExternalData would
  MD5_FILE="${REPO_ROOT}/Tracking/test/inputFiles/SimpleGatrIDEAv3o1.onnx.md5"
  if [[ ! -f "${MD5_FILE}" ]]; then
    echo "ERROR: No MODEL_PATH provided and md5 file not found at ${MD5_FILE}"
    echo "Usage: $0 <MODEL_PATH.onnx>"
    exit 2
  fi
  HASH="$(cat "${MD5_FILE}")"
  MODEL_PATH="${BUILD_DIR}/SimpleGatrIDEAv3o1.onnx"
  if [[ ! -f "${MODEL_PATH}" ]]; then
    echo "Downloading ONNX model to ${MODEL_PATH} ..."
    curl -fsSL "https://key4hep.web.cern.ch/testFiles/k4RecTracker/${HASH}" -o "${MODEL_PATH}"
  fi
fi

# If someone passed DATA{...}, strip the wrapper (ctest expands it, manual run might pass it literally)
if [[ "${MODEL_PATH}" == DATA{*} ]]; then
  MODEL_PATH="${MODEL_PATH#DATA\{}"; MODEL_PATH="${MODEL_PATH%\}}"
fi

if [[ ! -f "${MODEL_PATH}" ]]; then
  echo "ERROR: MODEL_PATH does not exist: ${MODEL_PATH}"
  exit 3
fi
echo "Using ONNX model: ${MODEL_PATH}"

# Now call the steering script from the repo root
# Choose ONE of these lines depending on how runTestTrackFinder.py parses it:

# If runTestTrackFinder.py expects a *named* option:
#k4run "${REPO_ROOT}/Tracking/test/testTrackFinder/runTestTrackFinder.py" --modelPath "${MODEL_PATH}" "$@"

# If it expects the model as *positional* arg, use this instead:
# k4run "${REPO_ROOT}/Tracking/test/testTrackFinder/runTestTrackFinder.py" "${MODEL_PATH}" "$@"


XML_FILE=$K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml
TBETA=0.6
TD=0.3

    
k4run /afs/cern.ch/user/c/cglenn/FCCWork/b0721d0/k4RecTracker/Tracking/test/testTrackFinder/runTestTrackFinder.py \
      --inputFile /afs/cern.ch/user/c/cglenn/FCCWork/b0721d0/k4RecTracker/Tracking/test/testTrackFinder/out_sim_edm4hep.root \
      --outputFile /afs/cern.ch/user/c/cglenn/FCCWork/b0721d0/k4RecTracker/Tracking/test/testTrackFinder/out_tracks.root \
      --modelPath $MODEL_PATH \
      --tbeta $TBETA \
      --td $TD 
