#!/usr/bin/env bash
set -euo pipefail

# --- edit these if needed ---
BUILD=${BUILD:-/afs/cern.ch/user/c/cglenn/FCCWork/k4RecTracker/build}
INSTALL=${INSTALL:-/afs/cern.ch/user/c/cglenn/FCCWork/k4RecTracker/install}
GENFIT_DIR=${GENFIT_DIR:-/afs/cern.ch/user/c/cglenn/FCCWork/genfit2/install}  # optional
OUT_TGZ=${OUT_TGZ:-k4rt_runtime.tar.gz}

# Optional geometry bundling (recommended)
GEOM_DIR=${GEOM_DIR:-/eos/user/c/cglenn/FCCWork/GithubRepos/k4geoMax/FCCee/IDEA/compact/IDEA_o1_v03/}
GEOM_COMPACT=${GEOM_COMPACT:-IDEA_o1_v03CF_2umAu.xml}

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP"/{Tracking,DCHdigi,ARCdigi,VTXdigi,python,meta}

echo "[pack] BUILD=$BUILD"
echo "[pack] INSTALL=$INSTALL"
echo "[pack] GENFIT_DIR=${GENFIT_DIR:-<none>}"
echo "[pack] GEOM_DIR=${GEOM_DIR:-<none>}"
echo "[pack] GEOM_COMPACT=${GEOM_COMPACT:-<none>}"

have() { command -v "$1" >/dev/null 2>&1; }

# --- libs (with fallbacks) ---
copy_one() {
  local srcA="$1" srcB="$2" dest="$3" name="$4"
  if [[ -f "$srcA" ]]; then cp -v "$srcA" "$dest"
  elif [[ -f "$srcB" ]]; then cp -v "$srcB" "$dest"
  else
    echo "[ERROR] Missing $name at '$srcA' or '$srcB'"; return 1
  fi
}

# Tracking plugin
copy_one "$BUILD/Tracking/libTracking.so" \
         "$INSTALL/lib/libTracking.so" \
         "$TMP/Tracking/" "libTracking.so"

# Copy any Tracking dicts if present
for c in "$BUILD/Tracking" "$INSTALL/lib"; do
  [[ -d "$c" ]] || continue
  shopt -s nullglob
  for f in "$c"/*Tracking*Dict*.so "$c"/*Tracking*Dict*.pcm "$c"/*Tracking*.rootmap; do
    cp -v "$f" "$TMP/Tracking/" || true
  done
  shopt -u nullglob
done

# DCHdigi libs + dicts
copy_one "$BUILD/DCHdigi/libDCHdigi.so" "$INSTALL/lib/libDCHdigi.so" "$TMP/DCHdigi/" "libDCHdigi.so"
# Some builds place extension libs under DCHdigi; include generously
shopt -s nullglob
for c in "$BUILD/DCHdigi" "$INSTALL/lib"; do
  [[ -d "$c" ]] || continue
  for f in "$c"/libextension*.so "$c"/*DCHdigi*Dict*.so "$c"/*DCHdigi*Dict*.pcm "$c"/*DCHdigi*.rootmap; do
    cp -v "$f" "$TMP/DCHdigi/" || true
  done
done
shopt -u nullglob

# --- Gaudi registries/components ---
for where in "$BUILD" "$INSTALL"; do
  shopt -s nullglob
  for f in "$where"/*.components "$where"/*.confdb "$where"/*.confdb2; do
    cp -v "$f" "$TMP/" || true
  done
  for PKG in DCHdigi ARCdigi VTXdigi Tracking ; do
    for f in "$where/$PKG/$PKG.components" "$where/$PKG/$PKG.confdb" "$where/$PKG/$PKG.confdb2"; do
      [[ -f "$f" ]] && cp -v "$f" "$TMP/" || true
    done
  done
  shopt -u nullglob
done

# --- PYTHON: TrackingConf.py (+ __init__.py) ---
FOUND=""
for base in "$INSTALL" "$BUILD"; do
  c=$(find "$base" -maxdepth 5 -type f -name TrackingConf.py 2>/dev/null | head -n1 || true)
  if [[ -n "$c" ]]; then FOUND="$c"; break; fi
done
if [[ -n "$FOUND" ]]; then
  cp -v "$FOUND" "$TMP/python/TrackingConf.py"
else
  echo "[WARN] TrackingConf.py not found in INSTALL/BUILD; import may fall back to Configurables only."
fi
: > "$TMP/python/__init__.py"

# --- OPTIONAL: bundle GenFit2 (libs + dicts) ---
if [[ -n "${GENFIT_DIR:-}" && -d "$GENFIT_DIR" ]]; then
  echo "[pack] bundling GenFit2 from $GENFIT_DIR"
  for LIBSUB in lib64 lib ; do
    if [[ -d "$GENFIT_DIR/$LIBSUB" ]]; then
      mkdir -p "$TMP/genfit/$LIBSUB"
      cp -v "$GENFIT_DIR/$LIBSUB"/libgenfit2*.so            "$TMP/genfit/$LIBSUB/" 2>/dev/null || true
      cp -v "$GENFIT_DIR/$LIBSUB"/*RootDict*.so             "$TMP/genfit/$LIBSUB/" 2>/dev/null || true
      cp -v "$GENFIT_DIR/$LIBSUB"/*.pcm                     "$TMP/genfit/$LIBSUB/" 2>/dev/null || true
      cp -v "$GENFIT_DIR/$LIBSUB"/*.rootmap                 "$TMP/genfit/$LIBSUB/" 2>/dev/null || true
    fi
  done
else
  echo "[pack] GENFIT_DIR not set or missing; skipping GenFit2."
fi

# --- Bundle geometry tree (if GEOM_DIR set) ---
if [[ -n "${GEOM_DIR:-}" ]]; then
  echo "[pack] bundling geometry from: ${GEOM_DIR}"
  mkdir -p "$TMP/geom/src"
  if [[ "${GEOM_DIR}" == root://* ]]; then
    BASENAME="$(basename "${GEOM_DIR%/}")"
    xrdcp -r -f "${GEOM_DIR%/}" "$TMP/geom/"
    mv "$TMP/geom/${BASENAME}/"* "$TMP/geom/src/" || true
    rmdir "$TMP/geom/${BASENAME}" 2>/dev/null || true
  else
    if have rsync; then
      rsync -aL --delete "${GEOM_DIR%/}/" "$TMP/geom/src/"
    else
      cp -Lr "${GEOM_DIR%/}/." "$TMP/geom/src/"
    fi
  fi
  if [[ -n "${GEOM_COMPACT:-}" ]]; then
    echo "${GEOM_COMPACT}" > "$TMP/geom/compact.path"
    echo "[pack] wrote geom/compact.path = ${GEOM_COMPACT}"
    if [[ ! -f "$TMP/geom/src/${GEOM_COMPACT}" ]]; then
      echo "[ERROR] Compact not found after copy: $TMP/geom/src/${GEOM_COMPACT}"
      exit 2
    fi
  fi
fi

# --- Tiny runtime env for the worker ---
cat > "$TMP/setup_runtime_env.sh" <<'ENVSH'
#!/usr/bin/env bash
set -euo pipefail
echo "[setup] CWD=$(pwd)"
[[ -d Tracking ]] || { echo "[error] Tracking/ not present after unpack"; exit 2; }
export GAUDI_PLUGIN_PATH="$PWD:$PWD/Tracking:$PWD/DCHdigi:$PWD/ARCdigi:$PWD/VTXdigi:${GAUDI_PLUGIN_PATH:-}"
export LD_LIBRARY_PATH="$PWD/Tracking:$PWD/DCHdigi:$PWD/ARCdigi:$PWD/VTXdigi:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$PWD/python:${PYTHONPATH:-}"
echo "[env] GAUDI_PLUGIN_PATH=${GAUDI_PLUGIN_PATH}"
echo "[env] LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
echo "[env] PYTHONPATH=${PYTHONPATH}"
ENVSH
chmod +x "$TMP/setup_runtime_env.sh"

# --- Manifest / fingerprint (very helpful on WNs) ---
{
  echo "build:  $(date -Is)"
  echo "who:    $(whoami)@$(hostname -f 2>/dev/null || hostname)"
  echo "k4rt git:  $( (cd "$BUILD" 2>/dev/null && git rev-parse --short=12 HEAD) || echo '<n/a>')"
  echo "genfit git:$( (cd "$GENFIT_DIR" 2>/dev/null && git rev-parse --short=12 HEAD) || echo '<n/a>')"
  echo "tracking.so strings (subset):"
  if [[ -f "$BUILD/Tracking/libTracking.so" ]]; then
    strings "$BUILD/Tracking/libTracking.so" | egrep -i 'GGTF_tracking|maxHitsPerEvent|after-onnx|after-bucket|TOTAL' | sed 's/^/  /' || true
  fi
} > "$TMP/meta/manifest.txt"

# --- pack tarball ---
rm -f "$OUT_TGZ"
tar -czf "$OUT_TGZ" -C "$TMP" .
echo "Wrote $OUT_TGZ"
tar -tzf "$OUT_TGZ" | sed 's/^/  /'
