#!/usr/bin/env bash
#
# install_geant4.sh — download, build, and install Geant4 from source.
#
# Installs into $BASE_DIR (default ~/geant4):
#   ~/geant4/geant4-src/geant4-<version>/      source tree
#   ~/geant4/geant4-build/                     out-of-tree build directory
#   ~/geant4/geant4-v<version>-install/        install prefix + data files
#
# Usage:
#   ./install_geant4.sh [options]
# Options:
#   -v, --version <ver>   Geant4 version to install (default: 11.4.3,
#                         matching the version pinned in docker/Dockerfile)
#   -p, --prefix <dir>    base directory (default: ~/geant4)
#   -j, --jobs <N>        parallel build jobs (default: nproc)
#       --no-data         skip downloading the ~1.4 GB physics data files
#                         (NOT recommended: simulation runs need them)
#   -f, --force           remove any existing build/install for this version
#                         and rebuild from scratch
#   -h, --help            show this help
#
# After a successful install, load the environment with:
#   source ~/geant4/geant4-v<version>-install/bin/geant4.sh

set -euo pipefail

G4_VERSION="11.4.3"
BASE_DIR="$HOME/geant4"
JOBS="$(nproc 2>/dev/null || echo 4)"
INSTALL_DATA="ON"
FORCE="0"

log()  { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33mWARNING:\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31mERROR:\033[0m %s\n' "$*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        -v|--version) G4_VERSION="$2"; shift 2 ;;
        -p|--prefix)  BASE_DIR="$2"; shift 2 ;;
        -j|--jobs)    JOBS="$2"; shift 2 ;;
        --no-data)    INSTALL_DATA="OFF"; shift ;;
        -f|--force)   FORCE="1"; shift ;;
        -h|--help)    grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) die "Unknown option: $1 (see --help)" ;;
    esac
done

SRC_DIR="$BASE_DIR/geant4-src"
BUILD_DIR="$BASE_DIR/geant4-build"
INSTALL_DIR="$BASE_DIR/geant4-v${G4_VERSION}-install"
TARBALL="v${G4_VERSION}.tar.gz"
URL="https://github.com/geant4/geant4/archive/refs/tags/${TARBALL}"

[[ "$JOBS" =~ ^[0-9]+$ ]] || die "--jobs must be a number, got: $JOBS"

# ---------------------------------------------------------------- prerequisites
need_pkg=()
have() { command -v "$1" &>/dev/null; }
have g++        || need_pkg+=(g++)
have make       || need_pkg+=(make)
have wget       || need_pkg+=(wget)
have tar        || need_pkg+=(tar)
cmake_ok=0
if have cmake; then
    cmake_ver="$(cmake --version | head -1 | awk '{print $3}')"
    [[ "$(printf '%s\n3.16\n' "$cmake_ver" | sort -V | head -1)" == "3.16" ]] && cmake_ok=1
fi
[[ "$cmake_ok" == 1 ]] || need_pkg+=(cmake)
# Headers needed by Geant4's find_dependency(EXPAT) / find_dependency(ZLIB)
[[ -e /usr/include/zlib.h ]]      || need_pkg+=(zlib1g-dev)
[[ -e /usr/include/expat.h ]]     || need_pkg+=(libexpat1-dev)

if [[ ${#need_pkg[@]} -gt 0 ]]; then
    log "Installing missing prerequisites: ${need_pkg[*]}"
    if have apt-get; then
        sudo apt-get update && sudo apt-get install -y --no-install-recommends "${need_pkg[@]}"
    elif have dnf; then
        sudo dnf install -y "${need_pkg[@]/zlib1g-dev/zlib-devel}" "${need_pkg[@]/libexpat1-dev/expat-devel}"
    elif have pacman; then
        sudo pacman -S --needed "${need_pkg[@]/zlib1g-dev/zlib}" "${need_pkg[@]/libexpat1-dev/expat}"
    else
        die "Missing: ${need_pkg[*]} — install them manually and re-run."
    fi
fi

# ------------------------------------------------------- idempotency / --force
if [[ -f "$INSTALL_DIR/bin/geant4.sh" ]]; then
    if [[ "$FORCE" == 1 ]]; then
        log "Removing existing build and install (--force)..."
        rm -rf "$BUILD_DIR" "$INSTALL_DIR"
    else
        log "Geant4 $G4_VERSION is already installed at: $INSTALL_DIR"
        echo "  Load it with:  source \"$INSTALL_DIR/bin/geant4.sh\""
        echo "  Rebuild anyway with: $0 --force"
        exit 0
    fi
fi

mkdir -p "$SRC_DIR" "$BUILD_DIR"

# ------------------------------------------------------------------- download
if [[ -d "$SRC_DIR/geant4-${G4_VERSION}" ]]; then
    log "Source already extracted at: $SRC_DIR/geant4-${G4_VERSION}"
else
    log "Downloading Geant4 $G4_VERSION source (~90 MB)..."
    if have wget; then
        wget -O "$SRC_DIR/$TARBALL" "$URL"
    else
        curl -L -o "$SRC_DIR/$TARBALL" "$URL"
    fi
    log "Extracting..."
    tar -xzf "$SRC_DIR/$TARBALL" -C "$SRC_DIR"
fi

# ------------------------------------------------------------------ configure
log "Configuring (install prefix: $INSTALL_DIR, data: $INSTALL_DATA)..."
cmake -S "$SRC_DIR/geant4-${G4_VERSION}" -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
    -DGEANT4_INSTALL_DATA="$INSTALL_DATA" \
    -DGEANT4_BUILD_MULTITHREADED=ON

# ----------------------------------------------------------------- build+install
log "Building with $JOBS jobs (this takes ~15-25 min)..."
make -C "$BUILD_DIR" -j"$JOBS"
log "Installing (downloads physics data files if enabled)..."
make -C "$BUILD_DIR" install

# ---------------------------------------------------------------------- finish
log "Geant4 $G4_VERSION installed successfully!"
echo
echo "Load the environment (add this to your ~/.bashrc if you like):"
echo "  source \"$INSTALL_DIR/bin/geant4.sh\""
echo
echo "Then build and test the simulation natively:"
echo "  ./run_geant4_stopping.sh reactions/o18_si30_2n/macros/si30_backing.mac"
echo
echo "Disk usage note: keep geant4-build/ if you plan to rebuild; it can be"
echo "safely deleted to reclaim space (~3 GB): rm -rf $BUILD_DIR"
