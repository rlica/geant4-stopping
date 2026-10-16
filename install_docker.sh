#!/bin/bash
#===============================================================================
# install_docker.sh
#
# Installs Docker so that ./run_docker.sh works out of the box on:
#   - WSL (Ubuntu/Debian)
#   - native Linux: Debian/Ubuntu, Fedora/RHEL, and Arch based distros
#   - macOS (Intel & Apple Silicon), via Docker Desktop installed with Homebrew
#
# What it does:
#   1. Detects your platform and package manager.
#   2. Checks the environment and explains any problem in plain language.
#   3. Installs Docker Engine / Docker Desktop using the platform's best method.
#   4. Starts the Docker daemon and (on Linux) adds your user to the 'docker'
#      group so sudo is not needed at runtime.
#
# Usage:
#   ./install_docker.sh
#
# Safe to re-run: every step is skipped if it has already been done.
#===============================================================================
set -e

# Pretty output helpers
info()    { echo -e "\e[1;34m==> ${*}\e[0m"; }
ok()      { echo -e "\e[1;32m ✓  ${*}\e[0m"; }
warn()    { echo -e "\e[1;33m !  ${*}\e[0m"; }
fail()    { echo -e "\e[1;31m ✗  ${*}\e[0m"; }

# Crash handler: tell the user what to do when something goes wrong
on_error() {
    fail "Installation failed at line $1."
    echo "Common fixes:"
    echo "  - Check your internet connection, then re-run this script."
    echo "  - WSL: run 'wsl --shutdown' in Windows PowerShell, reopen Ubuntu, retry."
    echo "  - macOS: if Docker Desktop was just installed, open it once manually,"
    echo "    accept the security prompt and the service agreement, then re-run."
    exit 1
}
trap 'on_error $LINENO' ERR

#--- Platform detection ---------------------------------------------------------
OS="$(uname -s)"                       # "Linux" or "Darwin"
IS_WSL=false
DISTRO=""
DISTRO_PRETTY=""

if [ "$OS" == "Linux" ]; then
    grep -qiE '(microsoft|wsl)' /proc/version 2>/dev/null && IS_WSL=true
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        DISTRO_PRETTY="${PRETTY_NAME:-$ID}"
        case "$ID" in
            ubuntu|debian|linuxmint|pop|raspbian|kali)   DISTRO="debian" ;;
            fedora|rhel|rocky|almalinux|centos|amzn|ol)  DISTRO="fedora" ;;
            arch|manjaro|endeavouros|garuda|cachyos)     DISTRO="arch" ;;
            *)
                # Fall back on ID_LIKE (e.g. derivatives of the main families)
                case "$ID_LIKE" in
                    *debian*)        DISTRO="debian" ;;
                    *fedora*|*rhel*) DISTRO="fedora" ;;
                    *arch*)          DISTRO="arch" ;;
                esac
                ;;
        esac
    fi
fi

# The real user (Linux only; used for the 'docker' group). Works via sudo too.
REAL_USER="${SUDO_USER:-$USER}"

# Sudo prefix (empty when running as root). Not used on macOS (Homebrew).
SUDO=""
if [ "$(id -u)" -ne 0 ]; then SUDO="sudo"; fi

echo "=============================================================="
echo " Docker installer -- Geant4 Stopping Simulation"
echo "=============================================================="
if [ "$OS" == "Darwin" ]; then
    ok "Detected platform: macOS ($(uname -m))"
elif [ -n "$DISTRO" ]; then
    ok "Detected platform: ${DISTRO_PRETTY:-$DISTRO}$([ "$IS_WSL" == true ] && echo ' (WSL)')"
else
    warn "Detected platform: unknown Linux"
fi

#--- Step 1: skip if Docker is already installed and running --------------------
# Check via sudo too: the user may not be in the 'docker' group yet, in which
# case a plain `docker info` fails even though the daemon is healthy.
if command -v docker &>/dev/null && { docker info &>/dev/null || $SUDO docker info &>/dev/null; }; then
    ok "Docker $(docker --version | sed 's/,.*//') is already installed and running."
    # Even when Docker already works, the user may still be missing from the
    # 'docker' group (e.g. a previous run aborted before that step).
    if [ "$OS" == "Linux" ] && ! id -nG "$REAL_USER" | tr ' ' '\n' | grep -qx docker; then
        warn "You are not in the 'docker' group yet -- adding you now..."
        $SUDO usermod -aG docker "$REAL_USER"
        echo ""
        echo "  One last step: close this terminal and open a new one"
        echo "  (or run: newgrp docker), then start ./run_docker.sh"
    else
        warn "Nothing to do. If you still get permission errors with ./run_docker.sh,"
        warn "log out and back in (or run: newgrp docker), then re-run ./run_docker.sh"
    fi
    exit 0
fi

#--- Step 2 (Linux only): sanity checks -----------------------------------------
if [ "$OS" == "Linux" ]; then
    if [ -z "$DISTRO" ]; then
        fail "Unsupported Linux distribution: ${DISTRO_PRETTY:-unknown}."
        echo "This script supports Debian/Ubuntu, Fedora/RHEL and Arch based distros."
        echo "Please install Docker manually, then re-run ./run_docker.sh:"
        echo "  https://docs.docker.com/engine/install/"
        exit 1
    fi
    ok "Linux distribution family: ${DISTRO}"

    if [ "$IS_WSL" == true ] && ! command -v systemctl &>/dev/null; then
        fail "systemd is not running in this WSL distro."
        echo "Docker needs systemd to manage its service. To enable it:"
        echo "  1. Edit (or create) /etc/wsl.conf and add:"
        echo "       [boot]"
        echo "       systemd=true"
        echo "  2. From Windows PowerShell run:  wsl --shutdown"
        echo "  3. Reopen your Ubuntu terminal and re-run this script."
        exit 1
    fi

    if ! command -v sudo &>/dev/null && [ "$(id -u)" -ne 0 ]; then
        fail "'sudo' is not installed and you are not root."
        echo "Install sudo first, then re-run this script."
        exit 1
    fi

    # Nudge users away from the broken snap variant (Debian family / WSL)
    if [ "$DISTRO" == "debian" ] && command -v snap &>/dev/null && snap list docker &>/dev/null; then
        warn "A snap-installed Docker was detected. Snap Docker is broken under WSL."
        warn "Remove it first:  sudo snap remove docker  -- then re-run this script."
        exit 1
    fi
fi

#--- Install routines ------------------------------------------------------------
install_debian() {
    # Works on Ubuntu, Debian, Mint, Pop!_OS (uses Docker's official apt repo)
    local REPO_ID CODENAME
    if [ "${ID:-}" == "debian" ]; then REPO_ID="debian"; else REPO_ID="ubuntu"; fi
    CODENAME="${UBUNTU_CODENAME:-$VERSION_CODENAME}"

    info "Installing prerequisites (ca-certificates, curl)..."
    $SUDO apt-get update -qq
    $SUDO apt-get install -y -qq ca-certificates curl gnupg

    info "Adding Docker's official apt repository..."
    if [ ! -f /etc/apt/keyrings/docker.asc ]; then
        $SUDO install -m 0755 -d /etc/apt/keyrings
        $SUDO curl -fsSL "https://download.docker.com/linux/${REPO_ID}/gpg" -o /etc/apt/keyrings/docker.asc
        $SUDO chmod a+r /etc/apt/keyrings/docker.asc
        ok "Docker GPG key added"
    else
        ok "Docker GPG key already present"
    fi

    DOCKER_LIST="/etc/apt/sources.list.d/docker.list"
    if [ ! -f "$DOCKER_LIST" ]; then
        echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] \
https://download.docker.com/linux/${REPO_ID} ${CODENAME} stable" \
            | $SUDO tee "$DOCKER_LIST" > /dev/null
        ok "Docker apt repository added"
    else
        ok "Docker apt repository already configured"
    fi

    info "Installing Docker Engine (this may take a minute)..."
    $SUDO apt-get update -qq
    $SUDO apt-get install -y -qq docker-ce docker-ce-cli containerd.io docker-buildx-plugin
    ok "Docker packages installed: $(docker --version | sed 's/,.*//')"
}

install_fedora() {
    # Works on Fedora, RHEL, CentOS, Rocky, AlmaLinux (official Docker dnf repo)
    local REPO_ID
    case "${ID:-}" in
        fedora) REPO_ID="fedora" ;;
        rhel)   REPO_ID="rhel" ;;
        *)      REPO_ID="centos" ;;
    esac

    info "Adding Docker's official dnf repository..."
    if [ ! -f /etc/yum.repos.d/docker-ce.repo ]; then
        $SUDO curl -fsSL "https://download.docker.com/linux/${REPO_ID}/docker-ce.repo" \
            -o /etc/yum.repos.d/docker-ce.repo
        ok "Docker dnf repository added"
    else
        ok "Docker dnf repository already configured"
    fi

    info "Installing Docker Engine (this may take a minute)..."
    $SUDO dnf install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin
    ok "Docker packages installed: $(docker --version | sed 's/,.*//')"
}

install_arch() {
    # Works on Arch Linux, Manjaro, EndeavourOS
    info "Installing Docker Engine with pacman (this may take a minute)..."
    $SUDO pacman -Sy --noconfirm --needed docker docker-buildx
    ok "Docker packages installed: $(docker --version | sed 's/,.*//')"
}

install_macos() {
    # Installs Docker Desktop via Homebrew (no sudo needed for the install)
    if ! command -v brew &>/dev/null; then
        fail "Homebrew is required to install Docker Desktop automatically."
        echo "Either install Homebrew first:"
        echo '  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"'
        echo "and re-run this script, or download Docker Desktop manually:"
        echo "  https://www.docker.com/products/docker-desktop/"
        exit 1
    fi

    if [ ! -d "/Applications/Docker.app" ]; then
        info "Installing Docker Desktop via Homebrew (this may take a few minutes)..."
        brew install --cask docker
        ok "Docker Desktop installed"
    else
        ok "Docker Desktop is already installed"
    fi

    info "Launching Docker Desktop..."
    open -a Docker 2>/dev/null || true

    info "Waiting for the Docker daemon (up to ~2 minutes)..."
    echo "  If a macOS security prompt appears, click 'Allow' / enter your password."
    echo "  On first launch, accept the Docker service agreement."
    local i
    for i in $(seq 1 60); do
        docker info &>/dev/null && break
        sleep 2
    done
    if ! docker info &>/dev/null; then
        warn "Docker Desktop is installed but the daemon is not up yet."
        echo "  1. Open 'Docker' from your Applications folder."
        echo "  2. Complete the first-launch setup (accept the agreement, skip sign-in)."
        echo "  3. Wait until the whale icon shows 'Docker Desktop is running'."
        echo "  4. Then just run ./run_docker.sh -- no need to re-run this installer."
        exit 1
    fi
    ok "Docker daemon is up"
}

#--- Linux: start the daemon ------------------------------------------------------
start_docker_linux() {
    if command -v systemctl &>/dev/null && systemctl is-system-running &>/dev/null; then
        info "Enabling and starting the Docker service..."
        $SUDO systemctl enable --now docker
        ok "Docker service is active"
    else
        warn "systemd not available; trying 'service docker start'..."
        $SUDO service docker start || true
    fi

    # NOTE: check via $SUDO -- the user may not be in the 'docker' group yet
    # (that happens after this step), so a plain `docker info` would fail with
    # a permission error even though the daemon is perfectly healthy.
    info "Waiting for the Docker daemon..."
    local i
    for i in $(seq 1 15); do
        $SUDO docker info &>/dev/null && break
        sleep 2
    done
    if ! $SUDO docker info &>/dev/null; then
        fail "Docker daemon did not start."
        echo "Check the logs with:  journalctl -u docker --no-pager | tail -50"
        exit 1
    fi
    ok "Docker daemon is up"
}

#--- Linux: allow the user to run docker without sudo ----------------------------
add_user_to_docker_group() {
    info "Adding user '${REAL_USER}' to the 'docker' group..."
    if id -nG "$REAL_USER" | tr ' ' '\n' | grep -qx docker; then
        ok "User '${REAL_USER}' is already in the docker group"
    else
        $SUDO usermod -aG docker "$REAL_USER"
        ok "User '${REAL_USER}' added to the docker group"
    fi
}

#--- Dispatch: install, start, and verify per platform ----------------------------
if [ "$OS" == "Darwin" ]; then
    install_macos
else
    case "$DISTRO" in
        debian) install_debian ;;
        fedora) install_fedora ;;
        arch)   install_arch ;;
    esac
    start_docker_linux
    add_user_to_docker_group
fi

#--- Done ---------------------------------------------------------------------------
echo ""
echo "=============================================================="
ok  "Docker installation complete!"
echo "=============================================================="
echo ""
if [ "$OS" == "Darwin" ]; then
    echo "  Run the simulation:"
    echo "    cd $(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    echo "    ./run_docker.sh"
else
    echo "  IMPORTANT -- one last step (group changes need a new session):"
    echo "    Close this terminal and open a new one, or run:  newgrp docker"
    echo ""
    echo "  Then run the simulation:"
    echo "    cd $(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    echo "    ./run_docker.sh"
fi
echo ""
echo "  The first run builds the Geant4 Docker image (several GB,"
echo "  ~10-20 min). Afterwards everything is cached and fast."
echo ""

