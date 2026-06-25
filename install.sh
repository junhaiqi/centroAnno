#!/usr/bin/env bash
# =============================================================================
# install.sh — One-command deployment for centroAnno
# =============================================================================
# Automatically detects the OS, installs missing system dependencies,
# installs Python packages (pip or conda), compiles the C++ binary,
# and verifies the installation.
#
# Usage:
#   ./install.sh
#   ./install.sh --conda          # Force conda for Python packages
#   ./install.sh --no-sudo        # Skip system package installation (manual)
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

info()  { echo -e "${GREEN}[INFO]${NC}  $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error() { echo -e "${RED}[ERROR]${NC} $*" >&2; }
step()  { echo -e "${BLUE}[STEP]${NC}  $*"; }

# ---------------------------------------------------------------------------
# Parse CLI flags
# ---------------------------------------------------------------------------
FORCE_CONDA=0
NO_SUDO=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --conda)    FORCE_CONDA=1; shift ;;
        --no-sudo)  NO_SUDO=1; shift ;;
        -h|--help)
            echo "Usage: $0 [--conda] [--no-sudo]"
            echo "  --conda     Force conda for Python package installation"
            echo "  --no-sudo   Skip system package installation (assume already installed)"
            exit 0
            ;;
        *) error "Unknown option: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Detect OS
# ---------------------------------------------------------------------------
OS=""
if [[ -f /etc/os-release ]]; then
    . /etc/os-release
    case "$ID" in
        ubuntu|debian) OS="debian" ;;
        centos|rhel|fedora|rocky|almalinux) OS="redhat" ;;
        *) OS="unknown" ;;
    esac
elif [[ "$OSTYPE" == "darwin"* ]]; then
    OS="macos"
else
    OS="unknown"
fi

info "Detected OS family: ${OS}"

# ---------------------------------------------------------------------------
# Check sudo availability
# ---------------------------------------------------------------------------
HAS_SUDO=0
if command -v sudo &>/dev/null && sudo -n true 2>/dev/null; then
    HAS_SUDO=1
elif command -v sudo &>/dev/null; then
    info "sudo is available but may require a password prompt."
fi

# ---------------------------------------------------------------------------
# Helper: install system packages
# ---------------------------------------------------------------------------
install_system_deps() {
    if [[ ${NO_SUDO} -eq 1 ]]; then
        warn "--no-sudo set: skipping system package installation."
        warn "Please ensure the following are installed manually:"
        warn "  - g++ (>= 9.3.0), make, git"
        warn "  - zlib development headers (zlib1g-dev or zlib-devel)"
        warn "  - OpenMP runtime (libgomp1 or libomp)"
        return
    fi

    if [[ ${HAS_SUDO} -eq 0 ]]; then
        warn "sudo not available or not passwordless. Skipping automatic system install."
        warn "Please install the following manually:"
        warn "  - g++ (>= 9.3.0), make, git"
        warn "  - zlib development headers (zlib1g-dev or zlib-devel)"
        warn "  - OpenMP runtime (libgomp1 or libomp)"
        return
    fi

    step "Installing system dependencies..."
    case "$OS" in
        debian)
            sudo apt-get update -qq
            sudo apt-get install -y -qq \
                build-essential g++ make git \
                zlib1g-dev libgomp1 \
                python3 python3-pip python3-venv \
                2>/dev/null || true
            ;;
        redhat)
            if command -v dnf &>/dev/null; then
                sudo dnf install -y \
                    gcc-c++ make git \
                    zlib-devel libgomp \
                    python3 python3-pip \
                    2>/dev/null || true
            else
                sudo yum install -y \
                    gcc-c++ make git \
                    zlib-devel libgomp \
                    python3 python3-pip \
                    2>/dev/null || true
            fi
            ;;
        macos)
            if ! command -v brew &>/dev/null; then
                error "Homebrew not found. Please install Homebrew first: https://brew.sh"
                exit 1
            fi
            brew install gcc make git zlib libomp python3 2>/dev/null || true
            ;;
        *)
            warn "Unknown OS. Please install the following manually:"
            warn "  - g++ (>= 9.3.0), make, git"
            warn "  - zlib development headers"
            warn "  - OpenMP runtime library"
            ;;
    esac
}

# ---------------------------------------------------------------------------
# Helper: check individual system binaries
# ---------------------------------------------------------------------------
check_binary() {
    local name="$1"
    if command -v "$name" &>/dev/null; then
        info "  ✓ ${name}: $(command -v "$name")"
        return 0
    else
        error "  ✗ ${name}: not found"
        return 1
    fi
}

# ---------------------------------------------------------------------------
# Step 1: System dependency check
# ---------------------------------------------------------------------------
step "Step 1/5: Checking system dependencies..."

MISSING_SYS=0
check_binary g++  || MISSING_SYS=1
check_binary make || MISSING_SYS=1
check_binary git  || MISSING_SYS=1

# Check g++ version
gpp_version=$(g++ --version 2>/dev/null | head -n1 | grep -oP '\d+\.\d+\.\d+' | head -n1 || echo "0")
if [[ -n "$gpp_version" ]]; then
    major_ver=$(echo "$gpp_version" | cut -d. -f1)
    if [[ "$major_ver" -lt 9 ]]; then
        warn "g++ version ${gpp_version} detected. centroAnno requires >= 9.3.0."
        MISSING_SYS=1
    else
        info "  ✓ g++ version: ${gpp_version} (>= 9.3.0 OK)"
    fi
fi

# Check zlib header
if [[ -f /usr/include/zlib.h || -f /usr/local/include/zlib.h ]]; then
    info "  ✓ zlib headers found"
else
    warn "  ✗ zlib headers not found (zlib1g-dev / zlib-devel)"
    MISSING_SYS=1
fi

# Check OpenMP
if g++ -fopenmp -E - </dev/null &>/dev/null; then
    info "  ✓ OpenMP support available"
else
    warn "  ✗ OpenMP support missing (libgomp1 / libomp)"
    MISSING_SYS=1
fi

if [[ ${MISSING_SYS} -eq 1 ]]; then
    install_system_deps
fi

# Re-check after install
MISSING_SYS=0
check_binary g++  || MISSING_SYS=1
check_binary make || MISSING_SYS=1
if [[ ${MISSING_SYS} -eq 1 ]]; then
    error "Required system tools still missing after installation attempt."
    error "Please install them manually and rerun this script."
    exit 1
fi

# ---------------------------------------------------------------------------
# Step 2: Python environment & packages
# ---------------------------------------------------------------------------
step "Step 2/5: Checking Python environment..."

# Determine whether to use conda or pip
USE_CONDA=0
if [[ ${FORCE_CONDA} -eq 1 ]]; then
    USE_CONDA=1
    info "--conda flag set: forcing conda installation."
else
    # Auto-detect conda
    if command -v conda &>/dev/null; then
        USE_CONDA=1
        info "conda detected at $(command -v conda)"
    elif command -v mamba &>/dev/null; then
        USE_CONDA=1
        info "mamba detected at $(command -v mamba)"
    fi
fi

# Check python3
if ! command -v python3 &>/dev/null; then
    error "python3 not found. Please install Python 3."
    exit 1
fi
info "  ✓ python3: $(python3 --version)"

# Required Python packages (pip_name:import_name)
declare -A PY_IMPORT_MAP=(
    [biopython]=Bio
    [numpy]=numpy
    [scipy]=scipy
    [matplotlib]=matplotlib
    [pandas]=pandas
)
PYTHON_DEPS=(biopython numpy scipy matplotlib pandas)

install_python_deps() {
    if [[ ${USE_CONDA} -eq 1 ]]; then
        step "Installing Python packages via conda..."
        local conda_cmd="conda"
        command -v mamba &>/dev/null && conda_cmd="mamba"

        # Check if we're in a conda environment
        if [[ -z "${CONDA_DEFAULT_ENV:-}" ]]; then
            warn "No active conda environment detected."
            warn "Packages will be installed into the base environment."
            warn "To use a specific environment, activate it first and rerun."
        fi

        $conda_cmd install -y -c conda-forge -c bioconda "${PYTHON_DEPS[@]}" 2>/dev/null || {
            error "conda install failed. Falling back to pip..."
            USE_CONDA=0
            install_python_deps
        }
    else
        step "Installing Python packages via pip..."
        # Upgrade pip first
        python3 -m pip install --upgrade pip setuptools wheel -q 2>/dev/null || true
        python3 -m pip install "${PYTHON_DEPS[@]}" -q
    fi
}

# Check which packages are missing
MISSING_PY=()
for pkg in "${PYTHON_DEPS[@]}"; do
    import_name="${PY_IMPORT_MAP[$pkg]}"
    if python3 -c "import ${import_name}" 2>/dev/null; then
        info "  ✓ ${pkg}: $(python3 -c "import ${import_name}; print(${import_name}.__version__)" 2>/dev/null || echo 'installed')"
    else
        warn "  ✗ ${pkg}: not installed"
        MISSING_PY+=("$pkg")
    fi
done

if [[ ${#MISSING_PY[@]} -gt 0 ]]; then
    install_python_deps
    # Re-verify
    for pkg in "${MISSING_PY[@]}"; do
        import_name="${PY_IMPORT_MAP[$pkg]}"
        if python3 -c "import ${import_name}" 2>/dev/null; then
            info "  ✓ ${pkg}: installed successfully"
        else
            error "  ✗ ${pkg}: still missing after installation"
            exit 1
        fi
    done
fi

# ---------------------------------------------------------------------------
# Step 3: Initialize git submodule (spoa)
# ---------------------------------------------------------------------------
step "Step 3/5: Initializing git submodule (spoa)..."

if [[ -d .git ]]; then
    if git submodule status lib/spoa &>/dev/null && [[ -n "$(git submodule status lib/spoa)" ]]; then
        info "  ✓ spoa submodule already initialized"
    else
        git submodule update --init --recursive lib/spoa
        info "  ✓ spoa submodule initialized"
    fi
else
    warn "  Not a git repository. Cloning spoa manually..."
    rm -rf lib/spoa
    git clone https://github.com/rvaser/spoa.git lib/spoa
    cd lib/spoa && git checkout 4.1.5 && cd ../..
    info "  ✓ spoa cloned manually"
fi

# ---------------------------------------------------------------------------
# Step 4: Compile centroAnno
# ---------------------------------------------------------------------------
step "Step 4/5: Compiling centroAnno..."

make clean 2>/dev/null || true
make -j$(nproc 2>/dev/null || echo 4)

if [[ ! -x ./centroAnno ]]; then
    error "Compilation failed. Binary not found: ./centroAnno"
    exit 1
fi
info "  ✓ centroAnno binary compiled successfully"

# ---------------------------------------------------------------------------
# Step 5: Verify installation
# ---------------------------------------------------------------------------
step "Step 5/5: Verifying installation..."

./centroAnno 2>&1 | head -n1 || true

info ""
info "========================================================================"
info "Installation complete!"
info "========================================================================"
info "Binary:   $(pwd)/centroAnno"
info "Scripts:  $(pwd)/script/"
info "Pipeline: $(pwd)/run_centroAnno.sh"
info ""
info "Quick test:"
info "  ./centroAnno -o test_output -x anno-asm example/cen21.fa"
info ""
info "Full pipeline:"
info "  ./run_centroAnno.sh -i example/cen21.fa -o results/cen21"
info "========================================================================"
