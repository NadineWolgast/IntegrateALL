#!/bin/bash
# IntegrateALL Singularity Wrapper Script
# Simplifies running the pipeline with Singularity/Apptainer

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Function to print colored messages
print_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
print_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Default values
CORES=4
CONTAINER="integrateall.sif"

# Usage information
usage() {
    cat << EOF
Usage: $0 [COMMAND] [OPTIONS]

Commands:
    build       Build Singularity container from definition file
    setup       Run one-time setup (downloads ~21GB of reference data)
    run         Run the analysis pipeline
    interactive Start interactive shell in container
    test        Run a dry-run to check configuration

Options:
    -c, --cores NUM     Number of CPU cores to use (default: 4)
    -s, --sif FILE      Container file to use (default: integrateall.sif)
    -h, --help          Show this help message

Examples:
    $0 build                        # Build container (requires sudo/fakeroot)
    $0 setup                        # Run setup with default 4 cores
    $0 run -c 8                     # Run analysis with 8 cores
    $0 interactive                  # Start interactive shell
    $0 test                         # Dry-run to test configuration

EOF
    exit 1
}

# Check if Singularity/Apptainer is installed
check_singularity() {
    if command -v singularity &> /dev/null; then
        SINGULARITY_CMD="singularity"
    elif command -v apptainer &> /dev/null; then
        SINGULARITY_CMD="apptainer"
    else
        print_error "Neither Singularity nor Apptainer is installed."
        print_info "Install with: sudo apt-get install singularity-container"
        print_info "Or: https://apptainer.org/docs/admin/main/installation.html"
        exit 1
    fi
}

# Check if required files exist
check_files() {
    if [ ! -f "config.yaml" ]; then
        print_error "config.yaml not found in current directory"
        exit 1
    fi

    if [ "$1" != "build" ] && [ ! -f "$CONTAINER" ]; then
        print_error "Container $CONTAINER not found"
        print_info "Run '$0 build' to create the container first"
        exit 1
    fi

    if [ ! -f "samples.csv" ] && [ "$1" = "run" ]; then
        print_warn "samples.csv not found - required for running analysis"
    fi
}

# Build container
build_container() {
    print_info "Building Singularity container..."

    if [ ! -f "integrateall.def" ]; then
        print_error "integrateall.def not found"
        exit 1
    fi

    # Check if we can use fakeroot or need sudo
    if $SINGULARITY_CMD --version | grep -q "apptainer\|singularity"; then
        print_info "Trying to build with fakeroot (no sudo needed)..."
        if $SINGULARITY_CMD build --fakeroot "$CONTAINER" integrateall.def 2>/dev/null; then
            print_info "Container built successfully with fakeroot!"
            return 0
        else
            print_warn "Fakeroot build failed, trying with sudo..."
        fi
    fi

    # Fall back to sudo
    if ! sudo -v &>/dev/null; then
        print_error "sudo access required to build container"
        exit 1
    fi

    sudo $SINGULARITY_CMD build "$CONTAINER" integrateall.def

    print_info "Container built successfully: $CONTAINER"
    print_info "Container size: $(du -h $CONTAINER | cut -f1)"
}

# Run setup
run_setup() {
    print_info "Starting IntegrateALL setup (this may take 1-3 hours)..."
    print_info "Downloading ~21GB of reference data..."

    mkdir -p refs

    $SINGULARITY_CMD exec \
        --bind "$(pwd)/refs:/IntegrateALL/refs" \
        --bind "$(pwd)/config.yaml:/IntegrateALL/config.yaml" \
        --bind "$(pwd):/IntegrateALL" \
        --pwd /IntegrateALL \
        "$CONTAINER" \
        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                 conda activate integrateall && \
                 snakemake --snakefile setup.smk --cores $CORES --use-conda"

    print_info "Setup completed successfully!"
}

# Run analysis
run_analysis() {
    print_info "Starting IntegrateALL analysis with $CORES cores..."

    # Check if refs exist
    if [ ! -d "refs/GATK" ]; then
        print_error "Reference data not found. Please run '$0 setup' first."
        exit 1
    fi

    mkdir -p data/samples output

    $SINGULARITY_CMD exec \
        --bind "$(pwd)/data:/IntegrateALL/data" \
        --bind "$(pwd)/refs:/IntegrateALL/refs" \
        --bind "$(pwd):/IntegrateALL" \
        --pwd /IntegrateALL \
        "$CONTAINER" \
        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                 conda activate integrateall && \
                 snakemake --cores $CORES --use-conda"

    print_info "Analysis completed successfully!"
    print_info "Results are in the following directories:"
    print_info "  - Final_classification/: Classification results"
    print_info "  - interactive_output/: HTML reports"
    print_info "  - fusions/: Fusion detection results"
}

# Run interactive shell
run_interactive() {
    print_info "Starting interactive shell in IntegrateALL container..."
    print_info "Conda environment 'integrateall' will be activated."
    print_info "Type 'exit' to leave the container."

    $SINGULARITY_CMD shell \
        --bind "$(pwd):/IntegrateALL" \
        --pwd /IntegrateALL \
        "$CONTAINER"
}

# Run dry-run test
run_test() {
    print_info "Running dry-run to test configuration..."

    $SINGULARITY_CMD exec \
        --bind "$(pwd):/IntegrateALL" \
        --pwd /IntegrateALL \
        "$CONTAINER" \
        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                 conda activate integrateall && \
                 snakemake -n --cores $CORES"

    print_info "Dry-run completed. Configuration looks good!"
}

# Parse command line arguments
COMMAND=""
while [[ $# -gt 0 ]]; do
    case $1 in
        build|setup|run|interactive|test)
            COMMAND="$1"
            shift
            ;;
        -c|--cores)
            CORES="$2"
            shift 2
            ;;
        -s|--sif)
            CONTAINER="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            print_error "Unknown option: $1"
            usage
            ;;
    esac
done

# Main execution
if [ -z "$COMMAND" ]; then
    print_error "No command specified"
    usage
fi

check_singularity
check_files "$COMMAND"

case $COMMAND in
    build)
        build_container
        ;;
    setup)
        run_setup
        ;;
    run)
        run_analysis
        ;;
    interactive)
        run_interactive
        ;;
    test)
        run_test
        ;;
    *)
        print_error "Unknown command: $COMMAND"
        usage
        ;;
esac
