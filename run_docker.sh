#!/bin/bash
# IntegrateALL Docker Wrapper Script
# Simplifies running the pipeline in Docker

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
print_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Default values
CORES=4
DOCKER_IMAGE="integrateall:latest"

# Usage information
usage() {
    cat << EOF
Usage: $0 [COMMAND] [OPTIONS]

Commands:
    setup       Run one-time setup (downloads ~21GB of reference data)
    run         Run the analysis pipeline
    interactive Start interactive shell in container
    test        Run a dry-run to check configuration

Options:
    -c, --cores NUM     Number of CPU cores to use (default: 4)
    -h, --help          Show this help message

Examples:
    $0 setup                    # Run setup with default 4 cores
    $0 run -c 8                 # Run analysis with 8 cores
    $0 interactive              # Start interactive shell
    $0 test                     # Dry-run to test configuration

EOF
    exit 1
}

# Check if Docker is installed and running
check_docker() {
    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed. Please install Docker first."
        exit 1
    fi

    if ! docker info &> /dev/null; then
        print_error "Docker daemon is not running. Please start Docker."
        exit 1
    fi
}

# Check if required files exist
check_files() {
    if [ ! -f "config.yaml" ]; then
        print_error "config.yaml not found in current directory"
        exit 1
    fi

    if [ ! -f "samples.csv" ] && [ "$1" != "setup" ]; then
        print_warn "samples.csv not found - required for running analysis"
    fi
}

# Run setup
run_setup() {
    print_info "Starting IntegrateALL setup (this may take 1-3 hours)..."
    print_info "Downloading ~21GB of reference data..."

    mkdir -p refs

    docker run --rm \
        -v "$(pwd)/refs:/IntegrateALL/refs" \
        -v "$(pwd)/config.yaml:/IntegrateALL/config.yaml" \
        -v "$(pwd):/IntegrateALL" \
        -w /IntegrateALL \
        "$DOCKER_IMAGE" \
        snakemake --snakefile setup.smk --cores "$CORES" --use-conda

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

    docker run --rm \
        -v "$(pwd)/data:/IntegrateALL/data" \
        -v "$(pwd)/refs:/IntegrateALL/refs" \
        -v "$(pwd):/IntegrateALL" \
        -w /IntegrateALL \
        "$DOCKER_IMAGE" \
        snakemake --cores "$CORES" --use-conda

    print_info "Analysis completed successfully!"
    print_info "Results are in the following directories:"
    print_info "  - Final_classification/: Classification results"
    print_info "  - interactive_output/: HTML reports"
    print_info "  - fusions/: Fusion detection results"
}

# Run interactive shell
run_interactive() {
    print_info "Starting interactive shell in IntegrateALL container..."
    print_info "Conda environment 'IALL' is already activated."
    print_info "Type 'exit' to leave the container."

    docker run -it --rm \
        -v "$(pwd):/IntegrateALL" \
        -w /IntegrateALL \
        "$DOCKER_IMAGE" \
        bash
}

# Run dry-run test
run_test() {
    print_info "Running dry-run to test configuration..."

    docker run --rm \
        -v "$(pwd):/IntegrateALL" \
        -w /IntegrateALL \
        "$DOCKER_IMAGE" \
        snakemake -n --cores "$CORES"

    print_info "Dry-run completed. Configuration looks good!"
}

# Parse command line arguments
COMMAND=""
while [[ $# -gt 0 ]]; do
    case $1 in
        setup|run|interactive|test)
            COMMAND="$1"
            shift
            ;;
        -c|--cores)
            CORES="$2"
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

check_docker
check_files "$COMMAND"

case $COMMAND in
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
