#!/bin/bash
#
# CMS Model — Complete Pipeline Script
#
# Usage:
#   ./runit.sh build     - Build only (cmake + make, incremental)
#   ./runit.sh clean     - Clean rebuild (preserves data files)
#   ./runit.sh run       - Build + batch simulation (run.mac)
#   ./runit.sh digi      - Run digitisation on simulation output
#   ./runit.sh track     - Run track finder (no noise)
#   ./runit.sh track 13  - Run track finder with 13 noise hits
#   ./runit.sh plot      - Generate simulation + digitisation plots
#   ./runit.sh all       - Build + Simulate + Digitise + Plot
#   ./runit.sh scratch   - FULL pipeline: Build + Sim + Digi + Plot + Track(0) + Track(13)
#   ./runit.sh           - Default: build + interactive mode
#

set -e  # Exit on error

PROJECT_DIR="$HOME/Desktop/cmsModel"
BUILD_DIR="$PROJECT_DIR/build"
GEANT4_DIR="$HOME/software/geant4/geant4-v11.4.0"
CONDA_PREFIX="$HOME/miniforge3/envs/hep"

# Source Geant4 environment
source "$GEANT4_DIR/bin/geant4.sh"

ACTION="${1:-default}"
EXTRA_ARG="${2:-}"

#------------------------------------------------------------
# BUILD function
#------------------------------------------------------------
do_build() {
    echo ""
    echo "╔══════════════════════════════════════════╗"
    echo "║  STEP 1: BUILD                           ║"
    echo "╚══════════════════════════════════════════╝"
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"

    # Only run cmake if CMakeCache doesn't exist yet
    if [ ! -f CMakeCache.txt ]; then
        echo "--- Running cmake ---"
        cmake -DGeant4_DIR="$GEANT4_DIR/lib/cmake/Geant4" \
              -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
              ..
    fi

    echo "--- Running make ---"
    make -j$(sysctl -n hw.ncpu)
    echo "✅ Build complete"
    echo ""
}

#------------------------------------------------------------
# CLEAN BUILD function (preserves data files)
#------------------------------------------------------------
do_clean() {
    echo "=== Clean build ==="

    # SAFETY: Back up output data files before deleting build dir
    if [ -d "$BUILD_DIR" ]; then
        BACKUP_DIR="/tmp/cmsmodel_backup_$$"
        mkdir -p "$BACKUP_DIR"
        echo "--- Backing up data files to $BACKUP_DIR ---"

        # Preserve ROOT files, PDFs, log files
        for ext in root pdf log; do
            for f in "$BUILD_DIR"/*.$ext; do
                [ -f "$f" ] && cp "$f" "$BACKUP_DIR/" && echo "  Saved: $(basename $f)"
            done
        done
    fi

    rm -rf "$BUILD_DIR"
    do_build

    # Restore backed-up files
    if [ -d "$BACKUP_DIR" ]; then
        echo "--- Restoring data files ---"
        for f in "$BACKUP_DIR"/*; do
            [ -f "$f" ] && cp "$f" "$BUILD_DIR/" && echo "  Restored: $(basename $f)"
        done
        rm -rf "$BACKUP_DIR"
    fi

    echo "=== Clean build complete (data files preserved) ==="
}

#------------------------------------------------------------
# RUN SIMULATION function (batch mode)
#------------------------------------------------------------
do_run() {
    cd "$BUILD_DIR"
    echo ""
    echo "╔══════════════════════════════════════════╗"
    echo "║  STEP 2: SIMULATION (run.mac)            ║"
    echo "╚══════════════════════════════════════════╝"
    echo "Particle: $(grep '/Serc19/gun/pid' run.mac | awk '{print $2}')"
    echo "Energy:   $(grep '/Serc19/gun/energy' run.mac | awk '{print $2}') GeV"
    echo "Events:   $(grep '/run/beamOn' run.mac | awk '{print $2}')"
    echo ""
    ./cmsmodel run.mac
    echo ""
    echo "✅ Simulation complete"
    echo "ROOT output:"
    ls -lh *_run*.root 2>/dev/null || echo "  (no ROOT files found)"
    echo ""
}

#------------------------------------------------------------
# DIGITISATION function
#------------------------------------------------------------
do_digi() {
    cd "$BUILD_DIR"
    echo ""
    echo "╔══════════════════════════════════════════╗"
    echo "║  STEP 3: DIGITISATION                    ║"
    echo "╚══════════════════════════════════════════╝"

    # Find the latest simulation ROOT file
    ROOT_FILE=$(ls -t *_run*.root 2>/dev/null | grep -v 'xwoecl\|skeleton\|output' | head -1)
    if [ -z "$ROOT_FILE" ]; then
        echo "ERROR: No simulation ROOT file found in $BUILD_DIR"
        echo "Run the simulation first: ./runit.sh run"
        exit 1
    fi

    echo "Input: $ROOT_FILE"
    echo "ECAL noise=40 MeV, threshold=200 MeV"
    echo "HCAL noise=25 MeV, threshold=100 MeV"

    # Create input file list
    echo "$ROOT_FILE 10000" > test_pion_klong.log

    # Run digitisation
    ./ecal_hcal_digitisation 40 200 25 100

    echo ""
    echo "✅ Digitisation complete"
    echo "Output:"
    ls -lh *xwoecl*.root 2>/dev/null
    echo ""
}

#------------------------------------------------------------
# PLOT function
#------------------------------------------------------------
do_plot() {
    cd "$BUILD_DIR"
    echo ""
    echo "╔══════════════════════════════════════════╗"
    echo "║  STEP 4: PLOTS                           ║"
    echo "╚══════════════════════════════════════════╝"

    # Find simulation file
    SIM_FILE=$(ls -t *_run*.root 2>/dev/null | grep -v 'xwoecl\|skeleton\|output' | head -1)
    DIGI_FILE=$(ls -t *xwoecl*.root 2>/dev/null | head -1)

    if [ -z "$SIM_FILE" ]; then
        echo "ERROR: No simulation ROOT file found"
        exit 1
    fi

    if [ -n "$DIGI_FILE" ]; then
        echo "Simulation: $SIM_FILE"
        echo "Digitisation: $DIGI_FILE"
        root -l -b -q "../plot_results.C(\"$SIM_FILE\", \"$DIGI_FILE\")"
    else
        echo "Simulation: $SIM_FILE (no digi file)"
        root -l -b -q "../plot_results.C(\"$SIM_FILE\")"
    fi

    echo ""
    echo "✅ Plots complete"
    ls -lh plots_*.pdf 2>/dev/null
    echo ""
}

#------------------------------------------------------------
# TRACK FINDER function
#------------------------------------------------------------
do_track() {
    cd "$BUILD_DIR"
    NOISE="${1:-0}"
    echo ""
    echo "╔══════════════════════════════════════════╗"
    echo "║  STEP 5: TRACK FINDER (noise=$NOISE)     ║"
    echo "╚══════════════════════════════════════════╝"

    # Find simulation file for track finder
    SIM_FILE=$(ls -t *_run*.root 2>/dev/null | grep -v 'xwoecl\|skeleton\|output' | head -1)
    if [ -z "$SIM_FILE" ]; then
        echo "ERROR: No simulation ROOT file found"
        exit 1
    fi

    echo "$SIM_FILE 10000" > skeleton_track.log
    echo "Input: $SIM_FILE"
    echo "Noise: $NOISE hits/event"
    echo ""

    ./track_finder $NOISE

    # Save with noise label
    if [ "$NOISE" -gt 0 ] 2>/dev/null; then
        cp track_resolution.pdf "track_resolution_noise${NOISE}.pdf" 2>/dev/null
        cp track_event_display.pdf "track_event_display_noise${NOISE}.pdf" 2>/dev/null
        echo ""
        echo "✅ Track finder complete (with noise=$NOISE)"
        echo "Saved as: track_resolution_noise${NOISE}.pdf"
        echo "          track_event_display_noise${NOISE}.pdf"
    else
        cp track_resolution.pdf "track_resolution_no_noise.pdf" 2>/dev/null
        cp track_event_display.pdf "track_event_display_no_noise.pdf" 2>/dev/null
        echo ""
        echo "✅ Track finder complete (no noise)"
        echo "Saved as: track_resolution_no_noise.pdf"
        echo "          track_event_display_no_noise.pdf"
    fi
    echo ""
}

#------------------------------------------------------------
# Main dispatch
#------------------------------------------------------------
case "$ACTION" in
    build)
        do_build
        ;;
    clean)
        do_clean
        ;;
    run)
        do_build
        do_run
        ;;
    digi)
        do_digi
        ;;
    track)
        do_track "$EXTRA_ARG"
        ;;
    plot)
        do_plot
        ;;
    all)
        do_build
        do_run
        do_digi
        do_plot
        ;;
    scratch)
        echo ""
        echo "╔══════════════════════════════════════════════════╗"
        echo "║  FULL PIPELINE: Build → Sim → Digi → Plot → Track ║"
        echo "╚══════════════════════════════════════════════════╝"
        echo ""
        do_build
        do_run
        do_digi
        do_plot
        do_track 0          # Track finder without noise
        do_track 13         # Track finder with 1% noise
        echo ""
        echo "╔══════════════════════════════════════════════════╗"
        echo "║  ✅ ALL DONE! Output files in build/:            ║"
        echo "╠══════════════════════════════════════════════════╣"
        echo "║  Simulation:                                     ║"
        echo "║    test_run10.root                                ║"
        echo "║  Digitisation:                                    ║"
        echo "║    test_pion_klong_xwoecl_40_200_25_100.root      ║"
        echo "║  Plots:                                           ║"
        echo "║    plots_simulation.pdf                           ║"
        echo "║    plots_digitisation.pdf                         ║"
        echo "║  Track Finder:                                    ║"
        echo "║    output_skeleton.root                           ║"
        echo "║    track_resolution_no_noise.pdf                  ║"
        echo "║    track_event_display_no_noise.pdf               ║"
        echo "║    track_resolution_noise13.pdf                   ║"
        echo "║    track_event_display_noise13.pdf                ║"
        echo "╚══════════════════════════════════════════════════╝"
        ;;
    default)
        do_build
        cd "$BUILD_DIR"
        echo "=== Starting interactive mode ==="
        ./cmsmodel
        ;;
    *)
        echo ""
        echo "CMS Model — Complete Pipeline Script"
        echo ""
        echo "Usage: ./runit.sh <command> [args]"
        echo ""
        echo "Commands:"
        echo "  build        Build only (incremental)"
        echo "  clean        Clean rebuild (preserves data)"
        echo "  run          Build + simulate (batch, run.mac)"
        echo "  digi         Digitise simulation output"
        echo "  plot         Generate plots (sim + digi)"
        echo "  track [N]    Track finder (N = noise hits, default 0)"
        echo "  all          Build + Sim + Digi + Plot"
        echo "  scratch      FULL pipeline: all + track(0) + track(13)"
        echo "  (no args)    Build + interactive mode"
        echo ""
        echo "Examples:"
        echo "  ./runit.sh scratch       # Run EVERYTHING from scratch"
        echo "  ./runit.sh track 50      # Track finder with 50 noise hits"
        echo "  ./runit.sh plot          # Regenerate plots only"
        echo ""
        exit 1
        ;;
esac
