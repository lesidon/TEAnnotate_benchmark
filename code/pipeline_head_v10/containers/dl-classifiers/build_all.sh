#!/usr/bin/env bash
# Build the five per-tool DL-classifier Apptainer images on the HPC login node.
# Run from the REPOSITORY ROOT (the .def %files paths resolve from there):
#   bash containers/dl-classifiers/build_all.sh
#
# Each image builds independently, so a transient network stall only affects one and
# you can re-run just that tool. Set APPTAINER_TMPDIR/CACHEDIR to beegfs (not the
# read-only /trinity home). Override the build command via $BUILD if --fakeroot is
# not permitted (e.g. BUILD="apptainer build" or BUILD="apptainer build --remote").
set -u

export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-/beegfs/home/$USER/tmp/apptainer_tmp}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-/beegfs/home/$USER/tmp/apptainer_cache}"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

BUILD="${BUILD:-apptainer build --fakeroot}"
defdir="containers/dl-classifiers"
rc=0

for tool in teclass2 create deepte neuralte terrier; do
    echo "==================== building ${tool} ===================="
    if $BUILD "teannotate-${tool}.sif" "${defdir}/${tool}.def"; then
        echo "OK: teannotate-${tool}.sif"
    else
        echo "FAILED: ${tool} — re-run: ${BUILD} teannotate-${tool}.sif ${defdir}/${tool}.def"
        rc=1
    fi
done

echo "Done. Built .sif files (if any) are in $(pwd)."
exit $rc
