#!/bin/bash

EXEC=./main2d.gnu.MPI.ex

# Config Files
CONFIGS=(
    configs/toro1
    configs/toro2
    configs/toro3
    configs/toro4
    configs/toro5
)

# No AMR Resolutions
RESOLUTIONS=(104 208 416)

# Fixed Resolution for AMR study
BASE_R=104

# Root Output Directory
OUTROOT=../../../../../../CPP/data/one_dimensional_test

# Toro Study Do Loop
for INPUT in "${CONFIGS[@]}"
do
    TESTNAME=$(basename "$INPUT")
    echo "Running test: $TESTNAME"

    # [1] NO AMR CONVERGENCE STUDY
    OUTDIR_NO_AMR=$OUTROOT/No_AMR/$TESTNAME
    mkdir -p "$OUTDIR_NO_AMR"

    for R in "${RESOLUTIONS[@]}"
    do
        echo "Running No_AMR resolution $R"

        mpirun -np 1 $EXEC $INPUT \
            amr.n_cell="$R 8" \
            amr.max_level=0 \
            adv.do_reflux=0 \
            amr.plot_file=$OUTDIR_NO_AMR/plt
    done


    # [2] AMR STUDY (R fixed = BASE_R)
    OUTDIR_AMR=$OUTROOT/with_AMR/$TESTNAME
    mkdir -p "$OUTDIR_AMR"

    for ML in 0 1 2
    do
        echo "Running AMR base=$BASE_R max_level=$ML"

        mpirun -np 1 $EXEC $INPUT \
            amr.n_cell="$BASE_R 8" \
            amr.max_level=$ML \
            adv.do_reflux=1 \
            amr.regrid_int=2 \
            amr.n_error_buf=2 \
            amr.plot_file=$OUTDIR_AMR/plt
    done

done

