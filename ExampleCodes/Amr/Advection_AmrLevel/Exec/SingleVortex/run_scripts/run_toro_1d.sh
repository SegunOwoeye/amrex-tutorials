#!/bin/bash

EXEC=./main2d.gnu.MPI.ex

# Config Files
CONFIGS=(
    configs/toro_tests/toro1
    configs/toro_tests/toro2
    configs/toro_tests/toro3
    configs/toro_tests/toro4
    configs/toro_tests/toro5
)

# No AMR Resolutions
RESOLUTIONS=(100 200 400)

# Fixed Resolution for AMR study
BASE_R=100

# Root Output Directory
OUTROOT=../../../../../../CPP/data/one_dimensional_test

# 1D Toro Study Do Loop
for INPUT in "${CONFIGS[@]}"
do
    TESTNAME=$(basename "$INPUT")
    echo "Running 1D validation for $TESTNAME"

    # [1] NO AMR CONVERGENCE STUDY
    OUTDIR_NO_AMR=$OUTROOT/No_AMR/$TESTNAME
    mkdir -p "$OUTDIR_NO_AMR"

    rm -rf "$OUTDIR_NO_AMR"/plt*
    rm -rf "$OUTDIR_NO_AMR"/*.txt
    rm -rf "$OUTDIR_NO_AMR"/*.png
    for R in "${RESOLUTIONS[@]}"
    do
        echo "Running No_AMR resolution $R"

      

        mpirun -np 1 $EXEC $INPUT \
            amr.n_cell="$R 4 1" \
            amr.max_level=0 \
            amr.blocking_factor=2 \
            adv.do_reflux=0 \
            adv.cfl=0.3 \
            amr.plot_file=$OUTDIR_NO_AMR/plt
    done


    # [2] AMR STUDY (R fixed = BASE_R)
    OUTDIR_AMR=$OUTROOT/with_AMR/$TESTNAME
    mkdir -p "$OUTDIR_AMR"

    rm -rf "$OUTDIR_AMR"/plt*
    rm -rf "$OUTDIR_AMR"/*.txt
    rm -rf "$OUTDIR_AMR"/*.png
    for ML in 0 1 2
    do
        echo "Running AMR base=$BASE_R max_level=$ML"

        

        mpirun -np 1 $EXEC $INPUT \
            amr.n_cell="$BASE_R 4 1" \
            amr.max_level=$ML \
            amr.blocking_factor=2 \
            adv.do_reflux=1 \
            amr.regrid_int=2 \
            adv.max_rhograd_lev=$ML \
            adv.cfl=0.2 \
            amr.n_error_buf=4 \
            amr.plot_file=$OUTDIR_AMR/plt
    done

done

