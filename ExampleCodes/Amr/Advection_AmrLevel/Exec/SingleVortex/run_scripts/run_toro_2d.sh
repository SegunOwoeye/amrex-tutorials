#!/bin/bash

EXEC=./main2d.gnu.MPI.ex

CONFIGS=(
    configs/toro_tests/toro1
    configs/toro_tests/toro2
    configs/toro_tests/toro3
    configs/toro_tests/toro4
    configs/toro_tests/toro5
)

PROB_TYPES=(0 1 2)

# Change to 400
RES=400

# CPU Cores
CCORES=6

OUTROOT=../../../../../../CPP/data/two_dimensional_test


for INPUT in "${CONFIGS[@]}"
do
    TESTNAME=$(basename "$INPUT")
    echo "Running 2D validation for $TESTNAME"

    for PTYPE in "${PROB_TYPES[@]}"
    do
        echo "Running prob.type = $PTYPE"

        # Set periodicity depending on split direction
        if [ "$PTYPE" -eq 0 ]; then
            GEOM_PERIODIC="0 1"
        elif [ "$PTYPE" -eq 1 ]; then
            GEOM_PERIODIC="1 0"
        else
            GEOM_PERIODIC="0 0"
        fi

        EXTRA=""
        if [ "$PTYPE" -eq 2 ]; then
            # Setting 45 degree angle
            EXTRA="prob.nx=0.707107 prob.ny=0.707107"
        fi

        # No AMR run
        OUTDIR_NO_AMR=$OUTROOT/${TESTNAME}/type${PTYPE}/No_AMR
        mkdir -p "$OUTDIR_NO_AMR"
        rm -rf "$OUTDIR_NO_AMR"/plt*
        rm -rf "$OUTDIR_NO_AMR"/*.txt
        rm -rf "$OUTDIR_NO_AMR"/*.png

        mpirun -np $CCORES $EXEC $INPUT \
            prob.type=$PTYPE \
            geometry.is_periodic="$GEOM_PERIODIC" \
            $EXTRA \
            amr.n_cell="$RES $RES" \
            amr.max_level=0 \
            adv.do_reflux=0 \
            adv.cfl=0.3 \
            amr.v=0 \
            adv.v=0 \
            amr.plot_file=$OUTDIR_NO_AMR/plt


        # With AMR run
        OUTDIR_AMR=$OUTROOT/${TESTNAME}/type${PTYPE}/With_AMR
        mkdir -p "$OUTDIR_AMR"
        rm -rf "$OUTDIR_AMR"/plt*
        rm -rf "$OUTDIR_AMR"/*.txt
        rm -rf "$OUTDIR_AMR"/*.png

        mpirun -np $CCORES $EXEC $INPUT \
            prob.type=$PTYPE \
            geometry.is_periodic="$GEOM_PERIODIC" \
            $EXTRA \
            amr.n_cell="$RES $RES" \
            amr.max_level=1 \
            adv.do_reflux=1 \
            amr.regrid_int=2 \
            amr.n_error_buf=2 \
            adv.max_rhograd_lev=1 \
            adv.cfl=0.3 \
            amr.v=0 \
            adv.v=0 \
            amr.plot_file=$OUTDIR_AMR/plt
    done
done


