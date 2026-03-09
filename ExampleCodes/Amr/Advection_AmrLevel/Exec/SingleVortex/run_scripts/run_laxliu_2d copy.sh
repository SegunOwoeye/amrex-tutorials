#!/bin/bash

EXEC=./main2d.gnu.MPI.ex
INPUT=configs/lax_and_liu_tests/laxliu_2d 

OUTDIR=../../../../../../CPP/data/timing_test
mkdir -p $OUTDIR

CSV=$OUTDIR/timing.csv
echo "mode,resolution,amr_level,cpus,time" > $CSV


# [1] Resolution list (doubling)
RESOLUTIONS=(256 512 1024)

# CPU scaling
CPUS=(1 2 4 8 12 16)


# [2] Uniform grid timing
for RES in "${RESOLUTIONS[@]}"
do
    echo "Uniform grid N=$RES"

    for C in "${CPUS[@]}"
    do
        echo "CPUs=$C"

        START=$(date +%s.%N)

        mpirun -np $C $EXEC $INPUT \
            prob.type=3 \
            geometry.is_periodic="0 0" \
            amr.n_cell="$RES $RES" \
            amr.max_level=0 \
            amr.plot_files_output=0 \
            adv.v=0 \
            amr.v=0

        END=$(date +%s.%N)

        TIME=$(echo "$END - $START" | bc)

        echo "uniform,$RES,0,$C,$TIME" >> $CSV
    done
done


# [3] AMR equivalent resolution tests
BASE_RES=256

for LEVEL in 1 2
do
    EFFECTIVE=$(( BASE_RES * (2**LEVEL) ))

    echo "AMR test: base=$BASE_RES levels=$LEVEL effective=$EFFECTIVE"

    for C in "${CPUS[@]}"
    do

        START=$(date +%s.%N)

        mpirun -np $C $EXEC $INPUT \
            prob.type=3 \
            geometry.is_periodic="0 0" \
            amr.n_cell="$BASE_RES $BASE_RES" \
            amr.max_level=$LEVEL \
            amr.regrid_int=2 \
            amr.n_error_buf=2 \
            adv.max_rhograd_lev=2 \
            amr.plot_files_output=0 \
            adv.v=0 \
            amr.v=0

        END=$(date +%s.%N)

        TIME=$(echo "$END - $START" | bc)

        echo "amr,$EFFECTIVE,$LEVEL,$C,$TIME" >> $CSV
    done
done


# [4] AMR parameter sensitivity
echo "Testing AMR parameter sensitivity"

THRESHOLDS=(0.02 0.01 0.005)

for T in "${THRESHOLDS[@]}"
do
    START=$(date +%s.%N)

    mpirun -np 8 $EXEC $INPUT \
        prob.type=3 \
        geometry.is_periodic="0 0" \
        amr.n_cell="$BASE_RES $BASE_RES" \
        amr.max_level=2 \
        adv.rho_grad="$T $T $T" \
        amr.regrid_int=2 \
        amr.plot_files_output=0 \
        adv.v=0 \
        amr.v=0

    END=$(date +%s.%N)

    TIME=$(echo "$END - $START" | bc)

    echo "amr_threshold,$BASE_RES,2,8,$TIME" >> $CSV
done


echo "All timing tests complete"
echo "Results saved to $CSV"