#!/bin/bash

# For last test maybe???
for np in 1 2 4 8 16
do
  mpirun -np $np ./main2d.gnu.MPI.ex inputs_laxliu \
      amr.plot_files_output=0 \
      adv.v=0 amr.v=0
done