#!/bin/sh
set -e

gfortran -c config.f90
gfortran -c si_declustering.f90
gfortran -c nn_cluster.f90
gfortran -c space_time_mag_count.f90
gfortran -c preprocess_step.f90

gfortran -O3 config.o si_declustering.o nn_cluster.o space_time_mag_count.o preprocess_step.o -o summary_stats.exe

./summary_stats.exe

rm -f *.o
rm -f *.mod
rm -f summary_stats.exe
