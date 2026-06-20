#!/bin/sh
set -e

gfortran -c config.f90
gfortran -c nn_cluster.f90
gfortran -c space_time_mag_count.f90
gfortran -c model_sims.f90
gfortran -c sbi_estimation.f90

gfortran -O3 config.o nn_cluster.o space_time_mag_count.o model_sims.o sbi_estimation.o -o summary_stats.exe

./summary_stats.exe

rm -f *.o
rm -f *.mod
rm -f summary_stats.exe
