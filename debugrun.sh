#!/bin/bash 

compiler="gfortran"
options="-Wall -Wno-unused-dummy-argument -Wunused-variable -Wno-unused-function -fimplicit-none -fcheck=all -g -pedantic -Wno-uninitialized"
source="ModVec.f90 Main.f90"
object="ModVec.o Main.o"
program="Main_gnu.exe"

mkdir -p obj_gnu;
for file in $source; do
  $compiler -J./obj_gnu -I./obj_gnu $options -c src/$file
  if [[ $? != 0 ]]; then
    exit
  fi
done
$compiler -I./obj_gnu -o $program $object -llapack -lblas
mv $object obj_gnu
./$program
