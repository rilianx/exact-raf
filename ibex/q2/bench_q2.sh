#!/bin/bash
# referencia de rendimiento de LinearizerQuad2: celdas y tiempo en instancias fijas (secuencial)
B=/home/user/jcgeraldo/ibex-lib; Q=/home/user/raa-affine-arithmetic/notas/raak-qplib/instancias_fqf; X=/home/user/ibex-local/ibexopt-raak/build
for a in ${VARIANTS:-q2 q2kn q2pc}; do
  for i in polynom/Pramanik non-polynom/Directkin; do
    r=$($X/ibexsolve-raak $B/benchs/solver/$i.bch --ctc=default-xn+affine --arith=$a -t 120 -q | grep -o "cells=[0-9]* time=[0-9.]*")
    echo "$a $(basename $i) $r"
  done
  for f in $B/benchs/optim/medium/ex6_2_8.bch $B/benchs/optim/medium/ex7_3_5bis.bch $B/raak/instancias/sinteticas/synq_n08_d50_pos_normal_s1.bch $Q/fqf_n20_d030_pos_s1.bch; do
    r=$($X/ibexopt-raak $f --lr=art --arith=$a -t 60 2>&1 | grep -E "number of cells|cpu time" | sed 's/\x1b\[[0-9;]*m//g; s/\t//g' | tr '\n' ' ')
    echo "$a $(basename $f .bch) $r"
  done
done
