#!/bin/bash
B=/home/user/jcgeraldo/ibex-lib; Q=/home/user/raa-affine-arithmetic/notas/raak-qplib/instancias_fqf
F="$(ls $B/benchs/optim/medium/*.bch $B/benchs/optim/hard/*.bch | tr '\n' ' ') $B/raak/instancias/reales/3pk-borne.bch $B/raak/instancias/reales/mistake.bch $B/raak/instancias/reales/s365mod.bch $(ls $B/raak/instancias/sinteticas/synq_*.bch | tr '\n' ' ') $(ls $Q/fqf_n10_*.bch $Q/fqf_n20_*.bch | tr '\n' ' ')"
S="$(sed "s#^#$B/#" $B/raak/results/20-benchs-solver/csp-900s/instancias.txt | tr '\n' ' ')"
X=/home/user/ibex-soplex/ibexopt-raak/build; R=/home/user/raa-affine-arithmetic/ibex/results
cd /home/user/raa-affine-arithmetic/ibex/q2
python3 barrido_opt.py $X/ibexopt-raak $R/opt_q2_soplex_t30.csv 30 4 af2,sh,q2,q2kn,q2pc $F > /tmp/claude-0/log_bs1.txt 2>&1
python3 barrido_solve.py $X/ibexsolve-raak $R/csp17_q2_soplex_t120.csv 120 4 af2,sh,q2,q2kn,q2pc $S > /tmp/claude-0/log_bs2.txt 2>&1
echo FIN > /tmp/claude-0/barrido_soplex.fin
