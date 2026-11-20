#!/usr/bin/env bash
# Barridos del solver con las cuatro aritmeticas (af2, raa, sh, exact; k = 0), con los scripts de
# JCGeraldo (raak/scripts/wsl_sweep_dir.sh y wsl_sweep_csp.sh).
#   ROOT=... RAAK=.../raak OUTROOT=... TLIM=300 P=4 bash ibex/run_experiments.sh
set -uo pipefail
ROOT="${ROOT:-$HOME/ibex-local}"
RAAK="${RAAK:-$ROOT/jcgeraldo-ibex-lib/raak}"
OUTROOT="${OUTROOT:-$ROOT/results}"
TLIM="${TLIM:-300}"; P="${P:-4}"
OPT="$ROOT/ibexopt-raak/build/ibexopt-raak"; SOLVE="$ROOT/ibexopt-raak/build/ibexsolve-raak"
for set in reales sinteticas sinteticas-log sinteticas-fun; do
  BIN="$OPT" DIR="$RAAK/instancias/$set" OUT="$OUTROOT/opt-$set" SET="$set" \
  VARIANTS="af2_k0 raa_k0 sh_k0 exact_k0" TLIM="$TLIM" P="$P" bash "$RAAK/scripts/wsl_sweep_dir.sh"
done
# CSP: las 17 instancias de benchs/solver con productos de factores no lineales (barrido 20 de JCGeraldo)
cd "$RAAK/.." && BIN="$SOLVE" LIST="$RAAK/results/20-benchs-solver/csp-900s/instancias.txt" OUT="$OUTROOT/csp" \
  CONFIGS="affine_af2 affine_raa affine_sh affine_exact default+affine_af2 default+affine_raa default+affine_sh default+affine_exact" \
  TLIM="$TLIM" P="$P" bash "$RAAK/scripts/wsl_sweep_csp.sh"
