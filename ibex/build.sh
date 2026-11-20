#!/usr/bin/env bash
# Reconstruye Ibex + ibex-affine (parches RAA-K 0001-0015 de JCGeraldo + 0016 de este repo) y el CLI
# ibexopt-raak / ibexsolve-raak con --arith=af2|raa|sh|exact (y ibexopt-raak con --arith=q2|q2kn|q2pc|q2full).
#   ROOT=/ruta/de/trabajo [LP_LIB=soplex*|clp] bash ibex/build.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${ROOT:-$HOME/ibex-local}"
LP_LIB="${LP_LIB:-soplex}"   # SoPlex 4.0.2 (incluido en ibex-lib, lp_lib_wrapper/soplex/3rd); clp como alternativa
JOBS="${JOBS:-4}"
PREFIX="$ROOT/prefix"
mkdir -p "$ROOT"
cd "$ROOT"

# 1. ibex-lib: fork de JCGeraldo (rama raak-experimentos), que trae el cambio 0003 y los parches RAA-K
[ -d jcgeraldo-ibex-lib ] || git clone -q -b raak-experimentos https://github.com/JCGeraldo/ibex-lib jcgeraldo-ibex-lib
cmake -S jcgeraldo-ibex-lib -B jcgeraldo-ibex-lib/build -DCMAKE_BUILD_TYPE=Release -DLP_LIB="$LP_LIB" \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" > ibexlib-cmake.log 2>&1
make -C jcgeraldo-ibex-lib/build -j"$JOBS" > ibexlib-make.log 2>&1
make -C jcgeraldo-ibex-lib/build install > ibexlib-install.log 2>&1

# 2. ibex-affine en c219009c + parche consolidado (0001-0015) + 0016 + 0017
RAAK="$ROOT/jcgeraldo-ibex-lib/raak"
if [ ! -d ibex-affine ]; then
  git clone -q https://github.com/ibex-team/ibex-affine.git ibex-affine
  git -C ibex-affine checkout -q -B raak c219009cc6425826b760939b9486776f744ded96
  git -C ibex-affine apply "$RAAK/patches/raak-ibex-affine-consolidado.patch"
  git -C ibex-affine apply "$HERE/patches/0016-ibex-affine-multiplicacion-exacta.patch"
  git -C ibex-affine apply "$HERE/patches/0017-ibex-affine-multiplicacion-rapida.patch"
fi
cmake -S ibex-affine -B ibex-affine/build -DCMAKE_BUILD_TYPE=Release -DIBEX_DIR="$PREFIX" \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" > affine-cmake.log 2>&1
make -C ibex-affine/build -j"$JOBS" > affine-make.log 2>&1
make -C ibex-affine/build install > affine-install.log 2>&1

# 3. CLI (fuentes de JCGeraldo con --arith=sh|exact, y en este repo --arith=q2|q2kn|q2pc|q2full:
#    LinearizerQuad2, formas afines de grado 2, ibex/q2/)
CLI="$ROOT/ibexopt-raak"
mkdir -p "$CLI/src/bin" "$CLI/src/q2"
cp "$HERE/cli/cli_CMakeLists.txt"               "$CLI/CMakeLists.txt"
cp "$HERE/cli/cli_Optimizer04Config_Affine.h"   "$CLI/src/ibex_Optimizer04Config_Affine.h"
cp "$HERE/cli/cli_Optimizer04Config_Affine.cpp" "$CLI/src/ibex_Optimizer04Config_Affine.cpp"
cp "$HERE"/q2/*.h                               "$CLI/src/q2/"
cp "$HERE/cli/cli_ibexopt_raak.cpp"             "$CLI/src/bin/ibexopt_raak.cpp"
cp "$HERE/cli/cli_ibexsolve_raak.cpp"           "$CLI/src/bin/ibexsolve_raak.cpp"
cp "$ROOT/jcgeraldo-ibex-lib/src/bin/parse_args.h" "$CLI/src/bin/parse_args.h"
cmake -S "$CLI" -B "$CLI/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$PREFIX" > cli-cmake.log 2>&1
make -C "$CLI/build" -j"$JOBS" > cli-make.log 2>&1

# 4. pruebas
export PKG_CONFIG_PATH="$PREFIX/share/pkgconfig:$PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
for t in affine_containment_check bench_mult; do
  g++ -std=gnu++11 -O2 -msse3 -frounding-math "$HERE/tests/$t.cpp" $(pkg-config --cflags ibex-affine ibex) \
      -o "$ROOT/$t" $(pkg-config --libs ibex-affine ibex)
done
echo "listo: $CLI/build/ibexopt-raak, $CLI/build/ibexsolve-raak, $ROOT/affine_containment_check, $ROOT/bench_mult"
