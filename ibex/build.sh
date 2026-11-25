#!/usr/bin/env bash
# Reconstruye Ibex + ibex-affine (parches RAA-K 0001-0015 de JCGeraldo + 0016 de este repo) y el CLI
# ibexopt-raak / ibexsolve-raak con --arith=af2|raa|sh|exact.
#   ROOT=/ruta/de/trabajo LP_LIB=clp|soplex bash ibex/build.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${ROOT:-$HOME/ibex-local}"
LP_LIB="${LP_LIB:-clp}"
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
  git -C ibex-affine apply "$HERE/patches/0018-ibex-affine-cuadrado-chebyshev-af2.patch"
  git -C ibex-affine apply "$HERE/patches/0019-ibex-affine-dos-esquinas-una-pasada.patch"
  git -C ibex-affine apply "$HERE/patches/0020-ibex-affine-cuadrado-exacto.patch"
  git -C ibex-affine apply "$HERE/patches/0021-ibex-affine-diagnostico-s-mas-ancho.patch"
fi
cmake -S ibex-affine -B ibex-affine/build -DCMAKE_BUILD_TYPE=Release -DIBEX_DIR="$PREFIX" \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" > affine-cmake.log 2>&1
make -C ibex-affine/build -j"$JOBS" > affine-make.log 2>&1
make -C ibex-affine/build install > affine-install.log 2>&1

# 3. CLI (fuentes de JCGeraldo con --arith=sh|exact)
CLI="$ROOT/ibexopt-raak"
mkdir -p "$CLI/src/bin"
cp "$RAAK/cli/cli_CMakeLists.txt"               "$CLI/CMakeLists.txt"
cp "$RAAK/cli/cli_Optimizer04Config_Affine.h"   "$CLI/src/ibex_Optimizer04Config_Affine.h"
cp "$RAAK/cli/cli_Optimizer04Config_Affine.cpp" "$CLI/src/ibex_Optimizer04Config_Affine.cpp"
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
