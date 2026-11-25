#!/usr/bin/env bash
# Tabla 5 con las variantes intercaladas por instancia en la misma sesion (P corridas en paralelo), reanudable.
#   TASKS=<archivo con lineas "set ruta.bch tlim variante"> OUT=<dir> P=4 bash run_interleaved.sh
set -uo pipefail
B="${B:-${ROOT:-$HOME/ibex-local}/ibexopt-raak/build/ibexopt-raak}"
: "${TASKS:?}"; : "${OUT:?}"; P="${P:-4}"
mkdir -p "$OUT/rows" "$OUT/logs"
export B OUT
one() {
  local set="$1" f="$2" tl="$3" v="$4"; local inst; inst=$(basename "$f" .bch)
  local row="$OUT/rows/${set}__${inst}__${v}.csv"; [ -s "$row" ] && return 0
  local env="" a=""
  case "$v" in af2) a=af2;; splus) a=af2; env="IBEX_AF2_SQR=cheb";; sh) a=sh;; exact) a=exact;; exactsq) a=exact; env="IBEX_SQR_EXACT=1";; shsq) a=sh; env="IBEX_SQR_EXACT=1";; esac
  local log="$OUT/logs/${set}__${inst}__${v}.log" t0; t0=$(date +%s.%N)
  env $env timeout $((tl * 2 + 120)) $B --lr=art --arith=$a --knomial_k=0 -t "$tl" "$f" > "$log" 2>&1; local code=$?
  local wall; wall=$(awk -v a="$t0" -v b="$(date +%s.%N)" 'BEGIN{printf "%.2f", b-a}')
  local clean; clean="$(sed -r 's/\x1b\[[0-9;]*m//g' "$log")"
  local st
  if grep -q "optimization successful" <<<"$clean"; then st=SUCCESS
  elif grep -q "time limit" <<<"$clean"; then st=TIME_OUT
  elif grep -q "no feasible point found" <<<"$clean"; then st=NO_FEASIBLE_FOUND
  elif grep -q "infeasible problem" <<<"$clean"; then st=INFEASIBLE
  elif grep -q "unreached precision" <<<"$clean"; then st=UNREACHED_PREC
  elif [ "$code" -eq 124 ]; then st=KILLED; else st=ERROR; fi
  local flb fub cells cpu
  flb="$(awk -F'[][,]' '/f\* in/ {print $2; exit}' <<<"$clean")"; fub="$(awk -F'[][,]' '/f\* in/ {print $3; exit}' <<<"$clean")"
  cells="$(awk -F'\t+' '/number of cells/ {print $2; exit}' <<<"$clean")"
  cpu="$(awk -F'\t+' '/cpu time used/ {v=$2; sub(/s.*$/,"",v); print v; exit}' <<<"$clean")"
  echo "$set/$inst,art_${v}_k0,$st,$cells,$flb,$fub,$cpu,$wall,$code,," > "$row"
  echo "[$st] $set/$inst $v celdas=$cells cpu=$cpu"
}
export -f one
REF=/home/user/ibex-local/benchs-mh/ex6_2_9.bch
ref() { for k in 1 2 3 4; do echo "ref_$1 $REF 60 splus"; done | xargs -P4 -L1 bash -c 'rm -f "$OUT/rows/$0__ex6_2_9__splus.csv"; one "$0" "$1" "$2" "$3"'; }
ref inicio
xargs -P"$P" -L1 bash -c 'one "$0" "$1" "$2" "$3"' < "$TASKS"
ref fin
echo done
