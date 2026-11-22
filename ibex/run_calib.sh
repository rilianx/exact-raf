#!/usr/bin/env bash
# Corre S, S+, SH y E en la misma sesion sobre un subconjunto, intercalados por instancia, con referencia al inicio y al final.
B=/home/user/ibex-local/ibexopt-raak/build/ibexopt-raak
OUT=/home/user/ibex-local/results/calib
REF=/home/user/ibex-local/benchs-mh/ex6_2_9.bch
one() { # tag variante archivo
  local tag="$1" v="$2" f="$3" env=""
  case "$v" in S) a=af2;; Splus) a=af2; env="IBEX_AF2_SQR=cheb";; SH) a=sh;; E) a=exact;; esac
  local res; res=$(env $env $B --lr=art --arith=$a --knomial_k=0 -t 120 "$f" 2>&1 | awk '/number of cells/{c=$NF} /cpu time used/{t=$NF} /f\* in/{st="ok"} END{print c","t}')
  echo "$tag,$(basename $f .bch),$v,$res"
}
export -f one; export B
ref() { for k in 1 2 3 4; do echo "ref$1 Splus $REF"; done | xargs -P4 -L1 bash -c 'one "$0" "$1" "$2"'; }
{
ref inicio
for rep in 1 2; do
  while read f; do for v in S Splus SH E; do echo "rep$rep $v $f"; done; done < /home/user/ibex-local/subset.txt
done | xargs -P4 -L1 bash -c 'one "$0" "$1" "$2"'
ref fin
} > $OUT/calib.csv
echo done
