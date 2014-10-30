#!/bin/bash

# plot precision based on epsilon, varying all other params
plot-prec-eps() {
    datafile=$1
    output=$2
    label=$3
    cat << END | gnuplot
set terminal post eps enhanced font "Helvetica,28"
set output "${output}"
#set yrange [-0.1:1.1]
set xrange [0:1.1]
set xlabel "{/Symbol e}"
set ylabel "${label}"
set key center right
plot "${datafile}" using 1:2 title "filtered" w lp ps 2 lt 1 pt 4,\
     "${datafile}" using 1:3 title "noisy" w lp ps 2 lt 1 pt 5
END
}

plot-rel-eps() {
    datafile=$1
    output=$2
    label=$3
    cat << END | gnuplot
set terminal post eps enhanced font "Helvetica,28"
set output "${output}"
#set yrange [-0.1:1.1]
set xrange [0:1.1]
set xlabel "{/Symbol e}"
set ylabel "${label}"
set key center right
plot "${datafile}" using 1:6 title "filtered" w lp ps 2 lt 1 pt 4,\
     "${datafile}" using 1:7 title "noisy" w lp ps 2 lt 1 pt 5
END
}

datasets="pairs_100nodes_10lmax_10ktmax_42seed pairs_100nodes_10lmax_10ktmax_6seed pairs_100nodes_10lmax_1ktmax_42seed pairs_5nodes_10lmax_10ktmax_42seed"

for f in ${datasets}; do
    plot-prec-eps "${f}.csv" "${f}-abs.eps" "Absolute errors (L2/Frob)"
    plot-rel-eps "${f}.csv" "${f}-rel.eps" "Relative errors"
done
