#!/bin/bash
for sp in arabidopsis_thaliana drosophila_melanogaster zymoseptoria_tritici septoria_linicola; do /beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python harmonised_any.py $sp; done
