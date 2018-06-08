#!/bin/bash

cat islet_meta.results.annotated.qvalue.postfreq.bestEqtl.tbl | awk '{$1=substr($1, 1, length($1)-4); print}' | tr ' ' '\t' > islet_meta.results.annotated.qvalue.postfreq.bestEqtl.mod.tbl

cat WithoutWindows_All_caQTLs.txt | cut -f1 | tail -n +2 > WithoutWindows_All_caQTLs_rsids.txt

grep -f WithoutWindows_All_caQTLs_rsids.txt islet_meta.results.annotated.qvalue.postfreq.bestEqtl.mod.tbl > islet_meta.results.annotated.qvalue.postfreq.bestEqtl.mod.leadCaQTLsnpsOnly.tbl
