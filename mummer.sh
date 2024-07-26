nucmer --maxgap=500 --mincluster=100 --prefix=out panax_chr.fasta panax_notoginseng.fasta
show-coords -r out.delta > out.coords
delta-filter -i 95 -l 200 -r -g  out.delta > out.flt.delta
mummerplot --fat --filter --layout --png -s large out.flt.delta
