#! /bin/sh

# Anne-Sophie Denommé-Pichon
# asdeno@hotmail.fr

./aso.py "./data/Early onset or syndromic epilepsy.tsv" \
	 "./data/2023-12-09_genemap2.txt" \
	 "./data/clinvar_20231126.vcf" \
	 "./data/gnomad.v4.0.constraint_metrics.tsv" \
	 "./data/gnomad.v2.1.1.lof_metrics.by_gene.txt" \
	 > ./data/panelapplist.annotated.clinvar.gnomad.tsv

