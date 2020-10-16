# Improperly paired reads that cover aav contigs

samtools view -h ${sample} bam aav_c_contig aav_n_contig| perl -lane 'if(/^@/ || $F[6] ne "=" ) {print }' | samtools view -q 30 -F1792 -S -bh -o $sample.aav_contig.bam 
 
# Extraction of mouse target reads:

samtools view -F256 ${sample}.aav_contig.bam | perl -lane 'print join("\t", @F[2, 3,6,7]) if $F[6] !~ /aav/' | sort | uniq -c | awk -v s=${sample} -v OFS="\t" '{print "chr"$4, $5, $5+150, s"_"$2"_"$3} ' | sort -k1,1 -k3,3n > ${sample}.target.bed
 
# Overlap with mouse genes (padded 2kb to the gene definition at two ends by this point)

bedtools intersect -wao -a ${sample}.target.bed -b gencode.vM24.all_genes_2kb_flank.bed.gz > target_${sample}.bed

awk -v s=${sample} '$6>0 {print s, $0}' target_${sample}.bed; done | perl -lane '$F[0] =~ s/target_improper_pair_//; print join("\t", @F); print join("\t", qw(Sample HitChr HitStart HitEnd AAV9 GeneChr GeneStart GeneEnd EnsmblID Symbol Strand Type Overlap))' > aav9_targets.xls
