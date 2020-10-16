
gatk AddOrReplaceReadGroups --CREATE_INDEX true --TMP_DIR tmp -I /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/bam/Tumor_1.bam -O 01_9544_A_T.tmp.bam \
  -ID 1 \
  -LB 01_9544_A_T \
  -PL ILLUMINA \
  -PU 01_9544_A_T \
  -SM 01_9544_A_T && touch 01_9544_A_T.rg.bam.done

if [[ -e 01_9544_A_T.rg.bam.done ]]; then
  mv 01_9544_A_T.tmp.bam 01_9544_A_T.rg.bam
  mv 01_9544_A_T.tmp.bai 01_9544_A_T.rg.bam.bai
  rm 01_9544_A_T.rg.bam.done
fi

#markduplicate

sambamba markdup -p -t 8 --tmpdir tmp /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_01_replace_read_group/result/01_9544_A_T.rg.bam 01_9544_A_T.tmp.bam && touch 01_9544_A_T.duplicates_marked.bam.done

  if [[ -e 01_9544_A_T.duplicates_marked.bam.done ]]; then
    mv 01_9544_A_T.tmp.bam 01_9544_A_T.duplicates_marked.bam
    mv 01_9544_A_T.tmp.bam.bai 01_9544_A_T.duplicates_marked.bam.bai
  fi

#BaseRecalibrator

gatk --java-options "-Xmx40G" \
  BaseRecalibrator \
  -R /scratch/jbrown_lab/references/genome_AAV1_AAV2/combined_mm_aav.fa \
  -I /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_02_markduplicates_sambamba/result/16_9548_A_L.duplicates_marked.bam \
  --use-original-qualities \
  -O 16_9548_A_L.49.recal_data.csv \
  -L /scratch/jbrown_lab/references/genome_AAV1_AAV2/combined_mm_aav_intervals/0049-scattered.interval_list \
  --known-sites /scratch/cqs_share/references/dbsnp/mouse_10090_b150_GRCm38.p4.vcf.gz --known-sites /scratch/cqs_share/references/dbsnp/mgp.v5.indels.nochr.vcf.gz

#GatherBQSRReports

gatk --java-options "-Xmx5G" \
  GatherBQSRReports \
  -I /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_03_BaseRecalibratorScatter/result/16_9548_A_L/16_9548_A_L.00.recal_data.csv \
  -I /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_03_BaseRecalibratorScatter/result/16_9548_A_L/16_9548_A_L.01.recal_data.csv \
  -I /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_03_BaseRecalibratorScatter/result/16_9548_A_L/16_9548_A_L.02.recal_data.csv \
  -O 16_9548_A_L.recal_data.csv
 
#ApplyBQSR 
 
gatk --java-options "-Xmx40G" \
  ApplyBQSR \
  --add-output-sam-program-record \
  -R /scratch/jbrown_lab/references/genome_AAV1_AAV2/combined_mm_aav.fa \
  -I /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_02_markduplicates_sambamba/result/16_9548_A_L.duplicates_marked.bam \
  --use-original-qualities \
  -O 16_9548_A_L.49.tmp.recalibrated.bam \
  -bqsr /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_04_GatherBQSRReports/result/16_9548_A_L.recal_data.csv \
  -L /scratch/jbrown_lab/references/genome_AAV1_AAV2/combined_mm_aav_intervals/0049-scattered.interval_list

if [[ -s 16_9548_A_L.49.tmp.recalibrated.bai ]]; then
  mv 16_9548_A_L.49.tmp.recalibrated.bam 16_9548_A_L.49.recalibrated.bam
  mv 16_9548_A_L.49.tmp.recalibrated.bai 16_9548_A_L.49.recalibrated.bai
fi

#merge recalibrated bam

sambamba merge  -t 8 16_9548_A_L.tmp.recalibrated.bam \
  /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_05_ApplyBQSRScatter/result/16_9548_A_L/16_9548_A_L.00.recalibrated.bam \
  /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_05_ApplyBQSRScatter/result/16_9548_A_L/16_9548_A_L.01.recalibrated.bam \
  /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_05_ApplyBQSRScatter/result/16_9548_A_L/16_9548_A_L.02.recalibrated.bam \
  /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_05_ApplyBQSRScatter/result/16_9548_A_L/16_9548_A_L.03.recalibrated.bam \
  /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_05_ApplyBQSRScatter/result/16_9548_A_L/16_9548_A_L.04.recalibrated.bam

if [[ -s 16_9548_A_L.tmp.recalibrated.bam.bai ]]; then
  mv 16_9548_A_L.tmp.recalibrated.bam 16_9548_A_L.recalibrated.bam
  mv 16_9548_A_L.tmp.recalibrated.bam.bai 16_9548_A_L.recalibrated.bam.bai
fi

#HaplotypeCaller

gatk --java-options "-Xmx40G" \
  HaplotypeCaller \
    -contamination 0 \
    -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
    -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 -ERC GVCF \
  -L /scratch/jbrown_lab/references/genome_AAV1_AAV2/combined_mm_aav_intervals/0049-scattered.interval_list -XL /scratch/cqs_share/references/blacklist_files/mm10-blacklist.v2.nochr.bed \
  --native-pair-hmm-threads 8 \
  -R /scratch/jbrown_lab/references/genome_AAV1_AAV2/combined_mm_aav.fa \
  -I /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_06_GatherSortedBamFiles/result/16_9548_A_L.recalibrated.bam \
  -O 16_9548_A_L.49.tmp.g.vcf.gz

if [[ -s 16_9548_A_L.49.tmp.g.vcf.gz.tbi ]]; then
  mv 16_9548_A_L.49.tmp.g.vcf.gz 16_9548_A_L.49.g.vcf.gz
  mv 16_9548_A_L.49.tmp.g.vcf.gz.tbi 16_9548_A_L.49.g.vcf.gz.tbi
fi

#MergeVcfs

java -Xms10G -jar /usr/gitc/picard.jar \
      MergeVcfs \
      INPUT=/scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_07_HaplotypeCallerScatter/result/16_9548_A_L/16_9548_A_L.00.g.vcf.gz \
      INPUT=/scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_07_HaplotypeCallerScatter/result/16_9548_A_L/16_9548_A_L.01.g.vcf.gz \
      OUTPUT=16_9548_A_L.tmp.g.vcf.gz

if [[ -s 16_9548_A_L.tmp.g.vcf.gz.tbi ]]; then
  mv 16_9548_A_L.tmp.g.vcf.gz 16_9548_A_L.g.vcf.gz
  mv 16_9548_A_L.tmp.g.vcf.gz.tbi 16_9548_A_L.g.vcf.gz.tbi
fi

#GenomicsDBImport

gatk --java-options -Xms8g \
    GenomicsDBImport \
    --genomicsdb-workspace-path lmna_paper.49 \
    --batch-size 50 \
    -L /scratch/jbrown_lab/references/genome_AAV1_AAV2/combined_mm_aav_intervals/0049-scattered.interval_list \
    --sample-name-map /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_09_GenomicsDBImportScatter/result/lmna_paper_gbi.sample_name_map.txt \
    --reader-threads 5 \
    --merge-input-intervals \
    --consolidate

#GenotypeGVCFs

gatk --java-options -Xms8g \
  GenotypeGVCFs \
  -R /scratch/jbrown_lab/references/genome_AAV1_AAV2/combined_mm_aav.fa \
  -O lmna_paper.49.tmp.g.vcf.gz \
  -D /scratch/cqs_share/references/dbsnp/mouse_10090_b150_GRCm38.p4.vcf.gz \
  -G StandardAnnotation -G AS_StandardAnnotation \
  --only-output-calls-starting-in-intervals \
  -V gendb:///scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_09_GenomicsDBImportScatter/result/lmna_paper.49 \
  -L /scratch/jbrown_lab/references/genome_AAV1_AAV2/combined_mm_aav_intervals/0049-scattered.interval_list \
  --merge-input-intervals

if [[ -s lmna_paper.49.tmp.g.vcf.gz.tbi ]]; then
  mv lmna_paper.49.tmp.g.vcf.gz lmna_paper.49.g.vcf.gz
  mv lmna_paper.49.tmp.g.vcf.gz.tbi lmna_paper.49.g.vcf.gz.tbi
fi

#filter SNP/indel

gatk SelectVariants \
  -V /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_10_GenotypeGVCFs/result/lmna_paper.49.g.vcf.gz \
  -select-type SNP \
  -O snps.vcf.gz

gatk SelectVariants \
  -V /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_10_GenotypeGVCFs/result/lmna_paper.49.g.vcf.gz \
  -select-type INDEL \
  -select-type MIXED \
  -O indels.vcf.gz

gatk VariantFiltration \
  -V snps.vcf.gz \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"   \
  -O snps_filtered.vcf.gz

gatk VariantFiltration \
  -V indels.vcf.gz \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" \
  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"   \
  -O indels_filtered.vcf.gz

gatk MergeVcfs \
  -I snps_filtered.vcf.gz \
  -I indels_filtered.vcf.gz \
  -O merged_filtered.vcf.gz

gatk SelectVariants \
  -O lmna_paper.49.indels.snp.hardfilter.pass.vcf.gz \
  -V merged_filtered.vcf.gz \
  --exclude-filtered

#LeftAlignAndNormalize

if [[ -s /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_11_VariantFilterHard/result/lmna_paper.49/lmna_paper.49.indels.snp.hardfilter.pass.vcf.gz && ! -s lmna_paper.49.norm.vcf ]]; then
  echo LeftAlignAndNorm=`date`
  bcftools norm -m- -o lmna_paper.49.split.vcf /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_11_VariantFilterHard/result/lmna_paper.49/lmna_paper.49.indels.snp.hardfilter.pass.vcf.gz
  bcftools norm -f /scratch/jbrown_lab/references/genome_AAV1_AAV2/combined_mm_aav.fa -o lmna_paper.49.norm.vcf lmna_paper.49.split.vcf
fi

if [[ -s lmna_paper.49.norm.vcf && ! -s lmna_paper.49.indels.snp.hardfilter.pass.norm.nospan.vcf.gz ]]; then
  echo noSpanDeletion=`date`
  python /scratch/cqs_share/softwares/ngsperl/lib/GATK4/fixLeftTrimDeletion.py -i lmna_paper.49.norm.vcf -o lmna_paper.49.indels.snp.hardfilter.pass.norm.nospan.vcf
  bgzip lmna_paper.49.indels.snp.hardfilter.pass.norm.nospan.vcf
  tabix -p vcf lmna_paper.49.indels.snp.hardfilter.pass.norm.nospan.vcf.gz
fi

  if [[ -s lmna_paper.49.indels.snp.hardfilter.pass.norm.nospan.vcf.gz ]]; then
    rm lmna_paper.49.split.vcf lmna_paper.49.norm.vcf
  fi

#Merge normalized vcfs 

gatk --java-options "-Xms10G -Xmx10G" \
      MergeVcfs \
      -I lmna_paper_mv.input.list \
      -O lmna_paper.tmp.vcf.gz

if [[ -s lmna_paper.tmp.vcf.gz.tbi ]]; then
  mv lmna_paper.tmp.vcf.gz lmna_paper.vcf.gz
  mv lmna_paper.tmp.vcf.gz.tbi lmna_paper.vcf.gz.tbi
fi

#Annotation by ANNOVAR

if [[ ! -s lmna_paper.vcf.gz.annovar.mm10_multianno.txt && ! -s lmna_paper.vcf.gz.annovar.final.tsv ]]; then
  convert2annovar.pl -format vcf4old /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_13_MergeVcfs/result/lmna_paper.vcf.gz | cut -f1-7 | awk '{gsub(",\\*", "", $0); print}'> lmna_paper.vcf.gz.avinput
  if [ -s lmna_paper.vcf.gz.avinput ]; then
    table_annovar.pl lmna_paper.vcf.gz.avinput /scratch/cqs_share/references/annovar/mousedb/ -buildver mm10 -protocol refGene -operation g --remove --outfile lmna_paper.vcf.gz.annovar --remove
  fi
fi

echo find_protein_position_for_splicing=`date`
if [[ -s lmna_paper.vcf.gz.annovar.mm10_multianno.txt && ! -s lmna_paper.vcf.gz.annovar.splicing.mm10_multianno.txt ]]; then
  python /scratch/cqs_share/softwares/ngsperl/lib/Annotation/annovarSplicing.py -i lmna_paper.vcf.gz.annovar.mm10_multianno.txt -d /scratch/cqs_share/references/annovar/mousedb/ -o lmna_paper.vcf.gz.annovar.splicing.mm10_multianno.txt -b mm10
fi

if [[ -s lmna_paper.vcf.gz.annovar.splicing.mm10_multianno.txt && ! -s lmna_paper.vcf.gz.annovar.final.tsv ]]; then
  sed -n '/^[^#]/q;p' /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_13_MergeVcfs/result/lmna_paper.vcf.gz|sed '$ d' > lmna_paper.vcf.gz.annovar.final.tsv.header
  zcat /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_13_MergeVcfs/result/lmna_paper.vcf.gz | grep -v "^##" | cut -f7- > /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_13_MergeVcfs/result/lmna_paper.vcf.gz.clean
  grep -v "^##" lmna_paper.vcf.gz.annovar.splicing.mm10_multianno.txt > lmna_paper.vcf.gz.annovar.splicing.mm10_multianno.txt.clean
  paste lmna_paper.vcf.gz.annovar.splicing.mm10_multianno.txt.clean /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_13_MergeVcfs/result/lmna_paper.vcf.gz.clean > lmna_paper.vcf.gz.annovar.final.tsv.data
  cat lmna_paper.vcf.gz.annovar.final.tsv.header lmna_paper.vcf.gz.annovar.final.tsv.data > lmna_paper.vcf.gz.annovar.final.tsv
  rm /scratch/jbrown_lab/shengq2/projects/20200919_paper_result_snv/nih_bwa/gatk4_13_MergeVcfs/result/lmna_paper.vcf.gz.clean lmna_paper.vcf.gz.annovar.splicing.mm10_multianno.txt.clean lmna_paper.vcf.gz.annovar.final.tsv.header lmna_paper.vcf.gz.annovar.final.tsv.data
fi
