# Do variant calling.
betsy_run.py --network_png test01.png --network_json test01.txt \
  --input FastqFolder \
  --input ReferenceGenome \
  --input SampleGroupFile \
  --output VCFFolder \
  --dattr VCFFolder.caller=gatk \
  --dattr VCFFolder.aligner=bowtie2 \
  --dattr VCFFolder.vartype=all >& test01.log
# Do Classification.
betsy_run.py --network_png test02.png --network_json test02.txt \
  --input UnprocessedSignalFile \
  --dattr UnprocessedSignalFile.contents=class0,class1 \
  --input ClassLabelFile \
  --dattr ClassLabelFile.contents=class0,class1 \
  --output ClassifyFile \
  --dattr ClassifyFile.classify_alg=random_forest \
  --dattr ClassifyFile.loocv=yes \
  --dattr SignalFile.gene_order=ttest_p \
  --dattr SignalFile.preprocess=unknown \
  --dattr SignalFile.num_features=yes >& test02.log
# Run FastQC.
betsy_run.py --network_png test03.png --network_json test03.txt \
  --input FastqFolder \
  --input SampleGroupFile \
  --output FastQCSummary >& test03.log
# Run RSeQC.
betsy_run.py --network_png test04.png --network_json test04.txt \
  --input FastqFolder \
  --input SampleGroupFile \
  --input ReferenceGenome \
  --input STARReferenceGenome \
  --input GTFGeneModel \
  --output RSeQCResults \
  --dattr BamFolder.aligner=star \
  --mattr gene_model=hg19 >& test04.log
# Subtract mouse reads.
betsy_run.py --network_png test05.png --network_json test05.txt \
  --input FastqFolder \
  --input SampleGroupFile \
  --input ReferenceGenome \
  --input STARReferenceGenome \
  --input GTFGeneModel \
  --output PerfectAlignmentSummary \
  --dattr PerfectAlignmentSummary.adapters_trimmed=yes \
  --dattr BamFolder.aligner=star >& test05.log
# See how many reads aligned perfectly in human genome.
# TODO: There's a shorter pipeline possible here that was not chosen
# by the pruner.
betsy_run.py --network_png test06.png --network_json test06.txt \
  --input FastqFolder \
  --input SampleGroupFile \
  --input ReferenceGenome \
  --input STARReferenceGenome \
  --input GTFGeneModel \
  --output PerfectAlignmentSummary \
  --dattr PerfectAlignmentSummary.adapters_trimmed=yes \
  --dattr BamFolder.aligner=star >& test06.log
# Calculate the RSEM values.
betsy_run.py --network_png test07.png --network_json test07.txt \
  --input FastqFolder \
  --dattr FastqFolder.adapters_trimmed=yes \
  --dattr FastqFolder.reads_merged=yes \
  --input SampleGroupFile \
  --input ReferenceGenome \
  --input RSEMReferenceGenome \
  --input GTFGeneModel \
  --output UnprocessedSignalFile \
  --dattr UnprocessedSignalFile.logged=no \
  --dattr UnprocessedSignalFile.preprocess=tpm >& test07.log
# Summarize the reads mapped for HTSeq-Count.
betsy_run.py --network_png test08.png --network_json test08.txt \
  --input FastqFolder \
  --dattr FastqFolder.reads_merged=yes \
  --input SampleGroupFile \
  --input STARReferenceGenome \
  --input GTFGeneModel \
  --output HTSeqCountSummary \
  --dattr BamFolder.aligner=star \
  --dattr BamFolder.adapters_trimmed=yes \
  --dattr BamFolder.mouse_reads_subtracted=no >& test08.log
# Run MACS21.
betsy_run.py --network_png test09.png --network_json test09.txt \
  --input BamFolder \
  --dattr BamFolder.aligner=bwa_mem \
  --dattr BamFolder.sorted=coordinate \
  --input SampleGroupFile \
  --output MACS21Results >& test09.log
# Preprocess BAM folder.
betsy_run.py --network_png test10.png --network_json test10.txt \
  --input BamFolder \
  --dattr BamFolder.aligner=bwa_mem \
  --input ReferenceGenome \
  --output BamFolder \
  --dattr BamFolder.sorted=coordinate \
  --dattr BamFolder.has_read_groups=yes \
  --dattr BamFolder.duplicates_marked=yes \
  --dattr BamFolder.indel_realigned=yes \
  --dattr BamFolder.base_quality_recalibrated=yes \
  --dattr BamFolder.indexed=yes \
  --dattr BamFolder.aligner=bwa_mem >& test10.log
# Check coverage of each of the alignments.
betsy_run.py --network_png test11.png --network_json test11.txt \
  --input BamFolder \
  --dattr BamFolder.sorted=coordinate \
  --input ReferenceGenome \
  --output DepthOfCoverage >& test11.log
# Variant calling with SomaticSniper
betsy_run.py --network_png test12.png --network_json test12.txt \
  --input BamFolder \
  --dattr BamFolder.duplicates_marked=yes \
  --dattr BamFolder.has_read_groups=yes \
  --dattr BamFolder.indel_realigned=yes \
  --dattr BamFolder.base_quality_recalibrated=yes \
  --dattr BamFolder.indexed=yes \
  --dattr BamFolder.sorted=coordinate \
  --dattr BamFolder.aligner=bwa_mem \
  --input ReferenceGenome \
  --input NormalCancerFile \
  --output VCFFolder \
  --dattr VCFFolder.caller=somaticsniper \
  --dattr VCFFolder.vartype=snp \
  --dattr VCFFolder.aligner=bwa_mem \
  --dattr VCFFolder.somatic=yes \
  >& test12.log
# Index for STAR.
betsy_run.py --network_png test13.png --network_json test13.txt \
  --input ReferenceGenome \
  --input GTFGeneModel \
  --output STARReferenceGenome >& test13.log
# Merge VCFFolder.
betsy_run.py --network_png test14.png --network_json test14.txt \
  --input ManyCallerVCFFolders \
  --input UnprocessedSignalFile \
  --dattr UnprocessedSignalFile.preprocess=tpm \
  --input BamFolder \
  --dattr BamFolder.aligner=bwa_mem \
  --input BamFolder \
  --dattr BamFolder.aligner=star \
  --input ReferenceGenome \
  --output SimpleVariantMatrix \
  --dattr SimpleVariantMatrix.vartype=snp \
  --dattr SimpleVariantMatrix.filtered_calls=yes \
  --dattr SimpleVariantMatrix.annotated=yes \
  --dattr SimpleVariantMatrix.with_coverage=yes \
  --dattr SimpleVariantMatrix.with_rna_coverage=yes \
  --dattr SimpleVariantMatrix.with_gxp=yes \
  --dattr SimpleVariantMatrix.with_cancer_genes=yes \
  --mattr filter_by_min_total_reads=20 >& test14.log
# Call germline SNPs.
betsy_run.py --network_png test15.png --network_json test15.txt \
  --input BamFolder \
  --dattr BamFolder.has_read_groups=yes \
  --dattr BamFolder.duplicates_marked=yes \
  --dattr BamFolder.indel_realigned=yes \
  --dattr BamFolder.base_quality_recalibrated=yes \
  --dattr BamFolder.sorted=coordinate \
  --dattr BamFolder.indexed=yes \
  --dattr BamFolder.aligner=bwa_mem \
  --input ReferenceGenome \
  --output SimpleVariantMatrix \
  --dattr SimpleVariantMatrix.vartype=snp \
  --dattr VCFFolder.aligner=bwa_mem \
  --mattr wgs_or_wes=wes >& test15.log
# RNA-Seq RSEM.
betsy_run.py --network_png test16.png --network_json test16.txt \
  --input FastqFolder \
  --input SampleGroupFile \
  --input STARReferenceGenome \
  --input RSEMReferenceGenome \
  --input GTFGeneModel \
  --output UnprocessedSignalFile \
  --dattr UnprocessedSignalFile.logged=no \
  --dattr UnprocessedSignalFile.preprocess=tpm \
  --dattr RNASeqUnprocessedSignalFile.adapters_trimmed=yes \
  --mattr adapters_fasta=adapters/TruSeq3-PE-2.fa >& test16.log
# RNA-Seq HTSeqCount.
betsy_run.py --network_png test17.png --network_json test17.txt \
  --input FastqFolder \
  --input SampleGroupFile \
  --input STARReferenceGenome \
  --input GTFGeneModel \
  --output UnprocessedSignalFile \
  --dattr UnprocessedSignalFile.logged=no \
  --dattr UnprocessedSignalFile.preprocess=counts \
  --dattr BamFolder.aligner=star \
  --dattr BamFolder.adapters_trimmed=yes \
  --mattr adapters_fasta=adapters/TruSeq3-PE-2.fa >& test17.log
# Call somatic SNPs, no Radia caller.
betsy_run.py --network_png test18.png --network_json test18.txt \
  --input BamFolder \
  --dattr BamFolder.has_read_groups=yes \
  --dattr BamFolder.duplicates_marked=yes \
  --dattr BamFolder.indel_realigned=yes \
  --dattr BamFolder.base_quality_recalibrated=yes \
  --dattr BamFolder.sorted=coordinate \
  --dattr BamFolder.indexed=yes \
  --dattr BamFolder.aligner=bwa_mem \
  --input ReferenceGenome \
  --input NormalCancerFile \
  --output SimpleVariantMatrix \
  --dattr SimpleVariantMatrix.vartype=snp \
  --dattr VCFFolder.aligner=bwa_mem \
  --dattr ManyCallerVCFFolders.somatic=yes \
  --mattr wgs_or_wes=wes >& test18.log
# Variant calling with Thunderbolts bam file.
# duplicates_marked=no
# base_quality_recalibrated=no
betsy_run.py --network_png test19.png --network_json test19.txt \
  --input BamFolder \
  --dattr BamFolder.sorted=coordinate \
  --dattr BamFolder.has_read_groups=yes \
  --dattr BamFolder.duplicates_marked=no \
  --dattr BamFolder.indel_realigned=yes \
  --dattr BamFolder.base_quality_recalibrated=no \
  --dattr BamFolder.indexed=yes \
  --dattr BamFolder.aligner=bwa_mem \
  --input ReferenceGenome \
  --input NormalCancerFile \
  --output SimpleVariantMatrix \
  --dattr ManyCallerVCFFolders.somatic=yes \
  --dattr BamFolder.duplicates_marked=no \
  --dattr BamFolder.base_quality_recalibrated=no \
  --dattr VCFFolder*.aligner=bwa_mem \
  --dattr SimpleVariantFile.filtered=yes \
  --dattr SimpleVariantMatrix.vartype=snp \
  --dattr SimpleVariantMatrix.annotated=yes \
  --dattr SimpleVariantMatrix.with_cancer_genes=yes \
  --dattr SimpleVariantMatrix.with_coverage=yes \
  --dattr SimpleVariantMatrix.filtered_calls=yes \
  --dattr SimpleVariantMatrix.filtered_variants=yes >& test19.log
# Call somatic SNPs, with Radia caller, add everything.
betsy_run.py --network_png test20.png --network_json test20.txt \
  --input BamFolder \
  --dattr BamFolder.aligner=bwa_mem \
  --input BamFolder \
  --dattr BamFolder.aligner=star \
  --input ReferenceGenome \
  --input NormalCancerFile \
  --input UnprocessedSignalFile \
  --dattr UnprocessedSignalFile.preprocess=tpm \
  --output SimpleVariantMatrix \
  --dattr ManyCallerVCFFolders.somatic=yes \
  --dattr ManyCallerVCFFolders.with_rna_callers=yes \
  --dattr SimpleVariantMatrix.vartype=snp \
  --dattr VCFFolder*.aligner=bwa_mem \
  --dattr VCFFolder*.aligner=star \
  --mattr mutect_dbsnp_vcf=MuTect/dbsnp_132_b37.leftAligned.vcf \
  --mattr mutect_cosmic_vcf=MuTect/b37_cosmic_v54_120711.vcf \
  --mattr muse_dbsnp_vcf=MuSE/dbsnp_132_b37.leftAligned.vcf.gz \
  --mattr radia_genome_assembly=hg19 \
  --mattr snp_eff_genome=GRCh37.75 \
  --dattr SimpleVariantMatrix.annotated=yes \
  --mattr annovar_buildver=hg19 \
  --dattr SimpleVariantFile.filtered=yes \
  --mattr wgs_or_wes=wes \
  --mattr remove_samples=PIM001_G \
  --mattr remove_radia_rna_samples=yes \
  --dattr SimpleVariantMatrix.filtered_calls=yes \
  --mattr filter_by_min_total_reads=20 \
  --dattr SimpleVariantMatrix.filtered_variants=yes \
  --mattr nonsynonymous_and_stopgain_only=yes \
  --dattr SimpleVariantMatrix.with_coverage=yes \
  --dattr SimpleVariantMatrix.with_cancer_genes=yes \
  --mattr cancer_genes_file="project/08 Cancer Genes/cancer_genes.txt" \
  --dattr SimpleVariantMatrix.with_gxp=yes \
  --dattr SimpleVariantMatrix.with_rna_coverage=yes >& test20.log
# Call somatic indels.
betsy_run.py --network_png test21.png --network_json test21.txt \
  --input BamFolder \
  --dattr BamFolder.aligner=bwa_mem \
  --input NormalCancerFile \
  --input ReferenceGenome \
  --input UnprocessedSignalFile \
  --dattr UnprocessedSignalFile.preprocess=tpm \
  --output SimpleVariantMatrix \
  --dattr SimpleVariantMatrix.vartype=indel \
  --dattr ManyCallerVCFFolders.somatic=yes \
  --dattr VCFFolder*.aligner=bwa_mem \
  --mattr wgs_or_wes=wes \
  --dattr SimpleVariantFile.filtered=yes \
  --mattr remove_samples=PIM001_G \
  --dattr SimpleVariantMatrix.filtered_calls=yes \
  --mattr filter_by_min_total_reads=20 \
  --dattr SimpleVariantMatrix.annotated=yes \
  --mattr annovar_buildver=hg19 \
  --dattr SimpleVariantMatrix.with_coverage=yes \
  --dattr SimpleVariantMatrix.with_cancer_genes=yes \
  --mattr cancer_genes_file="project/08 Cancer Genes/cancer_genes.txt" \
  --dattr SimpleVariantMatrix.with_gxp=yes >& test21.log
# RNA-Seq analysis.
betsy_run.py --network_png test22.png --network_json test22.txt \
  --input FastqFolder \
  --input SampleGroupFile \
  --input ReferenceGenome \
  --input RSEMReferenceGenome \
  --input STARReferenceGenome \
  --input GTFGeneModel \
  --output CompleteRNASeqAnalysis \
  --mattr gene_model=hg19 \
  --mattr adapters_fasta=adapters/TruSeq3-PE-2.fa >& test22.log
# Generate a heatmap.
betsy_run.py --network_png test23.png --network_json test23.txt \
  --input GEOSeries \
  --output Heatmap \
  --mattr GSEID=GSE14934 >& test23.log
# Preprocess Illumina microarray data.
betsy_run.py --network_png test24.png --network_json test24.txt \
  --input ExpressionFiles \
  --output CompleteExpressionPreprocessing \
  --dattr CompleteExpressionPreprocessing.preprocess=illumina >& test24.log
