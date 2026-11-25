# Lab 4: Short-Read Mapping, Variant Calling, and Consensus Generation

In the last lab, we quality-checked and trimmed our reads. Now we will map them to a reference genome, mark duplicates, call variants, filter them, and build a consensus sequence.

> **Working with your own data?** Keep following along with the reads you trimmed in Part 1. Wherever you see `SRR8883900` or `Heliconius_telesiphe`, substitute your own accession and sample name. Use your closest available reference wherever `../HelSar1/HelSar1.2_genomic.fna` appears. **Answer every question for your own dataset**. **What reference did you choose, what are the genome stats (i.e., N50) for it, and how far diverged is it from your short read samples?**

## Before you start

Work from your `short_read_assembly` folder with the conda environment active. These variables are used throughout this part; include them at the top of each job file:

```bash
cd short_read_assembly
source /home/ddebaun/mendel-nas1/miniconda3/bin/activate
conda activate /home/ddebaun/mendel-nas1/miniconda3/envs/shortreadassembly

name="Heliconius_telesiphe"
REF="../HelSar1/HelSar1.2_genomic"   # reference path without the .fna extension
```

If you don't have the reference genome in one of your folders, copy the inputs over from ```/home/ddebaun/mendel-nas1/Comparative_Genomics_Course/HelSar1```

## Mapping to a reference

BWA-MEM will map our cleaned reads to the reference genome, *Heliconius sara*. First, we index the reference so BWA can search it efficiently:

```bash
bwa index ${REF}.fna
```

We also need a read group identifier for our sequences. Copy this script to your folder:

```bash
cp /home/ddebaun/mendel-nas1/Comparative_Genomics_Course/short_read_assembly/getrg.py .
```

The script reads the read group from the FASTQ file. Store its output in the variable `RG`, which we will pass to `bwa mem`:

```bash
RG=$(python3 getrg.py "SRR8883900_1_0.10.fastq" "Heliconius_telesiphe")
echo $RG
```
**Q1. What is a read group?**

Create a job file that runs `bwa mem` **on the trimmed paired reads** with 30 threads, 10 GB of memory, and 3 hours. Use these flags:

- `-M` marks the shorter part of split hits as secondary.
- `-R "$RG"` sets the read group.

Instead of writing the SAM file with `-o` or with output redirection (`>`), we will pipe (`|`) the output to another command. The `-` marks where the input file would normally go in that command:

```bash
bwa mem [your command] | samtools sort | samtools view -bh - > "assembly/${name}.bam"
```

*Expected run time: ~14 min. For those of you working with real data, this will take much longer since we did not subsample the reads down, continue on since the next step sets up a dependency*

Piping BWA's SAM output to `samtools sort` converts it to BAM and sorts it by coordinate, which is required for downstream steps such as duplicate marking, merging, and variant calling.

## Marking duplicates

Next we mark duplicate reads.

**Q2. What is the purpose of this step?**

First we make two indices that are used for duplicate marking and later for variant calling: a sequence dictionary and a FASTA index.

```bash
picard -Xmx40g CreateSequenceDictionary R=${REF}.fna O=${REF}.dict
samtools faidx ${REF}.fna
```

`samtools faidx` produces a FASTA index file (`.fai`).

**Q3. What information is in this file?**

Make a job file that runs the following command with 10 GB of memory and 1 CPU per task:

```bash
picard MarkDuplicates \
    -I assembly/${name}.bam \
    -M assembly/markeddups.txt \
    -O assembly/${name}-mkdp.bam
```

*Expected run time: ~5 min*

This job should run only after BWA finishes. Instead of waiting for the completion email or checking the job status over and over, we can submit the job with a dependency, which tells the scheduler "run this job if and when something else happens." Here, we tell it to run after the BWA job completes with no errors (`afterok`):

```bash
sbatch --dependency=afterok:[bwajobid] picard.sh
```

Now run `squeue` to get details about the job.

**Q4. What is the REASON the job's status is pending (PD)?**

## Coverage and depth

When Picard finishes, we can check coverage and depth across the genome. This tells us about the quality of our data and helps us choose thresholds when we filter variants later.

Run these commands on a compute node:

```bash
samtools coverage assembly/${name}-mkdp.bam -m > cov
samtools depth assembly/${name}-mkdp.bam | awk '{sum+=$3} END {print sum/NR}' > assembly/${name}.depth
```

**Q5. What is the average depth across the genome?**

**Q6. What percent of the sex chromosome is covered by reads?**

> **Own data:** Use your reference's annotation or NCBI record to find which sequence is the sex chromosome. If your reference has no identified sex chromosome, say so and report coverage for one chromosome of your choice instead.

## Variant calling

Variant calling compares reads to a reference to identify sites where the sample differs from that reference (or doesn't). The standard best-practice pipeline uses GATK, specifically a tool called HaplotypeCaller, which is very computationally intensive. To run it in parallel, we can split the genome by scaffold and run the program on each piece. Today we will do this for one chromosome.

**Q7. What does GATK stand for?**

The examples below run GATK on chromosome 17. **Run yours on the smallest chromosome instead** (not the smallest scaffold), editing the commands to match.

**Q8. What is the smallest chromosome? How did you determine this?**

Index the duplicate-marked BAM, extract your chromosome into its own BAM, and index that BAM:

```bash
samtools index assembly/${name}-mkdp.bam
samtools view -b assembly/${name}-mkdp.bam 17 > assembly/chr17.bam
samtools index assembly/chr17.bam
```

Using this BAM as input, run HaplotypeCaller in a job with 10 CPUs and 10 GB of memory:

```bash
gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R ${REF}.fna \
    -I assembly/chr17.bam \
    -O var/chr17.g.vcf.gz \
    -ERC GVCF \
    -L 17 \
    --native-pair-hmm-threads ${SLURM_CPUS_PER_TASK}
```

*Expected run time: ~21 min*

**Q9. Why did we output a GVCF instead of a VCF?**

## Genotyping

HaplotypeCaller produced an intermediate file called a GVCF. Next we genotype the GVCF. This expands the summarized blocks of non-variant regions so that every site, including sites without variants, has a genotype call and associated quality metrics.

If you had one individual per species and only cared about variant sites, you could have run HaplotypeCaller without the `-ERC` flag and stopped there. But we want to compare many species to each other, and whether a site stays the same as the reference or varies is important evolutionary information.

> **Note:** Strictly speaking, you cannot genotype each region independently and expect correct all-sites output. Normally this is where you would use `gatk CombineGVCFs` to make one genome-wide file, because all the scaffolds together inform genotyping. For this lab we will continue with our one chromosome.

```bash
gatk GenotypeGVCFs \
    -R ${REF}.fna \
    -L 17 \
    -V var/chr17.g.vcf.gz \
    -O var/chr17.vcf.gz \
    --include-non-variant-sites
```

*Expected run time: ~7 min*

## Hard filtering

Next we filter our variant and non-variant sites to remove calls we have low confidence in. Today we will remove calls where there are too few reads to make a call, where call quality is low, or where there is strong strand bias.

**Q10. What is strand bias?**

Many tools can filter VCFs. GATK's VariantFiltration keeps every site in the file and labels the FILTER column to show whether each site passed, but it runs more slowly than `bcftools filter`, which we will use today. Apply the following filter in a job with 1 GB of memory and 1 CPU:

```bash
bcftools filter -O z -o var/chr17.filtered.vcf.gz \
    -e 'FMT/DP == "." || QD<2 || INFO/DP<2 || MQ<40 || FS>60' \
    var/chr17.vcf.gz
```

Note: `INFO/DP` is total depth across all samples, while `FORMAT/DP` is per sample. In our single-sample VCF they are the same.

This command excludes (`-e`) any site that meets any of these criteria.

**Q11. What does each of these filters do?**

You can quickly inspect a zipped VCF with bcftools. Look at the first few non-header lines (`-H`):

```bash
bcftools view -H var/chr17.filtered.vcf.gz | head
```
**Q12. What additional filters might you do?**

Now we have a VCF of calls we are confident in!

## Consensus

The final step is to turn the VCF into a FASTA sequence, using the reference sequence as a backbone onto which we stitch our allele calls.

First we must decide whether to keep the same coordinate system as the reference or to include insertions and deletions. If we don't want to include indels, we have to filter them out.

**Q13. How would you filter out indels? What are the pros and cons of filtering them out or keeping them?**

We could also filter out multi-allelic SNPs. `bcftools consensus` handles these by using the first allele listed and printing a message that it is skipping the others, since there is already a call at that site.

Extract the chromosome of interest from the reference. Again, this step would not be needed if we were working with the whole genome.

```bash
samtools faidx ${REF}.fna 17 > assembly/Hsara_17.fa
samtools faidx assembly/Hsara_17.fa
```

Generate the consensus. Every site not represented in the VCF (that is, filtered out) will have its reference base replaced with an N (`-a N`), sites where the reference allele is confidently called will have thta allele and sites where an alternative allele is called will have that allele.

```bash
bcftools index var/chr17.filtered.vcf.gz
bcftools consensus \
    -a N \
    -f assembly/Hsara_17.fa \
    -o consensus.fa \
    -c consensus.chain \
    -e 'FILTER!="PASS"' \
    var/chr17.filtered.vcf.gz
```

**Q14. Which line of the above command is unnecessary or redundant? Why?**

Let's see how much missing data we have. Remember that our mean depth was fairly low to begin with, so we probably filtered out many sites. There are many ways to do this; here's one:

```bash
samtools faidx consensus.fa
awk '{print $1 "\t0\t" $2}' consensus.fa.fai | bedtools nuc -fi consensus.fa -bed - > nucdist.txt
```

**Q15. What percent of the chromosome is made up of Ns?**

