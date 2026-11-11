# Lab 3b: Short-Read Quality Control

Today we start working with short-read data. In this lab, we will quality-check, trim, and screen our reads for contamination. In the next lab, we map the cleaned reads to a reference, call variants, and build a consensus sequence.

> **Working with your own data?** Follow along using your own short reads instead of *Heliconius telesiphe*. Wherever you see `SRR8883900` or `Heliconius_telesiphe`, substitute your own SRA accession and sample name. The variables at the top of each script are there to make this easy. **Answer every question for your own dataset**, and note in your answers which accession and species you used.

## Set up

Navigate to your `Comparative_Genomics_Course` folder and make a directory structure for our short read analysis:

```bash
mkdir short_read_assembly
cd short_read_assembly
mkdir -p {qc,trimmed,assembly,logs,var}
```

Activate the environment we will be using for lab today:

```bash
source /home/ddebaun/mendel-nas1/miniconda3/bin/activate
conda activate /home/ddebaun/mendel-nas1/miniconda3/envs/shortreadassembly
# conda install -c bioconda fastqc trimmomatic bwa samtools gatk4 bcftools bbmap picard seqtk bedtools # already installed
```

## Subsampling

In Lab 1, we downloaded short-read data for *Heliconius telesiphe* (SRR8883900), with forward and reverse reads. Variant calling is computationally heavy, especially following the GATK Best Practices pipeline we will use in Part 2. To learn these pipelines in a reasonable amount of time, we will subsample our dataset to 10% of the reads:

> **Working with your own data?** Skip this step, you'll want as much data as you can get to have a great final product. 

```bash
seqtk sample -s 100 ../SRR8883900_1.fastq 0.10 > SRR8883900_1_0.10.fastq
seqtk sample -s 100 ../SRR8883900_2.fastq 0.10 > SRR8883900_2_0.10.fastq
```

The seed (`-s 100`) is fixed, so everyone using the *H. telesiphe* data will be working with the same reads.

Paired-end reads must stay synchronized. Using `seqtk sample` with the same seed works **only if read order is identical** in both files. For complete safety we would interleave the read pairs before subsampling, but for this lab we will do it this way.

## Quality control with FastQC

FastQC reads a set of sequence files and produces a quality control report for each one. The report is made of several modules, each of which helps identify a different kind of potential problem in your data.

**fastqc.sh**

```bash
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=10gb
#SBATCH --reservation=CG2
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=[name]@amnh.org
#SBATCH --output=logs/slurm-%x_%j.out

source /home/ddebaun/mendel-nas1/miniconda3/bin/activate
conda activate /home/ddebaun/mendel-nas1/miniconda3/envs/shortreadassembly

FOR="SRR8883900_1_0.10.fastq"
REV="SRR8883900_2_0.10.fastq"

fastqc $FOR $REV -t ${SLURM_CPUS_PER_TASK} -o qc
```

*Expected run time: ~2 min*

Check that the HTML reports were written to the `qc` folder. A wildcard (`*`) matches any string in that position; here we are looking for files that start with `SRR` and end with `html`:

```bash
ls qc
ls qc/SRR*html
```

### Copying the reports to your laptop

We will use `scp`, run **from your local computer**. Your local shell will try to expand `SRR8883900_*_fastqc.html` before `scp` runs, but we want the wildcard sent to the cluster so it is expanded there. A backslash (`\`) tells the shell to treat `*` as a literal asterisk instead of a wildcard. This escape character is used for any character that a command would otherwise interpret in a special way. For example, `sed` treats a period as a wildcard, so a literal period must be written `\.`. Quotation marks are another way around this.

```bash
scp [username]@mendel-master.sdmz.amnh.org:/path/to/short_read_assembly/qc/SRR8883900_\*_fastqc.html .

# Quotation marks also work:
scp "[username]@mendel-master.sdmz.amnh.org:/path/to/short_read_assembly/qc/SRR8883900_*_fastqc.html" .
```

Open the HTML files and inspect them.

**Q1. At what position does read quality typically start to drop?**

**Q2. Do you see evidence of adapter contamination?**

**Q3. Why are paired-end reads often lower quality toward the 3′ end?**

## Trimming

Next we remove adapters and trim low-quality bases using Trimmomatic.

Java programs like Trimmomatic use a dedicated pool of memory called heap space, which the program uses for its entire run. The default is usually small (around 1 GB), so we often need to increase it. The Java flag for this is `-Xmx[size]`, where `[size]` is the amount of memory the heap can use. The maximum is typically about 85% of physical memory.

**trim.sh**

```bash
#!/bin/bash
#SBATCH --job-name=trim
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=50gb
#SBATCH --reservation=CG2
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=[name]@amnh.org
#SBATCH --output=logs/slurm-%x_%j.out

source /home/ddebaun/mendel-nas1/miniconda3/bin/activate
conda activate /home/ddebaun/mendel-nas1/miniconda3/envs/shortreadassembly

FOR="SRR8883900_1_0.10.fastq"
REV="SRR8883900_2_0.10.fastq"
name="Heliconius_telesiphe"

trimmomatic PE -Xmx6g -threads ${SLURM_CPUS_PER_TASK} \
    $FOR $REV \
    trimmed/${name}_for_paired.fq.gz trimmed/${name}_for_unpaired.fq.gz \
    trimmed/${name}_rev_paired.fq.gz trimmed/${name}_rev_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:false \
    HEADCROP:10 TRAILING:3 MINLEN:40 SLIDINGWINDOW:4:15
```

*Expected run time: ~3 min*

> **Own data:** The adapter file (`TruSeq3-PE-2.fa`) should match how your library was prepared. Check the "Adapter Content" and "Overrepresented sequences" modules of your FastQC report, and the library information on your SRA record, to decide which adapter file to use. Use your FastQC results (Q1) to decide whether the `HEADCROP` and quality settings make sense for your reads.

**Q4. What does each parameter choice in this command indicate?**

**Q5. How many paired reads were kept?** (Hint: check the slurm log.)

## Contamination

Now we will look for contamination in our reads. Copy the following file to your directory:

```bash
cp /home/ddebaun/mendel-nas1/Comparative_Genomics_Course/short_read_assembly/merged_ref_3032904399224293960.fa.gz .
```

This is a gzipped FASTA file, meaning it contains genetic sequences. Let's look at the sequence headers. `zcat` displays (cats) a gzipped file, and we pipe that to `grep` to print only the header lines:

```bash
zcat merged_ref_3032904399224293960.fa.gz | grep "^>" | head -100
```

**Q6. What kind of information is in this file?**

We can use BBSplit (part of BBMap) to identify contamination. BBSplit maps reads to multiple references at once and writes the reads that do not map (the uncontaminated reads) to a file:

```bash
bbsplit.sh -Xmx8g \
    in1=trimmed/${name}_for_paired.fq.gz in2=trimmed/${name}_rev_paired.fq.gz \
    ref=merged_ref_3032904399224293960.fa.gz minid=0.95 \
    outu=trimmed/${name}_for_paired_clean.fq.gz outu2=trimmed/${name}_rev_paired_clean.fq.gz
```

We will skip this step because reads are typically uploaded to NCBI with contaminants already removed. I ran it for you already; copy over the log file:

```bash
cp /home/ddebaun/mendel-nas1/Comparative_Genomics_Course/logs/slurm-bbmap_6623098.out .
```

> **Own data:** This log is for *H. telesiphe* only. To answer Q7 for your own reads, put the `bbsplit.sh` command above in a job file (8 GB memory, 1 CPU) and run it on your trimmed reads.

**Q7. What percent of reads matched contaminants (mapped to those references)?**

## Checking the QC'd reads

Run FastQC again on the trimmed reads by editing your `fastqc.sh` job file. Change `FOR` and `REV` to point to the trimmed paired files (`trimmed/${name}_for_paired.fq.gz` and `trimmed/${name}_rev_paired.fq.gz`, or the `_clean` files if you ran BBSplit).

*Expected run time: ~2 min*

**Q8. Compare the FastQC results from earlier with the new ones. What has changed?**

**Q9. What additional trimming might you do?**

Keep your trimmed reads. You will map them to a reference in the next lab.

> **Working with your own data?** Congrats you have homework! Find the closest related species with a decent de novo genome assembly to serve as your reference. Try to find something with >1Mb N50 and download it to your Mendel folder for the next class.
