# Lab 3: Repeat Annotation and Scaffolding

Oct 28, 2026 · @Dylan

## Repeat Annotation

We will look at the repeat profile of our genome and soft-mask repeats to lowercase, since many downstream steps require this.

```bash
mkdir Heliconius_doris/repeat_results
```

Activate the shared conda environment with RepeatMasker and RepeatModeler:

```bash
source /home/ddebaun/mendel-nas1/miniconda3/bin/activate
conda activate /home/ddebaun/mendel-nas1/miniconda3/envs/genome-assembly
# conda install -c bioconda repeatmasker repeatmodeler   # already installed
```

**DO NOT RUN REPEATMODELER OR REPEATMASKER. Copy the final files from my folders to answer the questions.**

### RepeatModeler

RepeatModeler is a de novo transposable element (TE) family identification and modeling package. It combines three de novo repeat finders (RECON, RepeatScout, and LtrHarvest/Ltr\_retriever), which use complementary methods to identify repeat boundaries and family relationships.

```bash
BuildDatabase -name H_doris Heliconius_doris.bp.p_ctg.fa
RepeatModeler -dir repeat_results -threads ${SLURM_CPUS_PER_TASK} -database H_doris -quick
```

**What does BuildDatabase do?**

This takes a while to run, so use the classified consensus library from my directory:

```bash
/home/ddebaun/mendel-nas1/Comparative_Genomics_Course/Heliconius_doris/repeat_results/consensi.fa.classified
```

### RepeatMasker

RepeatMasker screens DNA for interspersed repeats and low-complexity sequence. It outputs a detailed repeat annotation plus a copy of the genome with repeats masked.

The job would request 24 cores, 5 GB memory, and 5 hours. RMBlast uses 4 cores per search, so 24 cores supports 24 / 4 = 6 parallel searches (`-pa 6`). Change the paths to where these files are on your system.

```bash
RepeatMasker -qq -pa 6 -lib consensi.fa.classified -dir H_doris_mask2 Heliconius_doris.bp.p_ctg.fa -xsmall
```

**What does `-xsmall` do?**

This job takes about an hour, so copy the finished output folder instead:

```bash
/home/ddebaun/mendel-nas1/Comparative_Genomics_Course/heldorresults/repeat_results/H_doris_mask2/
```

**Q3. What percent of the genome is repetitive? What kinds of repeats are present, and in what proportions?**

## Hi-C Scaffolding (H. sara)

Hi-C is the standard way to build chromosome-level assemblies without a reference. Hi-C captures DNA that is physically close in the nucleus. Contacts are far more frequent within a chromosome than between chromosomes, and they decay with distance, so contact frequency tells us which contigs belong together and in what order.

We have no Hi-C data for *H. doris*, so we will scaffold *H. sara*, which has both PacBio HiFi and Hi-C data. I have assembled *H. sara* contigs with hifiasm (as you did for *H. doris* in Lab 1) and subsampled the Hi-C reads. Because a published chromosome-level *H. sara* assembly exists, we will be checking for the 3rd C of a good genome assembly- correctness.&#32;

### Set up

Activate my Hi-C conda environment (YaHS, bwa, samtools, samblaster, seqkit, Java):

```bash
source /home/ddebaun/mendel-nas1/miniconda3/bin/activate
conda activate /home/ddebaun/mendel-nas1/miniconda3/envs/hic-scaffolding
```

Make a working directory and copy the contigs and Hi-C reads:

```bash
mkdir HiC_sara
cd HiC_sara
cp /home/ddebaun/mendel-nas1/Comparative_Genomics_Course/HiC_sara/HelSar_hifiasm.p_ctg.fa .
cp /home/ddebaun/mendel-nas1/Comparative_Genomics_Course/HiC_sara/sara_HiC_R1.fq.gz .
cp /home/ddebaun/mendel-nas1/Comparative_Genomics_Course/HiC_sara/sara_HiC_R2.fq.gz .
```

Look at the reads and the contigs before you start:

```bash
zcat sara_HiC_R1.fq.gz | head -8
seqkit stats -a HelSar_hifiasm.p_ctg.fa
```

**Q14. How many contigs are there, and what is the contig N50?**

### How it works

1. **Map** each Hi-C read to the contigs. The two reads of a pair are mapped independently, because they may come from distant parts of the genome.
2. **Filter** PCR duplicates and low-quality alignments.
3. **Scaffold** with [YaHS](https://github.com/c-zhou/yahs), which uses the contact matrix to order and orient contigs, joining them with gaps of 200 Ns.
4. **Visualize** the contact map in Juicebox to check the result.

Because pair mates are mapped separately, we use `bwa mem -5SP`:

- `-5` marks the 5'-most alignment as primary (needed for chimeric Hi-C reads).
- `-S` skips mate rescue and `-P` skips pairing, so bwa does not force mates near each other.

**Q15. Why would normal paired-end mapping behavior be a problem for Hi-C reads?**

### hic.sh

```bash
#!/bin/bash
#SBATCH --job-name=hic_yahs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=40G
#SBATCH --reservation=ddebaun_37
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=[name]@amnh.org
#SBATCH --output=slurm-%x_%j.out

source /home/ddebaun/mendel-nas1/miniconda3/bin/activate
conda activate /home/ddebaun/mendel-nas1/miniconda3/envs/hic-scaffolding

CONTIGS=HelSar_hifiasm.p_ctg.fa
R1=sara_HiC_R1.fq.gz
R2=sara_HiC_R2.fq.gz

# 1. Index the contigs
samtools faidx $CONTIGS
bwa index $CONTIGS

# 2. Map Hi-C reads, mark and remove duplicates, write BAM (kept in read-name order for YaHS)
bwa mem -5SP -t ${SLURM_CPUS_PER_TASK} $CONTIGS $R1 $R2 \
    | samblaster -r \
    | samtools view -b -@ 4 -o sara_hic.bam -

# 3. Scaffold with YaHS
yahs $CONTIGS sara_hic.bam -o sara_yahs

# 4. Make a Juicebox contact map
juicer pre -a -o sara_JBAT sara_yahs.bin sara_yahs_scaffolds_final.agp ${CONTIGS}.fai > sara_JBAT.log 2>&1
java -Xmx32G -jar /home/ddebaun/mendel-nas1/Comparative_Genomics_Course/HiC_sara/juicer_tools.jar \
    pre sara_JBAT.txt sara_JBAT.hic.part \
    <(grep PRE_C_SIZE sara_JBAT.log | awk '{print $2" "$3}') \
    && mv sara_JBAT.hic.part sara_JBAT.hic
```

Check the log as before. Your key outputs are:

- `sara_yahs_scaffolds_final.fa`: the scaffolded genome.
- `sara_yahs_scaffolds_final.agp`: how contigs were ordered, oriented, and joined.
- `sara_JBAT.hic` and `sara_JBAT.assembly`: the contact map for Juicebox.

### Compare contigs, Hi-C scaffolds, and the published assembly

```bash
seqkit stats -a HelSar_hifiasm.p_ctg.fa sara_yahs_scaffolds_final.fa ../HelSar1/HelSar1.2_genomic.fna
seqkit fx2tab -n -l sara_yahs_scaffolds_final.fa | sort -k2,2nr | head -30
```

**Q16. What is the scaffold N50 compared with the contig N50? Did total assembly length change? Why or why not?**

**Q17. How many large scaffolds are there, and how does that compare with the number of *H. sara* chromosomes (Q5)?**

### Inspect the contact map

Download [Juicebox](https://github.com/aidenlab/Juicebox/wiki/Download) to your laptop. Copy over the map and assembly file:

```bash
scp ddebaun@mendel-master.sdmz.amnh.org:/home/ddebaun/mendel-nas1/Comparative_Genomics_Course/HiC_sara/sara_JBAT.* ~/Documents/examplefolder/
```

In Juicebox, open `sara_JBAT.hic` (**File > Open**), then load `sara_JBAT.assembly` (**Assembly > Import Map Assembly**). Each square along the diagonal should be one chromosome.

**Q18. Describe the contact map. What does strong signal along the diagonal mean? What would strong off-diagonal signal between two scaffolds suggest?**

**Q19. Do you see any likely misjoins or misoriented contigs (e.g., an "X" or bow-tie pattern)? Where?**
