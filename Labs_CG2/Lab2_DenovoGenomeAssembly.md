# Lab 2: Bioinformatics and Genome Assembly
# By: Dylan DeBaun

Please complete all bolded questions. Answers are due at the beginning of the next week's class.

Today we will be completing our Bioinformatics training by learning how to submit a job on the cluster and perform version control with git. We will use these new skills to begin our genomics journey- assembling a reference genome.


## Uploading your files to Github

The last lab is due at the start of each class. Let's upload your answers to those questions to a GitHub that will be shared with me via email.

First, create an account at GitHub.com, this will follow your career so best to choose a non-affiliated email.
Create a new repository on GitHub called `ComparativeGenomics2026`. Set it to Private (or Public if you want) and Check **Add a README file**. In your new repo, go to **Settings** → **Collaborators** → **Add people**. Enter my email address (dylandebaun623@gmail.com) and send the invite.

On your local computer, make sure git is installed:

```bash
git --version
```

If you get an error, install git from https://git-scm.com/downloads (Mac users can also run `xcode-select --install`).

Every git commit has an associated author, with a name and email address provided. Let's initialize this.

```bash
git config --global user.name "Your Name"
git config --global user.email "you@example.com"   # same email as your GitHub account
```

Now we need git to be able to push things to our GitHub repo, install the GitHub CLI with homebrew (https://cli.github.com).

Login to your github from your local computer
```bash
gh auth login

? Where do you use GitHub? GitHub.com
? How would you like to authenticate GitHub CLI? Login with a web browser
? What is your preferred protocol for Git operations on this host? HTTPS
? Authenticate Git with your GitHub credentials? Yes

! One-time code (3C33-9E74) copied to clipboard
Press Enter to open https://github.com/login/device in your browser... 
✓ Authentication complete.
- gh config set -h github.com git_protocol https
✓ Configured git protocol
✓ Logged in as dylandebaun

```

**Troubleshooting:** If you see `mkdir /Users/<you>/.config/gh: permission denied`, your config folder is owned by the admin account. Fix it with:
 ```bash
sudo chown -R $(whoami) ~/.config
```
Enter your computer password when prompted (nothing will appear as you type), then rerun `gh auth login`.


Back on the GitHub page for your repo, click the green **Code** button and copy the HTTPS URL. Then, in the terminal, :

```bash
cd ~/Desktop     # or wherever you want the folder to live
git clone https://github.com/YOUR-USERNAME/ComparativeGenomics2026.git
cd ComparativeGenomics2026
```

Copy your answers into this folder (use a consistent name, e.g. `Lab01_LastName.docx`):

```bash
cp ~/Downloads/Lab01_LastName.docx . #or wherever your lab file is
git add Lab01_LastName.docx
git commit -m "Add Lab 1 answers"
git push
```

Refresh the repo page on GitHub — your file should appear.

You will do this for all future labs, try to recall it more from memory each week.

**Annotate (add comments) the code in this section to explain what it is doing**

## Preparing a Job

Now let's log back into Mendel and navigate to the class folder. **How do you log back on?**

To submit a job, you tell the scheduler the resources you need and other information about the job. Let’s break down what each piece of a job script does. We are working in SLURM scheduler so we use ```#SBATCH`` lines that the scheduler can read to get that info. You often will not know exactly what kind of resources a job will use, but use your best estimate. With practice, you will start to know what sensible requests are, and programs will often mention how many resources they consume. You typically want to just slightly overestimate- the computer reserves all those resources for you for the allotted time, so if they go unused they are wasted and can't be used by others. In these labs you will see the resource requests to give you an idea of what these numbers look like. After a job runs successfully, you can probe it for important information about how much of the requested resources it ate up.

Important SBATCH commands ([https://slurm.schedmd.com/sbatch.html](https://slurm.schedmd.com/sbatch.html) [https://slurm.schedmd.com/rosetta.pdf](https://slurm.schedmd.com/rosetta.pdf) ):  
**\-J**, **\--job-name**\=\<*jobname*\> — This will be what is listed when you qstat to see the job  
**\-N**, **\--nodes**\=1 — Almost all programs work on only one node. It is very tough to get different computers to talk to each other.  
**\-n, –ntasks=1** — Similar to nodes, each computer can maximally perform one kind of task at a time, that task can have internal parallelization capabilities.  
**\-c**, **\--cpus-per-task**\=\<*ncpus*\> — How many cpus per task we want, this is the internal parallelization capability. On Mendel, 1 CPU ≈ 1 core for our purposes.  
**\--output=%x\_%j.out** , SBATCH allows us to call variables with the %, x is the jobid and j is the job name, it is helpful to output this info for referring back to jobs. I personally dont like having separate error files but you could also include a –error option  
\--partition=batch, this indicates which kind of node the scheduler should put this job on, look at what each of the different kinds of nodes can offer in terms of resources to choose this. Batch is the default if this is not included.  
**What are the 3 partitions in mendel? Hint: shown upon login**

Putting all that together gives us a job script like this.

```
#!/bin/bash
#SBATCH --job-name=ncbi_extract      # Name of the job
#SBATCH --nodes=1                        # Run on a single node
#SBATCH --partition=batch                 # What kind of node to run on
#SBATCH --reservation=CG2       #reservation tag
#SBATCH --ntasks=1                       # Single main task (process)
#SBATCH --cpus-per-task=1          # Number of threads for the task
#SBATCH --mem=1G                        # Total memory
#SBATCH --time=25:00:00                  # Max runtime (HH:MM:SS)
#SBATCH --mail-type=ALL                  # When to send email
#SBATCH --mail-user=xxx@amnh.org        # Email address
#SBATCH --output=%x_%j.out               # Stdout file (jobname_jobid.out)
#SBATCH --error=%x_%j.err                # Stderr file (jobname_jobid.err)


###commands

```
Sajesh has created what is called a reservation for this class, setting aside nodes for our use only. Please add this reservation tag to all of your jobs so that you don’t have to compete for resources with other AMNH affiliates.

In order to edit or create files inside the terminal, you will learn and make use of a text editor program. I use VIM, so that is what we will be learning in this class, but feel free to use others you are familiar with like ```nano```.

First, use vi to open a new file (or edit an old one).
```vi hifiasm.sh ```

You will see a blank screen, that is your document. You can’t type in it just yet. To enter the insert mode, type ```i```. Now you can add text to the document, using arrows to navigate lines. Create the job script like above. To escape insert mode, click ```Esc``` key. Now, to exit and write to (i.e., save) this file you have opened, type the following

```
:wq
```
w-write
q-quit

Other helpful notes for vim commands can be found here but this is the basics\! [https://vim.rtorr.com/](https://vim.rtorr.com/)  
dd- delete entire line  
```/word ```searches for 'word' in the document  
```:n``` goes to the nth line in the document  
```shift``` G goes to the bottom of the document  
```$``` go to end of line  
```0``` beginning of line  
```u``` undo the last thing you did

Remember, if you ever get stuck with vim or make a mistake you don’t know how to correct, click escape and then type ```:q!``` to exit without saving.


## De novo Genome Assembly

Let's take a look at the reads we downloaded last lab (ERR15606379), these are long reads that have been sequenced with HiFi technology. **If you did not download these reads last lab, complete this now following the steps done previously?**

After the FASTQ files are downloaded, you can remove the folders containing the intermediate files. Remove the directories from this run and the previous as these were just produced by the prefetch command to get the data but it is not in the fastq format we need. So, we now have unnecessary and large files. To remove a directory, we use rm with the \-r flag:
```
ls ERR15606379 #check inside the dir to ensure you will be deleting what you want to delete
rm -r ERR15606379
```

We will view 4 lines of the FASTQ file to understand the format. Some lines have lots of characters (b/c they are sequences) and we don’t need all that so we want to both show the first 4 lines AND truncate the lines. To perform two commands in one line we can use a **pipe (|).** A pipe allows you to take the output of one command (stdout- standard output) and use it directly as the input for another command (stdin- standard input).

Let’s look at how to use the cut command. **What does “cut \-c 1-100” do?**

```
cut --help
```
Ok now,

```
head -4 [FASTQFILE] | cut -c 1-100
```

**What does the line beginning with “IIDIII” tell us? What do I and D specifically communicate?**

It is good practice to keep your directories tidy and clean of any intermediate files you won’t be needing. But, beware, once you remove a file or directory it cannot be recovered.

We will be [Assembling HiFi reads without additional data types](https://github.com/chhylp123/hifiasm?tab=readme-ov-file#hifionly), which conveniently has a section in the ReadMe file

HiFiasm is really convenient in that it includes an internal adapter detection & trimming step for PacBio HiFi reads.  
 It automatically:

* detects PacBio adapter sequences  
* trims them  
* discards severely short reads  
* handles chimeric reads

So we will just be running the command listed there to get our genome assembly.

### Hifiasm.sh 

```
#!/bin/bash
#SBATCH --job-name=hifi_doris
#SBATCH --nodes=1
#SBATCH --reservation=CG2
#SBATCH --cpus-per-task=10
#SBATCH --mem=35G
#SBATCH --time=30:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ddebaun@amnh.org
#SBATCH --output=%x_%j.out

export PATH=/home/ddebaun/mendel-nas1/Comparative_Genomics_Course/hifiasm:$PATH
name="Heliconius_doris"
mkdir $name \#let's put results in their own directory
hifiasm -o ${name}_mine/${name} -t 10 ERR15606379.fastq 
```

When we run **hifiasm**, lots of parameters are run with their default values since we did not specify here. This can be investigated with the \-h option (default in brackets). **How many rounds of read error correction do we perform?**

**What is the PATH variable doing in this script?**

Before submitting a job that requests substantial resources and may spend time waiting in the queue, it is good practice to first test that the script runs without immediate errors.

You can do this by executing the script directly with bash. When run this way, the shell treats any \#SBATCH lines as comments and ignores them—these directives are only interpreted by the job scheduler during submission.

The first line of the script, \#\!/bin/bash, is called a *shebang* and tells the operating system that the file should be executed using the Bash shell.

```
bash hifiasm.sh
```
It should start printing lines that indicate the program is running. If it is running without issue, we can exit with the Control-C option. This does not harm files; it only stops the running process. Now submit the job.

```
sbatch hifiasm.sh
```
This should take about 3 hours to run, if we increased the cpus per task it would run faster. Let’s check the status of our job. **What status is it in?**

```
qstat \-u \[yourusername\] 
```
When it is completed, check your email for the completed job file. **How did we do with resource usage? Is there any parameter you would change next time?**  (come back to this question when job finishes)

## De novo Genome Assembly Part 2

Now let's pretend that finished so we can continue today. Copy the directory to your folder: ```/home/ddebaun/mendel-nas1/Comparative_Genomics_Course/Heliconius_doris``` to your folder. 

**Command for that?**

Take a look at a completed log file for you run.  **What is the file?**

We see our 3 rounds of kmer error correction. There are helpful resources available for interpreting this log output. https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output

**What is the homozygous read coverage threshold in our final dataset? Heterozygous read coverage threshold?**

We can see that we had considerable improvement in the number of single k-mers, which often indicate sequencing errors (true biological sequences usually repeat). **What was this number in round 1? In the last round?**


We see the fasta for haplotype 1 and haplotype 2: ```Heliconius_doris.bp.hap1.p_ctg.gfa Heliconius_doris.bp.hap2.p_ctg.gfa```

We also see one single whole genome assembly ```Heliconius_doris.bp.p_ctg.gfa``` We will be using this one for most of our analyses.
**What is this file? Define primary contig.**W

Lastly, **what format are these files in? What is a GFA?** Let’s convert it to a fasta file.
Again, info on how to do this can be found in the docs. https://hifiasm.readthedocs.io/en/latest/faq.html 

```
awk '/^S/{print ">"$2;print $3}' Heliconius_doris.bp.p_ctg.gfa > Heliconius_doris.bp.p_ctg.fa
```

In GFA files, sequence lines begin with S, followed by the contig name and sequence.

## Quality Control

It is imperative to check the quality of this genome. **What are the three Cs of a good genome?**

Let’s first get the summary statistics. For these labs, we will be using my conda environments directly. **But how would you create a yaml of this and genome-assembly directory and create your own environment from it? What is a YAML file?**

```
source /home/ddebaun/mendel-nas1/miniconda3/bin/activate
conda activate /home/ddebaun/mendel-nas1/miniconda3/envs/genome-assembly
```

For the summary stats, we can run a bunch of different programs, this one is easy:
```
assembly-stats Heliconius_doris.bp.p_ctg.fa
```

**What is our contig N50 in Mbps?**

To check completeness, we can look at the BUSCO scores. We need to check our fasta against a scaffold of known conserved orthologs. What database should we use (a few correct options)? Check out available options with ```busco --list-datasets```

You should always select the dataset that is most closely related to the assembly or gene set you are assessing. If you are unsure, you can use the --auto-lineage option to automatically select the most appropriate dataset. BUSCO will automatically download the requested dataset if it is not already present in the download folder. 

Let’s run BUSCO. Make a job file, following the format we did earlier. When we request number of threads, we can instead of specifying the number of cpus per task requested we can point to internal SLURM variables that will know how to interpret the value. See these variables here: https://docs.hpc.shef.ac.uk/en/latest/referenceinfo/scheduler/SLURM/SLURM-environment-variables.html#gsc.tab=0
These are helpful when you need to change resource demands, you now only need to change them in one place (the #SLURM lines) and they will carry throught the rest of your job script.

```
#INCLUDE THE FOLLOWING LINES IN YOUR JOB FILE
source /home/ddebaun/mendel-nas1/miniconda3/bin/activate
conda activate /home/ddebaun/mendel-nas1/miniconda3/envs/genome-assembly # we have to activate the conda env on the node we are given

busco -i [genome] -c [cpuspertask] -m [??] -l [??] -o Heliconius_doris/busco_results -r
```

**What is the correct command (fill in the brackets)?**

This BUSCO run uses 6GB of memory, make sure you adjust this accordingly. BUSCO should complete in 7 minutes.
Inspect the .txt file that was produced (for example with cat). **What are the BUSCO scores for this genome?**





