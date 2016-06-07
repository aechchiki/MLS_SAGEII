# SAGE_2 - ex_1
# A. Echchiki; M. Borer
# group 8 - SA4

# login to server
ssh aechchik@prd.vital-it.ch
# copy from server: scp username@prd.vital-it.ch:/serverpath/serverfile /localpath/
# copy to server: scp /localpath/localfile username@prd.vital-it.ch/serverpath

# Check your working folder
cd /scratch/cluster/monthly/mls_2015/
# create wdir (if has been erased)
mkdir /scratch/cluster/monthly/mls_2015/aechchik

# Get the genome, locate the reads
# 1. move to wdir
cd /scratch/cluster/monthly/mls_2015/aechchik
# 2. get .fasta, .gbk, .gtf in wdir
cp -r /scratch/cluster/monthly/mls_2015/SAGE/genome .
# 3. check data in copied folder, get size content
cd /scratch/cluster/monthly/mls_2015/SAGE/genome; du -h *
# 4. we have 4 conditions,
# 5. we have 4 replicates per condition
# DON'T copy reads, just get a look at the files 
ls /scratch/cluster/monthly/mls_2015/SAGE/RNAseq/SA # we'll be using: SA4.fastq.gz, 4.1G
less SA4.fastq.gz

# Explore those files
# answers 
# 1. gbk contains annotations for coding regiosn in the genome
# 2. nucleotides in the genome: wc -c in .fasta (rm header)
cd /scratch/cluster/monthly/mls_2015/SAGE/genome
cat Pseud_protegens_S5_final.fasta | sed '1d' | wc -c
# same as: consider only the 2nd line, no newline
cat Pseud_protegens_S5_final.fasta | tail -n 1| tr -d '\n' | wc -c
# 3. genome is 6.8 Mb

# Submit a quality control job
# load fastqc 
module add UHTS/Quality_control/fastqc/0.11.2;
# create output directory for quality control 
mkdir /scratch/cluster/monthly/mls_2015/aechchik/SA4_qc
# submit job 
bsub -n 2 -J SA4_qc fastqc -t 2 -o /scratch/cluster/monthly/mls_2015/aechchik/SA4_qc /scratch/cluster/monthly/mls_2015/SAGE/RNAseq/SA/SA4.fastq.gz
#control running jobs 
bjobs
# output is .html 

# Extract CDs 
# load cufflinks
module add UHTS/Assembler/cufflinks/2.2.1; 
# prepare output directory 
mkdir /scratch/cluster/monthly/mls_2015/aechchik/SA4_cds
cd /scratch/cluster/monthly/mls_2015/aechchik/SA4_cds
# run job
gffread -g /scratch/cluster/monthly/mls_2015/SAGE/genome/Pseud_protegens_S5_final.fasta -x Pseud_protegens_S5_cds.fasta /scratch/cluster/monthly/mls_2015/SAGE/genome/Pseudomonas_protegens_S5_genome.gtf 
# how many transcripts? 
awk '/>/{ print }' Pseud_protegens_S5_cds.fasta | wc -l
# we have 6088 cds 

# Read a quality control output
scp aechchik@prd.vital-it.ch:/scratch/cluster/monthly/mls_2015/aechchik/SA4_qc/SA4_fastqc.html /home/aechchiki/Documents/SAGEII
# there are some red flags, but no problem, see Kamil's comments

# Adapter trimming and quality filtering
# load trimmomatic
module add UHTS/Analysis/trimmomatic/0.33;
# submit job 
bsub -n 8 -M 8388608 -J trim_adapt_SA4 trimmomatic SE -threads 8 /scratch/cluster/monthly/SAGE/RNAseq/SA/SA4.fastq.gz /scratch/cluster/monthly/mls_2015/aechchik/SA4_qc/SA4t.fastq.gz ILLUMINACLIP:<adapters.fasta>:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:60

# Fastqc of trimmed reads
# load fastqc 
module add UHTS/Quality_control/fastqc/0.11.2;
# create output directory for quality control 
mkdir /scratch/cluster/monthly/mls_2015/aechchik/SA4t_qc
# submit job 
bsub -n 2 -J SA4t_qc fastqc -t 2 -o /scratch/cluster/monthly/mls_2015/aechchik/SA4t_qc /scratch/cluster/monthly/mls_2015/aechchik/SA4_qc/SA4t.fastq.gz
#control running jobs 
bjobs
# output is .html 
