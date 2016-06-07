# SAGE_2 - ex_2
# A. Echchiki; M. Borer
# group 8 - SA4

# Indexing the reference
mkdir /scratch/cluster/monthly/mls_2015/aechchik/kallisto 

# Activate modules kallisto and samtools
module add UHTS/Analysis/kallisto/0.42.4;
module add UHTS/Analysis/samtools/latest;

# Use kallisto to index the file with coding regions.
bsub -q priority -J kallisto_idxSA4 "kallisto index -i Pseud_protegens_S5_cds.idx /scratch/cluster/monthly/mls_2015/aechchik/SA4_cds/Pseud_protegens_S5_cds.fasta"

# Use kallisto for pseudoalignment
bsub -J kallisto_pseudobamSA4 "kallisto quant -i  Pseud_protegens_S5_cds.idx -o /scratch/cluster/monthly/mls_2015/aechchik/kallisto -b 100 --single -l 320 -s 60 --pseudobam /scratch/cluster/monthly/mls_2015/aechchik/SA4_qc/SA4t.fastq.gz > SA4t.sam"

# Conversion of pseudoaligment for IGV
mkdir /scratch/cluster/monthly/mls_2015/aechchik/scripts
cp /scratch/cluster/monthly/mls_2015/SAGE/scripts/ /scratch/cluster/monthly/mls_2015/aechchik/scripts

# modify the content of submit_conversion.sh with real filenames for SA4:
#
# python3 kallisto_sam_convertor.py /scratch/cluster/monthly/mls_2015/aechchik/SA4_cds/SA4t.sam /scratch/cluster/monthly/mls_2015/SAGE/genome/Pseudomonas_protegens_S5_genome.gtf | samtools view -bS - | samtools sort - -o /scratch/cluster/monthly/mls_2015/aechchik/SA4_cds/SA4t.bam

# Submit a job using edited submit_conversion.sh script and bsub command
bsub < submit_conversion.sh

# Indexing of bam file
cd /scratch/cluster/monthly/mls_2015/aechchik/SA4_cds/
samtools index SA4t.bam

# Download files 

# Download IGV
# load .fasta; .gtf; .bam; .bam.bai to IGV
