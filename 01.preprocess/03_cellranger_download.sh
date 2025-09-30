### I. CellRanger Installation
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

## 1. Cell Ranger (tar.gz compression) - 9.0.1 (December 7, 2022)
# 1.1 Download
cd /mnt/rdisk/gliomas/data/cellranger
curl -o cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1756483256&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=l9R6usCDeFDO3CyfZeU9Kx9YmcF3Vn5u7OybO72WmZiIqq-PRPVIxwewvCU1dmgOtu2PUCGNSt4YteR5DpDO8QHyDNw4eQeBCGrJNGqTd~KQswvxUkkgsYSA4qcCkmycSgwl75ji~tgY~2SrEYxIOON4atS1vFwxSwJudfKTLM0nLa2yimHE01uc0RRt9kJredZ8ICFSByvDAs56kqAYH3VEXEk6M9Da5O0sWghL-F~mf5iXT6K72hyrDIbBk0WoAN-uaVeIV4cfQq5YcRZPz1nnACLcsR4kpw07-vCja9jGOR7UE~fYYaSVXnBqYtc0rVJ~HoMXkpr6ti8dp1X2Gw__"ls /mnt/rdisk/gliomas/data/cellranger
tar -xzvf cellranger-9.0.1.tar.gz

# 1.3 Prepend the Cell Ranger directory to your $PATH. 
export PATH=/mnt/12T/chibao/data/cellranger/cellranger-9.0.1/:$PATH

# 1.4 Verify installation
export PATH=/mnt/12T/chibao/data/cellranger/cellranger-9.0.1/:$PATH
cellranger testrun --id=tiny

## 2. Download Reference data
#References - 2020-A (July 7, 2020)
cd /mnt/12T/chibao/data/cellranger/ref_genome   #directory contain Reference genome - optional
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz

tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

# ## 3. Download dataset (FASTQ files from publicly-available data sets on the 10x Genomics)
# #https://www.10xgenomics.com/resources/datasets/200-sorted-cells-from-human-glioblastoma-multiforme-3-lt-v-3-1-3-1-low-6-0-0
# curl -O https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Brain_Tumor_3p_LT/Brain_Tumor_3p_LT_fastqs.tar
# tar -xvf Brain_Tumor_3p_LT_fastqs.tar

# cd /media/bio03/DATA/quyen/THESIS/trial/cellranger/Brain_Tumor_3p  #the output directory - optional
# export PATH=/opt/cellranger/7.1.0/:$PATH
# cellranger count --id=Brain_Tumor_3p_L001 \
#                  --fastqs=/media/bio03/DATA/quyen/THESIS/trial/raw/Brain_Tumor_3p/L001 \
#                  --sample=Brain_Tumor_3p_LT \
#                  --transcriptome=/media/bio03/DATA/annotations/genome/refdata-gex-GRCh38-2020-A \
#                  --localcores=1 \
#                  --localmem=8
 
# #--id: is the output directory name (optional)
# #--fastqs: a path to the directory containing the FASTQ files
# #--sample: the sample id at the beginning of the FASTQ file name
#       #example: we have 4 file fastq from dataset 10x genomics
#       Brain_Tumor_3p_LT_S1_L001_I1_001.fastq.gz  
#       Brain_Tumor_3p_LT_S1_L001_I2_001.fastq.gz  
#       Brain_Tumor_3p_LT_S1_L001_R1_001.fastq.gz  
#       Brain_Tumor_3p_LT_S1_L001_R2_001.fastq.gz  

#       --> --sample=Brain_Tumor_3p_LT
# #--transcriptome: a path to the directory containing the tramscriptome
# #--localcores: number core
# #--localmem: 

# ## 5. Output running CellRanger


# cellranger count --id="SAMN37643007" \
#                  --transcriptome="/mnt/rdisk/gliomas/data/cellranger/ref_genome/refdata-gex-GRCh38-2020-A" \
#                  --fastqs="/mnt/rdisk/gliomas/data/official_data/fastq_cellranger/PRJNA1023164/SAMN37643007" \
#                  --sample="SAMN37643007" \
#                  --create-bam=false