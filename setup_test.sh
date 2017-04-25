PROJECT_DIR=$(pwd)/test
N_LINES=10000

mkdir -p ${PROJECT_DIR}/genome
cd ${PROJECT_DIR}/genome
wget http://hgdownload-test.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes
wget http://hgdownload-test.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz
{ tar -O -xf chromFa.tar.gz && rm chromFa.tar.gz; } | bgzip -c > sacCer3.fa.gz
bwa index sacCer3.fa.gz

mkdir -p ${PROJECT_DIR}/fastq/MATalpha_R1/lane1/
mkdir -p ${PROJECT_DIR}/fastq/MATalpha_R1/lane2/
mkdir -p ${PROJECT_DIR}/fastq/MATalpha_R2/lane1/
mkdir -p ${PROJECT_DIR}/fastq/MATa_R1/lane1/
mkdir -p ${PROJECT_DIR}/fastq/MATa_R2/lane1/


cd ${PROJECT_DIR}/fastq/MATalpha_R1/lane1/
fastq-dump -F -X ${N_LINES} SRR2601842 --split-files --gzip
cd ${PROJECT_DIR}/fastq/MATalpha_R1/lane2/
fastq-dump -F -X ${N_LINES} SRR2601843 --split-files --gzip
cd ${PROJECT_DIR}/fastq/MATalpha_R2/lane1/
fastq-dump -F -X ${N_LINES} SRR2601845 --split-files --gzip
cd ${PROJECT_DIR}/fastq/MATa_R1/lane1/
fastq-dump -F -X ${N_LINES} SRR2601848 --split-files --gzip
cd ${PROJECT_DIR}/fastq/MATa_R2/lane1/
fastq-dump -F -X ${N_LINES} SRR2601851 --split-files --gzip

