# This code requires:
   
# Trimmomatic (V. 0.36)
 
# Read files should be formated as such:
# prefix_R#.fastq.gz
 
# The only argument needed for this script is the unique sample ID
prefix="$1"
 
# create an initial fastqc report of paired end reads
# This will be used to determine the effectiveness of the trimming
echo "initial fastqc report"
fastqc -o . -t 2 ../rawdata/"$prefix"_*.fastq.gz
mv "$prefix"_1_fastqc.html "$prefix"_1_unfiltered_fastqc.html
mv "$prefix"_2_fastqc.html "$prefix"_2_unfiltered_fastqc.html
mv "$prefix"_1_fastqc.zip "$prefix"_1_unfiltered_fastqc.zip
mv "$prefix"_2_fastqc.zip "$prefix"_2_unfiltered_fastqc.zip
 
#Remove adapter squences and stretches of low quality base calls
#It also sets the minimum read size to 32
 
echo "filtering reads"
 
java -jar $TRIMMOMATIC_HOME/trimmomatic-0.39.jar PE -threads 16 -phred33 \
        ../rawdata/"$prefix"_1.fastq.gz ../rawdata/"$prefix"_2.fastq.gz \
        ../data.trimmed/"$prefix"_1_P.fastq.gz ../data.trimmed/"$prefix"_1_UP.fastq.gz \
        ../data.trimmed/"$prefix"_2_P.fastq.gz ../data.trimmed/"$prefix"_2_UP.fastq.gz \
        ILLUMINACLIP:../resources/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:32
 
rm ../data.trimmed/"$prefix"_*_UP.fastq.gz
 
# concatenate all data for R1 inlcuding reads not found in R2
#cat data_trimmedFH/"$prefix"_R1_P.fastq.gz data_trimmedFH/"$prefix"_R1_UP.fastq.gz > \
#       data_trimmedFH/"$prefix"_R1_A.fastq.gz
 
# create a fastqc report of filtered reads
# you should see: high phred scores across the entire
# read length, a read length distribuition that reflects
# expected length for RNA-seq data and no overrepresented sequences
# looking at the GC distribution can also give you a
# sense of bacterial contamination (curvibacter GC% is around 60)
 
echo "post-filtering fastqc report"
 
fastqc -o . -t 2 ../data.trimmed/"$prefix"_*_P.fastq.gz
 
 
# Maps reads to the reference H. oligactis transcriptome
# and count the number of aligned reads per transcript
 
echo "mapping reads"
 
# Rsem will output a genes.results file that will be used downstream to
# generate the read count matrix

rsem-calculate-expression --num-threads 16 --no-bam-output \
        --bowtie2 --paired-end \
        ../data.trimmed/"$prefix"_1_P.fastq.gz ../data.trimmed/"$prefix"_2_P.fastq.gz \
        reducedRef/HoCR "$prefix"_RNA_A

rm "$prefix"*isoforms.results
rm -r "$prefix"*.stat

echo "done"

