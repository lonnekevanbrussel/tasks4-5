## 4. EN-TEx ATAC-seq data: downstream analyses
## TASK 4.1
## Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files. 
## Make sure the files are organized in a consistent way as done for ChIP-seq.
cd ..
mkdir ATAC-seq
cd ATAC-seq
mkdir data
mkdir data/bigBed.files data/bigWig.files
mkdir annotation
mkdir analyses
mkdir analyses/peaks.analysis

## TASK 4.2
## Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon 
## for the same donor used in the previous sections. 
## Hint: have a look at what we did here. Make sure your md5sum values coincide with the ones provided by ENCODE.

## From encode: select epigenomical data from the same individual (ENCDO451RUA) => experiment search
## assay type DNA accessibility, biosample  term name stomach and sigmoid colon, assay title ATAC-seq
## Download selected files=> download metadata file 

../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assay_slims=DNA+accessibility&assay_title=ATAC-seq&type=Experiment"

## bigBed peak calling files (bigBed narrow, pseudoreplicated peaks, assembly GRCh38, most recent file for each tissue):

grep -F "bigBed_narrowPeak" metadata.tsv |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.txt

## ENCFF287UHP	sigmoid_colon	
## ENCFF762IFP	stomach	

cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done

bigWig FC files (bigWig files, fold-change, assembly GRCh38, most recent file for each tissue):

grep -F "bigWig"  metadata.tsv |\
grep -F "fold_change_over_control" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigWig.FC.ids.txt

## ==>
## ENCFF997HHO sigmoid_colon
## ENCFF841ZHA stomach

cut -f1 analyses/bigWig.FC.ids.txt |\
while read filename; do
  wget -P data/bigWig.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigWig"
done

## Check the integrity of the downloaded files:

for file_type in bigBed bigWig; do
  # retrieve original MD5 hash from the metadata
  ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt

  # compute MD5 hash on the downloaded files 
  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  # make sure there are no files for which original and computed MD5 hashes differ
  awk '$2!=$3' data/"$file_type".files/md5sum.txt
done

## TASK 4.3
## For each tissue, run an intersection analysis using BEDTools: report 
## 1) the number of peaks that intersect promoter regions, 
## 2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions)

## 1) the number of peaks that intersect promoter regions
## Because we want to get the number of peaks that intersect promoter regions I put the data/bed.files/"$filename".bed first, because the unique peak names are called here

cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b ../ChIP-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -u |\
  cut -f1-4 |\
  sort -u > analyses/peaks.analysis/peaks.intersect.TSS."$tissue".ATACseq.bed
done	

## counts
wc -l analyses/peaks.analysis/peaks.intersect.TSS.*.bed
## 47871 analyses/peaks.analysis/peaks.intersect.TSS.sigmoid_colon.ATACseq.bed
## 44749 analyses/peaks.analysis/peaks.intersect.TSS.stomach.ATACseq.bed

## 2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions)
## bedtools -v Only report those entries in A that have no overlap in B. So the peaks that are nót in the gene body. And I want the unique peaks (column 4)

cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed -v |\
  cut -f1-4 |\
  sort -u > analyses/peaks.analysis/peaks.outside.gene.body."$tissue".ATACseq.bed
done

## counts
wc -l analyses/peaks.analysis/peaks.outside.gene.body.*.bed
## 37035 analyses/peaks.analysis/peaks.outside.gene.body.sigmoid_colon.ATACseq.bed
## 34537 analyses/peaks.analysis/peaks.outside.gene.body.stomach.ATACseq.bed

## TASK 5
## Task 5.1: Create a folder regulatory_elements inside epigenomics_uvic. This will be the folder where you store all your subsequent results.
cd ..
mkdir regulatory_elements

## Task 5.2: Distal regulatory regions are usually found to be flanked by both H3K27ac and H3K4me1. 
## From your starting catalogue of open regions in each tissue, select those that overlap peaks of H3K27ac AND H3K4me1 in the corresponding tissue. 
## You will get a list of candidate distal regulatory elements for each tissue. How many are they?

## Get metadata file from ChIPseq folder
mkdir data/bigBed.files data/bigWig.files
## bigBed peak calling files for H3K27ac and H3K4me1

grep -F H3K27ac ../ChIP-seq/metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u  > analyses/bigBed.H3K27ac.peaks.ids.txt

cut -f1 analyses/bigBed.H3K27ac.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done

grep -F H3K4me1 ../ChIP-seq/metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u  > analyses/bigBed.H3K4me1.peaks.ids.txt

cut -f1 analyses/bigBed.H3K4me1.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done

## bigbed to bed files 
mkdir data/bed.files

cut -f1 analyses/bigBed.H3K27ac.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done

cut -f1 analyses/bigBed.H3K4me1.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done

## From your starting catalogue of open regions in each tissue (ATACseq peaks), select those that overlap peaks of H3K27ac AND H3K4me1 in the corresponding tissue.
## first get the peaks outside gene body (output of ATAC-seq outside gene body) and use bedtools to intersect the H3K27ac and H3K4me1 peaks.. 

cut -f-2 analyses/bigBed.H3K27ac.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a ../ATAC-seq/analyses/peaks.analysis/peaks.outside.gene.body.$tissue.ATACseq.bed -b data/bed.files/"$filename".bed -u |\
  cut -f1-4 |\
  sort -u > analyses/peaks.analysis/peaks.outside.gene.body."$tissue".H3K27ac.bed
done

cut -f-2 analyses/bigBed.H3K4me1.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a ../ATAC-seq/analyses/peaks.analysis/peaks.outside.gene.body.$tissue.ATACseq.bed -b data/bed.files/"$filename".bed -u |\
  cut -f1-4 |\
  sort -u > analyses/peaks.analysis/peaks.outside.gene.body."$tissue".H3K4me1.bed
done

## concatenate the files for H3K27ac and H3K4me1 per tissue and get the unique peaks
cat analyses/peaks.analysis/peaks.outside.gene.body.sigmoid_colon.H3K27ac.bed analyses/peaks.analysis/peaks.outside.gene.body.sigmoid_colon.H3K4me1.bed | sort -u > analyses/peaks.analysis/peaks.outside.gene.body.sigmoid_colon.H3K27ac.H3K4me1.bed
cat analyses/peaks.analysis/peaks.outside.gene.body.stomach.H3K27ac.bed analyses/peaks.analysis/peaks.outside.gene.body.stomach.H3K4me1.bed | sort -u > analyses/peaks.analysis/peaks.outside.gene.body.stomach.H3K27ac.H3K4me1.bed

## counts
wc -l analyses/peaks.analysis/peaks.outside.gene.body.*.H3K27ac.H3K4me1.bed
## 23130 analyses/peaks.analysis/peaks.outside.gene.body.sigmoid_colon.H3K27ac.H3K4me1.bed
## 18028 analyses/peaks.analysis/peaks.outside.gene.body.stomach.H3K27ac.H3K4me1.bed


## Task 5.3: Focus on regulatory elements that are located on chromosome 1 (hint: to parse a file based on the value of a specific column, have a look at what we did here), 
## and generate a file regulatory.elements.starts.tsv that contains the name of the regulatory region (i.e. the name of the original ATAC-seq peak) 
## and the start (5') coordinate of the region.

#outside gene peak overlap files per tissue, filter on chr1 and get the unique peaks
for tissue in sigmoid_colon stomach; do
	cat analyses/peaks.analysis/peaks.outside.gene.body.$tissue.H3K27ac.H3K4me1.bed | awk 'BEGIN{FS="\t"}$1=="chr1"{print $4,$2}' | sort -u > regulatory.elements.starts.$tissue.tsv
done

## counts
wc -l regulatory.elements.starts*.tsv
## 2452 regulatory.elements.starts.sigmoid_colon.tsv
## 2036 regulatory.elements.starts.stomach.tsv

## Task 5.4: Focus on protein-coding genes located on chromosome 1. From the BED file of gene body coordinates that you generated here, 
## prepare a tab-separated file called gene.starts.tsv which will store the name of the gene in the first column, and the start coordinate of the gene on the second column 
## (REMEMBER: for genes located on the minus strand, the start coordinate will be at the 3'). Use the command below as a starting point:

awk 'BEGIN{FS=OFS="\t"}$1=="chr1"{if ($6=="+"){start=$2} else {start=$3}; print $4, start}' ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed  > gene.starts.tsv

## Task 5.5: Download or copy this python script inside the epigenomics_uvic/bin folder. Have a look at the help page of this script to understand how it works:

wget -P ../bin https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/get.distance.py

for line in open_input.readlines(): # for each line in the input file
        gene, y = line.strip().split('\t') # split the line into two columns based on a tab
        position = int(y) # define a variable called position that correspond to the integer of the start of the gene
        abs_value = abs(position - enhancer_start) # compute the absolute value of the difference between position and enhancer_start

        if abs_value < x: # if this absolute value is lower than x
                x = abs_value # this value will now be your current x
                selectedGene = gene # save gene as selectedGene
                selectedGeneStart = position # save position as selectedGeneStart

print "\t".join([selectedGene, str(selectedGeneStart), str(x)])

## Task 5.6. For each regulatory element contained in the file regulatory.elements.starts.tsv, 
## retrieve the closest gene and the distance to the closest gene using the python script you created above. Use the command below as a starting point:

for tissue in sigmoid_colon stomach; do
 cat regulatory.elements.starts.$tissue.tsv | while read element start; do 
  python ../bin/get.distance.py --input gene.starts.tsv --start $start 
 done > regulatoryElements.genes.distances.$tissue.tsv
done

## counts
wc -l regulatoryElements.genes.distances*.tsv
## 2452 regulatoryElements.genes.distances.sigmoid_colon.tsv
## 2036 regulatoryElements.genes.distances.stomach.tsv

## Task 7: Use R to compute the mean and the median of the distances stored in regulatoryElements.genes.distances.tsv

## RegEl.distances.R
suppressPackageStartupMessages(library("optparse"))
option_list <- list (
  make_option ( c("--input"),
                help = "input gene distance matrix" ),
  make_option ( c("--output"), help = "Output file name")
)
parser <- OptionParser(
  usage = "%prog [options]",
  option_list=option_list,
  description = "\nReports mean and median values of distance matrix."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

# 1. read distance matrix
distances.matrix <- read.table( file = opt$input, sep="\t", header=F)

# 2. calculate mean and median
mean <- round(mean(distances.matrix[,3]),2)
median <- round(median(distances.matrix[,3]),2)

# 3. print mean and media
print(paste("Mean: ",mean))
print(paste("Median: ",median))

#####

## get mean and median per tissue

Rscript ../bin/RegEl.distances.R --input regulatoryElements.genes.distances.stomach.tsv
##[1] "Mean:  54805.9"
##[1] "Median:  29458"

Rscript ../bin/RegEl.distances.R --input regulatoryElements.genes.distances.sigmoid_colon.tsv
[1] "Mean:  75054.64"
[1] "Median:  34426"