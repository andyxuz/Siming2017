#!/bin/bash
if [[ $# -eq 0 ]] ; then
    echo 'This bash script will take in an input from a tab-delimited file and run plot reads alignments for the normal and tumor samples'
	echo 'Run this as: ./plotsamples.bash filecontainingsampledata.txt'
	echo 'In order for this script to work, the header in the input file must be included and have the fields: '
	echo 'Tumor Type	Hugo_Symbol	Chromosome	Start_Position	Tumor_Sample_Barcode	Matched_Norm_Sample_Barcode	Mutation_Type	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2'
	echo 'The script also requires the gdc-client, a valid gdc token, a reference genome FASTA and FAI, liftOver and the hg19ToHg38 chain file, and the python plotreads script in the current directory'
	echo 'The python scripts called in the script require the packages matplotlib, pysam, numpy, and pyliftover'
	echo 'To apply the token, edit the token declaration in line 11, or run the script with a "."'
    exit 0
fi
token=$(<gdc-user-token.2017-07-17T03-00-34.828Z.txt)
echo "header" > uuids_TumorSampleBAM.txt
echo "header" > uuids_NormalSampleBAM.txt
echo "header" > startingpositions.txt
{
	read
	while IFS=$'\t', read -r type hugo chromosome start tumorbarcode normalbarcode classification type referenceallele tseqallele1 tseqallele2 protein moreinfo
	do
		uuid=$(curl -s 'https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B+%22field%22%3A%22cases.samples.portions.analytes.aliquots.submitter_id%22%2C%22value%22%3A%22'"$tumorbarcode"'%22%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%22BAM%22%7D%7D%5D%7D&fields=file_id' | grep -Po '(?<="file_id": ").*(?=", "i)')
		echo "$uuid" >> uuids_TumorSampleBAM.txt
		nuuid=$(curl -s 'https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B+%22field%22%3A%22cases.samples.portions.analytes.aliquots.submitter_id%22%2C%22value%22%3A%22'"$normalbarcode"'%22%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%22BAM%22%7D%7D%5D%7D&fields=file_id' | grep -Po '(?<="file_id": ").*(?=", "i)')
		echo "$nuuid" >> uuids_NormalSampleBAM.txt
	done
} < $1
paste uuids_TumorSampleBAM.txt $1 > $1.txt
paste uuids_NormalSampleBAM.txt $1.txt > $1.txt.txt && mv $1.txt.txt $1.txt
{
	read
	while IFS=$'\t', read -r nuuid uuid type hugo chromosome start tumorbarcode normalbarcode classification type referenceallele tseqallele1 tseqallele2 protein moreinfo
	do
		python hg19tohg38.py $chromosome $start > lifted.$tumorbarcode.txt
		{
			while IFS=$",", read -r chr startp
			do
				echo "$startp" >> startingpositions.txt
				curl --header "X-Auth-Token: "$token"" 'https://api.gdc.cancer.gov/slicing/view/'"$uuid"'?region='"$chr"':'"$((startp-1000))"'-'"$((startp+1000))"'' --output $tumorbarcode.bam
				curl --header "X-Auth-Token: "$token"" 'https://api.gdc.cancer.gov/slicing/view/'"$nuuid"'?region='"$chr"':'"$((startp-1000))"'-'"$((startp+1000))"'' --output $normalbarcode.bam
			done  < lifted.$tumorbarcode.txt 
		}
	done
} < $1.txt
paste startingpositions.txt $1.txt > $1.txt.txt && mv $1.txt.txt $1.txt
ls *.bam | xargs -n1 -P5 samtools index
{
	read
	while IFS=$'\t', read -r position nuuid uuid type hugo chromosome start tumorbarcode normalbarcode classification type referenceallele tseqallele1 tseqallele2 protein moreinfo
	do
		python p08plotreads_SZ.py $tumorbarcode.bam chr$chromosome $position "Tumor"$tumorbarcode
		python p08plotreads_SZ.py $normalbarcode.bam chr$chromosome $position "Normal"$normalbarcode
	done
} < $1.txt
tar -cvzf plotreadsTumor.tar.gz Tumor*.txt Tumor*.png --remove-files
tar -cvzf plotreadsNormal.tar.gz Normal*.txt Normal*.png --remove-files
#rm uuids_TumorSampleBAM.txt uuids_NormalSampleBAM.txt startingpositions.txt $1.txt liftedTCGA*.txt
#rm TCGA*.bam TCGA*.bam.bai
