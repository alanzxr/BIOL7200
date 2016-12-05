#!/bin/bash

#./week7.sh -a ./data/D2-DS3_paired1.fq -b ./data/D2-DS3_paired2.fq -r ./data/chr17.fa -e -o ./output/output.vcf.gz -m -c -z -v -i

outgz='v' #v is output as .vcf not zipped
hh=0
realign=0
index=0
verbose=0
merg=0
recali=0
#get opts
while getopts a:b:r:o:emczvih arg
do
    case $arg in
	a)
	    file1=$OPTARG
	    while [ ! -f $file1 ] #trying to receive first read file
	    do
		echo "pair 1 file - $file1 missing, input again"
		read file1
	    done;;
	b)
	    file2=$OPTARG #trying to receive second read file
	    while [ ! -f $file2 ]
	    do
		echo "pair 2 file - $file2 missing, input again "
		read file2
	    done;;
	r)
	    fileref=$OPTARG
	    while [ ! -f $fileref ]
	    do
		echo "reference file - $fileref missing, input again "
		read fileref
	    done;;
	o)
	    fileout=$OPTARG
	    i=0
	    while [ $i -eq 0 ]
	    do
		if [ ! -f $fileout ];then
		    i=1
		else
		    echo "output file - $fileout exists, input 1 to overwrite, input any other character to input again"
		    read del
		    if [ $del -eq 1 ];then
			rm $fileout
			i=1
		    else
			echo "input output file name : "
			read fileout
		    fi
		fi
	    done;;
	e)
	    realign=1;;
	m)
		merg=1;;
	c)
		recali=1;;
	z)
	    outgz='z';;
	v)
	    verbose=1
	    echo "SNP calling piprline - Xinrui Zhou";;
	i)
	    index=1;;
	h)
	    echo "SNP calling pipeline - Xinrui Zhou"
	    hh=1
	    echo "
-a Input reads file – pair 1
-b Input reads file – pair 2
-r Reference genome file
-o The output VCF file
-e Do reads re-alignment
-m Do reads merge
-c Do recalibration
-z If the output VCF file should be gunzipped (*.vcf.gz)
-v Verbose mode
-i Index your output BAM file (using samtools index)
-h Print usage information and exit ";;
	?)
	    echo "unknown argument $OPTARG";;
     esac
done

if [ $verbose == 0 ]
then  
    exec 1>>./output/log.txt #std output
    exec 2>>./output/log.txt #stad error output
fi

if [ $hh == 0 ]
then
    
    samtools faidx $fileref #prepare BWA indices
    
    #mapping
    bwa index $fileref
    
    if [ $verbose == 1 ]
    then
	echo "prepare the BWA indices"
    fi

    bwa mem -R '@RG\tID:foo\tSM:foo\tLB:library1' $fileref $file1 $file2 > "./data/lane.sam" #produce the alignments
    if [ $verbose == 1 ]
    then
	echo "produce the alignments"
    fi

    samtools fixmate -O bam "./data/lane.sam" "./data/lane_fixmate.bam" #clean up read pairing information and flags
    samtools sort -O bam -o "./data/lane_sorted.bam" -T "./tmp/lane_temp" "./data/lane_fixmate.bam" #sort them from name order into coordinate order

    
    java -jar ./lib/picard.jar CreateSequenceDictionary R= $fileref O="${fileref%.*}.dict" #create .dict
    samtools index ./data/lane_sorted.bam #index

    #realign
    if [ $realign == 1 ]
    then
	if [ $verbose == 1 ]
	then
	    echo "realigning procedure"
	fi
	
	java -Xmx2g -jar ./lib/GenomeAnalysisTK.jar -T RealignerTargetCreator -R "$fileref" -I "./data/lane_sorted.bam" -o "./data/lane.intervals" --known "./data/Mills1000G.b38.vcf"
	java -Xmx4g -jar ./lib/GenomeAnalysisTK.jar -T IndelRealigner -R "$fileref" -I "./data/lane_sorted.bam" -targetIntervals "./data/lane.intervals" -known "./data/Mills1000G.b38.vcf" -o "./data/lane_realigned.bam"
	samtools index ./data/lane_realigned.bam
    fi
    

    #recalibrate
    if [ $recali == 1 ]
    then
	if [ $verbose == 1 ]
	then
	    echo "recalibrating procedure"
	fi
	java -Xmx4g -jar ./lib/GenomeAnalysisTK.jar -T BaseRecalibrator -R "$fileref" -knownSites "./data/Mills1000G.b38.vcf" -I "./data/lane_sorted.bam" -o "./data/lane_recal.table"
	java -Xmx2g -jar ./lib/GenomeAnalysisTK.jar -T PrintReads -R "$fileref" -I "./data/lane_sorted.bam" --BSQR "./data/lane_recal.table" -o "./data/lane_recal.bam"
	samtools index "./data/lane_recal.bam"
    fi
    


    
    #merge reads
    if [ $merg == 1 ]
    then
	if [ $verbose == 1 ]
	then
	    echo "merging procedure"
	fi
	if [ $realign == 1 ]
	then
		samtools merge "./data/sample.bam" "./data/lane_sorted.bam" "./data/lane_realigned.bam"
	fi
	if [ $realign == 0 ]
	then
		samtools merge "./data/sample.bam" "./data/lane_sorted.bam" 
	fi
	
	samtools index "./data/sample.bam"
    fi
    


    
    #index
    if [ $index == 1 ]
    then
	if [ $verbose == 1 ]
	then
	    echo "indexing procedure"
	fi
	if [ $merg == 1 ]
	then
		samtools index "./data/sample.bam"
	fi
	if [ $merg == 0 ]
	then
		samtools index "./data/lane_sorted.bam"
	fi
	
    fi

    #variant calling
    #convert BAM file into genomic positions
    if [ $verbose == 1 ]
    then
	echo "convert BAM file to genomic positions"
    fi
	if [ $merg == 1 ]
	then 
		samtools mpileup -ugf $fileref ./data/sample.bam  | bcftools call -vmO v -o ./output/study.vcf
	fi
	if [ $merg == 0 ]
	then
    	samtools mpileup -ugf $fileref ./data/lane_sorted.bam  | bcftools call -vmO v -o ./output/study.vcf
    fi
	if [ $verbose == 1 ]
	then
	echo "filter data"
    fi
    bcftools filter -O $outgz -o $fileout  -s LOWQUAL -i'%QUAL>10' ./output/study.vcf


fi
