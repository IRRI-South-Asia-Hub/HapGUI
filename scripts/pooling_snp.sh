#!/bin/bash

helpFunction()
{
	echo ""
   	echo "Usage: bash $0 -h hapmap.txt -p ped/map -f fam -v vcf -t <hap2ped/ped2hap/fam2hap/fam2vcf>"
   	echo -e "\t-h trait"
	exit 1
}

# inputs
while getopts "h:" option
do
case "${option}"
in
h) infile=${OPTARG};;
esac
done

awk 'OFS="\t" {print $3,$4}' et_marker.afreq > $infile
sed -i "1s/.*/Ref\tAlt/" $infile
