
#####
#!/bin/bash

#$ -cwd
#$ -j yes

#Script that compares specified genomic region in query genotype to reference using Mummer


chromosome=$1
start1=$2 #query region
end1=$3
start2=$4 #subject region (Barke)
end2=$5

genotypequery=$6
genotypesubject=Barke
#genotypesubject=RGT_Planet
#genotypesubject=B1K-4-12
genome=$7

queryin=${genotypequery}_genes_in_gff_${chromosome}_${start1}_${end1}
subjectin=${genotypesubject}_genes_in_gff_${chromosome}_${start2}_${end2}


mummer_query=${genotypequery}_${chromosome}_${start1}_${end1}_sequence.fasta
mummer_subject=${genotypesubject}_${chromosome}_${start2}_${end2}_sequence.fasta

mummerprefix=${genotypequery}_${genotypesubject}

#reference_genome=RGT_Planet_pseudomolecule_v1.fasta
reference_genome=180903_Barke_Unfiltered_chloro_clean_pseudomolecules_v1.fasta
#reference_genome=B1K-04-12_pseudomolecules_v1.fasta


##### First get region of interest (.fasta). The script pulls out gff region of interest as well.


python /cluster/db/mecoulter/genomes/pan_genome_tools.py -g $genotypequery -c ${chromosome},${start1},${end1} -genome $genome -o $queryin

python /cluster/db/mecoulter/genomes/pan_genome_tools.py -g $genotypesubject -c ${chromosome},${start2},${end2} -genome $reference_genome -o $subjectin


#Now run mummer to align the two regions
source activate mummer

 
nucmer -p $mummerprefix $mummer_subject $mummer_query
dnadiff -d ${mummerprefix}.delta

delta-filter ${mummerprefix}.delta > ${mummerprefix}_filtered.delta -i 95 -l 1000 -g

mummerplot --filter --layout -p ${mummerprefix}_clean_${chromosome}_${start2}_${end2} -t png -R $mummer_subject -Q $mummer_query ${mummerprefix}_filtered.delta

python mummerplot_hack.py -i ${mummerprefix}_clean_${chromosome}_${start2}_${end2}.gp -o ${mummerprefix}_clean_${chromosome}_${start2}_${end2}_fix.gp #Fixes gnuplot input by commenting out lines
gnuplot ${mummerprefix}_clean_${chromosome}_${start2}_${end2}_fix.gp

conda deactivate