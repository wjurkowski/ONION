#workflow to integrate lipidomics and transcriptomics data
#
#run: bash adipogene.sh var1 var2 var3 var4
#
#set variables
#var1 - transcriptomics data
#var2 - lipidomics data
#var3 - lipids IDs mapping
#var4 - db to map genes: hgnc, ensp, entrez, refseq

#	STEPS
#	1. Normalize data
#	2. Find differentially expressed genes
#	3. Create subgroups of genes and lipids
#	4. Test associations within groups
#	5. Functional analysis 

#### step 1. Normalize data


#### step 2. Find differentially expressed genes 
# run differential expression analysis with siggenes
R --no-save --args trans.txt lipid.txt < SAM.R > output

# Important metabolites. 
# Provide data of detected metabolites 

#### step 3. Create subgroups of lipids and genes
#Build network comprising known molecular pathways and putative interactions from highthroughput experiments
#associate genes with lipids. Use generic lipid names to match them with reactome

#Get Reactome pathways
cut -f1 $2 | sed -n '2,$p' | sort -u > DELs.txt
while read L; do
  #map lipids and genes in reactome
  R --no-save --args $L < mapReactome.R 1> Reactome.out 2>Reactome.err
done < DELs.txt

#create groups by merging genes mapped to individual lipids/fatty acids/metabolites
cut -f2 $3 | sort -u > generics.txt
while read G; do
  #split lipid data into groups
  grep -w $G $3 | cut -f1 > gr 
  grep -w -f gr $2 > temp 
  head -1 $2 > headl
  cat headl temp > lipidomics-group_$G.txt
  rm -f temp 

  while read member; do
    if [ -f $member.$4 ]; then
      echo $member >> lipids_t
      cat $member.$4 >> group
    else
      echo 'group $G - no proteins associated with: $member'
    fi 
  done < gr
  sort -u group > g_$G.hgnc
  if [ -s "g_$G.hgnc" ]; then
    #for given set of genes iterate through the list
    #for each gene in question find first neighbors 
    #homo sapiens: 9606
    #mus musculus: 10090
    R --no-save --args g_$G.hgnc 9606  < getStringNeighbor.R 1> getN.out 2>getN.err
    #expected output file: Rextended-string.txt with list of String Ids
    #change String IDs to HGNC
    cut -d "." -f2 extended.string >  extended.ensp
    R --no-save --args extended.ensp hgnc < mapStringIDs.R 1> map.out 2>map.err
    cat g_$G.hgnc extended.ensp.hgnc | sort -u > string-g_$G.hgnc
    #select subset of transcriptomics data baed on string-extended list of gene IDs
    grep -w -f string-g_$G.hgnc $1 > temp
    #list of genes in groups
    cut -f1 temp >> genes_t
    head -1 $1 > headg
    cat headg temp > transcriptomics-group_$G.txt
    rm -f group temp gr
  else
    echo "group $G has no representation in data. Group illdefined in the input file?"
  fi
done < generics.txt
rm -f generics.txt Reactome.out Reactome.err getN.out getN.err map.out map.err

#group of ungrouped
sort -u genes_t > ingroup_g.txt
sort -u lipids_t > ingroup_l.txt
rm -f genes_t lipids_t
grep -w -v -f ingroup_l.txt $2 > lipidomics_ungrouped.txt
grep -w -v -f ingroup_g.txt $1 > transcriptomics_ungrouped.txt

#### step 4. Test associations within groups

#STATISTICS
# CCA
# Apply only for large number of observations 
# e.g. 5 diets, two genotypes, 4 samples = 40 observations
# 5 time points, 2 conditions, 4 samples
#R --no-save --args trans.txt lipid.txt 0.8 0.8 < YACCA.R 1> out1 2>err1
# for groups defined above with ensp mapping and network based enrichment
cut -f2 $1 | sort -u > generics.txt
while read gr; do
        R --no-save --args transcriptomics-group_$gr.txt lipidomics_$gr.txt 0.8 0.8 < YACCA.R 1> out1 2>err1
done < generics.txt

#run PLS
#R --no-save --args trans.txt lipid.txt < PLS.R 1> out2 2>err2
while read gr; do
        R --no-save --args transcriptomics-group_$gr.txt lipidomics_$gr.txt 0.8 0.8 < PLS.R 1> out2 2>err2
done < generics.txt

# for small number of observations
# e.g. 5 diets, two genotypes, 4 samples = 40 observations
# 5 time points, 2 conditions, 4 samples
#rCCA
R --no-save --args trans.txt lipid.txt < frcca.R 1> out3 2>err3
# for groups
while read gr; do
        R --no-save --args transcriptomics-group_$gr.txt lipidomics_$gr.txt < frcca.R 1> out3 2>err3
done < generics.txt
#sPLS
R --no-save --args trans.txt lipid.txt < sPLS.R 1> out4 2>err4
# for groups
while read gr; do
        R --no-save --args transcriptomics-group_$gr.txt lipidomics_$gr.txt < sPLS.R 1> out4 2>err4
done < generics.txt

rm -f out1 out2 out3 out4 err1 err2 err3 err4



#### step 5. Functional analysis
#GO enrichment in groups
#while read gr; do
#	R --no-save --args temp/$gr /projects/integromika/data/annotations/gene2go < letsGO.R 1> output 2>error
#done < entrez-groups_enr.txt


