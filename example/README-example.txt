To run the ONION workflow to integrate lipidomics and transcriptomics data follow steps below:

> bash onion.sh nm-transcriptomics.txt nm-lipidomics.txt fa_mapping.txt entrez

Both gene expression and metabolite levels should be provided as tab-delimited text files. 
The current version does not yet rely on automatized mapping of individual metabolites levels into IDs recognised by Reactome.

Coming soon:

1) clustering of small molecules based on 
	- molecular similarity
	- shared pathways
	- ChEBI ontology
2) improved strategies for defining functional groups

