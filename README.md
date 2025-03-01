# Paper_Liu2025
Codes for Paper by Jing Liu (2025)

Hi, 

Thank you for being interested in out exciting research!
In this repository you will find all the scripts and most data to generate all Figures in the Paper: "Engineering intrinsically disordered regions (IDRs) for guiding genome navigation" by Jing Liu et al. 2025.
Scripts for each figure is named by Fig#, # = 1-7. For similar plots showing up several times, one file is prepared for each. This includes abundance results, last subplots in fig.2&3, binding at individual promoters. And for figures in the supplementary materials, figS1-3 has a .mat file for each, others are included in the respective code for the main figures. 
Mat files / tables are processed sequencing data, necessary annotations for promoters, inititally processed data for sequence analysis from datasets.
What are they?
- TF : TF names
- TF_information: a table for all yeast TFs
- TF_reference_7mer, TF_reference, descTF, medianSumPromNewAll: truncation data from Kumar, 2023
- TFdescription: TF composition, DBD, nonDBD regions. And TF's known DBD-recoginized motifs
- allTFs: data from barkai lab for all TFs' summerized binding at all promtoers
- FirstBatch, denovoinfo: annotated de novo IDRs designed in this study
- info_ortho: orthologs information for selected TFs in yeast
- int_FL_ortho: residues for all orthologs above
- mer_occurpancy: sequences for in vivo 7mer motifs
- msn2_denovo_seq & nonDBD_seq: fasta files for designed mutants' sequences in this study
- opnScore: score for OPN/DPN characteristics about promoters
- promoterIDXvecFimp: identities for all sites in yeast genome, which means whether each site is inside promoters or not
- promoterLengthsORF: lengths of promoters
- seqLogoMin: to make seq logo plots by in vivo 7mer data
- seqOrdMeta, MetaPredDisO: disordered tendency predicted by MetaPredict
- Protype: promoter types
- GP: promoters annotations
- PromoterGCcontent: z score for GC% in promoters
- StrainSumProm: summerized binding at promoters for reproducible profiles
- strains: strains for above.

Sequencing data is availabel in GEO database, due to file size limitation here. 

We wanna thank you again for reading our paper and wish you enjoyed looking at our data. If anything is not clear or cannot run, please do not be hesitated to contact us. We will try our best to help you figure out as soon as we can. :) 

Jing Liu (florencejing.liu@gmail.com)
