## BinarizeBed/BinarizeBam
java -mx4000M -jar /path/ChromHMM.jar BinarizeBam /path/cs.genome /path/input_bam cell_mark.txt BinarizeBam_output

##LearnModel
unset DISPLAY
java -mx102400M -jar /path/ChromHMM.jar LearnModel -p 10 -lowmem -printposterior -printstatebyline -u /path/cs_genomefeature BinarizeBam_output LearnModel_output 15 /path/CS.fa

##overlap enrichment
java -mx4000M -jar /path/ChromHMM.jar OverlapEnrichment State_15_segments.bed /path/cs_genomefeature ./overlap_enrichment/KHL_overlap

##neighborhood enrichment
java -mx4000M -jar /path/ChromHMM.jar NeighborhoodEnrichment State_15_segments.bed /path/CS_TSS_position.bed ./anchor/KHL_TSS

##ano gene
bedtools intersect -wa -wb -a /path/cs_gene.bed -b KHL_15_segments.bed -F 0.5 > gene_ano.bed

####################################################################################################################################################################################
##reorder chromatin states for plots
unset DISPLAY
java -mx4000M -jar /path/ChromHMM.jar Reorder -o reorder.txt LearnModel15_output/model_15.txt ./Reorder_model

##make new segmentation
java -mx4000M -jar /path/ChromHMM.jar MakeSegmentation -printposterior -printstatebyline Reorder_model/model_15.txt BinarizeBam_output Reorder_seg2

##MakeBrowserFiles
java -mx4000M -jar /path/ChromHMM.jar MakeBrowserFiles ./State_15_segments.bed State_15 State_15


