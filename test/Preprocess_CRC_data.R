## Author: Sambhawa Priya
## priya030@umn.edu

## Oct 21, 2019
## Preprocess CRC gene expression and microbiome data 
## This code is adapted from Preprocess_CRC_gene_expr.R and Preprocess_CRC_taxa_abnd.R

## Edited in March 2020 -- remove low abundant contaminant taxa
## Edited in July 2020 -- output case+control taxa abundance without CLR transforming data 

rm(list=ls())

library(DESeq2) ## for normalizing gene expression data
library(data.table) ## for fread()


setwd("/Users/priya030/Documents/UMN_research/gene_taxa_interaction_diseases/") ## project directory

################# Functions for processing OTU table ####################
## Normalize taxa names
## Convert taxa names to form kingdom;phylum;...;species
normalize_taxaname_CRC <- function(taxa){
  ##debug
  # microbe_name <- "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__" ## for debugging
  ## debug
  microbe_name <- taxa
  split <- strsplit(microbe_name,split="\\; [a-z]__")
  microbe_name <- paste(split[[1]],collapse=";") # Create collapsed names
  microbe_name <- gsub("(;+$)","",microbe_name,perl=T) ## get rid of any trailing ;
  microbe_name <- gsub("[A-z]__","",microbe_name,perl=T) ## get rid of initial k__
  # microbe_name <- gsub("[[:digit:]]+","",microbe_name,perl=T) ## get rid of number added at end by rowsum in summarize_taxa() to create unique groups. 
  ## These are anyways unique taxa, so we can get rid of them during visualization after downstream analysis. 
  return(microbe_name)
}

## clr + imputation function
## Borrowed from Gabe (email) 
## Matrix transforms to achieve CLR transform (http://www.sediment.uni-goettingen.de/staff/tolosana/extra/CoDaNutshell.pdf)
## Note, this assumes taxa as rows and samples as columns
clr_transform <- function(taxa){
  clr.taxa <- taxa
  clr.taxa = t(clr.taxa); eps = 0.5 ## changing eps alters the number of columns filterd
  clr.taxa = clr.taxa*(1 - rowSums(clr.taxa==0)*eps/rowSums(clr.taxa))
  clr.taxa[clr.taxa==0]=eps
  clr.taxa = sweep(clr.taxa,1,rowSums(clr.taxa),'/');
  ls = log(clr.taxa)
  clr.taxa = t(ls - rowMeans(ls))
  clr.taxa = clr.taxa[,!is.nan(colSums(clr.taxa))]
  return(clr.taxa)
}

## collapse at different taxa levels
summarize_taxa <- function(otu,num_samples,sep){
  lsTaxa <- list()
  for(rank in 2:7){
    
    # ## debug block start
    # num_samples = dim(otu)[2]
    # sep = ";"
    # rank = 5
    # ## debug block end
    
    ## Split taxa names by seprator 
    split = strsplit(as.character(rownames(otu)),sep) ## list of lists
    ## generate taxa name at a given taxa rank
    ## TODO: Either modify this or remove taxa with trailing NAs below. 
    taxaStrings = sapply(split,function(x) paste(x[1:rank],collapse=";"))
    # select <- grep("(;NA)+$",taxaStrings)
    # taxaStrings <- taxaStrings[-select]
    ## get rid of k__ at begining of string
    # taxaStrings <- gsub('([[:alpha:]])__','', taxaStrings)
    ## Remove any trailing NAs
    # taxaStrings = gsub("(;NA)+$","",taxaStrings,perl=T)# Clean tips
    ## Collapse table grouped by taxa names
    taxa = rowsum(otu, group =  taxaStrings) 
    dim(taxa)
    
    ## remove taxa with trailing NAs
    select <- grep("(;NA)+$",rownames(taxa))
    if(length(select)!= 0){
      taxa <- taxa[-select,]
    }
    
    ## Prevalence-based filtering
    ## Tune the prevalence filter to find out which one works for this data. 
    # select <- rowSums(taxa > 0) >= num_samples*0.5 #Keep otus which are present in more than 50% of the samples.
    # taxa <- taxa[select,]
    
    ## Filter further by only keeping taxa with rel. abundance atleast 0.1% (or 0.001) in atleast half of the samples
    taxa.rel <- sweep(taxa,2,colSums(taxa),'/')
    ## Note that (taxa.rel >= 0.001) gives a boolean matrix with each cell marked as true or false indicating whether
    ## for a given taxa and a sample, the value in taxa.rel >= 0.001. rowSums on this boolean matrix results in 
    ## counting how many samples for a given taxa (row) qualify under this. Next, this is compared with 
    ## a certain percentage of samples. 
    select <- rowSums(taxa.rel >= 0.001) >= num_samples*0.1 ## this was 0.2 percent in old code
    # select <- rowSums(taxa.rel >= 0.0001) >= num_samples*0.1 ## from new paper on predicting metabolite from microbiome from Huttenhower lab
    taxa <- taxa[select,]
    dim(taxa)
    
    lsTaxa[[rank]] <- taxa
  }
  return(lsTaxa)
}

keep_top_variable_genes <- function(genes_table){
  
  ## TODO: Add a check here to make sure the samples are rows.
  
  ## compute std.dev. of expression of genes
  genes.sd <- transform(as.data.frame(genes_table), SD=apply(as.data.frame(genes_table),1, sd, na.rm = TRUE))
  ## select top genes with high SD (variability across samples)
  # SD_cutoff <- 0.5
  SD_quantile <- quantile(genes.sd$SD)
  # SD_quantile ## identical to summary(genes.sd$SD)
  # # 0%       25%       50%       75%      100% 
  # # 0.1854572 0.4944215 0.6869490 0.9722974 4.0741752
  SD_cutoff <- SD_quantile[2] 
  genes.sd <- genes.sd[order(genes.sd$SD, decreasing = T),]
  top.variable.genes <- rownames(genes.sd[genes.sd$SD > SD_cutoff,])
  length(top.variable.genes) #12513
  ## Subset gene expr table for these genes
  genes_table <- genes_table[top.variable.genes,]
  return(genes_table)
  
}
################# Input metadata #####################

dir <- "./data/raw/CRC/metadata/"
filename <- paste0(dir,"sample.info.subset.txt")
sample.info <- data.frame(fread(filename, sep="\t", head=T), stringsAsFactors = F, check.names = F, row.names = 1)
dim(sample.info) #[1] 88 21

## order samples by patient id, so that tumor/normal of same patient are together.
sample.info <- sample.info[order(sample.info$Patient_Blind_ID),]
normal.samples <- rownames(sample.info[sample.info$Description=="normal",]) 
length(normal.samples) #[1] 44
tumor.samples <- rownames(sample.info[sample.info$Description=="tumor",] )
length(tumor.samples) #[1] 44
## make sure tumor and normal samples in same order of patient id
all(sample.info[tumor.samples,]$Patient_Blind_ID == sample.info[normal.samples,]$Patient_Blind_ID)
# [1] TRUE

################# Preprocess gene expression data ########################

## read the protein-coding genes table
dir <- "./data/raw/CRC/gene_expr/"
filename <- paste0(dir,"all_samples_protein_coding_subread_counts.txt")
protein.coding.read.count <- data.frame(fread(filename, sep="\t", head=T), stringsAsFactors = F, check.names = F, row.names = 1)
dim(protein.coding.read.count) #[1] 19264    88

## Filter low count genes to keep genes that are expressed in at least half the number of samples
select <- rowSums(protein.coding.read.count > 0) >= dim(protein.coding.read.count)[2]/2
length(which(select)) #[1] 16685
protein.coding.subset <- protein.coding.read.count[select,]
dim(protein.coding.subset) #[1] 16685    88

################# OLD, SKIP: Preprocess microbiome data ######################

## Note, the Preprocess_CRC_taxa_abnd.R uses a different metadata file (filtered_sample_info.txt) than 
## the metadata (sample.info.subset.txt) used for gene expression data above. 

## Read in OTU table
## Edit on Feb 18, 2020 
## Read the otu table from which contaminants have been removed.
otu<- read.table("./data/raw/CRC/taxa_abnd/OTU_table.txt", sep='\t', comment='', head=T, row.names=1)
dim(otu) 
#[1] 1724   89, last column is taxa names -- OLD
# [1] 1225   88 -- after removal of potential contaminants

# convert sample IDs to same format as found in gene expression (i.e. Tissue.RNA.DNA_Tube_ID)
# convert sample id from "Sample_xx" to "syy" format to match gene expression table
sampleID <- sample.info[match(colnames(otu[,-ncol(otu)]),sample.info$SampleID),]$Tissue.RNA.DNA_Tube_ID
colnames(otu) <- c(as.character(sampleID),"taxonomy")

## Make taxonomy as the rownames
otu.taxonomy <- as.character(otu$taxonomy)
otu <- otu[,-ncol(otu)]
rownames(otu) <- otu.taxonomy

## get rid of Archea taxa
length(grep("k__Archaea",rownames(otu)))
select <- grep("k__Archaea",rownames(otu))
otu <- otu[-select,]
dim(otu) #[1] 1685   88

## Summarize at different taxa levels and concat them together
## Approach: summarize -> concat -> filter duplicates -> split into disease and Healthy -> CLR transform disease and Healthy
taxa.levels <- summarize_taxa(otu,dim(otu)[2],"; ") ## summarize taxa at different levels
## combine all taxa levels
taxa.combined <- do.call(rbind,taxa.levels) ## concat all levels
dim(taxa.combined) 
#[1] 366  88
# [1] 335  88 -- Oct 22, 2019
which(duplicated(rownames(taxa.combined))) #0
sum(duplicated(taxa.combined)) #[1] 107, length(which(duplicated(taxa.combined)))
sort(rownames(taxa.combined[duplicated(taxa.combined),]))
# This is due to same feature being collapsed at different levels, 
## so the value of the feature remains the same while the name of taxa changes
## Of the duplicates, keep the one with most characterized levels, so pick from the last of the duplicates
taxa.combined <- taxa.combined[sort(rownames(taxa.combined)),]
taxa.combined <- taxa.combined[!duplicated(taxa.combined, fromLast = T),] ## genus is last level
sum(duplicated(taxa.combined)) #0, all dups removed
dim(taxa.combined) 
#[1] 247  88 -- OLD
# [1] 228  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.3
# [1] 320  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.2
# [1] 480  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.1

## Normalize taxa names
## Convert taxa names to form kingdom;phylum;...;species 
## so that it's easy to compare across different datasets. 
norm.taxanames <- as.vector(sapply(rownames(taxa.combined),normalize_taxaname_CRC))
rownames(taxa.combined) <- norm.taxanames
dim(taxa.combined) 
#[1] 228  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.3
# [1] 320  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.2

################# Preprocess microbiome data ######################

## Note, the Preprocess_CRC_taxa_abnd.R uses a different metadata file (filtered_sample_info.txt) than 
## the metadata (sample.info.subset.txt) used for gene expression data above. 

## Edit on Feb 18, 2020 
## Read the otu table from which contaminants have been removed.
# otu<- read.table("./data/raw/CRC/taxa_abnd/OTU_table.txt", sep='\t', comment='', head=T, row.names=1)
otu<- read.table("./data/raw/CRC/taxa_abnd/OTU_table_no_contam_V5.txt", sep='\t', comment='', head=T, row.names=1)
dim(otu) 
# [1] 1180   88 -- Feb 17, 2020
# [1] 1244   88 -- Feb 25, 2020
# [1] 1048   88 -- March 9, 2020

## SKIP this, already done with contamination removal 
## normalize taxa labels to trim-off stuff at uncharacterized levels 
# otu.labels.fixed <- as.vector(sapply(rownames(otu),normalize_taxaname_CRC))
# rownames(otu) <- otu.labels.fixed

## Summarize at different taxa levels and concat them together
## Approach: summarize -> concat -> filter duplicates -> split into disease and Healthy -> CLR transform disease and Healthy
taxa.levels <- summarize_taxa(otu,dim(otu)[2],";") ## summarize taxa at different levels
## combine all taxa levels
taxa.combined <- do.call(rbind,taxa.levels) ## concat all levels
dim(taxa.combined) 
# [1] 335  88 -- Oct 22, 2019
# [1] 575  88 -- Feb 2020
# [1] 291  88 -- March 2020
## Any duplicated taxa names
which(duplicated(rownames(taxa.combined))) #0
## How many duplicated rows
sum(duplicated(taxa.combined))
## list duplicated rows
sort(rownames(taxa.combined[duplicated(taxa.combined),]))
# This duplication is due to same feature being collapsed at different levels, 
## so the value of the feature remains the same while the name of taxa changes
## Of the duplicates, keep the one with most characterized levels, 
## so pick from the last of the duplicates

taxa.combined <- taxa.combined[sort(rownames(taxa.combined)),]
taxa.combined <- taxa.combined[!duplicated(taxa.combined, fromLast = T),] ## genus is last level
sum(duplicated(taxa.combined)) #0, all dups removed
dim(taxa.combined) 
#[1] 247  88 -- OLD
# [1] 228  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.3
# [1] 320  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.2
# [1] 480  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.1
# [1] 369  88 -- contam removal + rowSums(taxa.rel >= 0.001) >= num_samples*0.1
# [1] 307  88 -- updated summarize_taxa + contam removal + rowSums(taxa.rel >= 0.001) >= num_samples*0.1 
# [1] 186  88 -- updated summarize_taxa + contam removal + rowSums(taxa.rel >= 0.001) >= num_samples*0.2
# [1] 182  88 -- Feb 25, 2020 updated summarize_taxa + contam removal + rowSums(taxa.rel >= 0.001) >= num_samples*0.2
# [1] 263  88 -- Mar 9, 2020 (contam removal + rowSums(taxa.rel >= 0.001) >= num_samples*0.1)

## taxa.combined at 0.001 abundance in 10% of samples
# write.table(taxa.combined, file="./data/clean/CRC/taxa_abnd/CRC_taxa_combined_0.001_0.1.txt",sep="\t", row.names = T, quote = F, col.names = NA)
## taxa.combined at 0.001 abundance in 10% of samples
# write.table(taxa.combined, file="./data/clean/CRC/taxa_abnd/CRC_taxa_combined_0.001_0.1_V5.txt",sep="\t", row.names = T, quote = F, col.names = NA)


## SKIP, completed normalization of taxa names above. 
## Normalize taxa names
## Convert taxa names to form kingdom;phylum;...;species 
## so that it's easy to compare across different datasets. 
# norm.taxanames <- as.vector(sapply(rownames(taxa.combined),normalize_taxaname_CRC))
# rownames(taxa.combined) <- norm.taxanames
# dim(taxa.combined) 
#[1] 228  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.3
# [1] 320  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.2


################# Filter out low abundant contaminants from combined taxa maxtrix #############
## Mar 24, 2020
## Dropping a few low abundant taxa from the taxa combined matrix. This shouldn't alter
## the overall abundance of bugs, but we don't want these to lurk around in our associations. 

## List of taxa to filter out from combined matrix
# Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Nocardioidaceae
# Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;[Paraprevotellaceae];YRC22
# Bacteria;Fibrobacteres
# Bacteria;BRC1
# Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Comamonadaceae;Methylibium
# Bacteria;Proteobacteria;Gammaproteobacteria;Xanthomonadales;Sinobacteraceae
# Bacteria;WS*
# Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Frankiaceae
# Bacteria;Chloroflexi;Anaerolineae;SBR1031
# Bacteria;Proteobacteria;Deltaproteobacteria;PB19
# Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Candidatus Azobacteroides
# Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Paludibacter
# Bacteria;Chloroflexi;SAR202
# Special casees -- all members of Synergistales are contam except Synergistaceae
# e.g. contam -- Bacteria;Synergistetes;Synergistia;Synergistales;Dethiosulfovibrionaceae, 
# Bacteria;Synergistetes;Synergistia;Synergistales;Anaerobaculum, etc.
## So filter out all Synergistales except Bacteria;Synergistetes;Synergistia;Synergistales;Synergistaceae 

contam_phylum <- c("Bacteria;Fibrobacteres",
                   "Bacteria;BRC1",
                   "Bacteria;WS*" ## WS1,2,4,6
                   )
contam_class <- c("Bacteria;Chloroflexi;SAR202")
contam_order <- c("Bacteria;Chloroflexi;Anaerolineae;SBR1031",
                  "Bacteria;Proteobacteria;Deltaproteobacteria;PB19",
                  "Bacteria;Synergistetes;Synergistia;Synergistales;*",
                  "Bacteria;Proteobacteria;Deltaproteobacteria;NB1-j",
                  "Bacteria;Proteobacteria;Gammaproteobacteria;34P16",
                  "Bacteria;Proteobacteria;Betaproteobacteria;Methylophilales")
contam_family <- c("Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Nocardioidaceae",
                   "Bacteria;Proteobacteria;Gammaproteobacteria;Xanthomonadales;Sinobacteraceae",
                   "Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Frankiaceae",
                   "Bacteria;Firmicutes;Clostridia;Clostridiales;\\[Acidaminobacteraceae\\]",
                   "Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Pseudonocardiaceae")
contam_genus <- c("Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;\\[Paraprevotellaceae\\];YRC22",
                  "Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Comamonadaceae;Methylibium",
                  "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Candidatus Azobacteroides",
                  "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Paludibacter",
                  "Bacteria;Proteobacteria;Deltaproteobacteria;Desulfobacterales;Desulfobulbaceae;Desulfocapsa",
                  "Bacteria;Spirochaetes;Spirochaetes;Sphaerochaetales;Sphaerochaetaceae;Sphaerochaeta")
contam_species <- c("Bacteria;Firmicutes;Clostridia;Clostridiales;Peptococcaceae;Desulfosporosinus;meridiei",
                    "Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Comamonadaceae;Acidovorax;caeni",
                    "Bacteria;Proteobacteria;Betaproteobacteria;Methylophilales;Methylophilaceae;Methylotenera;mobilis")

contam_list <- list(contam_phylum,contam_class,contam_order,contam_family,contam_genus, contam_species)

## Get average relative abundance of these contaminants from the summarized OTU table. 
print("Average relabnd:")
for(l in 2:6){
  level <- l
  dim(taxa.levels[[level]]) #[1] 24 88
  taxa.relabnd <- sweep(taxa.levels[[level]],2,colSums(taxa.levels[[level]]),'/')
  taxa.relabnd$avg_relabnd <- rowMeans(taxa.relabnd)
  
  contam <- contam_list[[l-1]] 
  for(i in 1:length(contam)){
    select <- grep(contam[i], rownames(taxa.relabnd))
    if(length(select) != 0){
      print(paste0(contam[i],":",taxa.relabnd[select,]$avg_relabnd));flush.console()
    }
  }
}

## Updated Mar 30, 2020
# [1] "Bacteria;Fibrobacteres:0.000559144655811304"
# [1] "Bacteria;BRC1:0.00066031585243497"
# [1] "Bacteria;WS*:0.000321251845002229"
# [1] "Bacteria;Chloroflexi;SAR202:0.000645604319303085"
# [1] "Bacteria;Chloroflexi;Anaerolineae;SBR1031:0.00355628500209869"
# [1] "Bacteria;Proteobacteria;Deltaproteobacteria;PB19:0.000493155122121636"
# [1] "Bacteria;Synergistetes;Synergistia;Synergistales;*:0.00173503387448014"
# [1] "Bacteria;Proteobacteria;Deltaproteobacteria;NB1-j:0.000482666143610193"
# [1] "Bacteria;Proteobacteria;Gammaproteobacteria;34P16:0.000427661021471147"
# [1] "Bacteria;Proteobacteria;Betaproteobacteria;Methylophilales:0.000363530730499408"
# [1] "Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Nocardioidaceae:0.000926479718165835"
# [1] "Bacteria;Proteobacteria;Gammaproteobacteria;Xanthomonadales;Sinobacteraceae:0.000546832672275227"
# [1] "Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Frankiaceae:0.000430491089681425"
# [1] "Bacteria;Firmicutes;Clostridia;Clostridiales;\\[Acidaminobacteraceae\\]:0.000375032778842584"
# [1] "Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Pseudonocardiaceae:0.000941355630767476"
# [1] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;\\[Paraprevotellaceae\\];YRC22:0.00304865656401272"
# [1] "Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Comamonadaceae;Methylibium:0.0197468090421391"
# [1] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Candidatus Azobacteroides:0.00360341839571059"
# [1] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Paludibacter:0.00149060257463893"
# [1] "Bacteria;Proteobacteria;Deltaproteobacteria;Desulfobacterales;Desulfobulbaceae;Desulfocapsa:0.00150487780333599"
# [1] "Bacteria;Spirochaetes;Spirochaetes;Sphaerochaetales;Sphaerochaetaceae;Sphaerochaeta:0.00123960731741972"

## All, except Methylibium, have less than 0.005 or 0.5% relative abundance across samples. 
## Methylibium has ~2% rel abnd on average across samples. 

## Filter these contaminants out from the combined taxa matrix
contam_to_filter <- c("Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Nocardioidaceae",
                      "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;\\[Paraprevotellaceae\\];YRC22",
                      "Bacteria;Fibrobacteres",
                      "Bacteria;BRC1;*",
                      "Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Comamonadaceae;Methylibium",
                      "Bacteria;Proteobacteria;Gammaproteobacteria;Xanthomonadales;Sinobacteraceae",
                      "Bacteria;WS*", ## WS1,2,4,6
                      "Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Frankiaceae",
                      "Bacteria;Chloroflexi;Anaerolineae;SBR1031",
                      "Bacteria;Proteobacteria;Deltaproteobacteria;PB19",
                      "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Candidatus Azobacteroides",
                      "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Paludibacter",
                      "Bacteria;Chloroflexi;SAR202",
                      "Bacteria;Synergistetes;Synergistia;Synergistales;Dethiosulfovibrionaceae",
                      ## Added Mar 30
                      "Bacteria;Firmicutes;Clostridia;Clostridiales;\\[Acidaminobacteraceae\\]",
                      "Bacteria;Firmicutes;Clostridia;Clostridiales;Peptococcaceae;Desulfosporosinus;meridiei",
                      "Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Comamonadaceae;Acidovorax;caeni",
                      "Bacteria;Proteobacteria;Betaproteobacteria;Methylophilales;Methylophilaceae;Methylotenera;mobilis",
                      "Bacteria;Proteobacteria;Deltaproteobacteria;Desulfobacterales;Desulfobulbaceae;Desulfocapsa",
                      "Bacteria;Proteobacteria;Deltaproteobacteria;NB1-j",
                      "Bacteria;Proteobacteria;Gammaproteobacteria;34P16",
                      "Bacteria;Spirochaetes;Spirochaetes;Sphaerochaetales;Sphaerochaetaceae;Sphaerochaeta",
                      "Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Pseudonocardiaceae",
                      "Bacteria;Proteobacteria;Betaproteobacteria;Methylophilales"
                      )

taxa.combined.filt <- taxa.combined
dim(taxa.combined.filt) #[1] 263  88
for(i in 1:length(contam_to_filter)){
  # i <- 2
  select <- grep(contam_to_filter[i], rownames(taxa.combined.filt))
  if(length(select) != 0){
    taxa.combined.filt <- taxa.combined.filt[-select,]
    print(paste0("Removed ",contam_to_filter[i],":",dim(taxa.combined.filt)[1]));flush.console()
  }
}
dim(taxa.combined.filt) 
#[1] 246  88 -- Mar 24
#[1] 235  88 -- Mar 20 

## No. of taxa dropped: 
dim(taxa.combined)[1] - dim(taxa.combined.filt)[1]
# [1] 28

taxa.combined <- taxa.combined.filt
dim(taxa.combined)
# [1] 235  88

################# Overlap gene expression and microbiome data ##############
## TODO: check out the IBS preprocessing code. 
## Identify overlapping samples. In the case of CRC, we know we have the
## same set of samples. 
common.samples <- intersect(colnames(taxa.combined),colnames(protein.coding.subset))
length(common.samples) #88 -- So all samples are common between the two datasets. 

## Make sure the samples are lined up in both the datasets for integration 
all(colnames(taxa.combined) == colnames(protein.coding.subset)) ## F
taxa.combined <- taxa.combined[,common.samples]
protein.coding.subset <- protein.coding.subset[,common.samples]
all(colnames(taxa.combined) == colnames(protein.coding.subset)) ## T

## July 2020 -- output unnormalized taxa table for all samples to file.
dim(taxa.combined)
# [1] 235  88
## transform to keep the output format uniform 
taxa.combined.t <- as.data.frame(t(taxa.combined))
dim(taxa.combined.t)
# [1]  88 235

## Write unnormalized taxa abundance to a file
filename <- "crc_unnormalized_all_samples_taxa_t_no_contam_0.001_0.1_V3.txt" # adding trailing "no_contam_0.001_0.1_V3" to comply with latest files post contam removal.
filepath <- paste0("./data/clean/CRC/taxa_abnd/",filename)
# write.table(taxa.combined.t, file=filepath, sep="\t", row.names = T, col.names = NA )


################# DEseq2 normalization for tumor and normal samples (separately) ################

## If we subset vst-tranformed matrix, dispersion estimates may differ for tumor and normal samples.  
## Hence we need to do this separately 

## subset read counts into tumor and normal
protein.coding.tumor <- protein.coding.subset[,tumor.samples]; dim(protein.coding.tumor) #[1] 16685    44
protein.coding.normal <- protein.coding.subset[,normal.samples]; dim(protein.coding.normal) #[1] 16685    44

## normalize for tumor samples 
dds <- DESeqDataSetFromMatrix(countData = protein.coding.tumor,
                              colData = sample.info[tumor.samples,],
                              design =  ~ 1) ## important to keep if combining tumor and normal samples, 

## vst transform
vsd <- vst(dds, blind = TRUE) ## Check "Blind dispersion estimation" in DESeq2 vignette: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html  
protein.coding.tumor.vsd <- assay(vsd)

## Ensure that no genes got variance stabilized to constant value across all samples
genes.sd <- transform(as.data.frame(protein.coding.tumor), SD=apply(as.data.frame(protein.coding.tumor),1, sd, na.rm = TRUE))
any(genes.sd$SD == 0) #[1] FALSE
## If any genes have SD = 0, check raw read count for them

## Filter by variablity of genes
tumor.genes <- protein.coding.tumor.vsd ## to pass to function
dim(tumor.genes) #[1] 16685    44
tumor.genes.filt <- keep_top_variable_genes(tumor.genes)
dim(tumor.genes.filt) #[1] 12513    44

## transform to make rows as samples and columns as features before writing to file
tumor.genes.filt.t <- t(tumor.genes.filt)
dim(tumor.genes.filt.t) # 44 12513

## write to file
# write.table(tumor.genes.filt.t, file="./data/clean/CRC/gene_expr/tumor_protein_coding_vsd_SDfilt_t.txt",sep="\t", row.names = T, col.names = NA )

## normalize for normal samples
dds <- DESeqDataSetFromMatrix(countData = protein.coding.normal,
                              colData = sample.info[normal.samples,],
                              design =  ~ 1) ## important to keep if combining tumor and normal samples, 

## vst transform
vsd <- vst(dds, blind = TRUE) ## Check "Blind dispersion estimation" in DESeq2 vignette: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html  
# head(assay(vsd), 3)
protein.coding.normal.vsd <- assay(vsd)

## Ensure that no genes got variance stabilized to constant value across all samples
genes.sd <- transform(as.data.frame(protein.coding.normal.vsd), SD=apply(as.data.frame(protein.coding.normal.vsd),1, sd, na.rm = TRUE))
any(genes.sd$SD == 0) #[1] FALSE
## If any genes have SD = 0, check raw read count for them

## Check overlaps between genes for both tumor and normal matrices
all(rownames(protein.coding.normal.vsd) == rownames(protein.coding.tumor.vsd)) #[1] TRUE
length(intersect(rownames(protein.coding.normal.vsd),rownames(protein.coding.tumor.vsd))) 

## Make sure that I didn't make tumor and normal vsd matrices identical by mistake.
identical(protein.coding.normal.vsd, protein.coding.tumor.vsd) #[1] FALSE


## Filter by variablity of genes
normal.genes <- protein.coding.normal.vsd ## to pass to function
dim(normal.genes) #[1] 16685    44
normal.genes.filt <- keep_top_variable_genes(normal.genes)
dim(normal.genes.filt) #[1] 12513    44 ## Interesting that same #genes remain after filtering by top variable genes 
                                        ## both for tumor and normal. 

## transform to make rows as samples and columns as features before writing to file
normal.genes.filt.t <- t(normal.genes.filt)
dim(normal.genes.filt.t) # 44 12513

## write to files
# write.table(normal.genes.filt.t, file="./data/clean/CRC/gene_expr/normal_protein_coding_vsd_SDfilt_t.txt",sep="\t", row.names = T, col.names = NA )

################# SKIP: Deseq2 normalization for all samples (tumor+normal) combined ###########
## Create dds dataset from read count matrix
## Note, normalization is not affected by the design formula; see Mike Love's response: https://support.bioconductor.org/p/106548/

## Make sure protein.coding.subset columns and sample.info rows are in same order before calling deseq2
protein.coding.subset <- protein.coding.subset[,rownames(sample.info)]
sample.info$Patient_Blind_ID <- as.factor(sample.info$Patient_Blind_ID)

dds <- DESeqDataSetFromMatrix(countData = protein.coding.subset,
                              colData = sample.info,
                              design =  ~ as.factor(Patient_Blind_ID) + Description) ## important to keep if combining tumor and normal samples, 
## Check "Blind dispersion estimation" in DESeq2 vignette: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html  

## vst transform
vsd <- vst(dds, blind = FALSE) ## Check "Blind dispersion estimation" in DESeq2 vignette: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html  
protein.coding.vsd <- assay(vsd)

## Ensure that no genes got variance stabilized to constant value across all samples
genes.sd <- transform(as.data.frame(protein.coding.vsd), SD=apply(as.data.frame(protein.coding.vsd),1, sd, na.rm = TRUE))
any(genes.sd$SD == 0) #[1] FALSE
## If any genes have SD = 0, check raw read count for them

## Filter by variablity of genes
all.sample.genes <- protein.coding.vsd ## to pass to function
dim(all.sample.genes) #[1] 16685    88
all.sample.genes.filt <- keep_top_variable_genes(all.sample.genes)
dim(all.sample.genes.filt) #[1] 12513    88

## Make sure samples in genes table line up with samples in taxa table
all.sample.genes.filt <- all.sample.genes.filt[,colnames(taxa.combined)]
all(colnames(all.sample.genes.filt) == colnames(taxa.combined)) #T

## transform to make rows as samples and columns as features before writing to file
all.sample.genes.filt.t <- t(all.sample.genes.filt)
dim(all.sample.genes.filt.t) # 88 12513


## Write to file
# write.table(all.sample.genes.filt.t, file="./data/clean/CRC/gene_expr/all_samples_protein_coding_vsd_SDfilt_t.txt",sep="\t", row.names = T, col.names = NA )

################# Normalize microbiome data -- CLR transform ##################
## 1. CLR transform the taxa matrix for all samples
all.taxa.clr <- as.data.frame(clr_transform(taxa.combined)); dim(all.taxa.clr) 
#[1] 247  88
# [1] 228  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.3
# [1] 320  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.2
# [1] 263  88 -- rowSums(taxa.rel >= 0.001) >= num_samples*0.1
# [1] 246  88 -- Mar 2020, rowSums(taxa.rel >= 0.001) >= num_samples*0.1 + contam dropped from taxa matrix
# [1] 237  88 -- Mar 30, 2020, additional contam dropped from combined taxa matrix

## 2. CLR transform tumor and normal separately 
## Split into tumor and normal samples
tumor.taxa <- taxa.combined[,tumor.samples];dim(tumor.taxa) 
normal.taxa <- taxa.combined[,normal.samples]; dim(normal.taxa) 

tumor.taxa.clr <- as.data.frame(clr_transform(tumor.taxa))
normal.taxa.clr <- as.data.frame(clr_transform(normal.taxa))

## transform to make samples as rows and taxa as columns
all.taxa.clr.t <- as.data.frame(t(all.taxa.clr));dim(all.taxa.clr.t)
tumor.taxa.clr.t <- as.data.frame(t(tumor.taxa.clr)); dim(tumor.taxa.clr.t) 
normal.taxa.clr.t <- as.data.frame(t(normal.taxa.clr)); dim(normal.taxa.clr.t) 
## Write to files: all.taxa.clr, tumor.taxa.clr, normal.taxa.clr
# write.table(all.taxa.clr.t, file="./data/clean/CRC/taxa_abnd/all_samples_clr_t_no_contam_0.001_0.1_V3.txt",sep="\t", row.names = T, col.names = NA )
# write.table(tumor.taxa.clr.t, file="./data/clean/CRC/taxa_abnd/tumor_taxa_clr_t_no_contam_0.001_0.1_V3.txt",sep="\t", row.names = T, col.names = NA )
# write.table(normal.taxa.clr.t, file="./data/clean/CRC/taxa_abnd/normal_taxa_clr_t_no_contam_0.001_0.1_V3.txt",sep="\t", row.names = T, col.names = NA )

################# SKIP FOR NOW: Edit: July 31, 2019: Write sample diagnosis to file for sparse sCCA##########
sample_condition_df <- data.frame( sampleID = as.character(sample.info$Tissue.RNA.DNA_Tube_ID), condition = sample.info$Description )
table(sample_condition_df$condition)
# normal  tumor 
# 44     44 
rownames(sample_condition_df) <- sample_condition_df$sampleID
## Make sure sample IDs aligned for sample_condition_df and colnames of taxa table
all(rownames(sample_condition_df) == colnames(all.taxa.clr)) #[1] TRUE
sample_condition_df$sampleID <- NULL

## write to file
# write.table(sample_condition_df, file="./data/clean/CRC/taxa_abnd/sample_type.txt",sep="\t", row.names = T, col.names = NA )


################# Create covariate matrix ################

## Create a subset of mapping file with study_id and relevant covariates to include with predictor matrix
covariates <- as.data.frame(cbind(sampleID = rownames(sample.info),
                                  gender = as.vector(sample.info$Sex), 
                                  age = as.vector(sample.info$Age),
                                  condition = as.vector(sample.info$Description)),
                            stringsAsFactors = FALSE)

## TODO: Add CRC_subtype by using stage of CRC.
covariates$gender <- ifelse(covariates$gender == "male", 0, 1) ## male=0
covariates$condition <- ifelse(covariates$condition == "normal", 0, 1) ## normal=0
covariates$age <- as.numeric(covariates$age)
## Let's check what is the range of age in our cohort: 
summary(sample.info$Age) ## note duplicates as paired samples
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.00   56.00   67.00   64.91   76.00   91.00 
median(sample.info$Age)
# [1] 67
sum(unique(sample.info$Age) < 30)
# [1] 2 -- one at 17 yrs, another at 26. 

## since age is a continuous variable, let's normalize age on the scale of 0 to 1.
x <- covariates$age
x_norm <- (x-min(x))/(max(x)-min(x))
covariates$age <- x_norm
covariates <- covariates[!duplicated(covariates),];dim(covariates) #[1] 88  4
rownames(covariates) <- covariates$sampleID

################# Add covariates separately to tumor and normal tables #################
## covariates to add: gender,age
## Oct 28: Ran said last week to drop age from covariates, since age at collection is not available for IBD
## And we want to keep the models consistent. 

## tumor
tumor.taxa.clr.t$gender <- covariates[rownames(tumor.taxa.clr.t),]$gender
head(tumor.taxa.clr.t$gender)
# tumor.taxa.clr.t$age <-  covariates[rownames(tumor.taxa.clr.t),]$age
# head(tumor.taxa.clr.t$age)
dim(tumor.taxa.clr.t) 
#[1]  44 229 -- 228 taxa + 1 covariates
# [1]  44 321 -- 320 taxa + 1 covariate
# [1]  44 238 -- 237 taxa + 1 covariate (Feb 2020)
# [1]  44 264 -- 264 taxa + 1 covariate (Mar 2020)
# [1]  44 247 -- 246 taxa + 1 covariate (Mar 24, 2020) post filtering of additional contam from taxa matrix

##normal
normal.taxa.clr.t$gender <- covariates[rownames(normal.taxa.clr.t),]$gender
# normal.taxa.clr.t$age <-  covariates[rownames(normal.taxa.clr.t),]$age
dim(normal.taxa.clr.t) 
#[1]  44 229 -- 228taxa + 1 covariates
# [1]  44 321 -- 320 taxa + 1 covariate
# [1]  44 238 -- 237 taxa + 1 covariate (Feb 2020)
# [1]  44 264 -- 264 taxa + 1 covariate (Mar 2020)
# [1]  44 247 -- 246 taxa + 1 covariate (Mar 24, 2020) post filtering of additional contam from taxa matrix
 

## double check that sample IDs (rownames) match between gene expr and microbiome data
all(rownames(normal.genes.filt.t) == rownames(normal.taxa.clr.t) ) #T
all(rownames(tumor.genes.filt.t) == rownames(tumor.taxa.clr.t)) #T

write.table(tumor.taxa.clr.t, file="./data/clean/CRC/taxa_abnd/tumor_taxa_clr_w_gender_no_contam_0.001_0.1_V3.txt",sep="\t", row.names = T, col.names = NA )
write.table(normal.taxa.clr.t, file="./data/clean/CRC/taxa_abnd/normal_taxa_clr_w_gender_no_contam_0.001_0.1_V3.txt",sep="\t", row.names = T, col.names = NA )

################# SKIP: Add covariates to all samples (tumor+normal) combined ###############
## covariates to add to all samples: gender, age, sample_type

## Add covariates to the all samples table
all.taxa.clr.t$gender <- covariates[rownames(all.taxa.clr.t),]$gender
# all.taxa.clr.t$age <-  covariates[rownames(all.taxa.clr.t),]$age
all.taxa.clr.t$condition <- covariates[rownames(all.taxa.clr.t),]$condition
dim(all.taxa.clr.t) #[1]  88 231 -- 228 taxa + 3 covariates

## double chack that sample IDs (rownames) are lined up between gene expr and microbiome table
all(rownames(all.sample.genes.filt.t) == rownames(all.taxa.clr.t))#T

# write.table(all.taxa.clr.t, file="./data/clean/CRC/taxa_abnd/all_samples_clr_w_covariates_0.001_0.1.txt",sep="\t", row.names = T, col.names = NA )

################# Split genes list into 5 subsets #################### 
## Are tumor genes different than normal genes? What about genes in all_samples?
all(colnames(tumor.genes.filt.t) == colnames(normal.genes.filt.t)) #F
length(intersect(colnames(tumor.genes.filt.t),colnames(normal.genes.filt.t))) #10548
length(intersect(colnames(tumor.genes.filt.t),colnames(all.sample.genes.filt.t))) #11705
## Yep, these sets are different. Need to have condition specific subset.

##tumor genes list
filename <- "./data/clean/CRC/gene_expr/tumor_protein_coding_vsd_SDfilt_t.txt"
genes <- data.frame(fread(filename, sep="\t", head=T), stringsAsFactors = F, check.names = F, row.names = 1)
genes <- t(genes)
dim(genes)
chunk.size <- dim(genes)[1]/5
chunk.size 
genes.split <- split(rownames(genes), ceiling(seq_along(rownames(genes))/chunk.size))
length(genes.split) #5
for(split in 1:length(genes.split)){
  write(unlist(genes.split[split]), file=paste0("./data/clean/CRC/gene_expr/tumor_genes_split_dir_5nodes/genes_split_",split,".txt"), sep = "\n")
}

##normal genes list
filename <- "./data/clean/CRC/gene_expr/normal_protein_coding_vsd_SDfilt_t.txt"
genes <- data.frame(fread(filename, sep="\t", head=T), stringsAsFactors = F, check.names = F, row.names = 1)
genes <- t(genes)
dim(genes)
chunk.size <- dim(genes)[1]/5
chunk.size 
genes.split <- split(rownames(genes), ceiling(seq_along(rownames(genes))/chunk.size))
length(genes.split) #5
for(split in 1:length(genes.split)){
  write(unlist(genes.split[split]), file=paste0("./data/clean/CRC/gene_expr/normal_genes_split_dir_5nodes/genes_split_",split,".txt"), sep = "\n")
}

## tumor + normal genes list
filename <- "./data/clean/CRC/gene_expr/all_samples_protein_coding_vsd_SDfilt_t.txt"
genes <- data.frame(fread(filename, sep="\t", head=T), stringsAsFactors = F, check.names = F, row.names = 1)
genes <- t(genes)
dim(genes)
chunk.size <- dim(genes)[1]/5
chunk.size 
genes.split <- split(rownames(genes), ceiling(seq_along(rownames(genes))/chunk.size))
length(genes.split) #5
for(split in 1:length(genes.split)){
  write(unlist(genes.split[split]), file=paste0("./data/clean/CRC/gene_expr/all_samples_genes_split_dir_5nodes/genes_split_",split,".txt"), sep = "\n")
}
