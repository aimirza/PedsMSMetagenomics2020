---
title: "Essential Workflow"
author: "Ali Mirza"
output: 
  html_document: 
    keep_md: yes
    toc: yes
---

The following document contains the essential workflow used to processes and analyze the metagenomic sequences from the study. The purpose of this document is to show which parameters were used for each software and make our analysis more understandable and reproducable. While the parmaters listed in the document are identitical to what was used in the study, other non-essential codes used in the study are not included in the workflow, such as submitting a job to a cluster, command loops, ploting, creating or moving files and directories. 

# Pre-processing Workflow (Linux Bash)
This is the workflow used to generate the data for the metagenomics Pediatric multiple sclerosis study. 

## Quality control with kneaddata 

Create kneaddata (version 0.7.4) environment via conda  
```
$ conda create --name kneaddata
$ source activate kneaddata
$ conda install -c cyclus java-jre
$ conda install -c bioconda fastqc
$ conda install -c bioconda trf
$ conda install -c bioconda kneaddata
$ conda install pip
$ pip install kneaddata

$ mkdir main # this is where all the files of each sample will be placed
$ mkdir final # This is where the concatenated files will be placed

```

Install Trimmomatic version 0.33. Not via conda  
```
$ curl -o Trimmomatic-0.33.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
$ unzip Trimmomatic-0.33.zip

#keep note of the install location. Will need to provide its path to keanddata
```

Install kneaddata human contamination databases.  
Everything within and including the greater/less than symsbols "< >" should be replaced.  
```
$ kneaddata_database --download human_genome bowtie2 <path/to/install/kneaddata_database>
```

Run Kneaddata with removing tandom repeats  
```
$ kneaddata --input </path/forward_sample> --input </path/reverse_sample> -db <path/to/install/kneaddata_database> --output main --output-prefix <prefix>  --run-trf --threads 8  --verbose  --trimmomatic <path/to/Trimmomatic-0.33> --cat-final-output

#Move concatinated files to the directory "final"
$ mv main/<prefix>.fastq final/
```

## Generate metabolic pathway and gene families abundance profiles using humann3
Install the new version of humann (version 3.0 alpha) via conda  
```
$ conda create --name biobakery3 python=3.7
$ source activate biobakery3
$ conda install humann -c biobakery
$ git clone https://github.com/biobakery/MetaPhlAn.git
$ conda install metaphlan=3.0=pyh5ca1d4c_2 --no-channel-priority
```

Install required databases  
We will installed the full UniRef90 database      
```
#pangenome database:
$ humann_databases --download chocophlan full <path/to/install/humann_database> --update-config yes

#protein database:
$ humann_databases --download uniref uniref90_diamond <path/to/install/humann_database> --update-config yes

3annotations database:
$ humann_databases --download utility_mapping full <path/to/install/humann_database>  --update-config yes
```

Run humann3
```
$ humann --input final/<prefix>.fastq --output main --threads 8
```

Join tables of each sample into one file.
Move joined files into the directory "final"  
```
$ humann_join_tables --input main --output humann3_genefamilies.tsv --file_name genefamilies && mv humann3_genefamilies.tsv final
$ humann_join_tables --input main --output humann3_pathcoverage.tsv --file_name pathcoverage && mv humann3_pathcoverage.tsv final
$ humann_join_tables --input main --output humann3_pathabundance.tsv --file_name pathabundance && mv humann3_pathabundance.tsv final
```

Split a tables into two files (one stratified and one unstratified)  
```
$ humann_split_stratified_table --input final/humann3_pathabundance.tsv --output final/
$ humann_split_stratified_table --input final/humann3_genefamilies.tsv --output final/
$ humann_split_stratified_table --input final/humann3_pathcoverage.tsv --output final/
```

Regroup genes info into level 4 enzyme comisions (EC), gene ontologies (GO), and KEGG orthologs (KO)
```
$ humann_regroup_table --input final/humann3_genefamilies_unstratified.tsv --output final/humann3_unstratified_level4ec.tsv --groups uniref90_level4ec

$ humann_regroup_table --input final/humann3_genefamilies_unstratified.tsv --output final/humann3_unstratified_uniref90_go.tsv --groups uniref90_go

$ humann_regroup_table --input final/humann3_genefamilies_unstratified.tsv --output final/humann3_unstratified_uniref90_ko.tsv --groups uniref90_ko
```


## Estimate microbial tax abundances using Kraken2 and Bracken

You may ask why not use the Biobakery program "metaphlan" to calculate taxonomy abundances. Our statistical methods (log ratio transformation via ALDex2) requires that relative abundances are in counts, not fractions. Metaphlan by default outputs relative abndances as fractions not counts. Unfortunately, the only way to get counts from metaphlan is to get the estimated number of reads comming from each clade, not from each species.  Thus we cannot use metaphlan using our statistical methods. Not a problem! There are other great tools like Kraken abd Bracken.  


Install Kraken2 and Bracken via conda and activate the environment. I had kraken2 vsersion 2.0.9 beta and bracken version 2.6 installed. Download the full kraken2 reference database.  
```
$ kraken2-build --standard --threads 32 --db <path/to/kraken2_reference_database/> --use-ftp
```

For each sample, combine the forward reads with pairs and orphan ("unmatched") reads into one file. Repeat for reverse reads.  
```
$ cat <prefix>.repeats.removed.1.fastq <prefix>.repeats.removed.unmatched.1.fastq
$ cat <prefix>.repeats.removed.2.fastq <prefix>.removed.unmatched.2.fastq`
```

Run Kraken2 and combine each report
```
$ kraken2 --db </path/to/kraken2_reference_database/> --threads 16  --paired final/<prefix>.repeats.removed.all.1.fastq final/<prefix>.repeats.removed.all.2.fastq  --output <prefix>.repeats.removed.all.classified --use-names --report <prefix>.repeats.removed.all.classified.report

$ python ~/kraken2_scripts/KrakenTools/combine_kreports.py -o kraken2_combined_report -r `ls *report | tr "\n" " " `
```

Run Bracken
```
#Re-estimate species abundance
$ bracken -d </path/to/kraken2_reference_database/> -i <prefix>.repeats.removed.all.classified.report -o <prefix>.S.bracken -r 150 -l 'S' -t 2

#Re-estimate genus abundance
$ bracken -d </path/to/kraken2_reference_database/> -i <prefix>.repeats.removed.all.classified.report -o <prefix>.G.bracken -r 150 -l 'G' -t 2
```

# Metagenomic analysis (R programming)


Load packages  
```
> library(dplyr)
> library(ALDEx2) # version 1.21.1
> library(phyloseq)
> library(ggplot2)
```

## Functional diversity analysis 
Calculate diversity of enzyme commisions (EC) and KEGG ortholog (KO) proteins

EC diversity
```
#Upload EC abundance profiles generated from humann3
> level4ec <- read.csv("/humann3_unstratified_level4ec.tsv", header=TRUE, check.names = FALSE, fill = TRUE, row.names = 1, stringsAsFactors = FALSE, sep="\t", comment.char = "")

# clean names
> names(level4ec) <- gsub(x = names(level4ec), pattern = "_Abundance-RPKs", replacement = "")

# Remove 1st two rows. 
> level4ec <- level4ec[-c(1,2),]

# Rename sample names
> row.names(metadata) <- metadata$sample.id

# Save data into phyloseq object
> ps_ec = phyloseq(sample_data(metadata), otu_table(level4ec, taxa_are_rows=TRUE))

> set.seed(12345)

# Rarefy the table
> ps_ec_r <- rarefy_even_depth(ps_ec)

#caculate and plot alpha-diversity 
> p_ec_r <- plot_richness(ps_ec_r, measures=c("Observed", "InvSimpson"), x="diagnosis_category")
> p_ec_r$layers[1] <- NULL
> p_ec_r + theme_bw() + labs(x = "")  + geom_boxplot(lwd=1) + 
  theme(axis.text=element_text(size=14,face="bold"))
```

Run Wilxocon rank sum test on each diversity measure 
```
# Observed EC
> wilcox.test(p_ec_r$data %>% filter(variable=="Observed") %>% filter(sample.id %in% MS) %>% pull(value), p_ec_r$data %>% filter(variable=="Observed") %>% filter(sample.id %in% CONTROL) %>% pull(value))
#W = 212.5, p-value = 0.7454

# InvSimpson EC
> wilcox.test(p_ec_r$data %>% filter(variable=="InvSimpson") %>% filter(sample.id %in% MS) %>% pull(value), p_ec_r$data %>% filter(variable=="InvSimpson") %>% filter(sample.id %in% CONTROL) %>% pull(value))
#W = 183, p-value = 0.6588
```

Repeat for DMD status
```
> p_ec_r_dmd <- plot_richness(ps_ec_r, measures=c("Observed", "InvSimpson"), x="MS_DMD_ever_pre_stool_sample3")
> p_ec_r_dmd$layers[1] <- NULL
> p_ec_r_dmd + theme_bw() + labs(x = "") + geom_boxplot(lwd=1) + 
  theme(axis.text=element_text(size=14,face="bold"))
```

Repeat for kEGG Orthologs (proteins)
```
#Upload KEGG protein data
> ko <- read.csv("humann3_unstratified_uniref90_ko.tsv", header=TRUE, check.names = FALSE, fill = TRUE, row.names = 1, stringsAsFactors = FALSE, sep="\t", comment.char = "")

> names(ko) <- gsub(x = names(ko), pattern = "_Abundance-RPKs", replacement = "")

> ko <- ko[-c(1,2),]

> ps_ko = phyloseq(sample_data(metadata), otu_table(ko, taxa_are_rows=TRUE))

> ps_ko_r <- rarefy_even_depth(ps_ko)
#178OTUs were removed because they had zero abundance of rarefying 

#plot alpha-diversity 
> p_ko_r <- plot_richness(ps_ko_r, measures=c("Observed", "InvSimpson"), x="diagnosis_category")
> p_ko_r$layers[1] <- NULL
> p_ko_r + theme_bw() + labs(x = "") + geom_boxplot(lwd=1) + 
  theme(axis.text=element_text(size=14,face="bold"))

> p_ko_r$data %>% filter(variable=="Observed") %>% pull(value) %>% density() %>%  plot()
> p_ko_r$data %>% filter(variable=="InvSimpson") %>% pull(value) %>% density() %>%  plot()
```
Run Wilcoxon rank sum tests
```
# observed KO
> wilcox.test(p_ko_r$data %>% filter(variable=="Observed") %>% filter(sample.id %in% MS) %>% pull(value), p_ko_r$data %>% filter(variable=="Observed") %>% filter(sample.id %in% CONTROL) %>% pull(value))
#W = 239, p-value = 0.3013

#InvSImpson KO
> wilcox.test(p_ko_r$data %>% filter(variable=="InvSimpson") %>% filter(sample.id %in% MS) %>% pull(value), p_ko_r$data %>% filter(variable=="InvSimpson") %>% filter(sample.id %in% CONTROL) %>% pull(value))
#W = 207, p-value = 0.862

> p_ko_r_dmd <- plot_richness(ps_ko_r, measures=c("Observed", "InvSimpson"), x="MS_DMD_ever_pre_stool_sample3")
> p_ko_r_dmd$layers[1] <- NULL
> p_ko_r_dmd + theme_bw() + labs(x = "") + geom_boxplot(lwd=1) + 
  theme(axis.text=element_text(size=14,face="bold"))
#W = 239, p-value = 0.3013
```
Create function to run Wilcoxon rank sum test between specified groups
```
> functional_diversity <- function(data, g1, g2, diversity){
   ## data: data genreated from phyloseq::plot_richness
   ## g1: 1st group
   ## g2: 2nd group
   ## diversity: diversity measure
   
  wilcox.test(data %>% filter(variable==diversity) %>% filter(sample.id %in% g1) %>% pull(value), data %>% filter(variable==diversity) %>% filter(sample.id %in% g2) %>% pull(value))
}
```

Test significance by DMD status
```
#KO DMD observed
> functional_diversity(data=p_ko_r_dmd$data, g1=Ever, g2=CONTROL, diversity="Observed")
# W = 139, p-value = 0.4772

> functional_diversity(data=p_ko_r_dmd$data, g1=Never, g2=CONTROL, diversity="Observed")
#W = 100, p-value = 0.3285

> functional_diversity(data=p_ko_r_dmd$data, g1=Never, g2=Ever, diversity="Observed")
#W = 53, p-value = 0.7345

# All p > 0.33

#KO DMD InvSimpson

> functional_diversity(data=p_ko_r_dmd$data, g1=Ever, g2=CONTROL, diversity="InvSimpson")
# W = 143, p-value = 0.3867

> functional_diversity(data=p_ko_r_dmd$data, g1=Never, g2=CONTROL, diversity="InvSimpson")
#W = 64, p-value = 0.4385

> functional_diversity(data=p_ko_r_dmd$data, g1=Never, g2=Ever, diversity="InvSimpson")
#W = 34, p-value = 0.3054

# all p > 0.31


#EC DMD observed
> functional_diversity(data=p_ec_r_dmd$data, g1=Ever, g2=CONTROL, diversity="Observed")
# W = 135.5, p-value = 0.5592

> functional_diversity(data=p_ec_r_dmd$data, g1=Never, g2=CONTROL, diversity="Observed")
# W = 77, p-value = 0.8988

> functional_diversity(data=p_ec_r_dmd$data, g1=Never, g2=Ever, diversity="Observed")
# W = 49, p-value = 0.9692

# All p > 0.56

#EC DMD observed

> functional_diversity(data=p_ec_r_dmd$data, g1=Ever, g2=CONTROL, diversity="InvSimpson")
# W = 95, p-value = 0.3456

> functional_diversity(data=p_ec_r_dmd$data, g1=Never, g2=CONTROL, diversity="InvSimpson")
# W = 88, p-value = 0.7086

> functional_diversity(data=p_ec_r_dmd$data, g1=Never, g2=Ever, diversity="InvSimpson")
# W = 60, p-value = 0.3837

# All p > 0.35

```
## Pathway analysis
```
#Load pathway abundances unstratified
> pathabundance_unstratified <- read.csv("humann3_pathabundance_unstratified.tsv", header=TRUE, check.names = FALSE, fill = TRUE, row.names = 1, stringsAsFactors = FALSE, sep="\t", comment.char = "")


> names(pathabundance_unstratified) <- gsub(x = names(pathabundance_unstratified), pattern = "_Abundance", replacement = "")
#row names contain symbol and full name of each pathway, seprated by ": ". Convert rownames to symbol but save a table with both the symbol and name to use in the future. 

> pathabundance_unstratified_names <- reshape2::colsplit(row.names(pathabundance_unstratified), ": ", c("symbol","name"))
> row.names(pathabundance_unstratified) <- pathabundance_unstratified_names$symbol

#Replace NAs with 0
> pathabundance_unstratified[is.na(pathabundance_unstratified)] <- 0
```

Some pathways are identical.  Want to remove to reduce the multiple comparisons problem.  
After statistical testing, you can check if any of the significant findings match with the removed pathways.  
```
> table.uniq <- unique(pathabundance_unstratified)  
```

Upload pathway abundance stratified
```
> pathabundance_stratified <- read.csv("humann3_pathabundance_stratified.tsv", header=TRUE, check.names = FALSE, fill = TRUE, row.names = 1, stringsAsFactors = FALSE, sep="\t", comment.char = "")

#clean up names
> names(pathabundance_stratified) <- gsub(x = names(pathabundance_stratified), pattern = "_Abundance", replacement = "") # clean names
```

Exclude pathways with coverage (present) in less than 5% of samples
```
> keep <- which(rowSums(table.uniq > 0) > 0.5*40)
> table.uniq.filtered <- table.uniq[keep,]
```

### Pathway differential analysis using ALDEx2

Based on the simulation, the metacyc pathway HEXITOLDEGSUPER-PWY (superpathway of hexitol degradation (bacteria)) had the smallest clr variance. Find the row number of this pathway to index when using ALDEx2. When giving ALDEx2 a single index number as the denominator (paramater "denom") ALDex2 will effectively perform additive-log ratio (alr) transformation.
```
> ref <- which(rownames(table.uniq.filtered)=="HEXITOLDEGSUPER-PWY")
```

Alr transform data and generate Monte Carlo samples of the Dirichlet distribution for each sample.  
ALDEx2 only accepts integers. We will have to convert from float to integer. We used the function ceiling() to do that.  
```
#We chose to create 1000 Monte Carlo samples for best accuracy. This will take a while to run. 
> x.table.uniq.filtered <- aldex.clr(ceiling(table.uniq.filtered), conds=as.vector(factor(metadata$diagnosis_category, levels=c("MS","CONTROL"))), mc.samples=1000, denom=as.integer(ref), verbose=TRUE)
```

Calculate the expected values of the Wilcoxon Rank Sum test on the data returned by aldex.clr.
```
> table.uniq.filtered.ttest <- aldex.ttest(x.table.uniq.filtered) # this function also performs the Welch's t-test
```
aldex.ttest calculates outputs the students t-test and WIlcoxon Rank Sum tests results. We are only interested in the later.  

Repeat for assessing differences by disease-modyifying drug (DMD) exposure status. Filter pathways to include only those initially significantly different by disease status (MS vs controls). The file metadata contains a column named "MS_DMD_ever_pre_stool_sample3" which includes DMD exposure status (ever, never, and control).  
```
> sig_pathways <- row.names(table.uniq.filtered.ttest)[which(table.uniq.filtered.ttest$wi.ep < 0.05)]

> table.uniq.filtered.sig <- table.uniq.filtered[c(sig_pathways, "HEXITOLDEGSUPER-PWY"),] # Include the reference pathway

#DMD ever exposed VS naive (never) exposed MS cases

> keep <- which(metadata$MS_DMD_ever_pre_stool_sample3 %in% c("Ever", "Never"))  

> MS_DMD_ever_vs_never <- as.vector(metadata$MS_DMD_ever_pre_stool_sample3[keep])

ref <- which(rownames(table.uniq.filtered.sig)=="HEXITOLDEGSUPER-PWY")

> x.table.uniq.filtered.dmd.ever_vs_never <- aldex.clr(ceiling(table.uniq.filtered.sig[,keep]), 
                                   conds=MS_DMD_ever_vs_never,
                                   mc.samples=1000, 
                                   denom=as.integer(ref), verbose=TRUE)


> x.table.uniq.filtered.dmd.ever_vs_never.ttest <- aldex.ttest(x.table.uniq.filtered.dmd.ever_vs_never)


#EVER vs CONTROL

> keep <- which(metadata$MS_DMD_ever_pre_stool_sample3 %in% c("Ever", "Control"))

> MS_DMD_ever_vs_control <- as.vector(metadata$MS_DMD_ever_pre_stool_sample3[keep])

> ref <- which(rownames(table.uniq.filtered.sig)=="HEXITOLDEGSUPER-PWY")

> x.table.uniq.filtered.dmd.ever_vs_control <- aldex.clr(ceiling(table.uniq.filtered.sig[,keep]), 
                                   conds=MS_DMD_ever_vs_control,
                                   mc.samples=1000, 
                                   denom=as.integer(ref), verbose=TRUE)


> x.table.uniq.filtered.dmd.ever_vs_control.ttest <- aldex.ttest(x.table.uniq.filtered.dmd.ever_vs_control)


#NEVER vs CONTROL

> keep <- which(metadata$MS_DMD_ever_pre_stool_sample3 %in% c("Never", "Control"))

> MS_DMD_never_vs_control <- as.vector(metadata$MS_DMD_ever_pre_stool_sample3[keep])

> ref <- which(rownames(table.uniq.filtered.sig)=="HEXITOLDEGSUPER-PWY")

> x.table.uniq.filtered.dmd.never_vs_control <- aldex.clr(ceiling(table.uniq.filtered.sig[,keep]), 
                                   conds=MS_DMD_never_vs_control,
                                   mc.samples=1000, 
                                   denom=as.integer(ref), verbose=TRUE)


> x.table.uniq.filtered.dmd.never_vs_control.ttest <- aldex.ttest(x.table.uniq.filtered.dmd.never_vs_control)
```

### Pathway differential proportional analysis using Fisher Exact Test
Create a function "FtestAll" that converts the abundance table to a contingency table and outputs the odds ratio, p value, and odds ratio 95% confidense intervals.  The contingency table is created by classifying the each row (pathway) in the relative abundace table into two variables: disease status and prevelance (presence/absence) of each pathway.  
```
#First save MS cases and controls stool sample IDs from metadata

> MS <- metadata %>% filter(diagnosis_category=="MS") %>% pull(sample.id)
> CONTROL <- metadata %>% filter(diagnosis_category=="CONTROL") %>% pull(sample.id)

# Create function
> FtestAll <- function(table.uniq.filtered, MS=MS, CONTROL=CONTROL){
   
   ## table.uniq.filtered : Any abundance table with rows as features and columns as samples
   ## MS : the first compared group. Default MS cases
   ## MS : the second compared group, Default controls
   
    fisher_names <- names(table.uniq.filtered)
    
    #Convert abundance table into binary variables: "Yes" if abundance > 0 and "No" if abundance = 0
    table.uniq.filtered.fisher <- table.uniq.filtered
    table.uniq.filtered.fisher[table.uniq.filtered.fisher > 0]  <- "yes"
    table.uniq.filtered.fisher[table.uniq.filtered.fisher == 0]  <- "no"
    
    Ftest <- lapply(1:nrow(table.uniq.filtered.fisher), function(i){
        MS_yes <-      sum( table.uniq.filtered.fisher[i,MS]=="yes" )         # number of MS cases with non zero counts in each pathway
        MS_no <-       sum( table.uniq.filtered.fisher[i,MS]=="no" )          # number of MS cases with zero counts in each pathway
        CONTROL_yes <- sum( table.uniq.filtered.fisher[i,CONTROL]=="yes" )    # number of controls with non zero counts in each pathway
        CONTROL_no <-  sum( table.uniq.filtered.fisher[i,CONTROL]=="no" )     # number of controls with zero counts in each pathway
        
        contingency_table <- data.frame(Present=c(MS_yes, CONTROL_yes),
                                        Not_present=c(MS_no,CONTROL_no),
                                        stringsAsFactors = FALSE)
        
        # Compute Fisher exact test
        Ftest <- fisher.test(contingency_table)
        
        odds_ratio <- Ftest$estimate %>% as.numeric()
        p.value <- Ftest$p.value
        conf.int.95 <- paste(round(Ftest$conf.int[1],2),round(Ftest$conf.int[2], 2), sep="-")
        
        cbind(odds_ratio=odds_ratio, p.value=p.value, conf.int.95=conf.int.95, ci.low=Ftest$conf.int[1], ci.high=Ftest$conf.int[2])
    })
    
    Ftest.all <- do.call("rbind", Ftest)
    rownames(Ftest.all) <- rownames(table.uniq.filtered.fisher)
    Ftest.all <- as.data.frame( cbind(Ftest.all, coverage=rowSums(table.uniq.filtered > 0)) )
}    
```

Run function to get the Fisher test results  
```
> Ftest.all = FtestAll(table.uniq.filtered)
```

Get coverage of each pathway and only include pathways without 100% coverage (not present in all samples). It makes no sense to perform a Fisher exact test when a pathway is present in 100% of samples. Then adjust P values using Benjamini-Hochberg.  
```
> MS_yes <- rowSums( table.uniq.filtered[,MS] > 0 ) 
> CONTROL_yes <- rowSums( table.uniq.filtered[,CONTROL] > 0 ) 

> Ftest.all$MS_yes <- MS_yes
> Ftest.all$MS_yes_prop <- MS_yes / length(MS)

> Ftest.all$CONTROL_yes <- CONTROL_yes
> Ftest.all$CONTROL_yes_prop <- CONTROL_yes / length(CONTROL)

> Ftest.all <- Ftest.all %>% filter(as.integer(coverage) < 40)

# Add pathway descriptions
> Ftest.all$description <- pathabundance_unstratified_names[rownames(Ftest.all),]$name

# Adjust P values
> p.adjust(as.numeric(Ftest.all$p.value), method = "BH", n = length(as.numeric(Ftest.all$p.value)))
```

Repeat for assessing differences by disease-modyifying drug (DMD) exposure status. Filter pathways to include only those initially significantly different by disease status (MS vs controls). Below is how we campared DMD ever exposed MS cases VS DMD naive MS cases

```
> Ever <- metadata$sample.id[which(metadata$MS_DMD_ever_pre_stool_sample3=="Ever")]   # DMD exposed MS cases
> Never <- metadata$sample.id[which(metadata$MS_DMD_ever_pre_stool_sample3=="Never")] # DMD naive MS cases

# Filter pathways to only include those previously significant
> dataPwy.fisher <- Ftest.all[which(Ftest.all$p.value < 0.05),] # Get significant pathways

> table.uniq.filtered.sig <- table.uniq.filtered[row.names(dataPwy.fisher), ]

# Create contingency table and run Fisher Exact test while changing the compared groups
> Ftest.all.dmd.never_ever = FtestAll(table.uniq.filtered.sig, MS=Ever, CONTROL=Never)

# Calculate coverages of each pathway by DMD exposure group
> never_yes <- rowSums( table.uniq.filtered.sig[,never] > 0 ) 
> ever_yes <- rowSums( table.uniq.filtered.sig[,ever] > 0 ) 

> Ftest.all.dmd.never_ever$never_yes <- never_yes
> Ftest.all.dmd.never_ever$never_yes_prop <- never_yes / length(never)

> Ftest.all.dmd.never_ever$ever_yes <- ever_yes
> Ftest.all.dmd.never_ever$ever_yes_prop <- ever_yes / length(ever)

# Remove any pathways with 100% coverage
> Ftest.all.dmd.never_ever <- Ftest.all.dmd.never_ever %>% filter(as.integer(coverage) < 40)

# Addd pathway description names
> Ftest.all.dmd.never_ever$description <- pathabundance_unstratified_names[rownames(Ftest.all.dmd.never_ever),]$name

# Adjust P values
> p.adjust(as.numeric(Ftest.all.dmd.never_ever$p.value), method = "BH", n = length(as.numeric(Ftest.all.dmd.never_ever$p.value)))
```



## Gene Ontology Analysis

Upload gene ontology (GO) data and get GO category info
```
# Upload GO abundance table
> GO <- read.csv("humann3_unstratified_uniref90_go.tsv", header=TRUE, check.names = FALSE, fill = TRUE, row.names = 1, stringsAsFactors = FALSE, sep="\t", comment.char = "")

# Clean sample names
> names(GO) <- gsub(x = names(GO), pattern = "_Abundance-RPKs", replacement = "")

# Upload humann3 GO map names. Can obtain from here: https://github.com/biobakery/humann/tree/master/humann/data/misc
> map_go_name <- read.csv("map_go_name.txt", header=FALSE, check.names = FALSE, fill = TRUE, row.names = 1, stringsAsFactors = FALSE, sep="\t", comment.char = "")

#Get all Biological Processes GOs
> map_go_name_BP <- grep("[BP]", map_go_name$V2, fixed = TRUE)
> map_go_name_BP <- rownames(map_go_name)[map_go_name_BP]

#Get all Molecular Functions GOs
> map_go_name_MF <- grep("[MF]", map_go_name$V2, fixed = TRUE)
> map_go_name_MF <- rownames(map_go_name)[map_go_name_MF]

#Get all Cellular Components GOs
> map_go_name_CC <- grep("[CC]", map_go_name$V2, fixed = TRUE)
> map_go_name_CC <- rownames(map_go_name)[map_go_name_CC]

#split GO table by category
> GO_BP <- GO[rownames(GO) %in% map_go_name_BP,]
> GO_MF <- GO[rownames(GO) %in% map_go_name_MF,]
> GO_CC <- GO[rownames(GO) %in% map_go_name_CC,]
```

### GO differential analysis using ALDEx2

Starting with Biological Processes GO  

Find GO with the lowest variance after centered-log ratio (clr) transforming the abundances  
```
> GO_BP_ubiq = apply(GO_BP, 1, function(row) all(row !=0 )) # Keep only non-zero rows

> GO_BP_in_all_samples <- row.names(GO_BP[GO_BP_ubiq,])  

# Geometric mean function for centered-log ratio
> g.mean.log_average <- function(x){ 
  x.ln <- log(x)
  x.ln[!is.finite(x.ln)] <- 0 # replace Inf with 0
  exp( mean(x.ln) )
}

#clr transform by column (sample)
> GO_BP.counts.clr <- apply(GO_BP[GO_BP_ubiq,], 2, function(x){ log(x/g.mean.log_average(x)) })
> GO_BP.counts.clr[!is.finite(GO_BP.counts.clr)] <- 0                                              # Take care of infinites
> GO_BP.counts.clr.var <- apply(GO_BP.counts.clr, 1, var)                                          # Calculate variance for each row (each GO)

> which(GO_BP.counts.clr.var==min(GO_BP.counts.clr.var))
# GO:0006094 

> map_go_name["GO:0006094",]
"[BP] gluconeogenesis"
```
Gluconeogenesis is the GO term with the lowest clr variance

RUN ALDEx2
```
# Remove GOs equal or less than 5% coverage
> GO_BP_coverage <- rowSums(GO_BP > 0)
> keep <- which(GO_BP_coverage > 40*0.05)
> GO_BP.filtered <- GO_BP[keep,]

# Find the row index of Gluconeogenesis
> BP_index <- which(rownames(GO_BP.filtered)=="GO:0006094")

#Generate Monte Carlo samples of the Dirichlet distribution for each sample. 
# Convert each instance using the centred log-ratio transform This is the input for all further analyses.
> x.GO_BP.filtered <- aldex.clr(ceiling(GO_BP.filtered), conds=as.vector(factor(metadata$diagnosis_category, levels=c("MS","CONTROL"))), mc.samples=1000, denom=as.integer(BP_index), verbose=TRUE)

# Run WIlcoxon test
> x.GO_BP.filtered.ttest <- aldex.ttest(x.GO_BP.filtered)
```

### GO differential proportional analysis using Fisher Exact Test

```
# Convert abundance table to contingency table and perform Fisher exact test using previously created function
> GO_BP.filtered.fisher <- FtestAll(GO_BP.filtered)

# Combine GO descriptions, fisher stats, and coverage info
> GO_BP.filtered.fisher.all <- cbind(map_go_name_BP, GO_BP.filtered.fisher, GO_BP_coverage)

# Calculate proportions of each GO by disease status
> MS_yes <- rowSums( GO_BP.filtered[,MS] > 0 ) 
> CONTROL_yes <- rowSums( GO_BP.filtered[,CONTROL] > 0 ) 

> GO_BP.filtered.fisher.all$MS_yes <- MS_yes
> GO_BP.filtered.fisher.all$MS_yes_prop <- MS_yes / 20

> GO_BP.filtered.fisher.all$CONTROL_yes <- CONTROL_yes
> GO_BP.filtered.fisher.all$CONTROL_yes_prop <- CONTROL_yes / 20

# Remove any GOs with 100% coverage
> GO_BP.filtered.fisher.all <- GO_BP.filtered.fisher.all %>% filter(as.integer(coverage) < 40)

# Adjust P values
> p.adjust(as.numeric(sort.int(GO_BP.filtered.fisher.all$p.value)), method = "BH", n = length(GO_BP.filtered.fisher.all$p.value))
```
Repeat for Molecular Functions GO and Cellular Components GO.  


### GO differential analysis by drug exposure
Will only show codes for Biological processes GO

MS cases ever exposed to DMDs VS  MS cases naive for DMDs
```
# Remove terms equal or less than 5% coverage
> GO_BP_coverage <- rowSums(GO_BP > 0)
> keep <- which(GO_BP_coverage > 40*0.05)
> GO_BP.filtered <- GO_BP[keep,]


# filter GOs to include only previously significant GOs while including the reference GO (gluconeogenesis)
> GO_BP.sig <- row.names(x.GO_BP.filtered.ttest)[which(x.GO_BP.filtered.ttest$wi.ep < 0.05)]

> GO_BP.filtered.sig <- GO_BP.filtered[c(GO_BP.sig,"GO:0006094"),]

# Find row index of gluconeogenesis (GO:0006094)
> BP_index <- which(rownames(GO_BP.filtered.sig)=="GO:0006094")

# Run ALDEx2 clr transformation
> x.GO_BP.filtered.sig.ever_never <- aldex.clr(ceiling(GO_BP.filtered.sig[,metadata$sample.id[ever_never]]), conds=MS_DMD_ever_vs_never, mc.samples=1000, denom=as.integer(BP_index), verbose=TRUE)

> x.GO_BP.filtered.sig.ever_never.ttest <- aldex.ttest(x.GO_BP.filtered.sig.ever_never)

# Add GO description
> map_go_name_BP.filtered <- map_go_name_BP[keep]
> map_go_name_BP.filtered.sig <- map_go_name_BP.filtered[which(x.GO_BP.filtered.ttest$wi.ep < 0.05)]
```

Repeat for: 
(1) MS cases ever exposed to DMDs VS controls 
(2) MS cases naive for DMDs VS controls


### GO differential proportional analysis by drug exposure

MS cases ever exposed to DMDs VS controls 
```
# Get the symbol and description of GOs previously significant by disease status
> GO_BP.filtered.fisher.sig.GO_symbol <- row.names(GO_BP.filtered.fisher.all)[which(GO_BP.filtered.fisher.all$p.value < 0.05)]
> GO_BP.filtered.fisher.sig.GO_name <- GO_BP.filtered.fisher.all$map_go_name_BP[which(GO_BP.filtered.fisher.all$p.value < 0.05)]

# Create contingency tables and rUn fisher test using previously created function
> GO_BP.filtered.fisher.sig.Ever_CONTROL <- FtestAll(GO_BP.filtered[GO_BP.filtered.fisher.sig,], MS=Ever, CONTROL=CONTROL)

# Combine GO descriptions and Fisher results
> GO_BP.filtered.fisher.sig.Ever_CONTROL.all <- cbind(GO_name=GO_BP.filtered.fisher.sig.GO_name, GO_BP.filtered.fisher.sig.Ever_CONTROL)

# Calculate coverage for each GO by group
> Ever_yes <- rowSums( GO_BP.filtered[GO_BP.filtered.fisher.sig.GO_symbol, Ever] > 0 ) 
> CONTROL_yes <- rowSums( GO_BP.filtered[GO_BP.filtered.fisher.sig.GO_symbol, CONTROL] > 0 ) 

> GO_BP.filtered.fisher.sig.Ever_CONTROL.all$Ever_yes <- Ever_yes
> GO_BP.filtered.fisher.sig.Ever_CONTROL.all$Ever_yes_prop <- Ever_yes / length(Ever)

> GO_BP.filtered.fisher.sig.Ever_CONTROL.all$CONTROL_yes <- CONTROL_yes
> GO_BP.filtered.fisher.sig.Ever_CONTROL.all$CONTROL_yes_prop <- CONTROL_yes / length(CONTROL)

# Adjust P values
> p.adjust(as.numeric(sort.int(GO_BP.filtered.fisher.sig.Ever_CONTROL.all$p.value)), method = "BH", n = length(GO_BP.filtered.fisher.sig.Ever_CONTROL.all$p.value))
```

Repeat for: 
(1) MS cases ever exposed to DMDs VS MS cases naive for DMDs 
(2) MS cases naive for DMDs VS controls



## Taxonomy analysis

Upload genus relative abundance profiles
```
#identify bracken files to read in R
> filesToProcess <- list.files(path = "bracken/contigs/", pattern = "*.G.bracken$", recursive=TRUE, full.names=TRUE)

#Iterate over each of those file names with lapply
> listOfFiles <- lapply(filesToProcess, function(x) read.table(x, header=TRUE, check.names = FALSE, fill = FALSE, row.names = NULL, stringsAsFactors = FALSE, sep="\t"))

#Select columns of choice
> listOfFiles_abund <- lapply(listOfFiles, function(z) z[c(1,2,3,6)])

> names(listOfFiles_abund[[1]])
#[1] "name"          "taxonomy_id"   "taxonomy_lvl"  "new_est_reads.CHP-00002-PT-06_S11"


#Merge all of the objects in the list together with Reduce. 
> bracken.G.contigs <- Reduce(function(x,y) {full_join(x,y, by = c("name","taxonomy_id","taxonomy_lvl"))}, listOfFiles_abund)

#clean up the column names
> names(bracken.G.contigs) <- gsub(x = names(bracken.G.contigs), pattern = "new_est_reads.", replacement = "")

#Replace NAs with 0s
> bracken.G.contigs[is.na(bracken.G.contigs)] <- 0

#Move taxonomy ID to row names
> bracken.G.contigs <- tibble::column_to_rownames(bracken.G.contigs, var="name")
> bracken.G.contigs.info <- bracken.G.contigs[,1:2]
> bracken.G.contigs <- bracken.G.contigs[,-c(1:2)]

#Original number of ganera
> bracken.G.contigs %>% dim()
# 1757   40
```

Repeat for species level
```
#identify bracken files to read in R
> filesToProcess <- list.files(path = "bracken/contigs/", pattern = "*.S.bracken$", recursive=TRUE, full.names=TRUE)

#Iterate over each of those file names with lapply
> listOfFiles <- lapply(filesToProcess, function(x) read.table(x, header=TRUE, check.names = FALSE, fill = FALSE, row.names = NULL, stringsAsFactors = FALSE, sep="\t"))

#Select columns of choice
> listOfFiles_abund <- lapply(listOfFiles, function(z) z[c(1,2,3,6)])

> names(listOfFiles_abund[[1]])
#[1] "name"          "taxonomy_id"   "taxonomy_lvl"  "new_est_reads.CHP-00002-PT-06_S11"


#Merge all of the objects in the list together with Reduce. 
> bracken.S.contigs <- Reduce(function(x,y) {full_join(x,y, by = c("name","taxonomy_id","taxonomy_lvl"))}, listOfFiles_abund)

#clean up the column names
> names(bracken.S.contigs) <- gsub(x = names(bracken.S.contigs), pattern = "new_est_reads.", replacement = "")

#Replace NAs with 0s
> bracken.S.contigs[is.na(bracken.S.contigs)] <- 0

> bracken.S.contigs <- tibble::column_to_rownames(bracken.S.contigs, var="name")
> bracken.S.contigs.info <- bracken.S.contigs[,1:2]
> bracken.S.contigs <- bracken.S.contigs[,-c(1:2)]

#Original number of species
> bracken.S.contigs %>% dim()
# 6492  40
```

Get list of bactera, archaea, and virus taxnonomy
```
#Upload kraken2 abundance file
> kraken2 <- read.csv("kraken2_combined_report-table.txt", header=TRUE, check.names = FALSE, fill = TRUE, row.names = 86, stringsAsFactors = FALSE, sep="\t")

#Upload the mapping file generated from kraken2
> kraken2_mapped_samples <- read.csv("kraken2_combined_report-sample-map.txt", header=FALSE, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE, sep="\t")

> kraken2_taxa <- row.names(kraken2)

> kraken2_names <- grep("all", names(kraken2), value=TRUE)[-1]

> row.names(kraken2_mapped_samples) <- kraken2_names

> kraken2_taxa <- trimws(kraken2_taxa, which = "left") # remove leading white space

> row.names(kraken2) <- kraken2_taxa


#Find index of each domain
> which(kraken2$lvl_type=="D")
#     4 10342 10980 11010
> row.names(kraken2)[which(kraken2$lvl_type=="D")]
# "    Bacteria"  "    Archaea"   "    Eukaryota" "  Viruses"

> bacteria_index <- 4:10341
> Archaea_index <- 10342:10979
> Eukaryota_index <- 10980:11009
> Viruses_index <- 11010:11624

> bacteria_taxa <- row.names(kraken2)[bacteria_index]
> Archaea_taxa <- row.names(kraken2)[Archaea_index]
> Viruses_taxa <- row.names(kraken2)[Viruses_index]


```


### Taxa differential analysis with Aldex2

Genus only  

Find bacteria with the lowest variance after clr transformation
```
# remove leading white space from bacterial taxa names
> bacteria_taxa.trimws <- trimws(bacteria_taxa, which = "left") 

#filter to include only bacteria ganera
> bracken.G.contigs.bacteria <- bracken.G.contigs[row.names(bracken.G.contigs) %in%  bacteria_taxa.trimws,]

# clr transform each column
> bracken.G.contigs.bacteria.clr <- apply(bracken.G.contigs.bacteria, 2, clr)

# Remove rows with 0 counts
> bracken.G.contigs.bacteria.clr.noZero <- bracken.G.contigs.bacteria.clr[row.names(bracken.G.contigs.bacteria.clr) %in% row.names(bracken.G.contigs.bacteria)[rowSums(bracken.G.contigs.bacteria==0)==0],]

# Calculate variance for each taxa
> bracken.G.contigs.bacteria.clr.noZero.var <- apply(bracken.G.contigs.bacteria.clr.noZero, 1, var)

# Find the name of the taxa with the lowest variance
> ref_genus <- names(bracken.G.contigs.bacteria.clr.noZero.var)[which(bracken.G.contigs.bacteria.clr.noZero.var==min(bracken.G.contigs.bacteria.clr.noZero.var))]
#Neisseria 
```

Run aldex2 on all genera including virus and archaea
```
# filter taxa with coverage less than 5% coverage across samples
> keep <- which(rowSums(bracken.G.contigs > 0) > 0.05*40)
> bracken.G.contigs.filtered <-  bracken.G.contigs[keep,]

# Find the index of the taxa with the lowest variance
> inti <- which(row.names(bracken.G.contigs.filtered)==ref_genus)

# Run aldex2
> x.bracken.G.contigs <- aldex.clr(ceiling(bracken.G.contigs.filtered[,metadata$sample.id]), conds=as.vector(factor(metadata$diagnosis_category, levels=c("MS","CONTROL"))), mc.samples=1000, denom=as.integer(inti), verbose=TRUE)

> x.bracken.G.contigs.ttest <- aldex.ttest(x.bracken.G.contigs)
```

Repeat for species
```
> bracken.S.contigs.bacteria <- bracken.S.contigs[row.names(bracken.S.contigs) %in%  bacteria_taxa.trimws,]

> bracken.S.contigs.bacteria.clr <- apply(bracken.S.contigs.bacteria, 2, clr)

> bracken.S.contigs.bacteria.clr.noZero <- bracken.S.contigs.bacteria.clr[row.names(bracken.S.contigs.bacteria.clr) %in% row.names(bracken.S.contigs.bacteria)[rowSums(bracken.S.contigs.bacteria==0)==0],]

> bracken.S.contigs.bacteria.clr.noZero.var <- apply(bracken.S.contigs.bacteria.clr.noZero, 1, var)

> ref_species <- names(bracken.S.contigs.bacteria.clr.noZero.var)[which(bracken.S.contigs.bacteria.clr.noZero.var==min(bracken.S.contigs.bacteria.clr.noZero.var))]
#"Desulfosporosinus youngiae"

> keep <- which(rowSums(bracken.S.contigs > 0) > 0.05*40)
> bracken.S.contigs.filtered <-  bracken.S.contigs[keep,]

> inti <- which(row.names(bracken.S.contigs.filtered)==ref_species)

> x.bracken.S.contigs <- aldex.clr(ceiling(bracken.S.contigs.filtered[,metadata$sample.id]), conds=as.vector(factor(metadata$diagnosis_category, levels=c("MS","CONTROL"))), mc.samples=1000, denom=as.integer(inti), verbose=TRUE)

> x.bracken.S.contigs.ttest <- aldex.ttest(x.bracken.S.contigs)

``` 

Repeat for DMD exposure status

DMD ever vs naive exposed MS cases
```
#get only DMD ever and naive exposed MS cases sample IDs
> ever_never <- which(metadata$MS_DMD_ever_pre_stool_sample3 %in% c("Ever", "Never"))
> MS_DMD_ever_vs_never <- as.vector(metadata$MS_DMD_ever_pre_stool_sample3[ever_never])

# Get reference species name
> ref_species <- names(bracken.S.contigs.bacteria.clr.noZero.var)[which(bracken.S.contigs.bacteria.clr.noZero.var==min(bracken.S.contigs.bacteria.clr.noZero.var))]
#"Desulfosporosinus youngiae"

# Filter species to include only those significantly different between MS cases and controls
> sig <- row.names(x.bracken.S.contigs.Never_Control.ttest)[which(x.bracken.S.contigs.Never_Control.ttest$wi.ep < 0.05)]
> bracken.S.contigs.filtered.sig <- bracken.S.contigs.filtered[sig,]


> inti <- which(row.names(bracken.S.contigs.filtered.sig)==ref_species)

> x.bracken.S.contigs.ever_never <- aldex.clr(ceiling(bracken.S.contigs.filtered.sig[,metadata$sample.id[ever_never] ]), conds=MS_DMD_ever_vs_never, mc.samples=1000, denom=as.integer(inti), verbose=TRUE)

> x.bracken.S.contigs.ever_never.ttest <- aldex.ttest(x.bracken.S.contigs.ever_never)
```

Repeat for all other comparisons and for genera

### Taxa differential proportional analysis

Genus
```
# Filter taxa with coverage less than 5% across samples and with 100% coverage
> keep <- which(rowSums(bracken.G.contigs > 0) < 40 & rowSums(bracken.G.contigs > 0) > 40*0.05)
> bracken.G.contigs.fisher_filter <- bracken.G.contigs[keep,]

# Create contingency table and run Fisher exact test using previously created function
> bracken.G.contigs.fisher_filter.fisher <- FtestAll(table.uniq.filtered=bracken.G.contigs.fisher_filter, MS=MS, CONTROL=CONTROL)

# Include coverage of each taxa by group
> MS_yes <- rowSums( bracken.G.contigs.fisher_filter[, MS] > 0 ) 
> CONTROL_yes <- rowSums( bracken.G.contigs.fisher_filter[, CONTROL] > 0 ) 

> bracken.G.contigs.fisher_filter.fisher$MS_yes <- MS_yes
> bracken.G.contigs.fisher_filter.fisher$MS_yes_prop <- MS_yes / length(MS)

> bracken.G.contigs.fisher_filter.fisher$CONTROL_yes <- CONTROL_yes
> bracken.G.contigs.fisher_filter.fisher$CONTROL_yes_prop <- CONTROL_yes / length(CONTROL)

# Adjust P values
> p.adjust(as.numeric(sort.int(bracken.G.contigs.fisher_filter.fisher$p.value)), method = "BH", n = length(bracken.G.contigs.fisher_filter.fisher$p.value))
```

Repeat for species


Repeat for DMD exposure status  
Naive MS case vs control 
```
> bracken.G.contigs.fisher_filter.fisher.Never_Control <- FtestAll(table.uniq.filtered=bracken.G.contigs.fisher_filter[,Never_Control], MS=Never, CONTROL=CONTROL)

> Never_yes <- rowSums( bracken.G.contigs.fisher_filter[, Never] > 0 ) 
> CONTROL_yes <- rowSums( bracken.G.contigs.fisher_filter[, CONTROL] > 0 ) 

> bracken.G.contigs.fisher_filter.fisher.Never_Control$Never_yes <- Never_yes
> bracken.G.contigs.fisher_filter.fisher.Never_Control$Never_yes_prop <- Never_yes / length(Never)

> bracken.G.contigs.fisher_filter.fisher.Never_Control$CONTROL_yes <- CONTROL_yes
> bracken.G.contigs.fisher_filter.fisher.Never_Control$CONTROL_yes_prop <- CONTROL_yes / length(CONTROL)


> p.adjust(as.numeric(sort.int(bracken.G.contigs.fisher_filter.fisher.Never_Control$p.value)), method = "BH", n = length(bracken.G.contigs.fisher_filter.fisher.Never_Control$p.value))
```

Repeat for all other comparisons
