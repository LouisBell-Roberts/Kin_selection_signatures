################Calculates pi and Tajima’s D################
#####Ooceraea biroi
####Standardises pi to account for gene length
###Loops through each of the contigs (14 chromosomes in this case) and saves the outputs as individual csv files
##Louis Bell-Roberts
#18/12/2023

# Load required packages
library(PopGenome)
library(tidyverse)

##########Sections that need to be manually adjusted at the start###########
#Assign species under study
SDP <-  "Ooceraea_biroi"

# Set the vector of contigs based on folder names
## List only immediate subdirectories (excluding subdirectories within those)
immediate_subdirectories <- list.dirs(file.path("/Volumes/Pop_Gen/PopGenomeWD", SDP), full.names = FALSE, recursive = FALSE)
## Filter subdirectories based on the specified criteria
contigs <- grep("^N.*\\.1$", immediate_subdirectories, value = TRUE, ignore.case = TRUE)

######
######Also the section later on that creates a vector based on file names?????
######

################################################################################


#####

#Identify the genomic regions (position) associated with each coding sequence (CDS) and the locus that that CDS is part of
##Create a csv file for the locus and associated genomic region that can be read back in later and assigned to the polymorphism statistics based on CDS region

#Create directory for the locus region identifier dataframes
dir.create(file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Locus_identifiers"))

for (contig in contigs) {
  tryCatch({
    #Read lines of the gtf file for each contig
    gtf_content <- readLines(file.path("/Volumes/Pop_Gen/PopGenomeWD", SDP,contig,"vcf_files/gff",paste0(contig,".gtf")))
    
    # Initialize vectors to store gene locus and genomic region information
    gene_locus <- character()
    gene_location <- character()
    
    # Process each line in the GTF file
    for (line in gtf_content) {
      # Skip comment lines
      if (startsWith(line, "#")) next
      
      # Split the line into fields
      fields <- strsplit(line, "\t")[[1]]
      
      # Extract gene locus and genomic region information
      gene_id_match <- regmatches(fields[9], regexec("gene_id \"(\\S+)\";", fields[9]))
  
      if (!is.na(gene_id_match[[1]][2])) {
        gene_locus <- c(gene_locus, gene_id_match[[1]][2])
        gene_location <- c(gene_location, paste(fields[4], fields[5], sep = " - "))
      }
    }
    
    # Create a dataframe
    gene_df <- data.frame(Gene_Locus = gene_locus, GeneLocation = gene_location)
  
    #Remove all duplicate rows, while retaining the first occurrence of each duplicate
    gene_df_dup_remove <- unique(gene_df)
    
    #Write the dataframe to csv
    write.csv(gene_df_dup_remove, file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Locus_identifiers",paste0(contig,"_locus_identifier.csv")), row.names = F)
    
    cat("csv for", contig, "completed.\n")
  }, error = function(e) {
    cat("Error in analysis for", contig, ": ", conditionMessage(e), "\n")
  })
}


#####

#Calculate polymorphism statistics and write to csv

#Create directory for the polymorphism csv files
dir.create(file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Polymorphism_contigs"))

# Iterate over contigs
##If the for loop encounters an error in one of the contigs (e.g. the non-chromosomal contigs), it moves onto the next contig in the sequence
for (contig in contigs) {
  tryCatch({
    # Set the working directory
    setwd(file.path("/Volumes/Pop_Gen/PopGenomeWD", SDP, contig))
    
    # Load vcf file and gff file into R
    GENOME.class <- readData("vcf_files/vcf", format = "VCF", gffpath = "vcf_files/gff")
    
    # Load the sequence of the reference genome
    GENOME.class <- set.synnonsyn(GENOME.class, ref.chr = paste0(contig, ".fna"))
    
    # Find the position as to where each gene starts and ends
    gene.positions <- get_gff_info(gff.file = paste0("vcf_files/gff/", contig, ".gtf"), chr = contig, feature = "CDS")
    
    # Split the VCF file into individual genes
    genes <- split_data_into_GFF_features(GENOME.class, paste0("vcf_files/gff/", contig, ".gtf"), contig, "CDS")
    
    # Preliminary calculations for the FST module
    genes <- F_ST.stats(genes, detail = TRUE)
    
    # Present diversity measures as a data frame
    df2 <- as.data.frame(get.diversity(genes)[[1]])
    
    # Select the column containing the polymorphism values
    pisite <- df2[, 3]
    
    #calculate pi for synonymous sites
    genes <- F_ST.stats(genes,subsites="syn", detail=TRUE)
    df3=as.data.frame(get.diversity(genes)[[1]])
    head(df3)
    pisitesyn=df3[,3]
    
    #calculate pi for nonsynonymous sites
    genes <- F_ST.stats(genes,subsites="nonsyn", detail=TRUE)
    df4=as.data.frame(get.diversity(genes)[[1]])
    head(df4)
    pisitenonsyn=df4[,3]
    
    #calculate Tajima's D
    genes <- F_ST.stats(genes, detail=TRUE)
    genes <- neutrality.stats(genes, detail=TRUE)
    df5=as.data.frame(get.neutrality(genes)[[1]])
    head(df5)
    tajd <- df5
    
    #create the output as a data frame
    ##pi not yet standardised by the number of sites - will need to be done later
    output <- data.frame("GeneLocation"=genes@region.names,
                         "Gene_Length"=genes@n.sites,
                         "pi"=pisite, 
                         "pi_syn"=pisitesyn, 
                         "pi_nonsyn"=pisitenonsyn,
                         "Tajimas_D"=tajd)
    
    #Remove duplicate rows, while retaining the first occurrence of the duplicate
    output_dup_remove <- unique(output)
    
    # Write output_dup_remove to CSV file
    write.csv(output_dup_remove, file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Polymorphism_contigs",paste0(contig,"_polymorphism.csv")), row.names = F)
    
    cat("Analysis for", contig, "completed.\n")
  }, error = function(e) {
    cat("Error in analysis for", contig, ": ", conditionMessage(e), "\n")
  })
}



#######

#Identify which locus each CDS belongs to and sum together polymorphism statics before normalising by gene length

#Create directory to save the csv files for the summed and standardised polymorphism stats
dir.create(file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Poly_summed_standard"))

# Set the vector of contigs based on folder names
## List only immediate subdirectories (excluding subdirectories within those)
poly_stat_files <- list.files(file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen", SDP,"PopGenome/Polymorphism_contigs"))
## Filter subdirectories based on the specified criteria
poly_stat_contigs <- sub("_polymorphism.csv", "", poly_stat_files)


for (contig in poly_stat_contigs) {
  #Read in the locus identifier and polymorphism stats files according to contig
  locus_identifier <- read.csv(file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Locus_identifiers",paste0(contig,"_locus_identifier.csv")))
  poly_stats <- read.csv(file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Polymorphism_contigs",paste0(contig,"_polymorphism.csv")))
  
  #Create a new column in poly_stats dataframe called "Gene_Locus"
  poly_stats$Gene_Locus <- NA
  
  #Assign values to Gene_Locus in poly_stats from from locus_identifier$Gene_Locus based on matching values within the GeneLocation columns
  # Find matching rows in locus_identifier based on GeneLocation
  matches <- match(poly_stats$GeneLocation, locus_identifier$GeneLocation)
  
  # Update Gene_Locus in poly_stats where matches are found
  poly_stats$Gene_Locus[!is.na(matches)] <- locus_identifier$Gene_Locus[matches[!is.na(matches)]]
  
  #Select the columns of interest
  poly_stats_select <- poly_stats %>% dplyr::select(Gene_Length, pi, pi_syn, pi_nonsyn, Gene_Locus)
  
  # Group by Gene_Locus and calculate the sum for each group
  poly_stats_summed <- poly_stats_select %>%
    group_by(Gene_Locus) %>%
    summarise(
      Gene_Length_Sum = sum(Gene_Length, na.rm = TRUE),
      pi_Sum = sum(pi, na.rm = TRUE),
      pi_syn_Sum = sum(pi_syn, na.rm = TRUE),
      pi_nonsyn_Sum = sum(pi_nonsyn, na.rm = TRUE)
    ) %>%
    ungroup()  # Remove grouping
  
  #Calculate polymorphism standardised by gene length
  poly_stats_standardised <- poly_stats_summed %>%
    mutate(
      pi_Sum_site = pi_Sum / Gene_Length_Sum,
      pi_syn_Sum_site = pi_syn_Sum / Gene_Length_Sum,
      pi_nonsyn_Sum_site = pi_nonsyn_Sum / Gene_Length_Sum
    )
  
  # Write poly_stats_standardised to csv file
  write.csv(poly_stats_standardised, file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Poly_summed_standard",paste0(contig,"_poly_sum_stand.csv")), row.names = F)
  
  
  cat("Processing for", contig, "completed.\n")
}




################
#Combine the csv files created for each contig into a single csv file
################

# Specify the directory containing your CSV files
csv_path <- file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Poly_summed_standard")

# Get a list of CSV files in the directory
csv_files <- list.files(csv_path, pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine with existing data
for (csv_file in csv_files) {
  current_data <- read.csv(csv_file)
  combined_data <- rbind(combined_data, current_data)
}

#Make a new directory for the combined csv file
dir.create(file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Combined_csv"))

# Write the combined data to a new CSV file
write.csv(combined_data, file = file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Combined_csv",paste0(SDP,"_polymorphism.csv")), row.names = F)



########

#Identify loci in the combined csv that are queen/worker-biased

#Read in the csv from the DEG analysis and rename the first column to Gene_Locus
DEG_analysis <- read.csv(file.path("/Volumes/Pop_Gen/Differential_gene_expression/DGE_csv_outputs",paste0(SDP,".csv")))
DEG_analysis <- DEG_analysis %>% rename(Gene_Locus = X)

#Read in polymorphism data, per locus, for the whole genome
poly_genome <- read.csv(file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Combined_csv",paste0(SDP,"_polymorphism.csv")))

#Filter rows in poly_genome by the Loci present in the DEG_analysis
filtered_poly_genome <- poly_genome[poly_genome$Gene_Locus %in% DEG_analysis$Gene_Locus, ]
table(DEG_analysis$Gene_Locus %in% poly_genome$Gene_Locus)
setdiff(DEG_analysis$Gene_Locus, filtered_poly_genome$Gene_Locus) #Identify the loci where pomorphism estimates are mising for genes with differential expression

#Assign queen/worker status to each loci in new column in filtered_poly_genome based on DEG_analysis
#Create new column in filtered_poly_genome called biased_gene
filtered_poly_genome$biased_gene <- NA

# Find matching rows in DEG_analysis based on Gene_Locus
Gene_Locus_matches <- match(filtered_poly_genome$Gene_Locus, DEG_analysis$Gene_Locus)

# Update biased_gene in filtered_poly_genome where matches are found
filtered_poly_genome$biased_gene[!is.na(Gene_Locus_matches)] <- DEG_analysis$biased_gene[Gene_Locus_matches[!is.na(Gene_Locus_matches)]]

#Make a new directory for the final output csv file
dir.create(file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/OUTPUT"))

#Write filtered_poly_genome to csv
write.csv(filtered_poly_genome, file = file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/OUTPUT",paste0(SDP,"_poly_OUTPUT.csv")), row.names = F)


#END
################


#Why are some loci missing from our calculations of polymorphism?
##There are five loci in total: "LOC105279263" "LOC105280165" "LOC105283506" "LOC105284284" "LOC113563616" (4 out of 5 are queen genes)
###All 5 are present in the gtf file for the individual contig and for the whole genome gtf, they're also present in locus_identifiers - therefore, it's not being missed in my scripts that are finding text in the gtf files
###However, it seems that the genomic region associated with this locus is not included in the calculation of polymorphism. This suggests that PopGenome is not analysing it under the setting of CDS. Maybe it's because this locus is a non-coding region? - having looked into the gtf file it suggests that these regions are long non-coding RNA sequences and therefore not counting as CDS
#Does there appear to be any bias towards higher polymorphism in worker/queen genes? - no strong patterns




# 
# ################
# #Combine the csv files created for each contig into a single csv file - before pivoting to Laurie's approach where we assign loci to each CDS
# ################
# 
# # Specify the directory containing your CSV files
# csv_path <- file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome")
# 
# # Get a list of CSV files in the directory
# csv_files <- list.files(csv_path, pattern = "\\.csv$", full.names = TRUE)
# 
# # Initialize an empty data frame to store the combined data
# combined_data <- data.frame()
# 
# # Loop through each CSV file, read it, and combine with existing data
# for (csv_file in csv_files) {
#   current_data <- read.csv(csv_file)
#   combined_data <- rbind(combined_data, current_data)
# }
# 
# #Make a new directory for the combined csv file
# dir.create(file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Combined_csv"))
# 
# # Write the combined data to a new CSV file
# write.csv(combined_data, file = file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Combined_csv",SDP,"_polymorphism.csv"), row.names = F)



##########
#For loop for calculating polymorphism and Tajima D stats before pivoting to Laurie's approach where we assign loci to each CDS
# for (contig in contigs) {
#   tryCatch({
#     # Set the working directory
#     setwd(file.path("/Volumes/Pop_Gen/PopGenomeWD", SDP, contig))
#     
#     # Load vcf file and gff file into R
#     GENOME.class <- readData("vcf_files/vcf", format = "VCF", gffpath = "vcf_files/gff")
#     
#     # Load the sequence of the reference genome
#     GENOME.class <- set.synnonsyn(GENOME.class, ref.chr = paste0(contig, ".fna"))
#     
#     # Find the position as to where each gene starts and ends
#     gene.positions <- get_gff_info(gff.file = paste0("vcf_files/gff/", contig, ".gtf"), chr = contig, feature = "CDS")
#     
#     # Split the VCF file into individual genes
#     genes <- split_data_into_GFF_features(GENOME.class, paste0("vcf_files/gff/", contig, ".gtf"), contig, "CDS")
#     
#     # Preliminary calculations for the FST module
#     genes <- F_ST.stats(genes, detail = TRUE)
#     
#     # Present diversity measures as a data frame
#     df2 <- as.data.frame(get.diversity(genes)[[1]])
#     
#     # Select the column containing the polymorphism values
#     pisite <- df2[, 3]
#     
#     #calculate pi for synonymous sites
#     genes <- F_ST.stats(genes,subsites="syn", detail=TRUE)
#     df3=as.data.frame(get.diversity(genes)[[1]])
#     head(df3)
#     pisitesyn=df3[,3]
#     
#     #calculate pi for nonsynonymous sites
#     genes <- F_ST.stats(genes,subsites="nonsyn", detail=TRUE)
#     df4=as.data.frame(get.diversity(genes)[[1]])
#     head(df4)
#     pisitenonsyn=df4[,3]
#     
#     #calculate Tajima's D
#     genes <- F_ST.stats(genes, detail=TRUE)
#     genes <- neutrality.stats(genes, detail=TRUE)
#     df5=as.data.frame(get.neutrality(genes)[[1]])
#     head(df5)
#     tajd <- df5
#     
#     #create the output as a data frame
#     ##standardise pi by the number of sites, so that we have estimates that are in relation to gene length
#     output <- data.frame("GeneLocation"=genes@region.names,
#                          "Gene_Length"=genes@n.sites, # standardise pi for gene length
#                          "pi_site"=pisite/genes@n.sites, 
#                          "pi_syn_site"=pisitesyn/genes@n.sites, 
#                          "pi_nonsyn_site"=pisitenonsyn/genes@n.sites,
#                          "Tajimas_D"=tajd)
#     
#     # Write output to CSV file
#     write.csv(output, file.path("/Users/louis.bell-roberts/Documents/Github/Pop_Gen",SDP,"PopGenome/Polymorphism_contigs",paste0(contig,"_polymorphism.csv")), row.names = F)
#     
#     cat("Analysis for", contig, "completed.\n")
#   }, error = function(e) {
#     cat("Error in analysis for", contig, ": ", conditionMessage(e), "\n")
#   })
# }







