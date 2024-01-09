########################Prepare files (vcf, gtf and fna) for analysis in PopGenome package in R######################
###PopGenome runs into an error when the reference genome is made up of multiple contigs. Therefore all files must be split up into individual contigs
##Of the total sum of the unique contigs present in the reference genome, not all of these are always present in either the SNP calling file (vcf) or the gene annotation file (gtf). Therefore, these contigs are not carried forward for the PopGenome analysis - check whether this is appropriate with Laurie
#Louis Bell-Roberts
#17/12/2023

#Set working directory
setwd("/Volumes/Pop_Gen/PopGenomeWD/")

###########
#IMPORTANT#
###########
#Prior to running the script, create a folder using the species name and place the fna, gtf and vcf files in it
#Rename fna, gtf and vcf files so that they are just Genus_species.xxx

##########
#Create the folder structuring
##########

SDP <- "Ooceraea_biroi"
dir.create(file.path(SDP,"fastas"), showWarnings = FALSE) #Create the folder to house the fna contigs
dir.create(file.path(SDP,"gtf_files"), showWarnings = FALSE) #Create the folder to house the gtf contigs
dir.create(file.path(SDP,"vcf_files"), showWarnings = FALSE) #Create the folder to house the vcf contigs




###############################

#Split up the vcf, gtf and fna files based on the separate contigs

##########
#Create the contigs.txt file
##Extract all contigs where gene annotation is available in the gtf file
##########


# Read lines of the gtf file for each contig
gtf_content <- readLines(file.path(SDP, paste0(SDP, ".gtf")))
contigs <- character()

# Process each line in the GTF file - takes a few minutes to run
for (line in gtf_content) {
  # Skip comment lines
  if (startsWith(line, "#")) next
  
  # Split the line into fields
  fields <- strsplit(line, "\t")[[1]]
  
  #Extract the contig names
  contigs <- c(contigs, fields[1])
}

# Sort and get unique contigs
sorted_contigs <- unique(sort(contigs)) #63 unique contigs

# Write sorted contigs to the output file
writeLines(sorted_contigs, file.path(SDP, "contigs.txt"))



##########
#Convert the reference genome into its constituent contigs and save each as a separate fasta file
##########

# Read contig names from the contigs.txt file - used just to check what the names of these contigs are
contig_names <- readLines(file.path(SDP, "contigs.txt"))

#Read in the lines of the reference genome
genome_lines <- readLines(file.path(SDP,paste0(SDP,".fna")))

# Initialize variables
current_contigs <- NULL
current_sequence <- NULL

# Iterate through each line in the sequence file
for (line in genome_lines) {
  if (startsWith(line, ">")) {
    # If the line starts with ">", it's a header
    # Check if we already have a current contig
    if (!is.null(current_contigs)) {
      # If yes, save the sequence to a new file
      output_filename <- file.path(SDP,"fastas", paste0(current_contigs, ".fna"))
      writeLines(c(paste0(">", current_contigs), current_sequence), output_filename)
    }
    
    # Extract the current contig name from the header
    current_contigs <- sub("^>(\\S+).+", "\\1", line)
    current_sequence <- NULL  # Reset the sequence for the new contigs
  } else {
    # If it's not a header line, it's part of the DNA sequence
    current_sequence <- c(current_sequence, line)
  }
}

# Save the last contig after the loop ends
if (!is.null(current_contigs)) {
  output_filename <- file.path(SDP,"fastas", paste0(current_contigs, ".fna"))
  writeLines(c(paste0(">", current_contigs), current_sequence), output_filename)
}



##########
#Create separate .gtf files based on the separate contigs
##########

# Set the file paths
contigs_file <- file.path(SDP, "contigs.txt")
gtf_file <- file.path(SDP, paste0(SDP, ".gtf"))

# Assign a directory to store the output GTF files
output_directory <- file.path(SDP, "gtf_files")

# Read contig names from contigs.txt
contigs <- readLines(contigs_file)

# #OLD LOOP # Loop through contig names and create separate GTF files - this loop doesn't identify the contigs that are missing in the main gtf file
# for (contig in contigs) {
#   # Set the output file path
#   output_file <- file.path(output_directory, paste0(contig, ".gtf"))
# 
#   # Read lines from the original GTF file
#   gtf_lines <- readLines(gtf_file)
# 
#   # Identify header lines (lines starting with "#")
#   header_lines <- gtf_lines[startsWith(gtf_lines, "#")]
# 
#   # Identify contig-specific lines
#   contig_lines <- gtf_lines[startsWith(gtf_lines, paste0(contig, "\t"))]
# 
#   # Write the header lines and contig-specific lines to the output file
#   writeLines(c(header_lines, contig_lines), output_file)
# 
#   cat("Created", output_file, "\n")
# }

# Loop through contig names and create separate GTF files
for (contig in contigs) {
  # Set the output file path
  output_file <- file.path(output_directory, paste0(contig, ".gtf"))
  
  # Read lines from the original GTF file
  gtf_lines <- readLines(gtf_file)
  
  # Check if the contig is present in the GTF lines
  if (any(startsWith(gtf_lines, paste0(contig, "\t")))) {
    # Identify header lines (lines starting with "#")
    header_lines <- gtf_lines[startsWith(gtf_lines, "#")]
    
    # Identify contig-specific lines
    contig_lines <- gtf_lines[startsWith(gtf_lines, paste0(contig, "\t"))]
    
    # Write the header lines and contig-specific lines to the output file
    writeLines(c(header_lines, contig_lines), output_file)
    
    cat("Created", output_file, "\n")
  } else {
    cat("contig", contig, "not found. Skipping.\n")
  }
}




##########
#Create separate .vcf files based on the separate contigs - does it matter that in each of the separate .vcf files, there are hashed out names of all of the contigs - not just the specific contig for that vcf file?
##########

# Set file paths
contigs_file <- file.path(SDP, "contigs.txt")
vcf_file <- file.path(SDP, paste0(SDP, ".vcf"))

# Assign a directory to store the output VCF files
output_directory <- file.path(SDP,"vcf_files")

# Read contig names from contigs.txt
contigs <- readLines(contigs_file)

# Loop through contig names and create separate VCF files
for (contig in contigs) {
  # Set the output file path
  output_file <- file.path(output_directory, paste0(contig, ".vcf"))
  
  # Read lines from the original VCF file
  vcf_lines <- readLines(vcf_file)
  
  # Identify header lines (lines starting with "#")
  header_lines <- grep("^#", vcf_lines, value = TRUE)
  
  # Identify contig-specific lines
  contig_lines <- grep(sprintf("^%s\\b", contig), vcf_lines, value = TRUE)
  
  # Write the header lines and contig-specific lines to the output file
  writeLines(c(header_lines, contig_lines), output_file)
  
  cat("Created", output_file, "\n")
}




###############################

#Move all of the files created into a PopGenome friendly format

# Set variables
base_directory <- "/Volumes/Pop_Gen/PopGenomeWD"

# Function to get unique contigs from file names
get_unique_names <- function(directory, extension) {
  files <- list.files(directory, pattern = extension, full.names = TRUE)
  unique_names <- tools::file_path_sans_ext(basename(files))
  return(unique_names)
}

# Get unique contigs
contigs_gtf <- get_unique_names(file.path(base_directory, SDP, "gtf_files"), ".gtf")

# Iterate over contigs
for (contig in contigs_gtf) {
  contig_directory <- file.path(base_directory, SDP, contig) #contig is the name of the contig in the contigs.txt file that is being looped through
  dir.create(contig_directory) #Creates new directory based on the name of the contig in the loop
  setwd(contig_directory) #Sets working directory in this contig in the loop
  
  # Create the necessary subdirectories within the contig folder
  dir.create("gtf_files", showWarnings = FALSE)
  dir.create("vcf_files", showWarnings = FALSE)
  dir.create(file.path("vcf_files", "gff"), showWarnings = FALSE)
  dir.create(file.path("vcf_files", "vcf"), showWarnings = FALSE)
  
  # Move files
  ##gtf files (gene annotation)
  file.rename(
    file.path(base_directory, SDP, "gtf_files", paste0(contig, ".gtf")),
    file.path(contig_directory, "vcf_files", "gff", paste0(contig, ".gtf"))
  )
  ##vcf files (from SNP calling)
  file.rename(
    file.path(base_directory, SDP, "vcf_files", paste0(contig, ".vcf")),
    file.path(contig_directory, "vcf_files", "vcf", paste0(contig, ".vcf"))
  )
  ##fna files (from the reference genome)
  file.rename(
    file.path(base_directory, SDP, "fastas", paste0(contig, ".fna")),
    file.path(contig_directory, paste0(contig, ".fna"))
  )
}













