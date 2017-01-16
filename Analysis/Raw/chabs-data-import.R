# Import mothur files and sample metadata
sharedfile = "mothur/chabs.shared"
taxfile = "mothur/chabs-fwdb-silva.taxonomy"
mapfile = "subdata-metadata.csv"

mothurdata = import_mothur(mothur_shared_file = sharedfile, 
                           mothur_constaxonomy_file = taxfile)

# Add the OTU number as a column in the taxonomy file
tax_table(mothurdata) <- cbind(tax_table(mothurdata), 
                               row.names(tax_table(mothurdata)))

# Rename the taxonomy columns
colnames(tax_table(mothurdata)) <- 
  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Rank7", "Rank8", "Species")

# Import sample metadata and transform into sample_data class
map <- read.csv(mapfile)

# Convert map to sample_data class
# Merge mothurdata object with sample_data
map <- sample_data(map)
rownames(map) <- map$SampleID
mothur.merge <- merge_phyloseq(mothurdata, map)


# Filter out non-samples (i.e. water, mock) and samples from intensive cruises.
erie.subset <-
  mothur.merge %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family != "mitochondria" &
      Class != "Chloroplast"
  ) %>%
  subset_samples(
    Type == "sample" & 
      !(Date %in% c("5/27", "6/10", "8/11", "11/3"))
  )

erie <- subset_samples(erie.subset, Fraction == "CNA")


# Also prune out taxa which were only present in removed samples
erie <- prune_taxa(taxa_sums(erie) > 0, erie)

save(erie, file = "erie-data.Rdata")
