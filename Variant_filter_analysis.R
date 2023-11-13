################################################################################
##                            Name
##                            Date  
##                          Institute
##          Verification of variant filter in Ion-Reporter
################################################################################

# 1. Load libraries ------------------------------------------------------------
library(xlsx)
library(tidyverse)
library(janitor)

# 3. Import the unfiltered variant files ------------------------------
path <- file.path("path/to/unfiltered_variants_file") #tsv file after selecting "no filter" from Ion Reporter

unfiltered_filenames.list <- list.files(path = path, pattern='*_unfiltered.tsv',
                                        full.names = TRUE) 

unfiltered_variants.list <- lapply(unfiltered_filenames.list, readr::read_tsv, sheet=1, skip=17) # this takes everything (CNVs and non-CNVs), regardless of applied filter in excel

# Import the filtered variants
path2 <- file.path("path/to/filtered_variants_file")

filtered_filenames.list <- list.files(path = path2, pattern='*.tsv',
                                      full.names = TRUE) 

filtered_variants.list <- lapply(filtered_filenames, readr::read_tsv, skip=17) #these are the filtered variants

## 3.1 name all dfs in the list -----------------------------------------------

### 3.1.1 name all dfs in the list ---------------------------------------------

# count the no of characters in the path & apply pattern matching code to extract filenames from each filenames list
nchar(path)
names(unfiltered_variants.list) <- gsub("^.{77}|_OCAplus.*", "", unfiltered_filenames.list)

# repeat for filtered
nchar(path2)
names(filtered_variants.list) <- gsub("^.{85}|_OCAplus.*", "", filtered_filenames.list)

### 3.1.2 remove columns that might interfere with row binding later on---------
cols_to_remove <- c("Subtype", "Call", "No Call Reason", "PhyloP", "SIFT",
                    "Grantham", "PolyPhen") #cols class conflict for manual and Alil lists

unfiltered_variants_shortened.list <- lapply(unfiltered_variants.list, 
                                             function(x){x[,!names(x) %in% cols_to_remove]})

filtered_variants_shortened.list <- lapply(filtered_variants.list, 
                                           function(x){x[,!names(x) %in% cols_to_remove]})

### 3.1.3 change the class of some cols preventing row-merging of list----------------
filtered_variants_shortened.list <- lapply(filtered_variants_shortened.list, 
                                           function(x) mutate_at(x, .vars = c("Ref", "MAF"), as.character))

### 3.1.4 flatten and merge all dfs in list-------------------------------------
unfiltered_variants <- bind_rows(unfiltered_variants_shortened.list, .id = "Sample.ID")
filtered_variants <- bind_rows(filtered_variants_shortened.list, .id = "Sample.ID")

# 4. Apply NGS filter protocol to unfiltered variants df and compare with filtered results --------------------------------------

# Start with non-CNVs
R_filtered_variants_non_CNVs <- unfiltered_variants %>%
  # you need to create a dummy numeric MAF column
  mutate(MAF2 = case_when(grepl("(ref)", MAF) ~ 1.00, #make values that contain "ref" higher than 0.001, since you need to exclude them anyway
                          grepl("0.0,0.001", MAF) ~ 0.00055, #double values that are <=0.001 should be given a unique value < 0.001
                          TRUE ~ as.numeric(MAF)), #NAs should be left as well as other values.
         .after = MAF) %>%
  replace_na(list(`Variant Effect` = "unknown")) %>% #replace NAs here, so they aren't lost by filter
  filter(Type %in% c("SNV", "MNV", "INDEL"),
         !`UCSC Common SNPs` %in% c("YES"),
         !(`Variant Effect` == "synonymous" | `Variant Effect` == "synonymous, synonymous"),
         MAF2 <= 0.001 | is.na(MAF), # keep NAs
         GMAF <= 0.001 | is.na(GMAF)) %>% # keep NAs
  replace_na(list(Location = "unknown")) %>% #replace NAs here with unknown, so NA rows aren't removed with the grepl function
  filter(grepl("unknown|exonic|utr|splicesite|downstream|upstream", Location), #keeps rows with =>2 gene locations (e.g., intronic and utr)
         `Phred QUAL Score` >=30,
         Coverage >= 100,
         `Allele Frequency %` >=3) 

# review total counts and non-CNVs for each sample.
R_filtered_variants_non_CNVs %>%
  group_by(Sample.ID) %>%
  count(n = n())

# 5. Find discordant somatic variants between Ion Torrent filter and R-filter above----------
discordant_variants_non_CNVS <- R_filtered_variants_non_CNVs %>%
  anti_join(filtered_variants, by = c("Sample.ID", "Locus", "Genes"))

# review discordant total counts for each sample.
discordant_variants_non_CNVS %>%
  group_by(Sample.ID) %>%
  count(n = n())

# 8. Save discordant df --------------------------------
write.xlsx2(as.data.frame(discordanct_variants_non_CNVs), "./Output/discordant_variants_non_CNVs.xlsx", 
            row.names=FALSE, col.names=TRUE)

#####################################################
# Count the pathways (genes) that are highly represented in samples
#####################################################

# 9. Count all genes found in samples.
#Gene count by sample
Gene_counts <- R_filtered_variants_non_CNVs %>%
  mutate(Genes2 = as.factor(Genes),
         .after = Genes) %>%
  group_by(Sample.ID) %>%
  count(Genes2) %>%
  arrange(desc(n), .by_group = TRUE) 
#janitor::adorn_totals("row") (4130 is the total genes)

#top15
Gene_counts_top15 <- Gene_counts %>%
  top_n(15) 


# 10. Create doughnut plot.
# Compute percentages
Gene_counts_top15$fraction = Gene_counts_top15$n / sum(Gene_counts_top15$Genes2)

# Compute the cumulative percentages (top of each rectangle)
Gene_counts_top15$ymax = cumsum(Gene_counts_top15$fraction)

# Compute the bottom of each rectangle
Gene_counts_top15$ymin = c(0, head(Gene_counts_top15$ymax, n=-1))

# Compute label position
Gene_counts_top15$labelPosition <- (Gene_counts_top15$ymax + Gene_counts_top15$ymin) / 2

# Compute a good label
Gene_counts_top15$label <- paste0(Gene_counts_top15$Genes2, ":", Gene_counts_top15$n)

# Make the plot
Gene_counts_doughnut_plot <- ggplot(Gene_counts_top15, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Genes2)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=2.8) +
  coord_polar(theta="y") + # this changes rectangular chart to a circle (pie) chart
  xlim(c(2, 4)) + # remove to make a pie chart
  theme_void() +
  scale_fill_manual(values = c("darkgreen", "brown", "cornflowerblue", "purple", "coral", "pink",
                                          "deeppink", "cornsilk3", "yellow", "aquamarine", "cadetblue",
                                          "skyblue", "grey", "violet", "cyan")) +
                                            ggtitle("Top 15 genes") +
  theme(legend.position = "none",
        plot.title.position = "plot") # align title to center

ggsave("./Plots/Top_15_genes.pdf", Top_15_genes_plot,
       width = 16, height = 22, units = "cm")

