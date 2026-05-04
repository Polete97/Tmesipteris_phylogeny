library(ape)
library(geiger)
library(phytools)
library(nlme)
library(dendextend)

####################Genome size reconstruction############################

setwd("")

tree <- ape::read.nexus(file="") #cronogram

plot(tree)

# Trim tree to specific species based on your new data
# Species codes to keep from your table
species_to_keep <- c("")


# Alternative: specify which species to remove
species_to_remove <- c("")

# Method 1: Keep only specified species using ape package
trimmed_tree <- ape::keep.tip(tree, species_to_keep)

# Method 2: Remove specified species (alternative approach)
# trimmed_tree <- ape::drop.tip(tree, species_to_remove)

# Ladderize the trimmed tree for better visualization
trimmed_tree <- ape::ladderize(trimmed_tree)
plot(trimmed_tree)
# Rename tip labels to species names based on your data table
# Load the data first to get the mapping
data <- read.csv("tmesepteris_data.csv")

# Create a mapping from code to species name
code_to_species <- setNames(data$Species, data$code)

# Rename the tips in the trimmed tree
for(i in 1:length(trimmed_tree$tip.label)) {
  code <- trimmed_tree$tip.label[i]
  if(code %in% names(code_to_species)) {
    trimmed_tree$tip.label[i] <- code_to_species[code]
  }
}

# Plot the trimmed tree
plot(trimmed_tree)
title("Trimmed Tree")

# Continue with your analysis using the trimmed tree
library(readr)
# Load the new Tmesipteris data
data <- read.csv("tmesepteris_data.csv")

# Handle duplicate T_horomaka entries by averaging or taking the first entry
# Option 1: Take mean of duplicates
library(dplyr)
data_clean <- data %>%
  group_by(Species) %>%
  summarise(
    code = first(code),
    GS = mean(GS, na.rm = TRUE),
    sd = mean(sd, na.rm = TRUE),
    ploidy = first(ploidy),
    cx = mean(cx, na.rm = TRUE),
    .groups = 'drop'
  )

# Create named vectors using species names that match tree tip labels
cvalues2 <- setNames(data_clean$GS, data_clean$Species)
cxvalues <- setNames(data_clean$cx, data_clean$Species)

# Remove any NA values
cvalues2 <- cvalues2[!is.na(cvalues2)]
cxvalues <- cxvalues[!is.na(cxvalues)]

# Check which species are in both tree and data
cat("Tree tips:", trimmed_tree$tip.label, "\n")
cat("Data species (GS):", names(cvalues2), "\n")
cat("Data species (cx):", names(cxvalues), "\n")

# Make sure your data matches the trimmed tree
obj_gs <- name.check(trimmed_tree, cvalues2)
obj_cx <- name.check(trimmed_tree, cxvalues)
cat("GS data check:\n")
print(obj_gs)
cat("cx data check:\n")
print(obj_cx)

# Subset data to only include species present in both tree and data
common_species_gs <- intersect(trimmed_tree$tip.label, names(cvalues2))
common_species_cx <- intersect(trimmed_tree$tip.label, names(cxvalues))

cvalues2_final <- cvalues2[common_species_gs]
cxvalues_final <- cxvalues[common_species_cx]

# If needed, also subset the tree to match available data (use the intersection of both datasets)
common_species_both <- intersect(common_species_gs, common_species_cx)
if(length(common_species_both) < length(trimmed_tree$tip.label)) {
  trimmed_tree <- ape::keep.tip(trimmed_tree, common_species_both)
  cvalues2_final <- cvalues2_final[common_species_both]
  cxvalues_final <- cxvalues_final[common_species_both]
  cat("Tree trimmed to", length(common_species_both), "species with both GS and cx data\n")
}

# Use the trimmed tree for analysis
# GS analysis
fitted <- fastAnc(trimmed_tree, log(cvalues2_final))
cMap <- contMap(trimmed_tree, cvalues2_final, plot = FALSE, method="user", anc.states=exp(fitted), res=1000)
library(viridis)
cMap <- setMap(cMap, rev(plasma(100)))

# cx analysis
fitted_cx <- fastAnc(trimmed_tree, log(cxvalues_final))
cMap_cx <- contMap(trimmed_tree, cxvalues_final, plot = FALSE, method="user", anc.states=exp(fitted_cx), res=1000)
cMap_cx <- setMap(cMap_cx, rev(viridis(100)))

# Plot the genome size reconstruction
plot(cMap, fsize=0.8, ftype="i")
title("Genome Size (GS) Ancestral State Reconstruction")

# Plot the cx reconstruction
plot(cMap_cx, fsize=0.8, ftype="i")
title("cx Value Ancestral State Reconstruction")

# Alternative barplot visualizations
obj_gs <- plotTree.barplot(trimmed_tree, cvalues2_final, args.plotTree=list(fsize=0.7))
title("GS Values")

obj_cx <- plotTree.barplot(trimmed_tree, cxvalues_final, args.plotTree=list(fsize=0.7))
title("cx Values")

# Plot tree with ancestral node values displayed for GS
plot(trimmed_tree, show.tip.label=TRUE, cex=0.8, main="Phylogeny with Ancestral GS Values")
nodelabels(round(exp(fitted), 2), cex=0.6, bg="lightblue", frame="circle")
tiplabels(round(cvalues2_final, 2), cex=0.5, bg="lightgreen", frame="circle")

# Plot tree with ancestral node values displayed for cx
plot(trimmed_tree, show.tip.label=TRUE, cex=0.8, main="Phylogeny with Ancestral cx Values")
nodelabels(round(exp(fitted_cx), 2), cex=0.6, bg="lightcoral", frame="circle")
tiplabels(round(cxvalues_final, 2), cex=0.5, bg="lightyellow", frame="circle")

# Export node value reconstructions to Excel
library(openxlsx)

# Create data frames with node information
node_data_gs <- data.frame(
  Node = 1:length(fitted),
  Ancestral_GS = round(exp(fitted), 3),
  Log_Ancestral_GS = round(fitted, 3)
)

node_data_cx <- data.frame(
  Node = 1:length(fitted_cx),
  Ancestral_cx = round(exp(fitted_cx), 3),
  Log_Ancestral_cx = round(fitted_cx, 3)
)

# Create data frames with tip information
tip_data_gs <- data.frame(
  Species = names(cvalues2_final),
  Observed_GS = round(cvalues2_final, 3),
  Log_Observed_GS = round(log(cvalues2_final), 3)
)

tip_data_cx <- data.frame(
  Species = names(cxvalues_final),
  Observed_cx = round(cxvalues_final, 3),
  Log_Observed_cx = round(log(cxvalues_final), 3)
)

# Create a workbook and add worksheets
wb <- createWorkbook()
addWorksheet(wb, "Ancestral_States_GS")
addWorksheet(wb, "Ancestral_States_cx")
addWorksheet(wb, "Tip_Values_GS")
addWorksheet(wb, "Tip_Values_cx")
addWorksheet(wb, "Summary")

# Write data to worksheets
writeData(wb, "Ancestral_States_GS", node_data_gs)
writeData(wb, "Ancestral_States_cx", node_data_cx)
writeData(wb, "Tip_Values_GS", tip_data_gs)
writeData(wb, "Tip_Values_cx", tip_data_cx)

# Create summary information
summary_data <- data.frame(
  Statistic = c("Number_of_Species", "Number_of_Nodes", "Min_GS", "Max_GS", "Mean_GS", "Root_Ancestral_GS",
                "Min_cx", "Max_cx", "Mean_cx", "Root_Ancestral_cx"),
  Value = c(
    length(cvalues2_final),
    length(fitted),
    round(min(cvalues2_final), 3),
    round(max(cvalues2_final), 3),
    round(mean(cvalues2_final), 3),
    round(exp(fitted[1]), 3),  # Root node GS
    round(min(cxvalues_final), 3),
    round(max(cxvalues_final), 3),
    round(mean(cxvalues_final), 3),
    round(exp(fitted_cx[1]), 3)  # Root node cx
  )
)
writeData(wb, "Summary", summary_data)

# Save the Excel file
saveWorkbook(wb, "GS_cx_Ancestral_Reconstruction.xlsx", overwrite = TRUE)
cat("Ancestral state reconstruction for both GS and cx exported to GS_cx_Ancestral_Reconstruction.xlsx\n")


####################Tanglegram############################

#Read your trees (Newick/Nexus are fine)
nuc <- ape::read.nexus(file="final_epa_result_dated.tog.tre.nex")
plast <-ape::read.nexus(file="tangle_tree.nex")


plot(plast)

# Optional synonym mapping (if needed)
# synonyms.csv: two columns: from,to (both already normalized)
if (file.exists("synonyms.csv")) {
  map <- read.csv("synonyms.csv", stringsAsFactors = FALSE)
  renamer <- function(lbls, map) {
    mm <- setNames(map$to, map$from)
    lbls <- ifelse(lbls %in% names(mm), mm[lbls], lbls)
    as.character(lbls)
  }
  nuc$tip.label   <- renamer(nuc$tip.label, map)
  plast$tip.label <- renamer(plast$tip.label, map)
}


common <- intersect(nuc$tip.label, plast$tip.label)

cat("Nuclear tips:", length(nuc$tip.label),
    "\nPlastid tips:", length(plast$tip.label),
    "\nCommon tips:", length(common), "\n")

# See what’s unique to each tree (useful to inspect)
setdiff(nuc$tip.label, common)
setdiff(plast$tip.label, common)

# Prune both trees to the same set
nuc_pruned   <- drop.tip(nuc,   setdiff(nuc$tip.label,   common))
plast_pruned <- drop.tip(plast, setdiff(plast$tip.label, common))


# Useful cleanups
nuc_pruned   <- ladderize(nuc_pruned, right = TRUE)
plast_pruned <- ladderize(plast_pruned, right = TRUE)
nuc_pruned <- force.ultrametric(nuc_pruned, method = "extend")
plast_pruned <- force.ultrametric(plast_pruned, method = "extend")

set.seed(123)  # for reproducible rotation search
co <- cophylo(nuc_pruned, plast_pruned, rotate = TRUE)  # tries to minimize crossings

# Inspect crossings
cat("Estimated crossings:", co$Ncross, "\n")

# Plot
plot(
  co,
  link.type = "curved",   # "curved" or "straight"
  link.lwd  = 1,
  fsize     = 0.6,
  col       = "black"
)

# Save to file
pdf("tanglegram_phytools.pdf", width = 1600, height = 1000, res = 180)
plot(co, link.type = "curved", link.lwd = 1, fsize = 0.6)
dev.off()











