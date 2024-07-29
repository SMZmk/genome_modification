#install.packages("rtracklayer")
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
################################################################################ [2024-02-22]
setwd("~/Desktop/Rhome/Mutation-Simulator/")
###################################### NAME INPUT
gff_file <- "./ref/Arabidopsis_thaliana.TAIR10.55.gff3"
load("~/Desktop/Rhome/moca_blue/0MOTIFS/He_etal_Models/AthAly.genome.snp.RData")
###################################### NAME OUTPUT
new_gff_file1 <- "Arabidopsis_thaliana.TAIR10.55_IGS.gff"
new_gff_file2 <- "Arabidopsis_thaliana.TAIR10.55_IGS-SNPs.gff"
####################################### IMPORT DATA
AthAly.genome.snp_df <- as.data.frame(AthAly.genome.snp)
#########################################
gff_data <- rtracklayer::import(gff_file)
gff_df=as.data.frame(gff_data)
genes <- gff_df %>%
  filter(type == "gene")
levels_seqnames <- levels(genes$seqnames)
(sum(genes$width))/27655
#### CREATE DUMMY ANNOTATIONS FOR START OF CHROMOSOMES ######## DATA WRANGLING
dummy_df <- data.frame()
for (level in levels_seqnames) {
  dummy <- data.frame(seqnames = level,
                      start = 0,
                      end = 1,
                      width = 1,
                      strand = "+",
                      source = "araport11",
                      type = NA,
                      score = NA,
                      phase = NA,
                      ID = NA,
                      Alias = NA,
                      Name = NA,
                      biotype = NA,
                      description = NA,
                      gene_id = NA,
                      logic_name = NA,
                      Parent = NA,
                      tag = NA,
                      transcript_id = NA,
                      constitutive = NA,
                      ensembl_end_phase = NA,
                      ensembl_phase = NA,
                      exon_id = NA,
                      rank = NA, 
                      protein_id = NA,
                      Is_circular = NA)
  dummy_df<-rbind(dummy_df, dummy)
}
################################################################################
######### CREATE DUMMY ANNOTATIONS FOR THE END OF CHROMOSOMES ##################
dummy_df2 <- data.frame()
for (level in levels_seqnames) {
  max_start <- max(gff_df$end[gff_df$seqnames == level])
  dummy <- data.frame(seqnames = level,
                      start = max_start,
                      end = max_start,
                      width = 1,
                      strand = "+",
                      source = "araport11",
                      type = NA,
                      score = NA,
                      phase = NA,
                      ID = NA,
                      Alias = NA,
                      Name = NA,
                      biotype = NA,
                      description = NA,
                      gene_id = NA,
                      logic_name = NA,
                      Parent = NA,
                      tag = NA,
                      transcript_id = NA,
                      constitutive = NA,
                      ensembl_end_phase = NA,
                      ensembl_phase = NA,
                      exon_id = NA,
                      rank = NA, 
                      protein_id = NA,
                      Is_circular = NA)
  dummy_df2 <- rbind(dummy_df2, dummy)
}
################################
dummy_df2
dummy_dfs <- rbind(dummy_df, dummy_df2)
genes0 <- rbind(genes, dummy_dfs)
genes0 <- genes0 %>% 
  arrange(seqnames, start)
genes0 <- genes0 %>%
  mutate(idx1 = row_number())
genes0 <- genes0 %>%
  mutate(idx2 = lag(idx1))
#
idx_start<- genes0 %>%
  select(idx2, start)
idx_end <- genes0 %>%
  select(idx1, end)
names(idx_start)[1] <- "idx"
names(idx_end)[1] <- "idx"
#
idx_start0 <- idx_start[-c(1),]
idx_end0 <- idx_end[-c(max(idx_end$idx)),]
igs_ranges <- cbind(idx_start0, idx_end0)
colnames(igs_ranges)[c(1:4)] <- c("idx1","end","idx2","start")
igs_ranges0 <- igs_ranges %>%
  select(-idx1)
############################################################################
genes1 <- genes0 %>%
  select(-start, -end, -idx1, -type, -width)
igs_df0<- merge(genes1, igs_ranges0, by = "idx2")
igs_df0$type <- c("igs")
igs_df0$width <- c(igs_df0$end- igs_df0$start)
igs_df0 <- igs_df0[igs_df0$width >= 2, ]
igs_gr <- with(igs_df0, GRanges(seqnames, IRanges(start, end)))
export(igs_gr, new_gff_file1, format = "GFF3")
################################################################################
AthAly.genome.snp_df0 <- AthAly.genome.snp_df %>%
  filter(tha != lyr)
AthAly.genome.snp_df0$chr <- gsub("Chr", "", AthAly.genome.snp_df0$chr)
head(AthAly.genome.snp_df0)
modifications <- data.frame(chr = AthAly.genome.snp_df0$chr,
                            start = as.numeric(AthAly.genome.snp_df0$Atpos),
                            end = as.numeric(AthAly.genome.snp_df0$Atpos),
                            REF = paste0(AthAly.genome.snp_df0$tha,">",AthAly.genome.snp_df0$lyr))
snp_gr <- GRanges(
  seqnames = modifications$chr,
  ranges = IRanges(start = modifications$start, end = modifications$end)
)
combined_gr <- c(igs_gr, snp_gr)
combined_gr<-sort(combined_gr)
head(combined_gr)
overlapping_rows <- findOverlaps(combined_gr, type = "any")
rows_to_remove <- queryHits(overlapping_rows)
combined_gr_filtered <- combined_gr[-rows_to_remove]
################################################################################
export(combined_gr_filtered, new_gff_file2, format = "GFF3")
################################################################################