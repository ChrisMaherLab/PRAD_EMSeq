library(data.table)
library(dplyr)
library(ggplot2)

rm(list = ls())

#********************************* ----
# <<<<< Genetic and methylation annotation for each sample>>>>> -----

# 1.0 Plot the genetic and CpG annotation for each sample ----
# 1.1 Genetic annotation----
# 1.1.1 genetic annotation for each sample ----
# This should be run by /storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/validation_cohort/methylkit/methylkitAnalysis/general_exploration/genetic_annotation_for_each_sample.R on Server
library(data.table)
library(dplyr)
library(GenomicRanges)
library(methylKit)
library(genomation)
library(DMRcate)
#BiocManager::install("DMRcate")

setwd("/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/validation_cohort/methylkit/methylkitAnalysis/methykit_input")
file_list = list.files()
file_list <- grep(".txt", file_list, value = T)
name_list = gsub('.methylkit.txt','',file_list)#remove all the suffix(additional words)
name = as.vector(unlist(name_list))# turn the list into vector

genetic_proportion_df <- data.frame()
genetic_count_df <- data.frame()

for (j in 1:length(file_list)) {
  data = fread(file_list[j])
  data <- data[, c("chr", "base")]
  data$end <- data$base
  data$start <- data$base
  #Annotate
  annot <- "/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/files/annot.hg38.ncbirefseq.12col.bed" #Annotation file
  gene.obj <- readTranscriptFeatures(annot,
                                     up.flank=1000,
                                     down.flank=1000,
                                     remove.unusual = TRUE)
  #Convert methylation data to GRanges object
  ann.obj <- annotateWithGeneParts(as(data,"GRanges"), gene.obj)
  genetic_count_df <- rbind(genetic_count_df, ann.obj@num.annotation)
  genetic_proportion_df <- rbind(genetic_proportion_df, genomation::getTargetAnnotationStats(ann.obj,percentage=TRUE,precedence=TRUE))
}
setwd("/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/validation_cohort/methylkit/methylkitAnalysis/general_exploration/output")
names(genetic_count_df) <- c("promoter", "exon", "intron", "intergenic")
genetic_count_df <- data.frame(genetic_count_df, name = name)
names(genetic_proportion_df) <- c("promoter", "exon", "intron", "intergenic")
genetic_proportion_df <- data.frame(genetic_proportion_df, name = name)
fwrite(genetic_count_df, "each_sample_genetic_annotation_count.csv")
fwrite(genetic_proportion_df, "each_sample_genetic_annotation_proportion.csv")

# 1.1.2 Count plot----
# setwd("/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/validation_cohort/methylkit/methylkitAnalysis/general_exploration/output")
setwd("D:\\OneDrive - Washington University in St. Louis\\WashU_Graduate_Study\\Maher_Lab\\ProstateCancer\\output\\Validation_120samples\\generall_exploration")
genetic_count_df <- fread("each_sample_genetic_annotation_count.csv")
genetic_proportion_df <- fread("each_sample_genetic_annotation_proportion.csv")

genetic_count_df_long <- reshape::melt(genetic_count_df, id.vars = "name", variable_name = "roles")
#genetic_count_df_long <- genetic_count_df_long[genetic_count_df_long$roles != "count_all",]
genetic_count_df_long$roles <- factor(genetic_count_df_long$roles, levels = c("promoter", "exon", "intron", "intergenic"))# to adjust the order in the barplot

genetic_count_plot = ggplot(genetic_count_df_long, aes(x = name, y = value, fill = roles)) +
  # ????????????????????????position?????????????????????dodge???????????????????????????????????????
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "annotation",
                    values = c("#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2"),
                    labels = c("Promoter", "Exon", "Intron", "Intergenic"))+ # adjust the legend
  
  theme_bw() + 
  theme(#legend.position="none", #???????????????
    axis.text.x=element_text(#amily="Times",
      colour="black",size=9, angle = 90, hjust = 1), #??????x??????????????????????????????
    axis.text.y=element_text(#family="Times",
      size=12,face="plain"), #??????x??????????????????????????????
    axis.title.y=element_text(#family="Times",
      size = 14,face="bold"), #??????y???????????????????????????
    axis.title.x=element_text(#family="Times",
      size = 14,face="bold"), #??????x???????????????????????????
    legend.title = element_text(#family="Times",
      size = 14,face="bold"), # Adjust legend title size here
    legend.text = element_text(#family="Times",
      size = 12,face="plain"), # Adjust legend keys (levels) text size here
    plot.title = element_text(#family="Times",
      size=10,face="bold",hjust = 0.5), #??????????????????????????????
    panel.grid.major = element_blank(), #??????????????????
    panel.grid.minor = element_blank(),
    plot.background = element_blank(), # No background
    panel.border = element_blank(), # No border
    axis.line = element_line(colour = "black") # Add axis lines
  )+#ylim(0,30000000)
  scale_y_continuous(labels = label_comma(), limits = c(0,30000000))+ # Cancel scientific notation on Y axis
  labs(x = "name", y = "Counts", title = " ") #??????x??????y????????????

# 1.1.3 Peoportion plot----
genetic_count_df
genetic_proportion_df

genetic_proportion_df_long <- reshape::melt(genetic_proportion_df, id.vars = "name", variable_name = "roles")
#genetic_proportion_df_long <- genetic_proportion_df_long[genetic_proportion_df_long$roles != "count_all",]
genetic_proportion_df_long$roles <- factor(genetic_proportion_df_long$roles, levels = c("promoter", "exon", "intron", "intergenic"))# to adjust the order in the barplot

genetic_proportion_plot = ggplot(genetic_proportion_df_long, aes(x = name, y = value, fill = roles)) +
  # ????????????????????????position?????????????????????dodge???????????????????????????????????????
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "annotation",
                    values = c("#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2"),
                    labels = c("Promoter", "Exon", "Intron", "Intergenic"))+ # adjust the legend
  theme_bw() + 
  theme(#legend.position="none", #???????????????
    axis.text.x=element_text(#amily="Times",
      colour="black",size=9, angle = 90, hjust = 1), #??????x??????????????????????????????
    axis.text.y=element_text(#family="Times",
      size=12,face="plain"), #??????x??????????????????????????????
    axis.title.y=element_text(#family="Times",
      size = 14,face="bold"), #??????y???????????????????????????
    axis.title.x=element_text(#family="Times",
      size = 14,face="bold"), #??????x???????????????????????????
    legend.title = element_text(#family="Times",
      size = 14,face="bold"), # Adjust legend title size here
    legend.text = element_text(#family="Times",
      size = 12,face="plain"), # Adjust legend keys (levels) text size here
    plot.title = element_text(#family="Times",
      size=10,face="bold",hjust = 0.5), #??????????????????????????????
    panel.grid.major = element_blank(), #??????????????????
    panel.grid.minor = element_blank(),
    plot.background = element_blank(), # No background
    panel.border = element_blank(), # No border
    axis.line = element_line(colour = "black") # Add axis lines
  )+
  scale_y_continuous(labels = label_comma())+ # Cancel scientific notation on Y axis
  labs(x = "name", y = "Proportion", title = " ") #??????x??????y????????????


library(patchwork)
genetic_combined_plot <- (genetic_count_plot / genetic_proportion_plot) + 
  plot_layout(guides = 'collect')  
#&theme(legend.position = "bottom") # Move legend to bottom
ggsave("each_sample_genetic_annotation_Counts&Proportion.pdf", genetic_combined_plot, height = 5, width = 16.5)



# 1.2 CpG annotation
# 1.2.1 CpG annotation for each sample ----
# This could be run by /storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/validation_cohort/methylkit/methylkitAnalysis/general_exploration/cpg_annotation_for_each_sample.R on Server or locally
# Access the CoG annotation information from USCS
library(data.table)
library(dplyr)
library(GenomicRanges)
library(methylKit)
library(genomation)

island <- "/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/validation_cohort/methylkit/methylkitAnalysis/general_exploration/CpG_annotation_hg38_from_UCSC.bed"
gr_cpg_islands <- rtracklayer::import(island)
gr_cpg_shores <- GRanges(seqnames = seqnames(gr_cpg_islands),
                         ranges = IRanges(start = pmax(start(gr_cpg_islands) - 2000, 1),
                                          end = end(gr_cpg_islands) + 2000))
gr_cpg_shores <- GenomicRanges::reduce(gr_cpg_shores)
gr_cpg_shores <- GenomicRanges::setdiff(gr_cpg_shores, gr_cpg_islands)
gr_cpg_shelves <- GRanges(seqnames = seqnames(gr_cpg_shores),
                          ranges = IRanges(start = pmax(start(gr_cpg_shores) - 2000, 1),
                                           end = end(gr_cpg_shores) + 2000))
gr_cpg_shelves <- GenomicRanges::reduce(gr_cpg_shelves)
gr_cpg_shelves <- GenomicRanges::setdiff(gr_cpg_shelves, gr_cpg_shores)
gr_cpg_shelves <- GenomicRanges::setdiff(gr_cpg_shelves, gr_cpg_islands)


CpG_anno_df <- data.frame()
setwd("/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/validation_cohort/methylkit/methylkitAnalysis/methykit_input")
file_list = list.files()
file_list <- grep(".txt", file_list, value = T)
name_list = gsub('.methylkit.txt','',file_list)#remove all the suffix(additional words)
name = as.vector(unlist(name_list))# turn the list into vector

CpG_anno_df <- data.frame()
for (j in 1:length(file_list)) {
  data = fread(file_list[j])
  data <- data[, c("chr", "base")]
  data$end <- data$base
  data$start <- data$base
  
  # Convert your data frame to a GRanges object
  # Make sure to convert the start and end columns to numeric if they're not already
  gr_data <- GRanges(
    seqnames = data$chr,
    ranges = IRanges(start = as.numeric(data$start), end = as.numeric(data$end))
  )
  
  # Calculate overlaps
  overlap_islands <- countOverlaps(gr_data, gr_cpg_islands)
  overlap_shores <- countOverlaps(gr_data, gr_cpg_shores)
  overlap_shelves <- countOverlaps(gr_data, gr_cpg_shelves)
  
  # Count CpG sites in each category
  count_islands <- sum(overlap_islands > 0)
  count_shores <- sum(overlap_shores > 0)
  count_shelves <- sum(overlap_shelves > 0)
  
  # Print the counts
  CpG_anno <- c(count_islands, count_shores, count_shelves, nrow(data))
  CpG_anno_df <- rbind(CpG_anno_df, CpG_anno)
  print(file_list[j])
}
names(CpG_anno_df) <- c("count_islands", "count_shores", "count_shelves", "count_all") # correct the variable names
CpG_anno_df$count_other <- CpG_anno_df$count_all - (CpG_anno_df$count_islands + CpG_anno_df$count_shores + CpG_anno_df$count_shelves)
CpG_anno_df <- data.frame(CpG_anno_df, name = name)
setwd("/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/validation_cohort/methylkit/methylkitAnalysis/general_exploration/output")
fwrite(CpG_anno_df, "CpG_annotation_count.csv")



# 1.2.2 count plot -----
# setwd("/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/validation_cohort/methylkit/methylkitAnalysis/general_exploration/output")
setwd("D:\\OneDrive - Washington University in St. Louis\\WashU_Graduate_Study\\Maher_Lab\\ProstateCancer\\output\\Validation_120samples\\generall_exploration")

CpG_anno_df = fread("CpG_annotation_count.csv")
CpG_anno_df_long <- reshape::melt(CpG_anno_df, id.vars = "name", variable_name = "roles")
count_df_long <- CpG_anno_df_long[CpG_anno_df_long$roles != "count_all",]
count_df_long$roles <- factor(count_df_long$roles, levels = c("count_islands", "count_shores", "count_shelves", "count_other"))# to adjust the order in the barplot


count_plot = ggplot(count_df_long, aes(x = name, y = value, fill = roles)) +
  # ????????????????????????position?????????????????????dodge???????????????????????????????????????
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "annotation",
                    values = c("#2878B5", "#9AC9DB", "#F8AC8C", "#C82423"),
                    labels = c("Islands","Shores", "Shelves", "Other"))+ # adjust the legend
  theme_bw() + 
  theme(#legend.position="none", #???????????????
    axis.text.x=element_text(#amily="Times",
      colour="black",size=9, angle = 90, hjust = 1), #??????x??????????????????????????????
    axis.text.y=element_text(#family="Times",
      size=12,face="plain"), #??????x??????????????????????????????
    axis.title.y=element_text(#family="Times",
      size = 14,face="bold"), #??????y???????????????????????????
    axis.title.x=element_text(#family="Times",
      size = 14,face="bold"), #??????x???????????????????????????
    legend.title = element_text(#family="Times",
      size = 14,face="bold"), # Adjust legend title size here
    legend.text = element_text(#family="Times",
      size = 12,face="plain"), # Adjust legend keys (levels) text size here
    plot.title = element_text(#family="Times",
      size=10,face="bold",hjust = 0.5), #??????????????????????????????
    panel.grid.major = element_blank(), #??????????????????
    panel.grid.minor = element_blank(),
    plot.background = element_blank(), # No background
    panel.border = element_blank(), # No border
    axis.line = element_line(colour = "black") # Add axis lines
  )+
  scale_y_continuous(labels = label_comma(), limits = c(0,30000000))+ # Cancel scientific notation on Y axis
  labs(x = "name", y = "Counts", title = " ") #??????x??????y????????????


# 1.2.3 Proportion plot ----
CpG_anno_df
proportion_df <- CpG_anno_df
proportion_df$count_other <- proportion_df$count_other/proportion_df$count_all *100
proportion_df$count_islands <- proportion_df$count_islands/proportion_df$count_all*100
proportion_df$count_shores <- proportion_df$count_shores/proportion_df$count_all*100
proportion_df$count_shelves<- proportion_df$count_shelves/proportion_df$count_all*100
proportion_df_long <- reshape::melt(proportion_df, id.vars = "name", variable_name = "roles")
proportion_df_long <- proportion_df_long[proportion_df_long$roles != "count_all",]
proportion_df_long$roles <- factor(proportion_df_long$roles, levels = c("count_islands", "count_shores", "count_shelves", "count_other"))# to adjust the order in the barplot
fwrite(proportion_df, "each_sample_CpG_annotation_proportion.csv")

proportion_plot = ggplot(proportion_df_long, aes(x = name, y = value, fill = roles)) +
  # ????????????????????????position?????????????????????dodge???????????????????????????????????????
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "annotation",
                    values = c("#2878B5", "#9AC9DB", "#F8AC8C", "#C82423"),
                    labels = c("Islands","Shores", "Shelves", "Other"))+ # adjust the legend
  theme_bw() + 
  theme(#legend.position="none", #???????????????
    axis.text.x=element_text(#amily="Times",
      colour="black",size=9, angle = 90, hjust = 1), #??????x??????????????????????????????
    axis.text.y=element_text(#family="Times",
      size=12,face="plain"), #??????x??????????????????????????????
    axis.title.y=element_text(#family="Times",
      size = 14,face="bold"), #??????y???????????????????????????
    axis.title.x=element_text(#family="Times",
      size = 14,face="bold"), #??????x???????????????????????????
    legend.title = element_text(#family="Times",
      size = 14,face="bold"), # Adjust legend title size here
    legend.text = element_text(#family="Times",
      size = 12,face="plain"), # Adjust legend keys (levels) text size here
    plot.title = element_text(#family="Times",
      size=10,face="bold",hjust = 0.5), #??????????????????????????????
    panel.grid.major = element_blank(), #??????????????????
    panel.grid.minor = element_blank(),
    plot.background = element_blank(), # No background
    panel.border = element_blank(), # No border
    axis.line = element_line(colour = "black") # Add axis lines
  )+
  scale_y_continuous(labels = label_comma())+ # Cancel scientific notation on Y axis
  labs(x = "name", y = "Proportion", title = " ") #??????x??????y????????????

# Method2
combined_plot <- (count_plot / proportion_plot) + 
  plot_layout(guides = 'collect')  
#&theme(legend.position = "bottom") # Move legend to bottom
ggsave("each_sample_CpG_annotation_Counts&Proportion.pdf", combined_plot, height = 5, width = 16.5)



#********************************* ----
# DMR_signature_in_120_samples <- fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/Validation_120samples/Validation_of_signature_from_34_samples/sianature_annotaion/DMR_signature_in_120_samples.csv")
# str(DMR_signature_in_120_samples)
# DMR_signature_in_120_samples$patient <- substr(DMR_signature_in_120_samples$patient, 1, 5)
# setdiff(clinicals_for_120$id, DMR_signature_in_120_samples$patient)[order(setdiff(clinicals_for_120$id, DMR_signature_in_120_samples$patient))]

# import clinical outcome infor----
clinicals_for_120 <- fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/data/120validationsamples/120sample_clinical.csv")
clinicals_for_120 <- clinicals_for_120[if_seq == 1,]
# good outcome
set.seed(19970604)
good_120_cohort <- clinicals_for_120[clinicals_for_120$dangerous == 0,]
good_120_cohort$id
# str(good_120_cohort)

# bad outcome
set.seed(19970604)
bad_120_cohort <- clinicals_for_120[clinicals_for_120$dangerous == 1,]
bad_120_cohort$id
# str(bad_120_cohort)

#define the cohort ----
discovery_cohort <- rbind(good_for_discov, bad_for_discov)
validation_cohort <- rbind(good_for_validation, bad_for_validation)
# discovery_cohort$id
#********************************* ----
library(dplyr)
clinicals_for_120 <- clinicals_for_120 %>%
  mutate(risk_level = case_when(
    dangerous == 1 ~ "high-risk",
    dangerous == 0 ~ "low-risk",
    TRUE ~ "intermediate-risk" # 或者 is.na(dangerous) ~ NA_character_, 来处理NA
  ))

library(dplyr)
clinicals_processed <- clinicals_for_120 %>%
  mutate(
    # 处理id列：将首字母"N"替换为"Patient "
    id_processed = gsub("^N", "Patient ", id),
    
    # 处理Bad_outcome列：将1/0替换为描述性文本
    outcome_processed = ifelse(Bad_outcome == 1, "tumor progression", "no tumor progression")
  ) %>%
  # 合并三列（处理后的id、risk-level和outcome），用逗号分隔
  tidyr::unite("combined_column", 
        c("id_processed", "risk_level", "outcome_processed"), 
        sep = ", ", 
        remove = FALSE)  # remove=FALSE保留原始列，设为TRUE则不保留
fwrite(clinicals_processed, "clinicals_processed.csv")

#********************************* ----
BiocManager::install(c("bsseq", "DMRcate"))



# <<<<< DMR annotation >>>>> -----

library(data.table)       # For data handling and file I/O
library(tidyr)            # For data manipulation (e.g., unite)
library(dplyr)            # For data manipulation (e.g., right_join)
library(GenomicRanges)    # For handling genomic ranges
library(IRanges)          # For handling IRanges (required for GenomicRanges)
library(bsseq)            # For handling bisulfite sequencing data
library(BiocParallel)     # For parallel processing in bsseq operations
library(methylKit)        # For handling methylation data
library(DMRcate)          # For DMR annotation and visualization
library(ggplot2)          # Optional: For additional visualizations
library(scales)
library(RColorBrewer)
#rm(list = ls())
# setwd("/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/discovery_cohort_of_120_for_dangerous/general_exploration/output/DMR_summary_by_different_strategy")
dir = "/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/discovery_cohort_of_120_for_dangerous/methylkit_input"
dbdir = "/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/discovery_cohort_of_120_for_dangerous/methylkit_DB"

setwd("H:/discovery_cohort_of_120_for_dangerous/general_exploration/output/DMR_summary_by_different_strategy")
dir = "H:/discovery_cohort_of_120_for_dangerous/methylkit_input"
dbdir = "H:/discovery_cohort_of_120_for_dangerous/methylkit_DB"


all_files <- list.files(path = dir, pattern = ".methylkit.txt", all.files = TRUE, full.names = TRUE)

good_outcome_sample <- c("N1088", "N0354", "N0525", "N0039", "N0133", "N1247", "N0225", "N0313", "N1167", "N0446",
                         "N0009", "N0031", "N0079", "N0095", "N0099", "N0257", "N0468", "N1032", "N0429", "N0326")
bad_outcome_sample <- c("N0974", "N1344", "N0381", "N1343", "N1340", "N0305", "N1796", "N1439", "N1011", "N1445",
                        "N1207", "N1679", "N0323", "N0357", "N1250", "N1132", "N1501", "N0641", "N0771", "N0676",
                        "N0788", "N1197", "N0366", "N0954", "N1920")

high_files <- all_files[grep(paste(bad_outcome_sample, collapse = "|"), all_files)]
high_names <- gsub(".methylkit.txt", "", basename(high_files))
low_files <- all_files[grep(paste(good_outcome_sample, collapse = "|"), all_files)]
low_names <- gsub(".methylkit.txt", "", basename(low_files))

bsseq_obj <- readRDS(paste0(dbdir, "/bsseq_obj.filtered.united5.rds"))
str(bsseq_obj)
# bsseq_obj_smoothed <- BSmooth(bsseq_obj,
#                              ns = 35,
#                              h = 500,
#                              maxGap = 10^8,
#                              verbose = TRUE,
#                              BPPARAM = SerialParam())

# 1.1 DSS DMR annotation----
DSS_DMS <- readRDS(paste0(dbdir,"/DSS_DMS.Pthresh0.01.smoothing500.notSmoothed.united5.rds"))
DSS_DMR <- readRDS(paste0(dbdir,"/DSS_DMR.Pthresh0.01.minlen50.minCG3.dismerge100.pctsig0.5.smoothing500.notSmoothed.united5.rds"))
DSS_DMR$abs_mean_diff <- abs(DSS_DMR$diff.Methy)
DSS_DMR_sig <- DSS_DMR[DSS_DMR$abs_mean_diff > 0.10, ]
DSS_DMR_sig <- DSS_DMR_sig %>% tidyr::unite("DMR", chr:end, remove = F)
DSS_DMR_sig_GR <- makeGRangesFromDataFrame(DSS_DMR_sig)

file_names <- gsub(".methylkit.txt", "", basename(all_files))
color = rep("red", length(file_names))
color[grep(paste(good_outcome_sample, collapse = "|"), file_names)] <- "green"
names(color) = ifelse(color == "red", "Aggressive", "Indolent")
# DMR.plot(DSS_DMR_sig_GR, 
#          dmr = which(DSS_DMR_sig$abs_mean_diff == max(DSS_DMR_sig$abs_mean_diff)), 
#          CpGs = bsseq_obj, 
#          phen.col = color, 
#          genome = "hg38", 
#          flank = 1000)

# 1.1.1 extract methy matrix for significant DSS DMR ----
DSS_DMR_sig_methy_matrix <- getMeth(bsseq_obj, regions = DSS_DMR_sig_GR, type = "raw", what = "perRegion")
sum(is.na(DSS_DMR_sig_methy_matrix))
DSS_DMR_sig_methy_matrix <- data.table(DSS_DMR_sig_methy_matrix)
DSS_DMR_sig_methy_matrix[] <- lapply(DSS_DMR_sig_methy_matrix, function(x) as.numeric(trimws(x)))
DSS_DMR_sig_methy_matrix <- cbind(DMR = DSS_DMR_sig$DMR, DSS_DMR_sig_methy_matrix)
DSS_DMR_sig_methy_matrix[DSS_DMR_sig_methy_matrix == NaN] = NA
fwrite(DSS_DMR_sig_methy_matrix, "DSS_DMR_sig_methy_matrix.csv")
str(DSS_DMR_sig_methy_matrix)

rowSums(is.na(DSS_DMR_sig_methy_matrix))
DSS_DMR_sig_methy_matrix_noNA <- DSS_DMR_sig_methy_matrix[rowSums(is.na(DSS_DMR_sig_methy_matrix)) == 0,]
DSS_DMR_sig_methy_matrix_noNA <- right_join(DSS_DMR_sig, DSS_DMR_sig_methy_matrix_noNA,)
DSS_DMR_sig_methy_matrix_noNA_GR <- makeGRangesFromDataFrame(DSS_DMR_sig_methy_matrix_noNA)
DMR.plot(DSS_DMR_sig_methy_matrix_noNA_GR, 
         dmr = which(DSS_DMR_sig_methy_matrix_noNA$abs_mean_diff == max(DSS_DMR_sig_methy_matrix_noNA$abs_mean_diff)), 
         CpGs = bsseq_obj, 
         phen.col = color, 
         genome = "hg38"
         #, flank = 1000
         )

# 1.2 DMRcate DMR annotation----
DMRcate_DMR_lambda500 <- readRDS(paste0(dbdir,"/DMRcate_DMR.lambda1000.C2.smoothing500.notSmoothed.united5.rds"))
DMRcate_DMR_lambda500_df <- data.frame(DMRcate_DMR_lambda500)
DMRcate_DMR_lambda500_df$abs_mean_diff <- abs(DMRcate_DMR_lambda500_df$meandiff)
DMRcate_DMR_lambda500_sig <- DMRcate_DMR_lambda500_df[DMRcate_DMR_lambda500_df$Stouffer < 0.05 & 
                                                        DMRcate_DMR_lambda500_df$HMFDR < 0.05 & 
                                                        DMRcate_DMR_lambda500_df$abs_mean_diff > 0.10, ]
DMRcate_DMR_lambda500_sig <- DMRcate_DMR_lambda500_sig %>% tidyr::unite("DMR", seqnames:end, remove = F)
DMRcate_DMR_lambda500_sig_GR <- makeGRangesFromDataFrame(DMRcate_DMR_lambda500_sig)

file_names <- gsub(".methylkit.txt", "", basename(all_files))
color = rep("red", length(file_names))
color[grep(paste(good_outcome_sample, collapse = "|"), file_names)] <- "green"
names(color) = ifelse(color == "red", "Aggressive", "Indolent")

# DMR.plot(DMRcate_DMR_lambda500_sig_GR, 
#          dmr = which(DMRcate_DMR_lambda500_sig$abs_mean_diff == max(DMRcate_DMR_lambda500_sig$abs_mean_diff)), 
#          CpGs = bsseq_obj, 
#          phen.col = color, 
#          genome = "hg38", 
#          flank = 1000)

# 1.2.1 extract methy matrix for significant DMRcate DMR ----
DMRcate_DMR_lambda500_sig_methy_matrix <- getMeth(bsseq_obj, regions = DMRcate_DMR_lambda500_sig_GR, type = "raw", what = "perRegion")
sum(is.na(DSS_DMR_sig_methy_matrix))
DMRcate_DMR_lambda500_sig_methy_matrix <- data.table(DMRcate_DMR_lambda500_sig_methy_matrix)
DMRcate_DMR_lambda500_sig_methy_matrix[] <- lapply(DMRcate_DMR_lambda500_sig_methy_matrix, function(x) as.numeric(trimws(x)))
DMRcate_DMR_lambda500_sig_methy_matrix <- cbind(DMR = DMRcate_DMR_lambda500_sig$DMR, DMRcate_DMR_lambda500_sig_methy_matrix)
DMRcate_DMR_lambda500_sig_methy_matrix[DMRcate_DMR_lambda500_sig_methy_matrix == NaN] = NA
fwrite(DMRcate_DMR_lambda500_sig_methy_matrix, "DMRcate_DMR_lambda500_sig_methy_matrix.csv")
str(DMRcate_DMR_lambda500_sig_methy_matrix)

rowSums(is.na(DMRcate_DMR_lambda500_sig_methy_matrix))
DMRcate_DMR_lambda500_sig_methy_matrix_noNA <- DMRcate_DMR_lambda500_sig_methy_matrix[rowSums(is.na(DMRcate_DMR_lambda500_sig_methy_matrix)) == 0,]
DMRcate_DMR_lambda500_sig_methy_matrix_noNA <- right_join(DMRcate_DMR_lambda500_sig, DMRcate_DMR_lambda500_sig_methy_matrix_noNA,)
DMRcate_DMR_lambda500_sig_methy_matrix_noNA_GR <- makeGRangesFromDataFrame(DMRcate_DMR_lambda500_sig_methy_matrix_noNA)
DMR.plot(DMRcate_DMR_lambda500_sig_methy_matrix_noNA_GR, 
         dmr = which(DMRcate_DMR_lambda500_sig_methy_matrix_noNA$abs_mean_diff == max(DMRcate_DMR_lambda500_sig_methy_matrix_noNA$abs_mean_diff)), 
         CpGs = bsseq_obj, 
         phen.col = color, 
         genome = "hg38"
         #, flank = 1000
         )


# 1.3 find overlapping DMR between DSS and DMRcate ----

# 1.3.1 Method 1 of extract all DMR: merge all DMR first and then select the significant one ----

# From all DMR between DSS and DMRcate ---- 
# 1.3.1.1 check if overlap ----
DSS_DMR_GR <- makeGRangesFromDataFrame(DSS_DMR)
DMRcate_DMR_lambda500_GR <-  makeGRangesFromDataFrame(DMRcate_DMR_lambda500)
DSS_DMRcate500_overlaps <- findOverlaps(DSS_DMR_GR, DMRcate_DMR_lambda500_GR)
DSS_DMR$overlap <- ifelse(seq_along(DSS_DMR_GR) %in% queryHits(DSS_DMRcate500_overlaps), TRUE, FALSE)
DMRcate_DMR_lambda500$overlap <- ifelse(seq_along(DMRcate_DMR_lambda500_GR) %in% subjectHits(DSS_DMRcate500_overlaps), TRUE, FALSE)
sum(DSS_DMR$overlap)
sum(DMRcate_DMR_lambda500$overlap)

# 1.3.1.2 merge all DMR from DSS and DMRcate ----
# Separate overlapping and non-overlapping DMRs from each data frame
DSS_overlap <- DSS_DMR[DSS_DMR$overlap, ]
DSS_non_overlap <- DSS_DMR[!DSS_DMR$overlap, ]
DMRcate_overlap <- DMRcate_DMR_lambda500[DMRcate_DMR_lambda500$overlap, ]
DMRcate_non_overlap <- DMRcate_DMR_lambda500[!DMRcate_DMR_lambda500$overlap, ]
# Create unified data frame for overlapping DMRs
# Merge overlapping DMRs based on common regions, using the start and end of the overlapping regions
overlap_final <- data.frame(
  chr = seqnames(DMRcate_DMR_lambda500_GR[subjectHits(DSS_DMRcate500_overlaps)]),
  start = start(DMRcate_DMR_lambda500_GR[subjectHits(DSS_DMRcate500_overlaps)]),
  end = end(DMRcate_DMR_lambda500_GR[subjectHits(DSS_DMRcate500_overlaps)]),
  source = "Overlap"
)
# Create data frames for non-overlapping DMRs from both sources
DSS_non_overlap$source <- "DSS"
DMRcate_non_overlap$source <- "DMRcate"
# Merge all DMRs together: overlapping and non-overlapping
DMRcate_non_overlap_df = data.frame(DMRcate_non_overlap)
names(DMRcate_non_overlap_df)[names(DMRcate_non_overlap_df) == "seqnames"] <- "chr"
merged_DMRs <- rbind(overlap_final,
                     DSS_non_overlap[, c("chr", "start", "end", "source")],
                     DMRcate_non_overlap_df[, c("chr", "start", "end", "source")])
merged_DMRs$width <- merged_DMRs$end - merged_DMRs$start  # Adjust column name after merging
# Print summary of the merged data
table(merged_DMRs$source)
# Save the merged data frame to a file

str(merged_DMRs)

#1.3.1.3 Make a venn plot for DMR source----
source_counts <- table(merged_DMRs$source)
DSS_total <- source_counts["DSS"] + source_counts["Overlap"]
DMRcate_total <- source_counts["DMRcate"] + source_counts["Overlap"]
overlap_count <- source_counts["Overlap"]

library(ggVennDiagram)
library(ggplot2)

# 构造列表数据
venn_list <- list(
  DSS = which(merged_DMRs$source %in% c("DSS", "Overlap")),
  DMRcate = which(merged_DMRs$source %in% c("DMRcate", "Overlap"))
)

# 自定义颜色设置
light_blue <- "#ADD8E6"  # 浅蓝色
light_green <- "#E8D0A9"  # 浅绿色

# 绘制韦恩图
p <- ggVennDiagram(
  venn_list,
  label_alpha = 0,
  label_geom = c("label"),
  category.names = c("DSS", "DMRcate"),
  edge_lty = "solid"
  ,set_color = c(light_green, light_blue)  # 设置集合边框颜色
) +
  # 使用渐变填充
  scale_fill_gradient(
    low = light_blue,     # 低值区域用浅蓝色
    high = light_green,   # 高值区域用浅绿色
    guide = "none"        # 隐藏填充色图例
  ) +
  labs(title = "") +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    text = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 14)
  ) +
  coord_fixed(ratio = 1)

# 查看图形
print(p)

# 保存图形
ggsave("DSS_DMRcate_venn_colored.pdf", plot = p, height = 5, width = 5)

# 1.3.1.3 Calculate the new mean difference ----
# extract methy matrix
merged_DMRs <- merged_DMRs %>% tidyr::unite("DMR", chr:end, remove = F)
merged_DMRs_GR <- makeGRangesFromDataFrame(merged_DMRs)
merged_DMRs_methy_matrix <- getMeth(bsseq_obj, regions = merged_DMRs_GR, type = "raw", what = "perRegion")
sum(is.na(merged_DMRs_methy_matrix))
merged_DMRs_methy_matrix <- data.table(merged_DMRs_methy_matrix)
merged_DMRs_methy_matrix[] <- lapply(merged_DMRs_methy_matrix, function(x) as.numeric(trimws(x)))
merged_DMRs_methy_matrix <- cbind(DMR = merged_DMRs$DMR, merged_DMRs_methy_matrix)
merged_DMRs_methy_matrix[merged_DMRs_methy_matrix == NaN] = NA
str(merged_DMRs_methy_matrix)

# Calculate the difference in methylation between High and Low groups
merged_DMRs_methy_matrix$High_Avg <- rowMeans(merged_DMRs_methy_matrix[, high_names, with = FALSE], na.rm = TRUE)
merged_DMRs_methy_matrix$Low_Avg <- rowMeans(merged_DMRs_methy_matrix[, low_names, with = FALSE], na.rm = TRUE)
merged_DMRs_methy_matrix$Differential <- abs(merged_DMRs_methy_matrix$High_Avg - merged_DMRs_methy_matrix$Low_Avg)
merged_DMRs <- left_join(merged_DMRs, merged_DMRs_methy_matrix[,c("DMR", "High_Avg", "Low_Avg", "Differential")])
merged_DMRs <- unique(merged_DMRs)
str(merged_DMRs)
merged_DMRs_sig <- merged_DMRs[merged_DMRs$Differential > 0.10, ]
fwrite(merged_DMRs, "merged_DMRs.csv")
fwrite(merged_DMRs_sig, "merged_DMRs_sig.csv")

# 1.3.1.4 methy matrix of merged DMR in 45 discovory cohort ----
merged_DMRs_sig <- merged_DMRs_sig %>% tidyr::unite("DMR", chr:end, remove = F)
#bsseq_obj <- readRDS(paste0(dbdir, "/bsseq_obj.filtered.united5.rds"))
merged_DMRs_sig_GR <- makeGRangesFromDataFrame(merged_DMRs_sig)
merged_DMRs_sig_methy_matrix <- getMeth(bsseq_obj, regions = merged_DMRs_sig_GR, type = "raw", what = "perRegion")
merged_DMRs_sig_methy_matrix <- data.table(merged_DMRs_sig_methy_matrix)
merged_DMRs_sig_methy_matrix[] <- lapply(merged_DMRs_sig_methy_matrix, function(x) as.numeric(trimws(x)))
merged_DMRs_sig_methy_matrix <- cbind(DMR = merged_DMRs_sig$DMR, merged_DMRs_sig_methy_matrix)
merged_DMRs_sig_methy_matrix[merged_DMRs_sig_methy_matrix == NaN] = NA
str(merged_DMRs_sig_methy_matrix)
fwrite(merged_DMRs_sig_methy_matrix, "merged_DMRs_sig_methy_matrix.csv")

# *** ----
#generate methylation distributino in each chromosome ***----
# Calculate average methylation for high and low groups
setwd("/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/discovery_cohort_of_120_for_dangerous/general_exploration/output/DMR_summary_by_different_strategy")
dir = "/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/discovery_cohort_of_120_for_dangerous/methylkit_input"
all_files <- list.files(path = dir, pattern = ".methylkit.txt", all.files = TRUE, full.names = TRUE)
good_outcome_sample <- c("N1088", "N0354", "N0525", "N0039", "N0133", "N1247", "N0225", "N0313", "N1167", "N0446",
                         "N0009", "N0031", "N0079", "N0095", "N0099", "N0257", "N0468", "N1032", "N0429", "N0326")
bad_outcome_sample <- c("N0974", "N1344", "N0381", "N1343", "N1340", "N0305", "N1796", "N1439", "N1011", "N1445",
                        "N1207", "N1679", "N0323", "N0357", "N1250", "N1132", "N1501", "N0641", "N0771", "N0676",
                        "N0788", "N1197", "N0366", "N0954", "N1920")
high_files <- all_files[grep(paste(bad_outcome_sample, collapse = "|"), all_files)]
high_names <- gsub(".methylkit.txt", "", basename(high_files))
low_files <- all_files[grep(paste(good_outcome_sample, collapse = "|"), all_files)]
low_names <- gsub(".methylkit.txt", "", basename(low_files))

merged_DMRs_sig_methy_matrix <- fread("merged_DMRs_sig_methy_matrix.csv")
DMR <- merged_DMRs_sig_methy_matrix$DMR
merged_DMRs_sig_methy_matrix_methy <- merged_DMRs_sig_methy_matrix[,-1]
merged_DMRs_sig_methy_matrix_methy <- t(scale(t(merged_DMRs_sig_methy_matrix_methy)))
merged_DMRs_sig_methy_matrix <- cbind(DMR, data.table(merged_DMRs_sig_methy_matrix_methy))
merged_DMRs_sig_methy_matrix[, High_Avg := rowMeans(.SD, na.rm = TRUE), .SDcols = high_names]
merged_DMRs_sig_methy_matrix[, Low_Avg := rowMeans(.SD, na.rm = TRUE), .SDcols = low_names]

# Extract chromosome name and position if not already extracted
if (!"Chromosome" %in% colnames(merged_DMRs_sig_methy_matrix)) {
  merged_DMRs_sig_methy_matrix[, `:=`(
    Chromosome = sub("_.*", "", DMR),
    Start = as.numeric(sub(".*_(\\d+)_\\d+$", "\\1", DMR)),
    End = as.numeric(sub(".*_(\\d+)$", "\\1", DMR)) 
  )]
}
merged_DMRs_sig_methy_matrix <- merged_DMRs_sig_methy_matrix[merged_DMRs_sig_methy_matrix$Chromosome != "chrM",]

# Create a summary data frame for chromosome lengths
chromosome_lengths <- merged_DMRs_sig_methy_matrix[, .(Max_Position = max(End)), by = Chromosome]# Y chromosome's length is 57.8Mb  

# Define approximate centromere positions for demonstration purposes (replace with actual values if available)
centromere_positions <- data.frame(
  Chromosome = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                 "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"),
  Centromere_Position = c(124.2, 93.4, 91.7, 50.9, 47.7, 60.5, 58.9, 45.2, 50.6, 40.3,  #chr1-10
                          52.9, 35.4, 16, 15.6, 17, 38.2, 22.2, 16.1, 28.5, 27.1, 12.3, 11.8, #chr11-22
                          59.4, 12) #chrx-y
)

# Merge centromere positions with chromosome lengths
chromosome_lengths <- merge(chromosome_lengths, centromere_positions, by = "Chromosome")

# Ensure factor levels are set correctly
merged_DMRs_sig_methy_matrix$Chromosome <- factor(merged_DMRs_sig_methy_matrix$Chromosome, levels = chromosome_lengths$Chromosome)
chromosome_lengths$Chromosome <- factor(chromosome_lengths$Chromosome, levels = chromosome_lengths$Chromosome)

# Prepare the data for plotting by melting (reshaping)
merged_DMRs_sig_methy_matrix_long <- melt(
  merged_DMRs_sig_methy_matrix, 
  id.vars = c("Chromosome", "Start", "End"),  # Only keep Chromosome, Start, and End as id variables
  measure.vars = c("High_Avg", "Low_Avg"),  # Variables to melt
  variable.name = "Group",  # Rename melted variable column to Group
  value.name = "Methylation"  # Rename melted value column to Methylation
)
str(merged_DMRs_sig_methy_matrix_long)

# Create a factor with the desired order
correct_order <- c(paste0("chr", 1:22), "chrX", "chrY")

# Apply this order to the Chromosome variable in both merged_DMRs_sig_methy_matrix and chromosome_lengths
merged_DMRs_sig_methy_matrix$Chromosome <- factor(merged_DMRs_sig_methy_matrix$Chromosome, levels = correct_order)
merged_DMRs_sig_methy_matrix_long$Chromosome <- factor(merged_DMRs_sig_methy_matrix_long$Chromosome, levels = correct_order)
chromosome_lengths$Chromosome <- factor(chromosome_lengths$Chromosome, levels = correct_order)

# Assign x-axis positions for each group using the correct column name ("variable" instead of "Group")
merged_DMRs_sig_methy_matrix_long$Chromosome_Group <- as.numeric(merged_DMRs_sig_methy_matrix_long$Chromosome) + 
  ifelse(merged_DMRs_sig_methy_matrix_long$Group == "High_Avg", -0.2, 0.2)


# Create the plot
# Updated plot using geom_segment() to make the lines thicker instead of just wider
merged_DMRs_sig_methy_matrix_long$Start <- merged_DMRs_sig_methy_matrix_long$Start - 1000
merged_DMRs_sig_methy_matrix_long$end <- merged_DMRs_sig_methy_matrix_long$end + 1000

COLS <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
ggplot() +
  # Draw chromosome shapes as background rectangles using chromosome_lengths data
  geom_rect(data = chromosome_lengths, aes(xmin = as.numeric(Chromosome) - 0.4,
                                           xmax = as.numeric(Chromosome) + 0.4,
                                           ymin = 0, ymax = Max_Position / 10^6), fill = NA, color = "grey", alpha = 0.3) +
  
  # Plot methylation values for high and low groups using segments for thicker lines
  geom_segment(data = merged_DMRs_sig_methy_matrix_long, 
               aes(x = Chromosome_Group, xend = Chromosome_Group, 
                   y = Start / 10^6, yend = End / 10^6, color = Methylation), linewidth = 5, alpha = 1) +  # Use linewidth to control thickness
  
  # Draw centromere representation (Yellow lines)
  # geom_point(data = chromosome_lengths,
  #              aes(x = as.numeric(Chromosome), y = Centromere_Position ),
  #              color = "grey", size = 2, ) +  # Increase linewidth for visibility
  geom_segment(data = chromosome_lengths,
               aes(x = as.numeric(Chromosome), xend = as.numeric(Chromosome),
                   y = Centromere_Position - 2, yend = Centromere_Position + 2),
               color = "grey", linewidth = 10) +  # Increase linewidth for visibility
  
  
  # Correct Chromosome Order
  scale_x_continuous(
    breaks = 1:length(correct_order), 
    labels = sub("chr", "", correct_order)  # Remove 'chr' prefix from chromosome labels
  ) +
  scale_y_continuous(name = "Genomic Position (Mb)") +
  # Manually define legend range and color scale for methylation values
  scale_colour_gradientn(colours = rev(COLS(100)[c(1:40, 60:100)]),
                         limits = c(-1, 1),# Set the range of the legend from -1 to 1
                         breaks = seq(-1, 1, by = 0.5)) +
  #scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Methylation Z-score") +  # Color scale for methylation
  labs(
    title = "",
    x = "Chromosome", y = "Genomic Position (Mb)", color = "Methylation Z-score"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, color = "black"),  # Increase font size and change color to black
    axis.text.y = element_text(size = 12, color = "black"),  # Increase font size and change color to black
    axis.title = element_text(size = 14, color = "black"),  # Increase axis title size and color
    panel.border = element_blank(),
    axis.line = element_line(size = 0.2, color = "black"),  # Make axes lines thicker and black
    panel.grid = element_blank(),  # Remove grid lines
    # axis.text.x=element_text(colour="black",family="Times",size=14), #设置x轴刻度标签的字体属性
    # axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
    # axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴的标题的字体属性
    # axis.title.x=element_text(family="Times",size = 20,face="plain"), #设置x轴的标题的字体属性
    
    plot.title = element_text(hjust = 0.5, size = 20, color = "black")  # Increase title size and color
  )
ggsave("DMR_distribution_on_each_Chromosome.pdf", width = 15, height = 4)


# generate a cover percentage of each chromosome ----
# Calculate the total DMR length for each chromosome
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
merged_DMRs_sig_methy_matrix <- fread("merged_DMRs_sig_methy_matrix.csv")
DMR <- merged_DMRs_sig_methy_matrix$DMR
merged_DMRs_sig_methy_matrix_methy <- merged_DMRs_sig_methy_matrix[,-1]
merged_DMRs_sig_methy_matrix_methy <- t(scale(t(merged_DMRs_sig_methy_matrix_methy)))
merged_DMRs_sig_methy_matrix <- cbind(DMR, data.table(merged_DMRs_sig_methy_matrix_methy))

# Extract chromosome name and position if not already extracted
if (!"Chromosome" %in% colnames(merged_DMRs_sig_methy_matrix)) {
  merged_DMRs_sig_methy_matrix[, `:=`(
    Chromosome = sub("_.*", "", DMR),
    Start = as.numeric(sub(".*_(\\d+)_\\d+$", "\\1", DMR)),
    End = as.numeric(sub(".*_(\\d+)$", "\\1", DMR)) 
  )]
}

dmr_lengths <- merged_DMRs_sig_methy_matrix %>%
  group_by(Chromosome) %>%
  summarise(total_length = sum(End - Start, na.rm = TRUE)) %>%
  arrange(Chromosome)

dmr_lengths <- full_join(dmr_lengths, centromere_positions)
dmr_lengths$percent <- dmr_lengths$total_length / 10^6 / dmr_lengths$Centromere_Position * 100

correct_order <- c(paste0("chr", 1:22), "chrX", "chrY")
dmr_lengths$Chromosome <- factor(dmr_lengths$Chromosome, levels = correct_order)
ggplot(data = dmr_lengths) +
  geom_col(aes(x = Chromosome, y = percent))+
  geom_col(aes(x = Chromosome, y = total_length))

# Normalize 'percent' to match the range of 'total_length'
dmr_lengths <- dmr_lengths %>%
  mutate(
    percent_scaled = percent * max(total_length) / max(percent)  # Scale percent to match total_length
  )

# Create a combined data frame for both stats
dmr_long <- dmr_lengths %>%
  pivot_longer(cols = c(total_length, percent_scaled), names_to = "Statistic", values_to = "Value")

# Create the plot
ggplot(dmr_long, aes(x = Chromosome, y = Value, fill = Statistic)) +
  # Use geom_col() instead of geom_bar() for better control over position
  geom_col(position = position_dodge2(preserve = "single", width = 0.8), width = 0.7) +
  # Left Y-axis for total_length
  scale_y_continuous(
    name = "Total Length (bp)",
    labels = scales::comma,  # Format Y-axis labels without scientific notation
    sec.axis = sec_axis(~ . * max(dmr_lengths$percent) / max(dmr_lengths$total_length), name = "Percent (%)")  # Secondary Y-axis for percent
  ) +
  # Custom colors for the bars
  scale_fill_manual(values = c("total_length" = "skyblue", "percent_scaled" = "orange")) +
  # Revise the X-axis label
  scale_x_discrete(
    labels = sub("chr", "", levels(dmr_lengths$Chromosome))
  ) +
  # Labels and theme adjustments
  labs(title = "DMR Total Length and Percent by Chromosome", x = "Chromosome", fill = "Statistics") +
  theme_bw() +
  theme(
    legend.position="none", #不需要图例
    axis.title.y.left = element_text(color = "skyblue", size = 14),
    axis.text.y.left = element_text(color = "skyblue", size = 12),
    axis.line.y.left = element_line(color = "skyblue",size = 0.2),
    axis.title.y.right = element_text(color = "orange", size = 14),
    axis.text.y.right = element_text(color = "orange", size = 12),
    axis.line.y.right = element_line(color = "orange",size = 0.2),
    axis.line.x = element_line(color = "black",size = 0.2),
    axis.title.x = element_text(color = "black", size = 14),
    axis.text.x = element_text(color = "black", size = 12),
    panel.border = element_blank(),
    # axis.line = element_line(size = 0.2, color = "black"),  # Make axes lines thicker and black
    panel.grid = element_blank(),  # Remove grid lines
    plot.title = element_text(hjust = 0.5, size = 20, color = "black")  # Increase title size and color
  )
ggsave("DMR_total_length_and_pecentage_on_each_Chromosome.pdf", width = 15, height = 4)

# generate DMR length distribution ----
merged_DMRs_sig_methy_matrix$width <- merged_DMRs_sig_methy_matrix$End - merged_DMRs_sig_methy_matrix$Start
custom_log_trans <- trans_new(
  name = "custom_log",
  transform = function(x) log10(x + 1), # Add 1 to avoid log of zero issues
  inverse = function(x) 10^x - 1, # Subtract 1 to revert the transformation correctly
  breaks = function(x) c(0, 1, 5, 10, 20, 50, 100), # Custom breaks for desired spacing
  domain = c(0, Inf)
)
ggplot(merged_DMRs_sig_methy_matrix, aes(x = width)) +
  geom_histogram(binwidth = 1000, fill = "skyblue", color = "black", alpha = 0.7) +  # Adjust binwidth as needed
  labs(
    title = "",
    x = "DMR Length (bp)",
    y = "Frequency"
  ) +
  scale_y_continuous(
    trans = custom_log_trans,  # Apply custom log transformation
    breaks = c(0, 10, 100, 1000, 4000),  # Define specific breaks for Y-axis
    labels = c("0", "10", "100", "1000", "4000"),  # Custom labels for the Y-axis
    limits = c(0, 5000)
  ) +
  scale_x_continuous(
    # trans = custom_log_trans,  # Apply custom log transformation
    # breaks = c(0, 100, 5000, 30000),  # Define specific breaks for Y-axis
    # labels = c("0", "100", "5000", "30000"),  # Custom labels for the Y-axis
    n.break = 6,
    limits = c(0, 30000)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, color = "black"),  # Increase font size and change color to black
    axis.text.y = element_text(size = 12, color = "black"),  # Increase font size and change color to black
    axis.title = element_text(size = 14, color = "black"),  # Increase axis title size and color
    panel.border = element_blank(),
    axis.line = element_line(size = 0.2, color = "black"),  # Make axes lines thicker and black
    panel.grid = element_blank(),  # Remove grid lines
    # axis.text.x=element_text(colour="black",family="Times",size=14), #设置x轴刻度标签的字体属性
    # axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
    # axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴的标题的字体属性
    # axis.title.x=element_text(family="Times",size = 20,face="plain"), #设置x轴的标题的字体属性
    
    plot.title = element_text(hjust = 0.5, size = 20, color = "black")  # Increase title size and color
  )

ggsave("DMR_length_distribution.pdf", width = 5, height = 4)

# ***----




# # (Optional) 1.3.2 Method 2 of extract all DMR: select the significant DMR first and then merge ----
# # From significant DMR between DSS and DMRcate ---- 
# # 1.3.2.1 check if overlap ----
# DSS_DMRcate500_sig_overlaps <- findOverlaps(DSS_DMR_sig_GR, DMRcate_DMR_lambda500_sig_GR)
# DSS_DMR_sig$overlap <- ifelse(seq_along(DSS_DMR_sig_GR) %in% queryHits(DSS_DMRcate500_sig_overlaps), TRUE, FALSE)
# DMRcate_DMR_lambda500_sig$overlap <- ifelse(seq_along(DMRcate_DMR_lambda500_sig_GR) %in% subjectHits(DSS_DMRcate500_sig_overlaps), TRUE, FALSE)
# sum(DSS_DMR_sig$overlap)
# sum(DMRcate_DMR_lambda500_sig$overlap)
# 
# # 1.3.2.2 check overlapping percentage ----
# intersect_ranges <- pintersect(DSS_DMR_sig_GR[queryHits(DSS_DMRcate500_sig_overlaps)], 
#                                DMRcate_DMR_lambda500_sig_GR[subjectHits(DSS_DMRcate500_sig_overlaps)])
# overlap_length <- width(intersect_ranges)
# # Calculate the percentage overlap relative to each DMR's length
# DSS_overlap_percentage <- overlap_length / width(DSS_DMR_sig_GR[queryHits(DSS_DMRcate500_sig_overlaps)]) * 100
# DMRcate_overlap_percentage <- overlap_length / width(DMRcate_DMR_lambda500_sig_GR[subjectHits(DSS_DMRcate500_sig_overlaps)]) * 100
# # Add overlap percentage information to the original data frames
# DSS_DMR_sig$overlap_percentage <- NA
# DMRcate_DMR_lambda500_sig$overlap_percentage <- NA
# # Fill in the overlap percentages for the overlapping DMRs
# DSS_DMR_sig$overlap_percentage[queryHits(DSS_DMRcate500_sig_overlaps)] <- DSS_overlap_percentage
# DMRcate_DMR_lambda500_sig$overlap_percentage[subjectHits(DSS_DMRcate500_sig_overlaps)] <- DMRcate_overlap_percentage
# 
# # DMR.plot(GRanges(seqnames = "chr6",
# #                  ranges = IRanges(start = 99602312, end = 99606241)), 
# #          dmr = 1, 
# #          CpGs = bsseq_obj, 
# #          phen.col = color, 
# #          genome = "hg38", 
# #          flank = 3000)
# 
# 
# # 1.3.2.3 merge significant DMR from DSS and DMRcate ----
# # Separate overlapping and non-overlapping DMRs from each data frame
# DSS_sig_overlap <- DSS_DMR_sig[DSS_DMR_sig$overlap, ]
# DSS_sig_non_overlap <- DSS_DMR_sig[!DSS_DMR_sig$overlap, ]
# DMRcate_sig_overlap <- DMRcate_DMR_lambda500_sig[DMRcate_DMR_lambda500_sig$overlap, ]
# DMRcate_sig_non_overlap <- DMRcate_DMR_lambda500_sig[!DMRcate_DMR_lambda500_sig$overlap, ]
# # Create unified data frame for overlapping DMRs
# # The DMRcate version of overlapping DMRs will be retained and DSS version removed
# # Create a new data frame for the final set of merged DMRs, retaining:
# # - All non-overlapping DSS DMRs (if needed)
# # - All non-overlapping DMRcate DMRs
# # - Only DMRcate's version of overlapping DMRs
# 
# # sig_overlap_merged <- data.frame(
# #   chr = seqnames(DSS_DMR_sig_GR[queryHits(DSS_DMRcate500_sig_overlaps)]),
# #   start = pmin(start(DSS_DMR_sig_GR[queryHits(DSS_DMRcate500_sig_overlaps)]), start(DMRcate_DMR_lambda500_sig_GR[subjectHits(DSS_DMRcate500_sig_overlaps)])),
# #   end = pmax(end(DSS_DMR_sig_GR[queryHits(DSS_DMRcate500_sig_overlaps)]), end(DMRcate_DMR_lambda500_sig_GR[subjectHits(DSS_DMRcate500_sig_overlaps)])),
# #   source = "Overlap"
# # )
# sig_overlap_final <- data.frame(
#   chr = seqnames(DMRcate_DMR_lambda500_sig_GR[subjectHits(DSS_DMRcate500_sig_overlaps)]),
#   start = start(DMRcate_DMR_lambda500_sig_GR[subjectHits(DSS_DMRcate500_sig_overlaps)]),
#   end = end(DMRcate_DMR_lambda500_sig_GR[subjectHits(DSS_DMRcate500_sig_overlaps)]),
#   source = "Overlap"
# )
# # Create data frames for non-overlapping DMRs from both sources
# DSS_sig_non_overlap$source <- "DSS"
# DMRcate_sig_non_overlap$source <- "DMRcate"
# # Merge all DMRs together: overlapping and non-overlapping
# names(DMRcate_sig_non_overlap)[names(DMRcate_sig_non_overlap) == "seqnames"] <- "chr"
# sig_merged_DMRs <- rbind(sig_overlap_merged,
#                          DSS_sig_non_overlap[, c("chr", "start", "end", "source")],
#                          DMRcate_sig_non_overlap[, c("chr", "start", "end", "source")])
# sig_merged_DMRs$width <- sig_merged_DMRs$end - sig_merged_DMRs$start  # Adjust column name after merging
# # Print summary of the merged data
# table(sig_merged_DMRs$source)
# str(sig_merged_DMRs)
# 
# 
# 
# # 1.3.2.4 Calculate the new mean difference ----
# # extract methy matrix
# sig_merged_DMRs <- sig_merged_DMRs %>% tidyr::unite("DMR", chr:end, remove = F)
# sig_merged_DMRs_GR <- makeGRangesFromDataFrame(sig_merged_DMRs)
# sig_merged_DMRs_methy_matrix <- getMeth(bsseq_obj, regions = sig_merged_DMRs_GR, type = "raw", what = "perRegion")
# sum(is.na(sig_merged_DMRs_methy_matrix))
# sig_merged_DMRs_methy_matrix <- data.table(sig_merged_DMRs_methy_matrix)
# sig_merged_DMRs_methy_matrix[] <- lapply(sig_merged_DMRs_methy_matrix, function(x) as.numeric(trimws(x)))
# sig_merged_DMRs_methy_matrix <- cbind(DMR = sig_merged_DMRs$DMR, sig_merged_DMRs_methy_matrix)
# sig_merged_DMRs_methy_matrix[sig_merged_DMRs_methy_matrix == NaN] = NA
# str(sig_merged_DMRs_methy_matrix)
# 
# # Calculate the difference in methylation between High and Low groups
# sig_merged_DMRs_methy_matrix$High_Avg <- rowMeans(sig_merged_DMRs_methy_matrix[, high_names, with = FALSE], na.rm = TRUE)
# sig_merged_DMRs_methy_matrix$Low_Avg <- rowMeans(sig_merged_DMRs_methy_matrix[, low_names, with = FALSE], na.rm = TRUE)
# sig_merged_DMRs_methy_matrix$Differential <- abs(sig_merged_DMRs_methy_matrix$High_Avg - sig_merged_DMRs_methy_matrix$Low_Avg)
# sig_merged_DMRs <- left_join(sig_merged_DMRs, sig_merged_DMRs_methy_matrix[,c("DMR", "High_Avg", "Low_Avg", "Differential")])
# sig_merged_DMRs <- unique(sig_merged_DMRs)
# sig_merged_DMRs <- sig_merged_DMRs[sig_merged_DMRs$Differential > 0.10,]
# str(sig_merged_DMRs)
# fwrite(sig_merged_DMRs, "sig_merged_DMRs.csv")
# 
# DMR.plot(GRanges(seqnames = "chr1",
#                  ranges = IRanges(start = 50415084, end = 50420000)), 
#          dmr = 1, 
#          CpGs = bsseq_obj, 
#          phen.col = color, 
#          genome = "hg38", 
#          flank = 2000)



#***********************----
# 2.0 Annotation of merged signidicant-DMR on all 120 samples ----
# 2.1 extract the methy matrix of DMR ----
# 2.1.0 create_bsseq_obj_for_all120_sample.R ---- check /rdcw/fs2/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/discovery_cohort_of_120_for_dangerous/general_exploration/script/create_bsseq_obj_for_all120_sample.R
library(data.table)    # For reading and manipulating data
library(methylKit)     # For handling methylation data
library(bsseq)         # For BSseq object creation and smoothing
library(GenomicRanges) # For handling genomic ranges
library(IRanges)       # For handling IRanges operations
#BiocManager::install("Gviz")
library(Gviz)
library(DMRcate)
library(biomaRt)

setwd("/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/discovery_cohort_of_120_for_dangerous/general_exploration/output/DMR_summary_by_different_strategy")
# setwd("H:/discovery_cohort_of_120_for_dangerous/general_exploration/output/DMR_summary_by_different_strategy")



do_1_unite_cpg = TRUE
do_2_bsseq_obj = TRUE

dir = "/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/validation_cohort/methylkit/methylkitAnalysis/methykit_input"
dbdir = "/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/discovery_cohort_of_120_for_dangerous/general_exploration/output/DMR_summary_by_different_strategy"
all_files <- list.files(path = dir, pattern = ".methylkit.txt", all.files = TRUE, full.names = TRUE)

good_outcome_sample <- c("N1088", "N0354", "N0525", "N0039", "N0133", "N1247", "N0225", "N0313", "N1167", "N0446",
                         "N0009", "N0031", "N0079", "N0095", "N0099", "N0257", "N0468", "N1032", "N0429", "N0326")
bad_outcome_sample <- c("N0974", "N1344", "N0381", "N1343", "N1340", "N0305", "N1796", "N1439", "N1011", "N1445",
                        "N1207", "N1679", "N0323", "N0357", "N1250", "N1132", "N1501", "N0641", "N0771", "N0676",
                        "N0788", "N1197", "N0366", "N0954", "N1920")

high_files <- all_files[grep(paste(bad_outcome_sample, collapse = "|"), all_files)]
high_names <- gsub(".methylkit.txt", "", basename(high_files))
low_files <- all_files[grep(paste(good_outcome_sample, collapse = "|"), all_files)]
low_names <- gsub(".methylkit.txt", "", basename(low_files))


if (do_1_unite_cpg) {
  message("*** Extract coverage and methylated counts for all samples ***")
  DSS_datalist <- list()
  for (i in 1:length(all_files)) {
    print(paste0("Reformatting and storing sample: ", basename(all_files)[i]))
    # Read the methylation data
    methylkit_data <- fread(file = all_files[i])
    # Transform data to DSS-required format
    CpGdata <- transform(methylkit_data[, 2:3],
                         N = methylkit_data[, 5],  # Coverage
                         X = round((methylkit_data[, 6] / 100) * methylkit_data[, 5], digits = 0))  # Methylated counts
    colnames(CpGdata) <- c("chr", "pos", "N", "X")
    # Store each sample's data in DSS_datalist
    DSS_datalist[[i]] <- CpGdata
  }
  #saveRDS(DSS_datalist, paste0(dbdir, "/DSS_datalist.coverage.united1_all120.rds"))
} else {
  DSS_datalist <- readRDS(paste0(dbdir, "/DSS_datalist.coverage.united1_all120.rds"))
}

if (do_2_bsseq_obj) {
  message("# Modified function to filter CpG sites based on sample coverage")
  filter_CpG_sites <- function(DSS_datalist, sample_names, group1, group2, min_coverage = 5) {
    # Initialize combined_data with the first sample
    combined_data <- copy(DSS_datalist[[1]])
    setnames(combined_data, c("N", "X"), paste0(c("N", "X"), "_", 1))
    combined_data[, coverage_count := 1]  # Initialize coverage count with 1 for the first sample
    
    # Loop through the rest of the samples, rename columns, and merge
    for (i in 2:length(DSS_datalist)) {
      # Rename columns to avoid duplication, using underscore instead of dot
      temp <- copy(DSS_datalist[[i]])
      setnames(temp, c("N", "X"), paste0(c("N", "X"), "_", i))
      # Merge on chr and pos while keeping track of coverage count
      combined_data <- merge(combined_data, temp, by = c("chr", "pos"), all = TRUE)
      # Update coverage count for non-NA values in the current sample
      combined_data[, coverage_count := coverage_count + !is.na(get(paste0("N_", i)))]
      # Print combined_data column names for debugging
      print(paste("Columns after merging sample", i))
    }
    
    # Calculate coverage for each group using underscores
    group1_indices <- which(sample_names %in% group1)
    group2_indices <- which(sample_names %in% group2)
    # Generate the correct column names for .SDcols using underscores
    group1_colnames <- paste0("N_", group1_indices)
    group2_colnames <- paste0("N_", group2_indices)
    # Print expected column names for debugging
    print(paste("Group 1 columns:", paste(group1_colnames, collapse = ", ")))
    print(paste("Group 2 columns:", paste(group2_colnames, collapse = ", ")))
    
    # Check if all required columns exist in the combined_data
    missing_cols <- setdiff(c(group1_colnames, group2_colnames), names(combined_data))
    if (length(missing_cols) > 0) {
      stop("The following columns are missing from combined_data: ", paste(missing_cols, collapse = ", "))
    }
    
    # Create logical vectors for each group based on coverage
    combined_data[, group1_coverage := rowSums(.SD >= min_coverage, na.rm = TRUE), .SDcols = group1_colnames]
    combined_data[, group2_coverage := rowSums(.SD >= min_coverage, na.rm = TRUE), .SDcols = group2_colnames]
    # Keep only rows (CpG sites) where both groups have at least 5 samples with coverage
    filtered_data <- combined_data[group1_coverage >= 5 & group2_coverage >= 5]
    
    return(filtered_data)
  }
  # Extract sample names
  file_names <- gsub(".methylkit.txt", "", basename(all_files))
  
  # Filter CpG sites based on coverage in each group
  message("*** Filter CpG sites based on sample coverage ***")
  filtered_data <- filter_CpG_sites(DSS_datalist, file_names, high_names, low_names, min_coverage = 5)
  saveRDS(filtered_data, paste0(dbdir, "/filtered_data_filtered.filtered.united5_all120.rds"))
  
  message("*** Create a list of filtered DSS data for each sample ***")
  DSS_datalist_filtered <- lapply(seq_along(DSS_datalist), function(i) {
    print(paste0("filter DSS datalist: ", i))
    merge(filtered_data[, .(chr, pos)], DSS_datalist[[i]], by = c("chr", "pos"), all = FALSE)
  })
  #saveRDS(DSS_datalist_filtered, paste0(dbdir, "/DSS_datalist_filtered.filtered.united5_all120.rds"))
  
  # Create the BSseq object with filtered data
  message("*** Create the BSseq object ***")
  bsseq_obj <- makeBSseqData(DSS_datalist_filtered, sampleNames = file_names)
  saveRDS(bsseq_obj, paste0(dbdir, "/bsseq_obj.filtered.united5_all120.rds"))
  print(bsseq_obj)
  
  # bsseq_obj_smoothed <- BSmooth(bsseq_obj,
  #                               ns = 35,
  #                               h = 500,
  #                               maxGap = 10^8,
  #                               verbose = TRUE,
  #                               BPPARAM = MulticoreParam(workers = 6))
  # saveRDS(bsseq_obj_smoothed, paste0(dbdir,"/bsseq_obj.Smoothed.united1.rds"))
  # print(bsseq_obj_smoothed)
  rm(DSS_datalist, DSS_datalist_filtered, filtered_data)
  
} else {
  bsseq_obj <- readRDS(paste0(dbdir, "/bsseq_obj.filtered.united5_all120.rds"))
  file_names <- gsub(".methylkit.txt", "", basename(all_files))
}


# 2.1.1 significant DMRs ----
bsseq_obj_all120 <- readRDS("bsseq_obj.filtered.united5_all120.rds")
# All_sig_DMR <- fread("sig_merged_DMRs.csv") # select the sig onee first and then merge
All_sig_DMR <- fread("merged_DMRs_sig.csv") # merge first and then select the sig ones
str(All_sig_DMR)
table(All_sig_DMR$source)
All_sig_DMR <- All_sig_DMR %>% tidyr::unite("DMR", chr:end, remove = F)
All_sig_DMR_GR <- makeGRangesFromDataFrame(All_sig_DMR)
# with NA
All_sig_DMR_methy_matrix <- getMeth(bsseq_obj_all120, regions = All_sig_DMR_GR, type = "raw", what = "perRegion")
All_sig_DMR_methy_matrix <- data.table(All_sig_DMR_methy_matrix)
All_sig_DMR_methy_matrix[] <- lapply(All_sig_DMR_methy_matrix, function(x) as.numeric(trimws(x)))
All_sig_DMR_methy_matrix <- cbind(DMR = All_sig_DMR$DMR, All_sig_DMR_methy_matrix)
All_sig_DMR_methy_matrix[All_sig_DMR_methy_matrix == NaN] = NA
str(All_sig_DMR_methy_matrix)
fwrite(All_sig_DMR_methy_matrix, "All_merged_DMRs_sig_methy_matrix.csv")
# without NA
All_sig_DMR_methy_matrix_noNA <- All_sig_DMR_methy_matrix[rowSums(is.na(All_sig_DMR_methy_matrix)) == 0,]
All_sig_DMR_methy_matrix_noNA <- right_join(All_sig_DMR, All_sig_DMR_methy_matrix_noNA,)
str(All_sig_DMR_methy_matrix_noNA)
fwrite(All_sig_DMR_methy_matrix_noNA, "All_merged_DMRs_sig_methy_matrix_noNA.csv")


# # (Optional) 2.1.2 simplely all DMRs ----
# # All_sig_DMR <- fread("sig_merged_DMRs.csv")
# All_DMR <- fread("merged_DMRs.csv")
# str(All_DMR)
# table(All_DMR$source)
# All_DMR <- All_DMR %>% tidyr::unite("DMR", chr:end, remove = F)
# All_DMR_GR <- makeGRangesFromDataFrame(All_DMR)
# # with NA
# All_DMR_methy_matrix <- getMeth(bsseq_obj_all120, regions = All_DMR_GR, type = "raw", what = "perRegion")
# All_DMR_methy_matrix <- data.table(All_DMR_methy_matrix)
# All_DMR_methy_matrix[] <- lapply(All_DMR_methy_matrix, function(x) as.numeric(trimws(x)))
# All_DMR_methy_matrix <- cbind(DMR = All_DMR$DMR, All_DMR_methy_matrix)
# All_DMR_methy_matrix[All_DMR_methy_matrix == NaN] = NA
# str(All_DMR_methy_matrix)
# fwrite(All_DMR_methy_matrix, "All_DMR_methy_matrix.csv")
# # without NA
# All_DMR_methy_matrix_noNA <- All_DMR_methy_matrix[rowSums(is.na(All_DMR_methy_matrix)) == 0,]
# All_DMR_methy_matrix_noNA <- right_join(All_DMR, All_DMR_methy_matrix_noNA,)
# str(All_DMR_methy_matrix_noNA)
# fwrite(All_DMR_methy_matrix_noNA, "All_DMR_methy_matrix_noNA.csv")


# quick annotation 1 ----
# Genome annotation of typical DMR ----
file_names_all120 <- gsub(".methylkit.txt", "", 
                          list.files(path = "/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/validation_cohort/methylkit/methylkitAnalysis/methykit_input", #""
                                     pattern = ".methylkit.txt", all.files = TRUE, full.names = FALSE))
file_names_all120

good_outcome_sample <- c("N1088", "N0354", "N0525", "N0039", "N0133", "N1247", "N0225", "N0313", "N1167", "N0446",
                         "N0009", "N0031", "N0079", "N0095", "N0099", "N0257", "N0468", "N1032", "N0429", "N0326")
bad_outcome_sample <- c("N0974", "N1344", "N0381", "N1343", "N1340", "N0305", "N1796", "N1439", "N1011", "N1445",
                        "N1207", "N1679", "N0323", "N0357", "N1250", "N1132", "N1501", "N0641", "N0771", "N0676",
                        "N0788", "N1197", "N0366", "N0954", "N1920")

color = rep("goldenrod1", length(file_names_all120))
color[grep(paste(good_outcome_sample, collapse = "|"), file_names_all120)] <- "#5b8efd"
color[grep(paste(bad_outcome_sample, collapse = "|"), file_names_all120)] <- "#dd217d"
names(color) = ifelse(color == "#dd217d", "High-risky", "Low-risky")
names(color)[color == "goldenrod1"] = "Middle-risky"

table(color)




library(GenomicRanges)
library(IRanges)
library(AnnotationHub)
library(biomaRt)
library(minfi)
library(bsseq)
library(Gviz)
library(RColorBrewer)
library(S4Vectors)
library(GenomeInfoDb)
setwd("/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/discovery_cohort_of_120_for_dangerous/general_exploration/output/DMR_summary_by_different_strategy")
bsseq_obj_all120 <- readRDS("bsseq_obj.filtered.united5_all120.rds")

DMR_plot <- function(ranges, 
                     dmr, 
                     CpGs, 
                     what = c("Beta", "M"), 
                     arraytype = c("EPICv2", "EPICv1", "450K"),
                     phen.col,
                     genome = c("hg19", "hg38", "mm10"),
                     labels = names(ranges),
                     flank = 5000,
                     heatmap = TRUE,
                     extra.ranges = NULL, 
                     extra.title = names(extra.ranges),
                     # 修改轨道大小设置：热图轨道可以按组设置
                     track_sizes = list(
                       ideogram = 2,      # 染色体意象图轨道高度
                       genome_axis = 2,   # 基因组坐标轴轨道高度
                       gene = 7,          # 基因注释轨道高度
                       cpgs = 2,          # CpG位点轨道高度
                       dmrs = 3,          # DMR区域轨道高度
                       extras = 3,        # 额外区域轨道高度
                       heatmap = NULL,    # 热图轨道高度（可以是向量，每个元素对应一个组）
                       mean = 6           # 均值轨道高度
                     ),
                     # 标题背景不透明度
                     title_alpha = 1,
                     # 标题颜色参数
                     title_color = "black",
                     # 新增：标题和图例字体大小
                     title_cex = 0.8,     # 标题字体大小
                     legend_cex = 0.7,    # 图例字体大小
                     # 是否显示热图图例及其设置
                     show_heatmap_legend = TRUE,
                     legend_title = "甲基化水平",
                     # 自定义热图颜色
                     heatmap_colors = c("cornflowerblue", "white", "darksalmon")) 
{
  # 参数检查与初始化
  what <- match.arg(what)
  arraytype <- match.arg(arraytype)
  genome <- match.arg(genome)
  stopifnot(class(CpGs)[1] %in% c("matrix", "BSseq", "GenomicRatioSet"))
  
  # 初始化 trackstoplot 作为空列表
  trackstoplot <- list()
  
  if(arraytype == "EPICv2" & genome == "hg19"){
    stop("Error: genome must be hg38 for EPICv2 data.")
  }
  
  if(flank < 10 | flank > 10000){
    stop("Error: DMR flanking region needs to be between 10bp and 10000bp")
  }
  
  # 检查title_alpha参数
  if(title_alpha < 0 | title_alpha > 1){
    stop("Error: title_alpha must be between 0 and 1")
  }
  
  stopifnot(dmr %in% 1:length(ranges))
  IDs <- unique(names(phen.col))
  
  # 处理热图轨道大小设置
  if (heatmap) {
    # 如果heatmap为NULL或不是向量，设置默认值
    if (is.null(track_sizes$heatmap)) {
      track_sizes$heatmap <- rep(10, length(IDs))
    } else if (length(track_sizes$heatmap) == 1) {
      # 如果只提供一个值，则对所有组使用相同的大小
      track_sizes$heatmap <- rep(track_sizes$heatmap, length(IDs))
    } else if (length(track_sizes$heatmap) != length(IDs)) {
      # 如果提供的大小与组数不匹配，发出警告并使用默认值
      warning("Length of heatmap track sizes does not match number of groups. Using default sizes.")
      track_sizes$heatmap <- rep(10, length(IDs))
    }
  }
  
  # 处理 CpGs 数据
  if (is(CpGs, "matrix") | is(CpGs, "GenomicRatioSet")) {
    if (is(CpGs, "matrix")) {
      if (arraytype == "450K") {
        grset <- makeGenomicRatioSetFromMatrix(CpGs, 
                                               array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", 
                                               mergeManifest = TRUE, what = what)
      }
      if (arraytype == "EPICv1") {
        grset <- makeGenomicRatioSetFromMatrix(CpGs, 
                                               array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19", 
                                               mergeManifest = TRUE, what = what)
      }
      if (arraytype == "EPICv2") {
        if (!is(CpGs, "matrix")){
          stop("Error, CpGs argument must be a matrix for EPICv2")
        }
      }
    } else {
      grset <- CpGs
    }
    
    # 加载 EPICv2 注释信息
    if (arraytype == "EPICv2") {
      message("Loading EPICv2 manifest...")
      ah <- AnnotationHub()
      EPICv2manifest <- ah[["AH116484"]]
      EPICv2manifest <- EPICv2manifest[EPICv2manifest$CHR != "chr0",]
      RSanno <- GRanges(EPICv2manifest$CHR, IRanges(EPICv2manifest$MAPINFO, EPICv2manifest$MAPINFO))
      names(RSanno) <- rownames(EPICv2manifest)
      RSanno <- sort(sortSeqlevels(RSanno))
      RSanno <- RSanno[names(RSanno) %in% rownames(CpGs)]
      CpGs <- CpGs[names(RSanno),]
      cpgs.ranges <- RSanno
    } else {
      CpGs <- getBeta(grset)
      RSanno <- getAnnotation(grset)
      RSanno <- RSanno[order(RSanno$chr, RSanno$pos), ]
      CpGs <- CpGs[rownames(RSanno), ]
      cpgs.ranges <- GRanges(RSanno$chr, IRanges(RSanno$pos, RSanno$pos))
    }
    values(cpgs.ranges) <- CpGs
    isbsseq <- FALSE
  } else {
    #处理 BSseq 对象
    if (any(width(CpGs) > 1)) {
      stop("Error: all ranges in the BSseq object must be single nucleotides with width 1.")
    }
    if (is.null(rownames(colData(CpGs)))) {
      stop("Error: BSseq object must be annotated with colData with sample IDs as rownames of the data.frame.")
    }
    stopifnot(ncol(CpGs) == length(phen.col))
    isbsseq <- TRUE
  }
  
  # 计算 DMR 区域
  ranges$ID <- rep("", length(ranges))
  ranges.reduce <- reduce(ranges + flank)
  dmrs.inplot <- ranges[queryHits(findOverlaps(ranges, ranges.reduce[subjectHits(findOverlaps(ranges[dmr], 
                                                                                              ranges.reduce))]))]
  ranges.inplot <- ranges.reduce[queryHits(findOverlaps(ranges.reduce, dmrs.inplot))]
  
  # 筛选 CpG 数据
  if (is(CpGs, "matrix")) {
    cpgs.ranges <- subsetByOverlaps(cpgs.ranges, ranges.inplot)
  } else {
    cpgs.ranges <- subsetByOverlaps(CpGs, ranges.inplot)
  }
  
  genome(cpgs.ranges) <- genome
  
  # 计算甲基化比例
  if (isbsseq) {
    methRatios <- GRanges(seqnames(cpgs.ranges), ranges(cpgs.ranges), 
                          mcols = as.matrix(getCoverage(cpgs.ranges, type = "M"))/as.matrix(getCoverage(cpgs.ranges, 
                                                                                                        type = "Cov")))
  } else {
    methRatios <- cpgs.ranges
  }
  values(methRatios) <- as.matrix(values(methRatios))
  colnames(values(methRatios)) <- gsub("mcols.", "", colnames(values(methRatios)))
  
  # 添加基因组注释轨道
  suppressWarnings(switch(genome, hg19 = {
    ensembl <- useEnsembl(host = "https://grch37.ensembl.org", 
                          biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
    grt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = as.character(seqnames(methRatios)[1]), 
                                  start = min(start(ranges.inplot)), end = max(start(ranges.inplot)), 
                                  name = "ENSEMBL", biomart = ensembl)
  }, hg38 = {
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    grt <- BiomartGeneRegionTrack(genome = "hg38", chromosome = as.character(seqnames(methRatios)[1]), 
                                  start = min(start(ranges.inplot)), end = max(start(ranges.inplot)), 
                                  name = "ENSEMBL", biomart = ensembl)
  }, mm10 = {
    ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    grt <- BiomartGeneRegionTrack(genome = "mm10", chromosome = as.character(seqnames(methRatios)[1]), 
                                  start = min(start(ranges.inplot)), end = max(start(ranges.inplot)), 
                                  name = "ENSEMBL", biomart = ensembl)
  }))
  
  # 设置基因组注释轨道的标题和文本参数
  suppressMessages(grt <- setPar(grt, "collapseTranscripts", "meta"))
  suppressMessages(grt <- setPar(grt, "exonAnnotation", "symbol"))
  suppressMessages(grt <- setPar(grt, "fontcolor.item", "black"))
  suppressMessages(grt <- setPar(grt, "cex", 0.6))
  suppressMessages(grt <- setPar(grt, "cex.title", title_cex)) # 设置标题字体大小
  suppressMessages(grt <- setPar(grt, "rotation.item", 45))
  suppressMessages(grt <- setPar(grt, "rotation.title", 0))
  suppressMessages(grt <- setPar(grt, "background.title", "white"))
  suppressMessages(grt <- setPar(grt, "fontcolor.title", title_color)) # 设置标题字体颜色
  
  # 添加 CpG 轨道
  cpgs.track <- AnnotationTrack(
    GRanges(seqnames(cpgs.ranges), ranges(cpgs.ranges)), 
    name = "CpGs",
    fill = "green", 
    stacking = "dense", 
    fontcolor = "black",
    background.title = "white",
    rotation.title = 0
  )
  suppressMessages(cpgs.track <- setPar(cpgs.track, "lty", "blank"))
  suppressMessages(cpgs.track <- setPar(cpgs.track, "fontcolor.title", title_color)) 
  suppressMessages(cpgs.track <- setPar(cpgs.track, "cex.title", title_cex)) # 设置标题字体大小
  
  # 添加 DMR 轨道
  if (all(is.null(names(dmrs.inplot)))) {
    names(dmrs.inplot) <- rep("DMR", length(dmrs.inplot))
  }
  
  dmrs.track <- AnnotationTrack(
    dmrs.inplot, 
    name = "DMRs", 
    showFeatureId = TRUE, 
    fill = "purple", 
    id = names(dmrs.inplot), 
    fontcolor = "black", 
    background.title = "white",
    rotation.title = 0
  )
  suppressMessages(dmrs.track <- setPar(dmrs.track, "fontcolor.title", title_color))
  suppressMessages(dmrs.track <- setPar(dmrs.track, "cex.title", title_cex)) # 设置标题字体大小
  
  # 创建意象图轨道
  ideogram_track <- IdeogramTrack(genome = genome, chromosome = as.character(seqnames(ranges.inplot)))
  suppressMessages(ideogram_track <- setPar(ideogram_track, "fontcolor.title", title_color))
  suppressMessages(ideogram_track <- setPar(ideogram_track, "col.axis", title_color))
  suppressMessages(ideogram_track <- setPar(ideogram_track, "cex.title", title_cex)) # 设置标题字体大小
  
  # 创建基因组坐标轴轨道
  genome_axis_track <- GenomeAxisTrack()
  suppressMessages(genome_axis_track <- setPar(genome_axis_track, "fontcolor.title", title_color))
  suppressMessages(genome_axis_track <- setPar(genome_axis_track, "col.axis", title_color))
  suppressMessages(genome_axis_track <- setPar(genome_axis_track, "cex.title", title_cex)) # 设置标题字体大小
  
  # 添加基本轨道到trackstoplot
  trackstoplot <- c(trackstoplot, ideogram_track, genome_axis_track, grt, cpgs.track, dmrs.track)
  
  # 添加额外区域轨道
  if (!is.null(extra.ranges)) {
    extra.ranges <- extra.ranges[extra.ranges %over% dmrs.inplot]
    extras.track <- AnnotationTrack(extra.ranges, showFeatureId = TRUE, 
                                    name = extra.title, fill = "pink", id = names(extra.ranges), 
                                    rotation.title = 0, background.title = "white")
    suppressMessages(extras.track <- setPar(extras.track, "fontcolor.title", title_color))                            
    suppressMessages(extras.track <- setPar(extras.track, "cex.title", title_cex)) # 设置标题字体大小
    trackstoplot <- c(trackstoplot, extras.track)
  }
  
  # 处理颜色不透明度
  # 将十六进制颜色转换为RGBA颜色带透明度的函数
  add_alpha <- function(color, alpha = 1) {
    if (startsWith(color, "#")) {
      # 如果是十六进制颜色
      rgb <- col2rgb(color)
      return(rgb(rgb[1]/255, rgb[2]/255, rgb[3]/255, alpha))
    } else {
      # 如果是命名颜色
      rgb <- col2rgb(color)
      return(rgb(rgb[1]/255, rgb[2]/255, rgb[3]/255, alpha))
    }
  }
  
  # 绘制热图,自定义热图颜色
  if (heatmap) {
    custom_colors <- colorRampPalette(heatmap_colors)(100)  # 自定义颜色方案
    dt.group <- list()
    for (i in seq_along(IDs)) {
      group_id <- IDs[i]
      # 为每个组应用透明度
      bg_color <- add_alpha(phen.col[group_id], title_alpha)
      dt <- DataTrack(methRatios[,names(phen.col) %in% group_id], 
                      name = group_id, 
                      background.title = bg_color,  # 使用带透明度的颜色
                      type = "heatmap", 
                      showSampleNames = FALSE, 
                      ylim = c(0,1), 
                      genome = genome, 
                      gradient = custom_colors,
                      rotation.title = 0,
                      # 新增：设置热图y轴刻度值仅在0, 0.5, 1.0处显示
                      yTicksAt = c(0, 0.25, 0.5, 0.75, 1),
                      yTicksLabels = c("0", "0.25", "0.5", "0.75", "1"),
                      # yAxisSide = "right",
                      legend = FALSE)  # 禁用内置的legend
      # 设置轨道标题字体颜色以及轴的颜色
      suppressMessages(dt <- setPar(dt, "fontcolor.title", title_color))
      suppressMessages(dt <- setPar(dt, "col.axis", title_color))
      suppressMessages(dt <- setPar(dt, "cex.title", title_cex)) # 设置标题字体大小
      dt.group[[i]] <- dt
    }
    trackstoplot <- c(trackstoplot, dt.group)  # 添加热图轨道
  }
  
  # 添加均值甲基化轨道
  meanmeth <- DataTrack(methRatios, groups = names(phen.col), 
                        type = "smooth",         # 改用平滑曲线
                        span = 0.1,            # 控制平滑程度（值越小越波动）
                        col = phen.col[sort(unique(names(phen.col)))], 
                        ylim = c(0, 1), 
                        name = "Smoothed mean\nmethylation",
                        background.title = "white",
                        rotation.title = 0,
                        na.rm = TRUE,
                        # 新增：设置y轴刻度值仅在0, 0.5, 1.0处显示
                        yTicksAt = c(0, 0.5, 1),
                        yTicksLabels = c("0", "0.5", "1"))
  
  # 设置轨道标题和图例的文本颜色和大小
  suppressMessages(meanmeth <- setPar(meanmeth, "fontcolor.title", title_color))
  suppressMessages(meanmeth <- setPar(meanmeth, "fontcolor.legend", title_color))
  suppressMessages(meanmeth <- setPar(meanmeth, "col.axis", title_color))
  suppressMessages(meanmeth <- setPar(meanmeth, "cex.title", title_cex)) # 设置标题字体大小
  suppressMessages(meanmeth <- setPar(meanmeth, "cex.legend", legend_cex)) # 设置图例字体大小
  suppressMessages(setPar(meanmeth, "groupAnnotation", "feature"))
  trackstoplot <- c(trackstoplot, meanmeth)  # 添加均值轨道
  
  # 设置轨道大小
  sizes <- c(
    track_sizes$ideogram,           # 意象图轨道
    track_sizes$genome_axis,        # 基因组坐标轨道
    track_sizes$gene,               # 基因注释轨道
    track_sizes$cpgs,               # CpG位点轨道
    track_sizes$dmrs                # DMR区域轨道
  )
  
  # 如果有额外区域，添加其大小
  if (!is.null(extra.ranges)) {
    sizes <- c(sizes, track_sizes$extras)
  }
  
  # 如果有热图，根据每个组单独设置大小
  if (heatmap) {
    # 为每个组分别设置大小
    sizes <- c(sizes, track_sizes$heatmap)
  }
  
  # 添加均值轨道大小
  sizes <- c(sizes, track_sizes$mean)
  
  # 绘制轨道
  plot_result <- plotTracks(trackstoplot, 
                            from = min(start(ranges.inplot)), 
                            to = max(end(ranges.inplot)),
                            title.width = 1.5,  # 增加标题宽度，确保能够显示完整标题
                            sizes = sizes,      # 使用自定义轨道大小
                            main = "",          # 确保没有主标题干扰
                            showTitle = TRUE)   # 明确指定显示标题
  
  return(invisible(plot_result))
}



pdf("/home/wg-liao/Genome_annotation_of_DMR.pdf", width = 6, height = 3)
DMR_plot(
  ranges = GRanges(seqnames = "chr9", ranges = IRanges(start = 101433398, end = 101433510)),
  dmr = 1, 
  CpGs = bsseq_obj_all120,  
  phen.col = color,
  genome = "hg38",
  flank = 800,
  track_sizes = list(
    ideogram = 2,
    genome_axis = 3,
    gene = 1,
    cpgs = 1,
    dmrs = 1,
    extras = 3,
    heatmap = c(75/20,20/10,25/10),
    mean = 8
  ),
  title_alpha = 0.7,
  title_color = "black",
  title_cex = 0.7,    # 设置标题字体大小
  legend_cex = 0.7,   # 设置图例字体大小
  heatmap_colors = c("cornflowerblue", "white", "darksalmon")
)
dev.off()



# quick annotation 2 ----
# DMR coverage of EM-seq and 450K/850K ----
setwd("/storage1/fs1/christophermaher/Active/maherlab/l.muheng/methylation/SeverityPrediction/discovery_cohort_of_120_for_dangerous/general_exploration/output/DMR_summary_by_different_strategy")

# 加载必要的包
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) # 450K array注释
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # 850K (EPIC) array注释
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(ggVennDiagram) # 使用ggVennDiagram绘制韦恩图
library(ggplot2)
library(ggforce)

# 载入DMR数据（假设已有All_sig_DMR对象）
# 如果没有，可以这样创建示例
All_sig_DMR <- fread("All_merged_DMRs_sig_methy_matrix_noNA.csv") # merge first and then select the sig ones

# 创建DMR的GRanges对象(hg38)
dmr_gr <- GRanges(
  seqnames = All_sig_DMR$chr,
  ranges = IRanges(start = All_sig_DMR$start, end = All_sig_DMR$end),
  DMR_id = All_sig_DMR$DMR
)

# 下载并导入chain文件进行从hg19到hg38的坐标转换
# 我们将把芯片探针从hg19转换到hg38，而不是把DMR从hg38转换到hg19
chain_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
chain_file_gz <- "hg19ToHg38.over.chain.gz"
chain_file <- "hg19ToHg38.over.chain"

# 下载chain文件（如果本地不存在）
if (!file.exists(chain_file)) {
  cat("下载chain文件...\n")
  download.file(chain_url, chain_file_gz)
  system(paste("gunzip", chain_file_gz))
  cat("chain文件下载完成\n")
}

# 导入chain文件
chain <- import.chain(chain_file)
cat("chain文件导入完成\n")

# 获取450K array探针位置
data(list = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
annot450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probes450k_gr <- GRanges(
  seqnames = annot450k$chr,
  ranges = IRanges(start = annot450k$pos, end = annot450k$pos),
  probe_id = rownames(annot450k)
)

# 获取850K (EPIC) array探针位置
data(list = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
annotEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
probesEPIC_gr <- GRanges(
  seqnames = annotEPIC$chr,
  ranges = IRanges(start = annotEPIC$pos, end = annotEPIC$pos),
  probe_id = rownames(annotEPIC)
)

# 查找DMR与450K和850K array探针的重叠情况
overlap_450k <- findOverlaps(dmr_gr, probes450k_gr)
overlap_EPIC <- findOverlaps(dmr_gr, probesEPIC_gr)

# 统计每个DMR中包含多少个450K和EPIC探针
dmr_450k_counts <- table(queryHits(overlap_450k))
dmr_EPIC_counts <- table(queryHits(overlap_EPIC))

# 初始化覆盖状态
dmr_coverage <- data.frame(
  DMR_id = All_sig_DMR$DMR,
  DMR_index = 1:nrow(All_sig_DMR),
  covered_450k = FALSE,
  covered_EPIC = FALSE
)

# 标记被450K覆盖的DMR
dmr_coverage$covered_450k[as.numeric(names(dmr_450k_counts))] <- TRUE

# 标记被EPIC覆盖的DMR
dmr_coverage$covered_EPIC[as.numeric(names(dmr_EPIC_counts))] <- TRUE

# 统计覆盖情况
not_covered <- sum(!dmr_coverage$covered_450k & !dmr_coverage$covered_EPIC)
only_450k <- sum(dmr_coverage$covered_450k & !dmr_coverage$covered_EPIC)
only_EPIC <- sum(!dmr_coverage$covered_450k & dmr_coverage$covered_EPIC)
both_covered <- sum(dmr_coverage$covered_450k & dmr_coverage$covered_EPIC)

# 打印结果
cat("DMR覆盖情况统计:\n")
cat("完全未覆盖: ", not_covered, " (", round(not_covered/nrow(All_sig_DMR)*100, 2), "%)\n", sep="")
cat("仅被450K覆盖: ", only_450k, " (", round(only_450k/nrow(All_sig_DMR)*100, 2), "%)\n", sep="")
cat("仅被EPIC覆盖: ", only_EPIC, " (", round(only_EPIC/nrow(All_sig_DMR)*100, 2), "%)\n", sep="")
cat("同时被两者覆盖: ", both_covered, " (", round(both_covered/nrow(All_sig_DMR)*100, 2), "%)\n", sep="")
cat("总DMR数量: ", nrow(All_sig_DMR), "\n", sep="")


# 使用ggplot2创建包含大圆的韦恩图


total_dmr = nrow(All_sig_DMR)
# 创建数据框用于绘制韦恩图
# 首先计算各圆的中心位置
# 大圆(EM-seq)居中
center_x_em <- 0
center_y_em <- 0
radius_em <- 10  # EM-seq的圆半径

# 450K和850K圆的位置
center_x_450k <- -3
center_y_450k <- 0
radius_450k <- 6  # 450K的圆半径

center_x_850k <- 3
center_y_850k <- 0
radius_850k <- 6  # 850K的圆半径

# 创建圆的数据框
circles <- data.frame(
  x0 = c(center_x_em, center_x_450k, center_x_850k),
  y0 = c(center_y_em, center_y_450k, center_y_850k),
  r = c(radius_em, radius_450k, radius_850k),
  name = c("EM-seq", "450K array", "850K array"),
  fill = c("#DDDDDD", "#3366CC", "#FF9900"),
  alpha = c(0.2, 0.5, 0.5)
)

# 创建文本标签的数据框
# 计算四个区域的位置
text_x_none <- 0  # 未被芯片覆盖的位置(EM-seq区域但不在450K和850K内)
text_y_none <- -7
text_x_450k <- -6  # 仅被450K覆盖的位置
text_y_450k <- 0
text_x_850k <- 6  # 仅被850K覆盖的位置
text_y_850k <- 0
text_x_both <- 0  # 同时被两者覆盖的位置
text_y_both <- 0

labels <- data.frame(
  x = c(text_x_none, text_x_450k, text_x_850k, text_x_both),
  y = c(text_y_none, text_y_450k, text_y_850k, text_y_both),
  label = c(
    paste0("Covered by\nEM-seq only: ", not_covered, "\n(", round(not_covered/total_dmr*100, 2), "%)"),
    paste0("Covered by\n450K array only: ", only_450k, "\n(", round(only_450k/total_dmr*100, 2), "%)"),
    paste0("Covered by\n850K array only: ", only_EPIC, "\n(", round(only_EPIC/total_dmr*100, 2), "%)"),
    paste0("Covered by\nboth arrays: ", both_covered, "\n(", round(both_covered/total_dmr*100, 2), "%)")
  ),
  size = c(13, 13, 13, 13),
  fontface = c("bold", "bold", "bold", "bold")
)

# 绘制韦恩图
venn_plot_ggplot <- ggplot() +
  # 添加圆
  geom_circle(aes(x0 = x0, y0 = y0, r = r, fill = name, alpha = alpha), data = circles) +
  # 添加文本标签
  geom_text(aes(x = x, y = y, label = label, size = size, fontface = fontface), data = labels) +
  # 设置颜色
  scale_fill_manual(values = c("EM-seq" = "#DDDDDD", "450K array" = "#3366CC", "850K array" = "#FF9900")) +
  scale_alpha_identity() +
  # 添加标题和主题
  labs(
    title = "Corverage analysis of EM-seq and array analysis",
    subtitle = paste0("Total DMRs: ", total_dmr)
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
  ) +
  coord_fixed() # 保持纵横比例一致

# 显示韦恩图
print(venn_plot_ggplot)

# 保存韦恩图
ggsave("/home/wg-liao/dmr_coverage_venn_ggplot_hg38.pdf", venn_plot_ggplot, width = 7.5, height = 7.5)





# 2.2 Make heatmap for top 200 different DMR ----
# 2.2.1 Add the clinics informaiton ----
library("pheatmap")
library("RColorBrewer")
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
# setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy")
setwd("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy")


All_sig_DMR_methy_matrix_noNA <- data.frame(fread("All_merged_DMRs_sig_methy_matrix_noNA.csv"))
# All_sig_DMR_methy_matrix_noNA <- data.frame(fread("All_sig_merged_DMRs_methy_matrix_noNA.csv"))
All_sig_DMR_methy_group_infor <- t(All_sig_DMR_methy_matrix_noNA[, c("DMR", grep("N", names(All_sig_DMR_methy_matrix_noNA), value = T))])
colnames(All_sig_DMR_methy_group_infor) <- All_sig_DMR_methy_group_infor[1,]
All_sig_DMR_methy_group_infor <- data.frame(All_sig_DMR_methy_group_infor[-1,])
patient <- substr(rownames(All_sig_DMR_methy_group_infor), 1, 5)
All_sig_DMR_methy_group_infor <- data.frame(lapply(All_sig_DMR_methy_group_infor, function(x) as.numeric(trimws(x))))
All_sig_DMR_methy_group_infor <- cbind(patient, All_sig_DMR_methy_group_infor)
low_risky_sample <- c("N1088", "N0354", "N0525", "N0039", "N0133", "N1247", "N0225", "N0313", "N1167", "N0446",
                         "N0009", "N0031", "N0079", "N0095", "N0099", "N0257", "N0468", "N1032", "N0429", "N0326")
high_risky_sample <- c("N0974", "N1344", "N0381", "N1343", "N1340", "N0305", "N1796", "N1439", "N1011", "N1445",
                        "N1207", "N1679", "N0323", "N0357", "N1250", "N1132", "N1501", "N0641", "N0771", "N0676",
                        "N0788", "N1197", "N0366", "N0954", "N1920")
risky <- rep("Middle", nrow(All_sig_DMR_methy_group_infor))
risky[All_sig_DMR_methy_group_infor$patient %in% low_risky_sample] <- "Low"
risky[All_sig_DMR_methy_group_infor$patient %in% high_risky_sample] <- "High"
All_sig_DMR_methy_group_infor <- cbind(risky, All_sig_DMR_methy_group_infor)
str(All_sig_DMR_methy_group_infor)

# All120_clinics <- fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/data/120validationsamples/120sample_clinical.csv")
All120_clinics <- fread("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/data/120validationsamples/120sample_clinical.csv")


names(All120_clinics)[1] = "patient"
All120_clinics$clinicalt[All120_clinics$clinicalt == "1C"] = "1c"
All120_clinics$clinicalt[All120_clinics$clinicalt == "2"] = "2a"
All120_clinics$clinicalt[All120_clinics$clinicalt == "2A"] = "2a"
All120_clinics$clinicalt[All120_clinics$clinicalt == "2C"] = "2c"

All120_clinics$surg_pstagt[All120_clinics$surg_pstagt == "2B"] = "2b"
All120_clinics$surg_pstagt[All120_clinics$surg_pstagt == "2C"] = "2c"
All120_clinics$surg_pstagt[All120_clinics$surg_pstagt == "3A"] = "3a"
All120_clinics$surg_pstagt[All120_clinics$surg_pstagt == "3B"] = "3b"

All120_clinics$Death[All120_clinics$Death == 0] = "Alive"
All120_clinics$Death[All120_clinics$Death == 1] = "Dead of disease"
All120_clinics$Death[All120_clinics$Death == 2] = "Dead of other cause"
#All120_clinics = All120_clinics[All120_clinics$if_seq == 1,]

summary(All120_clinics$surg_age)

All_sig_DMR_methy_group_tumor_infor <- right_join(All120_clinics[, c("patient","surg_pgleasnt", "surg_pstagt", "clinicalt", "recurrence", 
                                                                     "Death", "Bad_outcome", "PSA_recurrence", "surg_lnmets", "surg_age"
                                                                     )],
                                            All_sig_DMR_methy_group_infor, by = "patient")
str(All_sig_DMR_methy_group_tumor_infor)
table(All_sig_DMR_methy_group_tumor_infor[All_sig_DMR_methy_group_tumor_infor$risky == "Low", ]$Bad_outcome)

fwrite(All_sig_DMR_methy_group_tumor_infor, "All_merged_DMRs_sig_methy_group_clinics_infor.csv")
# fwrite(All_sig_DMR_methy_group_tumor_infor, "All_sig_merged_DMRs_methy_group_clinics_infor.csv")

# 2.2.2 Extract the top 200 different DMR ----
All_sig_DMR_methy_matrix_noNA <- All_sig_DMR_methy_matrix_noNA[order(All_sig_DMR_methy_matrix_noNA$Differential,decreasing = T), ]
top200_sig_DMR_methy_matrix_noNA <- All_sig_DMR_methy_matrix_noNA[1:200, ]
top200_sig_DMR_tumor_infor <- subset(All_sig_DMR_methy_group_tumor_infor, 
                                       select = c("patient", 
                                                  "surg_pgleasnt", "surg_pstagt", "clinicalt", "risky", top200_sig_DMR_methy_matrix_noNA$DMR))
top200_sig_DMR_clinics_infor <- subset(All_sig_DMR_methy_group_tumor_infor, 
                                     select = c("patient", "recurrence", "Death", "Bad_outcome", "PSA_recurrence", "surg_lnmets",
                                                "surg_pgleasnt", "surg_pstagt", "clinicalt", "risky", top200_sig_DMR_methy_matrix_noNA$DMR))
str(top200_sig_DMR_clinics_infor)


# 2.2.3 Make heatmap plot ----
library(circlize)
library(ComplexHeatmap)

heat_map_for_tumor_annotation_with_multi_label_certain_order <- function(raw_data, scale_or_not = T, show_row_names = T, show_column_names = T, 
                                                           clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value){
  # 新增：设置subtype的因子顺序，并重新排序数据
  raw_data$risky <- factor(raw_data$risky, 
                             levels = c("Low", "Middle", "High"))
  raw_data <- raw_data[order(raw_data$risky), ]  # 按新因子顺序排序数据
  
  rownames(raw_data) <- raw_data$patient
  #raw_data$surg_pgleasnt = factor(raw_data$surg_pgleasnt, levels = c(0,1), labels = c("Blood Healthy","Blood Tumor"))
  data_methy <- subset(raw_data, select = -c(patient, #recurrence, Death, Bad_outcome, PSA_recurrence,
                                             surg_pgleasnt, surg_pstagt, clinicalt, risky))
  if (scale_or_not == T) {data <- t(apply(data_methy, 2, scale))} else {data = t(data_methy)}
  colnames(data) <- raw_data$patient
  
  
  # set color
  #set main color
  if (scale_or_not == T) { 
    myBreaks <- seq(-5, 5, length.out = 100)
    myCol <- colorRampPalette(c('#712b8f', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  } else {
    myBreaks <- seq(0, 1, length.out = 100)
    myCol <- colorRampPalette(c('#712b8f', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  }
  
  #set annotation color
  # surg_pgleasnt <- raw_data$surg_pgleasnt
  # pick.col <- brewer.pal(9, 'Greens') # 可以设置一个主色（比如绿色），然后从中挑出对应的颜色.allowed maximum for palette Greens is 9
  # col.surg_pgleasnt <- colorRampPalette(pick.col)(length(unique(surg_pgleasnt)))
  
  # create annotation object
  ann <- data.frame(#Tumor_Recurrence = raw_data$recurrence, 
    #All_Death = raw_data$Death, 
    #Tumor_Progression = raw_data$Bad_outcome, 
    #PSA_recurrence = raw_data$PSA_recurrence, 
    
    Gleason_score = raw_data$surg_pgleasnt, 
    Surgical_Stage_T = raw_data$surg_pstagt, 
    Clinical_Stage_T = raw_data$clinicalt, 
    Risk = raw_data$risky)
  colors=list(#Tumor_Recurrence = c("1" = "black", "0" = "grey"),
    #All_Death = c("1" = "black", "0" = "grey"),
    #Tumor_Progression = c("1" = "black", "0" = "grey"),
    #PSA_recurrence = c("1" = "black", "0" = "grey"),
    
    Gleason_score = c("7" = "#FEDBC7", "8" = "#F6A482", "9" = "#D75F4C","10"= "#B31529"),
    Surgical_Stage_T = c("2a" = "#d9f1d5", "2b" = "#add4a0", "2c" = "#5cae63", "3a"="#1b7939", "3b"= "#295e11" ),
    Clinical_Stage_T = c("1c" = "#d1e5f0", "2a" = "#8ec4de", "2b" = "#3a93c3", "2c"="#1065ab"),
    Risk = c("Low" = "#5b8efd", "Middle" = "goldenrod1", "High"= "#dd217d")
  )
  colAnn <- HeatmapAnnotation(
    df = ann,
    col = colors,
    which = 'col', # set 'col' (samples) or 'row' (gene) annotation
    na_col = 'white', # NA颜色，默认白色
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'),
    # 下面这个都是类似的设置：占几行、标题、字体等
    annotation_legend_param = list(
      Gleason_score = list(
        nrow = 2, # 这个legend显示几行
        title = 'Gleason_score',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 10, fontface = 'bold'),
        labels_gp = gpar(fontsize = 10, fontface = 'plain')
      )
    )
  )
  
  # create bottom boxplot
  boxplotCol <- HeatmapAnnotation(
    boxplot = anno_boxplot(
      data,
      border = FALSE,
      gp = gpar(fill = '#CCCCCC'),
      pch = '.',
      size = unit(1, 'mm'),
      axis = TRUE,
      axis_param = list(
        gp = gpar(fontsize = 10),
        side = 'left')),
    annotation_width = unit(c(1.0), 'cm'),
    which = 'col')
  
  # Mark the rows we want
  genelabels <- rowAnnotation(
    Genes = anno_mark(
      at = seq(1, nrow(data), 1),
      labels = rownames(data)[seq(1, nrow(data), 1)],
      labels_gp = gpar(fontsize = 5, fontface = 'plain'),
      padding = 0.75),
    width = unit(2.0, 'cm') +
      
      max_text_width(
        rownames(data)[seq(1, nrow(data), 1)],
        gp = gpar(fontsize = 5,  fontface = 'plain')))
  
  
  # plot the heatmap
  hmap=Heatmap(data,
               # split the genes / rows according to the PAM clusters
               #split = NA, # kmeans_clus / hc_clus
               cluster_row_slices = F,
               column_split = raw_data$risky, # 按risky分组
               
               name = 'Methylation\nZ-score',
               
               col = colorRamp2(myBreaks, myCol),
               
               # parameters for the colour-bar that represents gradient of expression
               heatmap_legend_param = list(
                 color_bar = 'continuous',
                 legend_direction = 'vertical',
                 legend_width = unit(12, 'cm'),
                 legend_height = unit(5.0, 'cm'),
                 title_position = 'topleft',
                 title_gp=gpar(fontsize = 10, fontface = 'bold'),
                 labels_gp=gpar(fontsize = 10, fontface = 'bold')),
               
               
               # row (gene) parameters
               cluster_rows = T,
               show_row_dend = TRUE,
               row_title = '',
               row_title_side = 'left',
               row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
               row_title_rot = 90,
               show_row_names = show_row_names,
               row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
               row_names_side = 'right',
               row_dend_width = unit(20,'mm'),
               
               # column (sample) parameters
               cluster_columns = T,
               show_column_dend = TRUE,
               column_title = paste(clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value),
               column_title_side = 'top',
               column_title_gp = gpar(fontsize = 12, fontface = 'plain'),
               column_title_rot = 0,
               show_column_names = show_column_names,
               column_names_gp = gpar(fontsize = 10, fontface = 'plain'),
               column_names_side = "top",
               column_names_max_height = unit(20, 'cm'),
               column_dend_height = unit(2,'cm'),
               
               # cluster methods for rows and columns
               clustering_distance_columns = clustering_distance_columns_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_columns = clustering_method_columns_value, #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               clustering_distance_rows = clustering_distance_rows_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_rows = clustering_method_rows_value,#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               
               # specify top and bottom annotations
               top_annotation = colAnn,
               bottom_annotation = boxplotCol)
  
  return(hmap)
  
}



heat_map_for_tumor_annotation_with_multi_label <- function(raw_data, scale_or_not = T, show_row_names = T, show_column_names = T, 
                                                           clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value){
  rownames(raw_data) <- raw_data$patient
  #raw_data$surg_pgleasnt = factor(raw_data$surg_pgleasnt, levels = c(0,1), labels = c("Blood Healthy","Blood Tumor"))
  data_methy <- subset(raw_data, select = -c(patient, #recurrence, Death, Bad_outcome, PSA_recurrence,
                                             surg_pgleasnt, surg_pstagt, clinicalt, risky))
  if (scale_or_not == T) {data <- t(apply(data_methy, 2, scale))} else {data = t(data_methy)}
  colnames(data) <- raw_data$patient
  
  
  # set color
  #set main color
  if (scale_or_not == T) { 
    myBreaks <- seq(-5, 5, length.out = 100)
    myCol <- colorRampPalette(c('cornflowerblue', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  } else {
    myBreaks <- seq(0, 1, length.out = 100)
    myCol <- colorRampPalette(c('cornflowerblue', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  }
  
  #set annotation color
  # surg_pgleasnt <- raw_data$surg_pgleasnt
  # pick.col <- brewer.pal(9, 'Greens') # 可以设置一个主色（比如绿色），然后从中挑出对应的颜色.allowed maximum for palette Greens is 9
  # col.surg_pgleasnt <- colorRampPalette(pick.col)(length(unique(surg_pgleasnt)))
  
  # create annotation object
  ann <- data.frame(#Tumor_Recurrence = raw_data$recurrence, 
    #All_Death = raw_data$Death, 
    #Tumor_Progression = raw_data$Bad_outcome, 
    #PSA_recurrence = raw_data$PSA_recurrence, 
    
    Gleason_score = raw_data$surg_pgleasnt, 
    Surgical_Stage_T = raw_data$surg_pstagt, 
    Clinical_Stage_T = raw_data$clinicalt, 
    Risk = raw_data$risky)
  colors=list(#Tumor_Recurrence = c("1" = "black", "0" = "grey"),
    #All_Death = c("1" = "black", "0" = "grey"),
    #Tumor_Progression = c("1" = "black", "0" = "grey"),
    #PSA_recurrence = c("1" = "black", "0" = "grey"),
    
    Gleason_score = c("7" = "darkorange", "8" = "darkorange2", "9" = "darkorange3","10"= "darkorange4"),
    Surgical_Stage_T = c("2a" = "deepskyblue", "2b" = "dodgerblue1", "2c" = "dodgerblue2", "3a"="dodgerblue3", "3b"= "dodgerblue4" ),
    Clinical_Stage_T = c("1c" = "darkorchid1", "2a" = "darkorchid2", "2b" = "darkorchid3", "2c"="darkorchid4"),
    Risk = c("Low" = "forestgreen", "Middle" = "goldenrod1", "High"= "firebrick3")
  )
  colAnn <- HeatmapAnnotation(
    df = ann,
    col = colors,
    which = 'col', # set 'col' (samples) or 'row' (gene) annotation
    na_col = 'white', # NA颜色，默认白色
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'),
    # 下面这个都是类似的设置：占几行、标题、字体等
    annotation_legend_param = list(
      Gleason_score = list(
        nrow = 2, # 这个legend显示几行
        title = 'Gleason_score',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 10, fontface = 'bold'),
        labels_gp = gpar(fontsize = 10, fontface = 'plain')
      )
    )
  )
  
  # create bottom boxplot
  boxplotCol <- HeatmapAnnotation(
    boxplot = anno_boxplot(
      data,
      border = FALSE,
      gp = gpar(fill = '#CCCCCC'),
      pch = '.',
      size = unit(1, 'mm'),
      axis = TRUE,
      axis_param = list(
        gp = gpar(fontsize = 10),
        side = 'left')),
    annotation_width = unit(c(1.0), 'cm'),
    which = 'col')
  
  # Mark the rows we want
  genelabels <- rowAnnotation(
    Genes = anno_mark(
      at = seq(1, nrow(data), 1),
      labels = rownames(data)[seq(1, nrow(data), 1)],
      labels_gp = gpar(fontsize = 5, fontface = 'plain'),
      padding = 0.75),
    width = unit(2.0, 'cm') +
      
      max_text_width(
        rownames(data)[seq(1, nrow(data), 1)],
        gp = gpar(fontsize = 5,  fontface = 'plain')))
  
  
  # plot the heatmap
  hmap=Heatmap(data,
               # split the genes / rows according to the PAM clusters
               #split = NA, # kmeans_clus / hc_clus
               cluster_row_slices = F,
               
               name = 'Expression\nZ-score',
               
               col = colorRamp2(myBreaks, myCol),
               
               # parameters for the colour-bar that represents gradient of expression
               heatmap_legend_param = list(
                 color_bar = 'continuous',
                 legend_direction = 'vertical',
                 legend_width = unit(12, 'cm'),
                 legend_height = unit(5.0, 'cm'),
                 title_position = 'topleft',
                 title_gp=gpar(fontsize = 10, fontface = 'bold'),
                 labels_gp=gpar(fontsize = 10, fontface = 'bold')),
               
               
               # row (gene) parameters
               cluster_rows = T,
               show_row_dend = TRUE,
               row_title = '',
               row_title_side = 'left',
               row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
               row_title_rot = 90,
               show_row_names = show_row_names,
               row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
               row_names_side = 'right',
               row_dend_width = unit(20,'mm'),
               
               # column (sample) parameters
               cluster_columns = TRUE,
               show_column_dend = TRUE,
               column_title = paste(clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value),
               column_title_side = 'top',
               column_title_gp = gpar(fontsize = 12, fontface = 'plain'),
               column_title_rot = 0,
               show_column_names = show_column_names,
               column_names_gp = gpar(fontsize = 10, fontface = 'plain'),
               column_names_side = "top",
               column_names_max_height = unit(20, 'cm'),
               column_dend_height = unit(2,'cm'),
               
               # cluster methods for rows and columns
               clustering_distance_columns = clustering_distance_columns_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_columns = clustering_method_columns_value, #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               clustering_distance_rows = clustering_distance_rows_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_rows = clustering_method_rows_value,#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               
               # specify top and bottom annotations
               top_annotation = colAnn,
               bottom_annotation = boxplotCol)
  
  return(hmap)
  
}

heat_map_for_clinics_annotation_with_multi_label <- function(raw_data, scale_or_not = T, show_row_names = T, show_column_names = T, 
                                                     clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value){
  rownames(raw_data) <- raw_data$patient
  #raw_data$surg_pgleasnt = factor(raw_data$surg_pgleasnt, levels = c(0,1), labels = c("Blood Healthy","Blood Tumor"))
  data_methy <- subset(raw_data, select = -c(patient, recurrence, Death, Bad_outcome, PSA_recurrence, surg_lnmets, 
                                             surg_pgleasnt, surg_pstagt, clinicalt, risky))
  if (scale_or_not == T) {data <- t(apply(data_methy, 2, scale))} else {data = t(data_methy)}
  colnames(data) <- raw_data$patient
  
  
  # set color
  #set main color
  if (scale_or_not == T) { 
    myBreaks <- seq(-5, 5, length.out = 100)
    myCol <- colorRampPalette(c('cornflowerblue', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  } else {
    myBreaks <- seq(0, 1, length.out = 100)
    myCol <- colorRampPalette(c('cornflowerblue', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  }
  
  #set annotation color
  # surg_pgleasnt <- raw_data$surg_pgleasnt
  # pick.col <- brewer.pal(9, 'Greens') # 可以设置一个主色（比如绿色），然后从中挑出对应的颜色.allowed maximum for palette Greens is 9
  # col.surg_pgleasnt <- colorRampPalette(pick.col)(length(unique(surg_pgleasnt)))
  
  # create annotation object
  ann <- data.frame(Tumor_Recurrence = raw_data$recurrence, 
                    Lymph_Node_Mets = raw_data$surg_lnmets,
                    PSA_recurrence = raw_data$PSA_recurrence, 
                    All_Death = raw_data$Death, 
                    Tumor_Progression = raw_data$Bad_outcome,
                   
                    Gleason_score = raw_data$surg_pgleasnt, 
                    Surgical_Stage_T = raw_data$surg_pstagt, 
                    Clinical_Stage_T = raw_data$clinicalt, 
                    Risk = raw_data$risky)
  colors=list(Tumor_Recurrence = c("1" = "black", "0" = "gray80"),
              Lymph_Node_Mets = c("1" = "black", "0" = "gray80"),
              PSA_recurrence = c("1" = "black", "0" = "gray80"),
              All_Death = c("Dead of disease" = "black",  "Dead of other cause" = "gray50", "Alive" = "gray80"),
              Tumor_Progression = c("1" = "black", "0" = "gray80"),
                                 
              Gleason_score = c("7" = "darkorange", "8" = "darkorange2", "9" = "darkorange3","10"= "darkorange4"),
              Surgical_Stage_T = c("2a" = "deepskyblue", "2b" = "dodgerblue1", "2c" = "dodgerblue2", "3a"="dodgerblue3", "3b"= "dodgerblue4" ),
              Clinical_Stage_T = c("1c" = "darkorchid1", "2a" = "darkorchid2", "2b" = "darkorchid3", "2c"="darkorchid4"),
              Risk = c("Low" = "forestgreen", "Middle" = "goldenrod1", "High"= "firebrick3")
              )
  colAnn <- HeatmapAnnotation(
    df = ann,
    col = colors,
    which = 'col', # set 'col' (samples) or 'row' (gene) annotation
    na_col = 'white', # NA颜色，默认白色
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'),
    # 下面这个都是类似的设置：占几行、标题、字体等
    annotation_legend_param = list(
      Gleason_score = list(
        nrow = 2, # 这个legend显示几行
        title = 'Gleason_score',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 10, fontface = 'bold'),
        labels_gp = gpar(fontsize = 10, fontface = 'plain')
      )
    )
  )
  
  # create bottom boxplot
  boxplotCol <- HeatmapAnnotation(
    boxplot = anno_boxplot(
      data,
      border = FALSE,
      gp = gpar(fill = '#CCCCCC'),
      pch = '.',
      size = unit(1, 'mm'),
      axis = TRUE,
      axis_param = list(
        gp = gpar(fontsize = 10),
        side = 'left')),
    annotation_width = unit(c(1.0), 'cm'),
    which = 'col')
  
  # Mark the rows we want
  genelabels <- rowAnnotation(
    Genes = anno_mark(
      at = seq(1, nrow(data), 1),
      labels = rownames(data)[seq(1, nrow(data), 1)],
      labels_gp = gpar(fontsize = 5, fontface = 'plain'),
      padding = 0.75),
    width = unit(2.0, 'cm') +
      
      max_text_width(
        rownames(data)[seq(1, nrow(data), 1)],
        gp = gpar(fontsize = 5,  fontface = 'plain')))
  
  
  # plot the heatmap
  hmap=Heatmap(data,
               # split the genes / rows according to the PAM clusters
               #split = NA, # kmeans_clus / hc_clus
               cluster_row_slices = F,
               
               name = 'Expression\nZ-score',
               
               col = colorRamp2(myBreaks, myCol),
               
               # parameters for the colour-bar that represents gradient of expression
               heatmap_legend_param = list(
                 color_bar = 'continuous',
                 legend_direction = 'vertical',
                 legend_width = unit(12, 'cm'),
                 legend_height = unit(5.0, 'cm'),
                 title_position = 'topleft',
                 title_gp=gpar(fontsize = 10, fontface = 'bold'),
                 labels_gp=gpar(fontsize = 10, fontface = 'bold')),
               
               
               # row (gene) parameters
               cluster_rows = T,
               show_row_dend = TRUE,
               row_title = '',
               row_title_side = 'left',
               row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
               row_title_rot = 90,
               show_row_names = show_row_names,
               row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
               row_names_side = 'right',
               row_dend_width = unit(20,'mm'),
               
               # column (sample) parameters
               cluster_columns = TRUE,
               show_column_dend = TRUE,
               column_title = paste(clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value),
               column_title_side = 'top',
               column_title_gp = gpar(fontsize = 12, fontface = 'plain'),
               column_title_rot = 0,
               show_column_names = show_column_names,
               column_names_gp = gpar(fontsize = 10, fontface = 'plain'),
               column_names_side = "top",
               column_names_max_height = unit(20, 'cm'),
               column_dend_height = unit(2,'cm'),
               
               # cluster methods for rows and columns
               clustering_distance_columns = clustering_distance_columns_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_columns = clustering_method_columns_value, #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               clustering_distance_rows = clustering_distance_rows_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_rows = clustering_method_rows_value,#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               
               # specify top and bottom annotations
               top_annotation = colAnn,
               bottom_annotation = boxplotCol)
  
  return(hmap)
  
}



# 2.2.3.1 Tumor infor only ----
top200_sig_DMR_heatmap_with_tumor_all120_certain_order = heat_map_for_tumor_annotation_with_multi_label_certain_order(raw_data = top200_sig_DMR_tumor_infor, scale_or_not = F,
                                                       show_row_names = F, show_column_names = F, 
                                                       clustering_distance_columns_value = "euclidean",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                       clustering_method_columns_value = 'ward.D2', #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                                       clustering_distance_rows_value = "euclidean",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                       clustering_method_rows_value = 'ward.D2'#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
)
pdf("top200_sig_DMR_heatmap_with_tumor_all120_certain_order.pdf", height = 11, width = 9)
top200_sig_DMR_heatmap_with_tumor_all120_certain_order
dev.off()


top200_sig_DMR_heatmap_with_tumor_all120 = heat_map_for_tumor_annotation_with_multi_label(raw_data = top200_sig_DMR_tumor_infor, scale_or_not = F,
                                         show_row_names = F, show_column_names = F, 
                                         clustering_distance_columns_value = "euclidean",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                         clustering_method_columns_value = 'ward.D2', #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                         clustering_distance_rows_value = "euclidean",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                         clustering_method_rows_value = 'ward.D2'#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                         )
pdf("top200_sig_DMR_heatmap_with_tumor_all120.pdf", height = 10, width = 8)
top200_sig_DMR_heatmap_with_tumor_all120
dev.off()



top200_sig_DMR_heatmap_with_tumor_HighLow = heat_map_for_tumor_annotation_with_multi_label(raw_data = top200_sig_DMR_tumor_infor[patient %in% c(low_risky_sample, high_risky_sample)], scale_or_not = F,
                                                                         show_row_names = F, show_column_names = F, 
                                                                         clustering_distance_columns_value = "manhattan",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                         clustering_method_columns_value = 'ward.D2', #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                                                         clustering_distance_rows_value = "euclidean",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                         clustering_method_rows_value = 'ward.D2'#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                                                         )
pdf("top200_sig_DMR_heatmap_with_tumor_HighLow.pdf", height = 10, width = 8)
top200_sig_DMR_heatmap_with_tumor_HighLow
dev.off()




# 2.2.3.1 Tumor and survival infor ----
top200_sig_DMR_heatmap_with_clinics_all120 = heat_map_for_clinics_annotation_with_multi_label(raw_data = top200_sig_DMR_clinics_infor, scale_or_not = F,
                                                                                          show_row_names = F, show_column_names = F, 
                                                                                          clustering_distance_columns_value = "euclidean",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                          clustering_method_columns_value = 'ward.D2', #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                                                                          clustering_distance_rows_value = "euclidean",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                          clustering_method_rows_value = 'ward.D2'#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
)
pdf("top200_sig_DMR_heatmap_with_clinics_all120.pdf", height = 10, width = 8)
top200_sig_DMR_heatmap_with_clinics_all120
dev.off()

top200_sig_DMR_heatmap_with_clinics_HighLow = heat_map_for_clinics_annotation_with_multi_label(raw_data = top200_sig_DMR_clinics_infor[patient %in% c(low_risky_sample, high_risky_sample)], scale_or_not = F,
                                                                                     show_row_names = F, show_column_names = F, 
                                                                                     clustering_distance_columns_value = "manhattan",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                     clustering_method_columns_value = 'ward.D2', #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                                                                     clustering_distance_rows_value = "euclidean",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                     clustering_method_rows_value = 'ward.D2'#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
)
pdf("top200_sig_DMR_heatmap_with_clinics_HighLow.pdf", height = 10, width = 8)
top200_sig_DMR_heatmap_with_clinics_HighLow
dev.off()


# 2.3 Annotation of Hyper and Hypo significant DMR ----
library(data.table)       
library(dplyr)
library(genomation)
library(scales)
library(patchwork)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(scales) 


# 2.3.1 Genetic annotation ----
# setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy")
setwd("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy")

All_sig_DMR_methy_matrix_noNA <- data.frame(fread("All_merged_DMRs_sig_methy_matrix_noNA.csv"))
# All_sig_DMR_methy_matrix_noNA <- data.frame(fread("All_sig_merged_DMRs_methy_matrix_noNA.csv"))
hyper_hypo <- ifelse(All_sig_DMR_methy_matrix_noNA$High_Avg > All_sig_DMR_methy_matrix_noNA$Low_Avg, "Hyper", "Hypo")
All_sig_DMR_methy_matrix_noNA <- cbind(hyper_hypo, All_sig_DMR_methy_matrix_noNA)
Hyper_sig_DMR_methy_matrix_noNA <- All_sig_DMR_methy_matrix_noNA[All_sig_DMR_methy_matrix_noNA$hyper_hypo == "Hyper",]
Hypo_sig_DMR_methy_matrix_noNA <- All_sig_DMR_methy_matrix_noNA[All_sig_DMR_methy_matrix_noNA$hyper_hypo == "Hypo",]

# annot <- "/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/share/annot/annot.hg38.ncbirefseq.12col.bed" #Annotation file
annot <- "D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/share/annot/annot.hg38.ncbirefseq.12col.bed" #Annotation file

gene.obj <- readTranscriptFeatures(annot,
                                   up.flank=1000,
                                   down.flank=1000,
                                   remove.unusual = TRUE)
genetic_proportion_df <- data.frame()
genetic_count_df <- data.frame()
for (i in 1:2) {
  if (i == 1){data = Hyper_sig_DMR_methy_matrix_noNA } else {data = Hypo_sig_DMR_methy_matrix_noNA}
  ann.obj <- annotateWithGeneParts(as(data,"GRanges"), gene.obj)
  genetic_count_df <- rbind(genetic_count_df, ann.obj@num.annotation)
  genetic_proportion_df <- rbind(genetic_proportion_df, genomation::getTargetAnnotationStats(ann.obj,percentage=TRUE,precedence=TRUE))
}
names(genetic_count_df) <- c("promoter", "exon", "intron", "intergenic")
genetic_count_df <- data.frame(genetic_count_df, DMR_type = c("Hyper DMR", "Hypo DMR"))
names(genetic_proportion_df) <- c("promoter", "exon", "intron", "intergenic")
genetic_proportion_df <- data.frame(genetic_proportion_df, DMR_type = c("Hyper DMR", "Hypo DMR"))
fwrite(genetic_count_df, "genetic_annotation_count.csv")
fwrite(genetic_proportion_df, "genetic_annotation_proportion.csv")


# 2.3.1.1 Count plot----
genetic_count_df_long <- reshape::melt(genetic_count_df, id.vars = "DMR_type", variable_name = "roles")
genetic_count_df_long$roles <- factor(genetic_count_df_long$roles, levels = c("promoter", "exon", "intron", "intergenic"))# to adjust the order in the barplot
# plot
custom_log_trans <- trans_new(
  name = "custom_log",
  transform = function(x) log10(x + 1), # Add 1 to avoid log of zero issues
  inverse = function(x) 10^x - 1, # Subtract 1 to revert the transformation correctly
  breaks = function(x) c(0, 1, 5, 10, 20, 50, 100), # Custom breaks for desired spacing
  domain = c(0, Inf)
)
genetic_count_plot = ggplot(genetic_count_df_long, aes(x = DMR_type, y = value, fill = roles)) +
  # 条形图函数：未将position参数显示设置为dodge，则绘制出的条形图为堆积型
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(name = "",
                    values = c("#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2"),
                    labels = c("promoter", "exon", "intron", "intergenic"))+ # adjust the legend
  theme_bw() + 
  theme(#legend.position="none", #不需要图例
    axis.text.x=element_text(#amily="Times",
      colour="black",size=14, angle = 0, hjust = 0.5), #设置x轴刻度标签的字体属性
    axis.text.y=element_text(#family="Times",
      size=14,face="plain"), #设置x轴刻度标签的字体属性
    axis.title.y=element_text(#family="Times",
      size = 20,face="plain"), #设置y轴的标题的字体属性
    axis.title.x=element_text(#family="Times",
      size = 20,face="plain"), #设置x轴的标题的字体属性
    legend.title = element_text(#family="Times",
      size = 14,face="plain"), # Adjust legend title size here
    legend.text = element_text(#family="Times",
      size = 14,face="plain"), # Adjust legend keys (levels) text size here
    plot.title = element_text(#family="Times",
      size=23,face="bold",hjust = 0.5), #设置总标题的字体属性
    panel.grid.major = element_blank(), #不显示网格线
    panel.grid.minor = element_blank(),
    plot.background = element_blank(), # No background
    panel.border = element_blank(), # No border
    axis.line = element_line(colour = "black") # Add axis lines
  )+
  scale_y_continuous(
    trans = custom_log_trans,  # Apply custom log transformation
    breaks = c(0, 10, 100, 1000, 10000),  # Define specific breaks for Y-axis
    labels = c("0", "10", "100", "1000", "10000"),  # Custom labels for the Y-axis
    limits = c(0, 12000)
  ) +
  labs(x = "", y = "Counts", title = " ") #设置x轴和y轴的标题


# 2.3.1.2 Proportion plot----
genetic_proportion_df_long <- reshape::melt(genetic_proportion_df, id.vars = "DMR_type", variable_name = "roles")
genetic_proportion_df_long$roles <- factor(genetic_proportion_df_long$roles, levels = c("promoter", "exon", "intron", "intergenic"))# to adjust the order in the barplot

custom_log_trans <- trans_new(
  name = "custom_log",
  transform = function(x) log10(x + 1), # Add 1 to avoid log of zero issues
  inverse = function(x) 10^x - 1, # Subtract 1 to revert the transformation correctly
  breaks = function(x) c(0, 1, 5, 10, 20, 100), # Custom breaks for desired spacing
  domain = c(0, Inf)
)
genetic_proportion_plot = ggplot(genetic_proportion_df_long, aes(x = DMR_type, y = value, fill = roles)) +
  # 条形图函数：未将position参数显示设置为dodge，则绘制出的条形图为堆积型
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(name = "",
                    values = c("#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2"),
                    labels = c("promoter", "exon", "intron", "intergenic"))+ # adjust the legend
  theme_bw() + 
  theme(#legend.position="none", #不需要图例
    axis.text.x=element_text(#amily="Times",
      colour="black",size=14, angle = 0, hjust = 0.5), #设置x轴刻度标签的字体属性
    axis.text.y=element_text(#family="Times",
      size=14,face="plain"), #设置x轴刻度标签的字体属性
    axis.title.y=element_text(#family="Times",
      size = 20,face="plain"), #设置y轴的标题的字体属性
    axis.title.x=element_text(#family="Times",
      size = 20,face="plain"), #设置x轴的标题的字体属性
    legend.title = element_text(#family="Times",
      size = 14,face="plain"), # Adjust legend title size here
    legend.text = element_text(#family="Times",
      size = 14,face="plain"), # Adjust legend keys (levels) text size here
    plot.title = element_text(#family="Times",
      size=23,face="bold",hjust = 0.5), #设置总标题的字体属性
    panel.grid.major = element_blank(), #不显示网格线
    panel.grid.minor = element_blank(),
    plot.background = element_blank(), # No background
    panel.border = element_blank(), # No border
    axis.line = element_line(colour = "black") # Add axis lines
  )+
  scale_y_continuous(
    trans = custom_log_trans,  # Apply custom log transformation
    breaks = c(0, 1, 5, 20, 100),  # Define specific breaks for Y-axis
    labels = c("0%", "1%", "5%", "20%", "100%"),  # Custom labels for the Y-axis
    limits = c(0, 100)
  ) +
  labs(x = "", y = "Proportion", title = " ") #设置x轴和y轴的标题


# 2.3.1.3 Output the plot ----
genetic_combined_plot <- (genetic_count_plot / genetic_proportion_plot) + 
  plot_layout(guides = 'collect')  
#&theme(legend.position = "bottom") # Move legend to bottom
ggsave("genetic_annotation_Counts&Proportion.pdf", 
       genetic_combined_plot, height = 7, width = 5)

# 2.3.2 Methylation annotation ----
# island <- "/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/share/annot/CpG_annotation_hg38_from_UCSC.bed"
island <- "D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/share/annot/CpG_annotation_hg38_from_UCSC.bed"


gr_cpg_islands <- rtracklayer::import(island)

# Define shores as 2kb up- and downstream of islands, excluding the islands
gr_cpg_shores <- GRanges(seqnames = seqnames(gr_cpg_islands),
                         ranges = IRanges(start = pmax(start(gr_cpg_islands) - 2000, 1),
                                          end = end(gr_cpg_islands) + 2000))
gr_cpg_shores <- GenomicRanges::reduce(gr_cpg_shores)
gr_cpg_shores <- GenomicRanges::setdiff(gr_cpg_shores, gr_cpg_islands)

# Define shelves as 2kb up- and downstream of shores, excluding shores and islands
gr_cpg_shelves <- GRanges(seqnames = seqnames(gr_cpg_shores),
                          ranges = IRanges(start = pmax(start(gr_cpg_shores) - 2000, 1),
                                           end = end(gr_cpg_shores) + 2000))
gr_cpg_shelves <- GenomicRanges::reduce(gr_cpg_shelves)
gr_cpg_shelves <- GenomicRanges::setdiff(gr_cpg_shelves, gr_cpg_shores)
gr_cpg_shelves <- GenomicRanges::setdiff(gr_cpg_shelves, gr_cpg_islands)

CpG_anno_df <- data.frame()
for (i in 1:2) {
  if (i == 1){data = Hyper_sig_DMR_methy_matrix_noNA } else {data = Hypo_sig_DMR_methy_matrix_noNA}
  gr_data <- makeGRangesFromDataFrame(data)
  # Calculate overlaps
  overlap_islands <- countOverlaps(gr_data, gr_cpg_islands)
  overlap_shores <- countOverlaps(gr_data, gr_cpg_shores)
  overlap_shelves <- countOverlaps(gr_data, gr_cpg_shelves)
  
  # Count CpG sites in each category
  count_islands <- sum(overlap_islands > 0)
  count_shores <- sum(overlap_shores > 0)
  count_shelves <- sum(overlap_shelves > 0)
  
  # Print the counts
  CpG_anno <- c(count_islands, count_shores, count_shelves, nrow(data))
  CpG_anno_df <- rbind(CpG_anno_df, CpG_anno)
}
names(CpG_anno_df) <- c("count_islands", "count_shores", "count_shelves", "count_all") # correct the variable names
CpG_anno_df$count_other <- CpG_anno_df$count_all - (CpG_anno_df$count_islands + CpG_anno_df$count_shores + CpG_anno_df$count_shelves)
CpG_anno_df <- data.frame(CpG_anno_df, DMR_type = c("Hyper DMR", "Hypo DMR"))

CpG_anno_df_long <- reshape::melt(CpG_anno_df, id.vars = "DMR_type", variable_name = "roles")
count_df_long <- CpG_anno_df_long[CpG_anno_df_long$roles != "count_all",]
count_df_long$roles <- factor(count_df_long$roles, levels = c("count_islands", "count_shores", "count_shelves", "count_other"))# to adjust the order in the barplot
fwrite(CpG_anno_df, "CpG_annotation_count.csv")


# 2.3.2.1 count plot -----
custom_log_trans <- trans_new(
  name = "custom_log",
  transform = function(x) log10(x + 1), # Add 1 to avoid log of zero issues
  inverse = function(x) 10^x - 1, # Subtract 1 to revert the transformation correctly
  breaks = function(x) c(0, 1, 5, 10, 20, 50, 100), # Custom breaks for desired spacing
  domain = c(0, Inf)
)
count_plot = ggplot(count_df_long, aes(x = DMR_type, y = value, fill = roles)) +
  # 条形图函数：未将position参数显示设置为dodge，则绘制出的条形图为堆积型
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(name = "CpG Overlap",
                    values = c("#2878B5", "#9AC9DB", "#F8AC8C", "#C82423"),
                    labels = c("Islands","Shores", "Shelves", "Other"))+ # adjust the legend
  theme_bw() + 
  theme(#legend.position="none", #不需要图例
    axis.text.x = element_text(colour = "black", size = 14, angle = 0, hjust = 0.5), # x-axis labels
    axis.text.y = element_text(size = 14, face = "plain"), # y-axis labels
    axis.title.y = element_text(size = 20, face = "plain"), # y-axis title
    axis.title.x = element_text(size = 20, face = "plain"), # x-axis title
    legend.title = element_text(size = 14, face = "plain"), # legend title
    legend.text = element_text(size = 14, face = "plain"), # legend text
    plot.title = element_text(size = 23, face = "bold", hjust = 0.5), # plot title
    panel.grid.major = element_blank(), # no grid lines
    panel.grid.minor = element_blank(),
    plot.background = element_blank(), # no background
    panel.border = element_blank(), # no border
    axis.line = element_line(colour = "black") # axis lines
  )+
  scale_y_continuous(
    trans = custom_log_trans,  # Apply custom log transformation
    breaks = c(0, 10, 100, 1000, 10000),  # Define specific breaks for Y-axis
    labels = c("0", "10", "100", "1000", "10000"),  # Custom labels for the Y-axis
    limits = c(0, 12000)
  ) +
  labs(x = "", y = "Counts", title = " ") #设置x轴和y轴的标题


# 2.3.2.2 Proportion plot ----
CpG_anno_df
proportion_df <- CpG_anno_df
proportion_df$count_other <- proportion_df$count_other/proportion_df$count_all *100
proportion_df$count_islands <- proportion_df$count_islands/proportion_df$count_all*100
proportion_df$count_shores <- proportion_df$count_shores/proportion_df$count_all*100
proportion_df$count_shelves<- proportion_df$count_shelves/proportion_df$count_all*100
proportion_df_long <- reshape::melt(proportion_df, id.vars = "DMR_type", variable_name = "roles")
proportion_df_long <- proportion_df_long[proportion_df_long$roles != "count_all",]
proportion_df_long$roles <- factor(proportion_df_long$roles, levels = c("count_islands", "count_shores", "count_shelves", "count_other"))# to adjust the order in the barplot
fwrite(proportion_df, "CpG_annotation_proportion.csv")

custom_log_trans <- trans_new(
  name = "custom_log",
  transform = function(x) log10(x + 1), # Add 1 to avoid log of zero issues
  inverse = function(x) 10^x - 1, # Subtract 1 to revert the transformation correctly
  breaks = function(x) c(0, 1, 5, 10, 20, 50, 100), # Custom breaks for desired spacing
  domain = c(0, Inf)
)
# Create the bar plot with the custom transformation
proportion_plot <- ggplot(proportion_df_long, aes(x = DMR_type, y = value, fill = roles)) +
  geom_bar(stat = "identity", position = position_dodge()) + # Create grouped bar plot
  scale_fill_manual(
    name = "CpG Overlap",
    values = c("#2878B5", "#9AC9DB", "#F8AC8C", "#C82423"),
    labels = c("Islands","Shores", "Shelves", "Other")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(colour = "black", size = 14, angle = 0, hjust = 0.5), # x-axis labels
    axis.text.y = element_text(size = 14, face = "plain"), # y-axis labels
    axis.title.y = element_text(size = 20, face = "plain"), # y-axis title
    axis.title.x = element_text(size = 20, face = "plain"), # x-axis title
    legend.title = element_text(size = 14, face = "plain"), # legend title
    legend.text = element_text(size = 14, face = "plain"), # legend text
    plot.title = element_text(size = 23, face = "bold", hjust = 0.5), # plot title
    panel.grid.major = element_blank(), # no grid lines
    panel.grid.minor = element_blank(),
    plot.background = element_blank(), # no background
    panel.border = element_blank(), # no border
    axis.line = element_line(colour = "black") # axis lines
  ) +
  # Use the custom transformation for the Y-axis
  scale_y_continuous(
    trans = custom_log_trans,  # Apply custom log transformation
    breaks = c(0, 1, 5, 20, 100),  # Define specific breaks for Y-axis
    labels = c("0%", "1%", "5%", "20%", "100%"),  # Custom labels for the Y-axis
    limits = c(0, 100)
  ) +
  labs(x = "", y = "Proportion", title = "") # Set titles


# 2.3.2.3 Output the plot ----
combined_plot <- (count_plot / proportion_plot) + 
  plot_layout(guides = 'collect')  
#&theme(legend.position = "bottom") # Move legend to bottom
#ggsave("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DMS_annotation_by_diff_methods/CpG_Annotation/CpG_annotation_Counts&Proportion.pdf", 
#       combined_plot, height = 10, width = 6)
ggsave("CpG_annotation_Counts&Proportion.pdf", combined_plot, height = 7, width = 5)


# 3.0 TSS distance relevant anlysis ----


# 4.0 Gene Ontology analysis ----
# 4.1 Position Annotation ---- 
# 4.1.1 import annotation file ---
setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/Overall_methylation_summary")

library("rtracklayer")
gtf_data = import('/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/data/Homo_sapiens.GRCh38.111.gtf') #gtf的路径
gtf_data = as.data.frame(gtf_data)
gtf_data[gtf_data == "NA"] <- NA
gtf_data$seqnames2 <- paste0("chr", gtf_data$seqnames)

# 4.1.2 create GRanges objects ----
All_DMRs <- data.frame(fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_matrix_noNA.csv"))
All_DMRs <- All_DMRs[order(All_DMRs$Differential, decreasing = T), ]
#All_DMRs <- fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DMS_annotation_for_34_and_120/DMS_summary_by_different_strategy/All_common_DM_raw_Methy_value.csv")
str(All_DMRs)

correlation_ranges <- GRanges( 
  seqnames = Rle(All_DMRs$chr),
  ranges = IRanges(start=All_DMRs$start, end=All_DMRs$end)
)

gtf_ranges <- GRanges(
  seqnames = Rle(as.character(gtf_data$seqnames2)),
  ranges = IRanges(start=gtf_data$start, end=gtf_data$end),
  gene_name = gtf_data$gene_name
)

# 4.1.3 find overlap gene region ----
overlaps <- findOverlaps(correlation_ranges, gtf_ranges)

# 4.1.4 extract gene name corresponding to each row index in All_DMRs ----
matched_gene_names <- data.frame(
  index = queryHits(overlaps),
  gene_name = gtf_ranges$gene_name[subjectHits(overlaps)]
)
matched_gene_names <- matched_gene_names[!is.na(matched_gene_names$gene_name), ]
matched_gene_names$Differential <- (All_DMRs$High_Avg - All_DMRs$Low_Avg)[matched_gene_names$index]
matched_gene_names <- matched_gene_names[!duplicated(matched_gene_names$gene_name),]
matched_gene_names <- matched_gene_names[order(matched_gene_names$Differential, decreasing = T),]
fwrite(data.frame(matched_gene_names[,2:3]),sep="\t", "matched_gene_names.rnk")

# 4.1.5 extract gene information corresponding to each row index in All_DMRs ----
All_DMR_gene_annotation = cbind(
  All_DMRs[queryHits(overlaps),1:9],
  gtf_data[subjectHits(overlaps),]
)
fwrite(All_DMR_gene_annotation, "All_DMR_gene_annotation.csv")



# 4.2 Analysis ----
library(clusterProfiler)#富集包
library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
# 4.2.1 prepare data frame for Gene Ontology analysis (no this part)----
# 4.4 extract the ENTREZID of each gene name/symbol ----
gene.df <- bitr(matched_gene_names$gene_name,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)#TCGA数据框如果没有进行基因注释，那么fromType应该是Ensembl，各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求                     
gene <- gene.df$ENTREZID
# 4.2.3 extract GO result ----
##CC表示细胞组分，MF表示分子功能，BP表示生物学过程，ALL表示同时富集三种过程，选自己需要的,我一般是做BP,MF,CC这3组再合并成一个数据框，方便后续摘取部分通路绘图。
ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)

ego_result_BP <- ego_result_BP[order(ego_result_BP$qvalue,decreasing = F), ]
ego_result_CC <- ego_result_CC[order(ego_result_CC$qvalue,decreasing = F), ]
ego_result_MF <- ego_result_MF[order(ego_result_MF$qvalue,decreasing = F), ]

# 4.2.4 select top 15 gene pathways ----
display_number = c(10, 10, 10)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- ego_result_BP[1:display_number[1], ]
ego_result_CC <- ego_result_CC[1:display_number[2], ]
ego_result_MF <- ego_result_MF[1:display_number[3], ]

go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),      
  qvalue = c(ego_result_BP$qvalue, ego_result_CC$qvalue, ego_result_MF$qvalue),
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("Biological process", display_number[1]), 
                rep("Cellular component", display_number[2]),
                rep("Molecular function", display_number[3])), 
              levels=c("Biological process", "Cellular component","Molecular function" )))

# 4.2.5 GO plot ----
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
go_enrich_df <- go_enrich_df[!is.na(go_enrich_df$ID), ]

library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))#设定颜色
# for all
ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=`qvalue`)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  facet_grid(type~., scale="free")+ # ggplot2的分面语法
  #scale_fill_manual(values = COLS) + ###颜色
  scale_fill_gradientn(colours = rev(myPalette(100))) +
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") + 
  ylab("Covered Gene Number") + 
  labs(title = "The Most Enriched GO Terms (Ranked by q-value)")+
  theme_bw()+ #背景变为白色  
  theme(#legend.position="none", #不需要图例
    axis.text.x=element_text(#amily="Times",
      colour="black",size=12, angle = 45, hjust = 1), #设置x轴刻度标签的字体属性
    axis.text.y=element_text(#family="Times",
      size=12,face="plain"), #设置x轴刻度标签的字体属性
    axis.title.y=element_text(#family="Times",
      size = 16,face="plain"), #设置y轴的标题的字体属性
    axis.title.x=element_text(#family="Times",
      size = 16,face="plain"), #设置x轴的标题的字体属性
    legend.title = element_text(#family="Times",
      size = 14,face="plain"), # 设置图例标题的字体属性
    legend.text = element_text(#family="Times",
      size = 12,face="plain"), # # 设置图例文本的字体属性
    plot.title = element_text(#family="Times",
      size=20,face="bold",hjust = 0.5), #设置总标题的字体属性
    panel.grid.major = element_blank(), #不显示网格线
    panel.grid.minor = element_blank(),
    plot.background = element_blank(), # No background
    panel.border = element_blank(), # No border
    axis.line = element_line(colour = "black")) # Add axis lines
    
ggsave("GO_for_all_DMS.pdf", height = 10, width = 9.5)


# 4.0.0 Gene Set Enrichment Analysis (GSEA) ----
# 4.0.1 import MSigDB database ----
# BiocManager::install("msigdbr")
library(msigdbr)
Hall <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
C1 <- msigdbr(species = "Homo sapiens", category = "C1") %>% dplyr::select(gs_name, entrez_gene)
C2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::select(gs_name, entrez_gene)
C3 <- msigdbr(species = "Homo sapiens", category = "C3") %>% dplyr::select(gs_name, entrez_gene)
C4 <- msigdbr(species = "Homo sapiens", category = "C4") %>% dplyr::select(gs_name, entrez_gene)
C5 <- msigdbr(species = "Homo sapiens", category = "C5") %>% dplyr::select(gs_name, entrez_gene)
C6 <- msigdbr(species = "Homo sapiens", category = "C6") %>% dplyr::select(gs_name, entrez_gene)
C7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% dplyr::select(gs_name, entrez_gene)
C8 <- msigdbr(species = "Homo sapiens", category = "C8") %>% dplyr::select(gs_name, entrez_gene)

# 4.0.2 do GSEA ----
# setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/Overall_methylation_summary")
setwd("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/Overall_methylation_summary")

library(clusterProfiler)
matched_gene_names = fread("matched_gene_names.rnk")
geneList = matched_gene_names$Differential
names(geneList) = matched_gene_names$gene_name
id = bitr(names(geneList), "SYMBOL", "ENTREZID", "org.Hs.eg.db")
geneList = geneList[names(geneList) %in% id[,1]] #keep the genes names that have been transformed
names(geneList) = id[match(names(geneList), id[,1]), 2]  

result_H <- GSEA(geneList, TERM2GENE = Hall,
                  minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "fdr",)
result_H_df <- cbind(group = "H", as.data.table(result_H))
result_c1 <- GSEA(geneList, TERM2GENE = C1,
                  minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "fdr",)
result_c1_df <- cbind(group = "c1", as.data.table(result_c1))
result_c2 <- GSEA(geneList, TERM2GENE = C2,
                  minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "fdr",)
result_c2_df <- cbind(group = "c2", as.data.table(result_c2))
result_c3 <- GSEA(geneList, TERM2GENE = C3,
                  minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "fdr",)
result_c3_df <- cbind(group = "c3", as.data.table(result_c3))
result_c4 <- GSEA(geneList, TERM2GENE = C4,
                  minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.25, pAdjustMethod = "fdr",)
result_c4_df <- cbind(group = "c4", as.data.table(result_c4))
result_c5 <- GSEA(geneList, TERM2GENE = C5,
                  minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "fdr",)
result_c5_df <- cbind(group = "c5", as.data.table(result_c5))
result_c6 <- GSEA(geneList, TERM2GENE = C6,
                  minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.25, pAdjustMethod = "fdr",)
result_c6_df <- cbind(group = "c6", as.data.table(result_c6))
result_c7 <- GSEA(geneList, TERM2GENE = C7,
                  minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.25, pAdjustMethod = "fdr",)
result_c7_df <- cbind(group = "c7", as.data.table(result_c7))
result_c8 <- GSEA(geneList, TERM2GENE = C8,
                  minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.25, pAdjustMethod = "fdr",)
result_c8_df <- cbind(group = "c8", as.data.table(result_c8))

result_summary <- rbind(#result_H_df,
                        result_c1_df,
                        result_c2_df,
                        result_c3_df,
                        result_c4_df,
                        result_c5_df
                        #result_c6_df,
                        #result_c7_df
                        #DataFrame(result_c8)
                        )
fwrite(as.data.table(result_summary), "gsea_result_summary.csv")


# 4.0.2 GSEA plots----
library(enrichplot)
gseaplot2(result_c5,
          result_c5$ID[1],#富集的ID编号
          title = result_c5$Description[1],#标题
          color = "red", #GSEA线条颜色
          base_size = 15,#基础字体大小
          rel_heights = c(1.5, 0.5, 1),#副图的相对高度
          subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
          ES_geom = "line", #enrichment score用线还是用点"dot"
          pvalue_table = F) #显示pvalue等信息
  ggsave(filename = 'GSEA_c5_sig_1.pdf', width =6, height =6)

gseaplot2(result_c2,
          result_c2$ID[1],#富集的ID编号
          title = result_c2$Description[1],#标题
          color = "red", #GSEA线条颜色
          base_size = 15,#基础字体大小
          rel_heights = c(1.2, 0.3, 0.5),#副图的相对高度
          subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
          ES_geom = "line", #enrichment score用线还是用点"dot"
          pvalue_table = F) #显示pvalue等信息
ggsave(filename = 'GSEA_c2_sig_1.pdf', width =6, height =6)


# 默认
p1 <- cnetplot(result_c3
               , showCategory = 5 # 也可以直接写条目名字
               , layout = "kk" #网络形状，’star’, ’circle’, ’gem’, ’dh’, ’graphopt’, ’grid’, ’mds’, ’randomly’, ’fr’, ’kk’, ’drl’ or ’lgl’
               , node_label = "all" #显示哪些节点的标签，’category’, ’gene’, ’all’(默认), ’none’
               , shadowtext = "all" #哪些节点标签需要添加阴影，’category’, ’gene’, ’all’(默认), ’none’
               # 控制节点和连线的颜色
               , color.params = list(foldChange = NULL 
                                     , edge = FALSE #根据富集的不同条目上色
                                     , category = "#E5C494" #条目节点颜色
                                     , gene = "#B3B3B3" #基因节点颜色
               )
               
               # 控制标签和节点的大小
               , cex.params = list(category_node = .1
                                   , gene_node = .1
                                   , category_label = 1
                                   , gene_label = .5)
               
               # 控制哪些节点和连线高亮显示
               , hilight.params = list(category = NULL
                                       , alpha_hilight = 1
                                       , alpha_no_hilight = 0.3)
               
)
p1


# 5.0 DMR selection for prognostic analysis ----
# All_sig_DMR_methy_group_tumor_infor = fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")
All_sig_DMR_methy_group_tumor_infor = fread("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")


table(All_sig_DMR_methy_group_tumor_infor$risky)
table(All_sig_DMR_methy_group_tumor_infor$surg_pstagt)
all_DMR_clean = subset(All_sig_DMR_methy_group_tumor_infor, select = c("Bad_outcome",
                                                                       grep("chr*", names(All_sig_DMR_methy_group_tumor_infor), value = T)))
names(all_DMR_clean)[1] = "severity"

# 5.1 Elastic Net Regression feature selection ----
# setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
setwd("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")

source("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/Rscript/Elastic_Net_Regression.R")
source("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/Rscript/Elastic_Net_Regression.R")

DMR_selection <- elastic_net_regression(data = all_DMR_clean, 
                                        target_var = "severity", 
                                        num_iterations = 100, 
                                        target_num_vars = 100, 
                                        potential_lambda = exp(seq(-50,10,length.out = 150)))
fwrite(DMR_selection, "ENR_sites_selection.csv")




# 5.2 XGBoost----
sites_selection <- fread("ENR_sites_selection.csv")
index <- seq(0,ncol(sites_selection),2)
sites_name <- unlist(sites_selection[, ..index], use.names = F)
dim(sites_selection)
frequncy <- table(sites_name)
frequncy <- frequncy[order(frequncy,decreasing = T)]
#top100_frequncy <- frequncy[frequncy > frequncy[150]]
top100_frequncy <- frequncy #keep all DMRs instead of selecting top N ones.
sites_from_ENR <- names(top100_frequncy)[-which(names(top100_frequncy) == "severity")]# exclude the "severity" which represent the times of iteration
fwrite(data.frame(sites_from_ENR), "ENR_top100.csv")

sites_from_ENR <- fread("ENR_top100.csv")
sites_from_ENR <- sites_from_ENR[-1,]
str(all_DMR_clean)
DMR_from_ENR <- subset(all_DMR_clean, select = c("severity", sites_from_ENR$sites_from_ENR))
str(DMR_from_ENR)

source("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/Rscript/XGBoost.r")
summary_xgboost_for_parameter_design <- xgboost_for_feature_selection(data = DMR_from_ENR,
                                                                      target_var = "severity",
                                                                      num_iterations = 1,
                                                                      max_depth_values = 2:6,
                                                                      eta_values = 10^(-3:2),
                                                                      nround_values = seq(10, 150, by = 10),
                                                                      gamma_values = 10^(-3:2),
                                                                      min_child_weight_values = 10^(-3:1))
fwrite(summary_xgboost_for_parameter_design[[2]], "summary_xgboost_for_parameter_design.csv")

XGB_result <- xgboost_for_feature_selection(data = DMR_from_ENR,
                                            target_var = "severity",
                                            num_iterations = 20,
                                            max_depth_values = 2,
                                            eta_values = 10^(-1:-1),
                                            nround_values = seq(10, 200, by = 10),
                                            gamma_values = 10^(-3:-1),
                                            min_child_weight_values = 10^(-3:-1))

Xgboost_import <- XGB_result[[1]]
summary_xgboost <- XGB_result[[2]]
fwrite(summary_xgboost, "summary_xgboost_for_different_trainset.csv")
fwrite(Xgboost_import, "Xgboost_import_gene.csv")


# 5.3 Random Forest----
source("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/Rscript/RandomForest.R")
RF_result1 <- fit_random_forest(data = DMR_from_ENR,
                                target_var = "severity",
                                num_iterations = 100, ntree_values = seq(10, 100, by = 5), mtry_values = seq(1, 25))


RF_import <- RF_result1[[1]]
summary_RF <- RF_result1[[2]]
fwrite(summary_RF, "summary_RF_for_for_different_trainset.csv")
fwrite(RF_import, "RF_import_gene.csv")


# 5.4 results summary ----
setwd("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
#setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")

# 5.4.1 Random Froest ----
RF_import2 <- fread("RF_import_gene.csv")

#RF_import2 <- RF_import[,1:200]
summary_RF2 <- fread("summary_RF_for_for_different_trainset.csv")
#summary_RF2 <- summary_RF
summary_RF2$times <-  as.numeric(rownames(summary_RF2))

pdf("HIST_of_RF_Accuracy.pdf", width = 5, height = 4)
hist(summary_RF2$accuracy_test, xlim = c(0,1), xlab = "Random Forest Test Accuracy", main = "For cohort of 120 samples")
dev.off()

RF_import2 <- RF_import2[,-1]
top_RF <- summary_RF2[summary_RF2$accuracy_test >= 0.7,] #select all the models with test accuracy higher than 0.5
index <- sort(c((top_RF$times)*5, (top_RF$times)*5-2, (top_RF$times)*5-1))
top_RF_impor <- subset(RF_import2, select = index)
index <- seq(0, ncol(top_RF_impor),3)
top_RF_impor_sites <- subset(top_RF_impor, select = index)
top_RF_impor_sites_summary <- sort(table(unlist(top_RF_impor_sites)), decreasing = T)
top_RF_impor_sites_summary <- top_RF_impor_sites_summary[-1]

pdf("HIST_of_RF_sites.pdf", width = 5, height = 4)
hist(top_RF_impor_sites_summary/(ncol(RF_import2)/5), xlab = "Selected rate of Sites from Random Forest", main = "For cohort of 120 samples")
dev.off()

top_RF_impor_sites_summary_DF <- data.frame(select_frequency = top_RF_impor_sites_summary, select_rate = top_RF_impor_sites_summary/(ncol(RF_import2)/5))
fwrite(top_RF_impor_sites_summary_DF, "top_RF_impor_sites_summary_DF.csv")

# 5.4.2 XGBoost ----
Xgboost_import2 <- fread("Xgboost_import_gene.csv")
str(Xgboost_import2)
Xgboost_import2[Xgboost_import2 == ""] <- NA
summary_xgboost2 <- fread("summary_xgboost_for_different_trainset.csv")
summary_xgboost2$times <-  as.numeric(rownames(summary_xgboost2))

pdf("HIST_of_xgboost_Accuracy.pdf", width = 5, height = 4)
hist(summary_xgboost2$accuracy_test, xlim = c(0,1), xlab = "XGBoost Test Accuracy", main = "For cohort of 120 samples")
dev.off()

Xgboost_import2 <- Xgboost_import2[,-1]
top_xgboost2 <- summary_xgboost2[summary_xgboost2$accuracy_test >= 0.7,] #select all the models with test accuracy higher than 0.9
index <- top_xgboost2$times
top_xgboost_impor <- Xgboost_import2[, ..index]
#index <- seq(0, ncol(top_xgboost_impor),3)
top_xgboost_impor_sites <- top_xgboost_impor
top_xgboost_impor_sites_summary <- sort(table(unlist(top_xgboost_impor_sites)), decreasing = T)
#top_xgboost_impor_sites_summary <- top_xgboost_impor_sites_summary[-1]

pdf("HIST_of_xgboost_sites.pdf", width = 5, height = 3)
hist(top_xgboost_impor_sites_summary/(ncol(Xgboost_import2)), xlab = "Selected rate of Sites for XGBoost", main = "For cohort of 120 samples")
dev.off()

top_xgboost_impor_sites_summary_DF <- data.frame(site_names = names(top_xgboost_impor_sites_summary), 
                                                 select_frequency = top_xgboost_impor_sites_summary,
                                                 select_rate = top_xgboost_impor_sites_summary / (ncol(Xgboost_import2) -1))
fwrite(top_xgboost_impor_sites_summary_DF, "top_xgboost_impor_sites_summary_DF.csv") 

# 5.5 checking the jonit target DMRs ----
#setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
setwd("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")


common_DMR_XGBoost_RF <- intersect(top_xgboost_impor_sites_summary_DF$site_names, top_RF_impor_sites_summary_DF$select_frequency.Var1)
common_DMR_methy_XGBoost_RF <- subset(all_DMR_clean, select = c("severity", common_DMR_XGBoost_RF))
common_DMR_methy_XGBoost_RF$patient <- rownames(common_DMR_methy_XGBoost_RF)
fwrite(common_DMR_methy_XGBoost_RF, "common_DMR_methy_XGBoost_RF.csv")

unique_DMR_XGBoost_RF <- union(top_xgboost_impor_sites_summary_DF$site_names, top_RF_impor_sites_summary_DF$select_frequency.Var1)
unique_DMR_methy_XGBoost_RF <- subset(all_DMR_clean, select = c("severity", unique_DMR_XGBoost_RF))
unique_DMR_methy_XGBoost_RF$patient <- rownames(unique_DMR_methy_XGBoost_RF)
fwrite(unique_DMR_methy_XGBoost_RF, "unique_DMR_methy_XGBoost_RF.csv")



# 5.5.1 Common Sites from XGBoost and RF ----
# 5.5.1.1 Fit a logistic model with lasso penalty to select non-colinear predictors ----
#https://juliastats.org/MixedModels.jl/dev/rankdeficiency/#Fixed-effects
common_DMR_methy_XGBoost_RF <- fread("common_DMR_methy_XGBoost_RF.csv")
rownames(common_DMR_methy_XGBoost_RF) <- common_DMR_methy_XGBoost_RF$patient

source("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/Rscript/logistic_regression.R")
source("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/Rscript/logistic_regression.R")

common_logistic_common_sites = lasso_logistic_regression_for_select_noncolinear_predictor(data = common_DMR_methy_XGBoost_RF[,-c("patient")], 
                                                                                          target_var = "severity", 
                                                                                          scale_or_not = FALSE,
                                                                                          num_iterations = 100, 
                                                                                          family = "binomial",
                                                                                          typemeasure = "class",
                                                                                          nfolds = 5,
                                                                                          custom_lambda = NULL)

# 5.5.1.1.1 Cross valisation plot ----
cv_fit_list <- common_logistic_common_sites[[4]]
length(cv_fit_list)
pdf(paste("common_XGBoost_RF_CV_fit_plot.pdf", sep = ""), width = 4, height = 4)
for (i in 1:length(cv_fit_list)) {
  plot(cv_fit_list[[i]])
}
dev.off()

# 5.5.1.1.2 AUC plot ----
AUC_summary <- common_logistic_common_sites[[2]]
fwrite(AUC_summary, "common_XGBoost_RF_AUC_summary.csv")

ggplot(AUC_summary, aes(x = "",y = roc_train)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "text", aes(label=sprintf("%.2f", after_stat(y))), vjust = -4, size = 3) + # Add mean labels
  theme(axis.text.x = element_text(angle = 0)) +
  ylim(0,1.3) +
  labs(x = "AUC", y = "AUC Value", title = "")
ggsave(filename = "common_XGBoost_RF_AUC_distribution.pdf",width = 2, height = 3)

# 5.5.1.1.3 predictor summary ----
Common_DMR_predictor <- common_logistic_common_sites[[3]]
Common_DMR_predictor2 <- Common_DMR_predictor
fwrite(Common_DMR_predictor, "common_XGBoost_RF_predictor_selected_by_logistic_regression.csv")

Common_DMR_predictor <- fread("common_XGBoost_RF_predictor_selected_by_logistic_regression.csv")
index <- seq(0,ncol(Common_DMR_predictor),2)
Common_DMR_predictor_name <- unlist(data.table(Common_DMR_predictor)[, ..index], use.names = F)
Common_DMR_predictor_frequncy <- table(Common_DMR_predictor_name)
Common_DMR_predictor_frequncy <- Common_DMR_predictor_frequncy[order(Common_DMR_predictor_frequncy,decreasing = T)]

pdf("HIST_of_logistic_with_lasso_sites.pdf", width = 5, height = 4)
hist(Common_DMR_predictor_frequncy[-1]/100, xlab = "Selected rate of Sites for LR with lasso", main = "DMR for 34 and 120 samples")
dev.off()

Common_DMR_predictor_frequncy_select <- Common_DMR_predictor_frequncy[Common_DMR_predictor_frequncy == 100]
Common_DMR_predictor_target <- names(Common_DMR_predictor_frequncy_select)[-which(names(Common_DMR_predictor_frequncy_select) %in% c("(Intercept)", ""))]
fwrite(data.frame(Common_DMR_predictor_target), "common_XGBoost_RF_predictor_target.csv")

# ************* ----
# 6.0 Annotate the target DMR (althoug a lot of the items were named as DMR ) ----
rm(list = ls())
# 6.1 heat_map ---- 
library("pheatmap")
library("RColorBrewer")
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

heat_map_for_survival_annotation_with_multi_label_and_certain_order <- function(raw_data, scale_or_not = T, show_row_names = T, show_column_names = T, 
                                                              clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value){
  #raw_data = common_XGBoost_RF_predictor_target
  #rownames(raw_data) <- raw_data$patient
  #raw_data$surg_pgleasnt = factor(raw_data$surg_pgleasnt, levels = c(0,1), labels = c("Blood Healthy","Blood Tumor"))
  raw_data$risky <- factor(raw_data$risky, 
                           levels = c("Low", "Middle", "High"))
  raw_data <- raw_data[order(raw_data$risky), ]  # 按新因子顺序排序数据
  
  data_methy <- subset(raw_data, select = -c(patient,recurrence, Death, Bad_outcome, PSA_recurrence, risky,surg_lnmets))
  if (scale_or_not == T) {data <- t(apply(data_methy, 2, scale))} else {data = t(data_methy)}
  colnames(data) <- raw_data$patient
  
  
  # set color
  #set main color
  if (scale_or_not == T) { 
    myBreaks <- seq(-5, 5, length.out = 100)
    myCol <- colorRampPalette(c('#712b8f', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  } else {
    myBreaks <- seq(0, 1, length.out = 100)
    myCol <- colorRampPalette(c('#712b8f', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  }
  
  #set annotation color
  # surg_pgleasnt <- raw_data$surg_pgleasnt
  # pick.col <- brewer.pal(9, 'Greens') # 可以设置一个主色（比如绿色），然后从中挑出对应的颜色.allowed maximum for palette Greens is 9
  # col.surg_pgleasnt <- colorRampPalette(pick.col)(length(unique(surg_pgleasnt)))
  
  # create annotation object
  ann <- data.frame(Tumor_Recurrence = raw_data$recurrence, 
                    Lymph_Node_Mets = raw_data$surg_lnmets,
                    PSA_recurrence = raw_data$PSA_recurrence, 
                    All_Death = raw_data$Death, 
                    Tumor_Progression = raw_data$Bad_outcome,
                    Risk = raw_data$risky)
  colors=list(Tumor_Recurrence = c("1" = "black", "0" = "gray80"),
              Lymph_Node_Mets = c("1" = "black", "0" = "gray80"),
              PSA_recurrence = c("1" = "black", "0" = "gray80"),
              All_Death = c("Dead of disease" = "black",  "Dead of other cause" = "gray50", "Alive" = "gray80"),
              Tumor_Progression = c("1" = "black", "0" = "gray80"),
              Risk = c("Low" = "#5b8efd", "Middle" = "goldenrod1", "High"= "#dd217d")
  )
  
  colAnn <- HeatmapAnnotation(
    df = ann,
    col = colors,
    which = 'col', # set 'col' (samples) or 'row' (gene) annotation
    na_col = 'white', # NA color, default is white
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'),
    # Correct the legend parameter key to match the column names in `ann`
    annotation_legend_param = list(
      Risk = list(
        nrow = 2, # Display legend in 2 rows
        title = 'Risk', # Set the title of the legend
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 10, fontface = 'bold'),
        labels_gp = gpar(fontsize = 10, fontface = 'plain')
      )
    )
  )
  
  # create bottom boxplot
  boxplotCol <- HeatmapAnnotation(
    boxplot = anno_boxplot(
      data,
      border = FALSE,
      gp = gpar(fill = '#CCCCCC'),
      pch = '.',
      size = unit(1, 'mm'),
      axis = TRUE,
      axis_param = list(
        gp = gpar(fontsize = 10),
        side = 'left')),
    annotation_width = unit(c(1.0), 'cm'),
    which = 'col')
  
  # Mark the rows we want
  genelabels <- rowAnnotation(
    Genes = anno_mark(
      at = seq(1, nrow(data), 1),
      labels = rownames(data)[seq(1, nrow(data), 1)],
      labels_gp = gpar(fontsize = 5, fontface = 'plain'),
      padding = 0.75),
    width = unit(2.0, 'cm') +
      
      max_text_width(
        rownames(data)[seq(1, nrow(data), 1)],
        gp = gpar(fontsize = 5,  fontface = 'plain')))
  
  
  # plot the heatmap
  hmap=Heatmap(data,
               # split the genes / rows according to the PAM clusters
               #split = NA, # kmeans_clus / hc_clus
               cluster_row_slices = F,
               column_split = raw_data$risky, # 按subtype分组
               
               name = 'Methylation\nZ-score',
               
               col = colorRamp2(myBreaks, myCol),
               
               # parameters for the colour-bar that represents gradient of expression
               heatmap_legend_param = list(
                 color_bar = 'continuous',
                 legend_direction = 'vertical',
                 legend_width = unit(12, 'cm'),
                 legend_height = unit(5.0, 'cm'),
                 title_position = 'topleft',
                 title_gp=gpar(fontsize = 10, fontface = 'bold'),
                 labels_gp=gpar(fontsize = 10, fontface = 'bold')),
               
               
               # row (gene) parameters
               cluster_rows = T,
               show_row_dend = TRUE,
               row_title = '',
               row_title_side = 'left',
               row_title_gp = gpar(fontsize = 12,  fontface = 'plain'),
               row_title_rot = 90,
               show_row_names = show_row_names,
               row_names_gp = gpar(fontsize = 10, fontface = 'plain'),
               row_names_side = 'right',
               row_dend_width = unit(20,'mm'),
               
               # column (sample) parameters
               cluster_columns = T,
               show_column_dend = TRUE,
               column_title = paste(clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value),
               column_title_side = 'top',
               column_title_gp = gpar(fontsize = 12, fontface = 'plain'),
               column_title_rot = 0,
               show_column_names = show_column_names,
               column_names_gp = gpar(fontsize = 10, fontface = 'plain'),
               column_names_side = "top",
               column_names_max_height = unit(20, 'cm'),
               column_dend_height = unit(2,'cm'),
               
               # cluster methods for rows and columns
               clustering_distance_columns = clustering_distance_columns_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_columns = clustering_method_columns_value, #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               clustering_distance_rows = clustering_distance_rows_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_rows = clustering_method_rows_value,#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               
               # specify top and bottom annotations
               top_annotation = colAnn,
               bottom_annotation = boxplotCol)
  
  return(hmap)
  
}


heat_map_for_annotation <- function(raw_data, scale_or_not = T, show_row_names = T, show_column_names = T, 
                                    clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value){
  
  # rownames(raw_data) <- raw_data$patient
  # raw_data$patient <- substr(rownames(raw_data), nchar(raw_data$patient) -3, nchar(raw_data$patient))
  raw_data$severity = factor(raw_data$severity, levels = c(0,1), labels = c("Indolent","Aggressive"))
  data_methy <- subset(raw_data, select = -c(patient,severity))
  if (scale_or_not == T) {data <- t(apply(data_methy, 2, scale))} else {data = t(data_methy)}
  colnames(data) <- raw_data$patient
  
  # set color
  #set main color
  if (scale_or_not == T) { 
    myBreaks <- seq(-3, 3, length.out = 1000)
    myCol <- colorRampPalette(c('blue', 'white', 'red'))(1000)
    col = colorRamp2(myBreaks, myCol)
  } else {
    myBreaks <- seq(0, 1, length.out = 1000)
    myCol <- colorRampPalette(c('blue', 'white', 'red'))(1000)
    col = colorRamp2(myBreaks, myCol)
  }
  
  
  
  #set annotation color
  Severity <- raw_data$severity
  pick.col <- brewer.pal(9, 'Greens') # 可以设置一个主色（比如绿色），然后从中挑出对应的颜色.allowed maximum for palette Greens is 9
  col.Severity <- colorRampPalette(pick.col)(length(unique(Severity)))
  
  # create annotation object
  ann <- data.frame(Severity = Severity)
  colors=list(Severity = c("Indolent" = "lightsalmon4", "Aggressive" = "mistyrose"))
  
  colAnn <- HeatmapAnnotation(
    df = ann,
    col = colors,
    which = 'col', # set 'col' (samples) or 'row' (gene) annotation
    na_col = 'white', # NA颜色，默认白色
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'),
    # 下面这个都是类似的设置：占几行、标题、字体等
    annotation_legend_param = list(
      Severity = list(
        nrow = 2, # 这个legend显示几行
        title = 'Severity',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold'),
        labels_gp = gpar(fontsize = 10, fontface = 'plain')
      )
    )
  )
  
  # create bottom boxplot
  boxplotCol <- HeatmapAnnotation(
    boxplot = anno_boxplot(
      data,
      border = FALSE,
      gp = gpar(fill = '#CCCCCC'),
      pch = '.',
      size = unit(1, 'mm'),
      axis = TRUE,
      axis_param = list(
        gp = gpar(fontsize = 10),
        side = 'left')),
    annotation_width = unit(c(1.0), 'cm'),
    which = 'col')
  
  # Mark the rows we want
  genelabels <- rowAnnotation(
    Genes = anno_mark(
      at = seq(1, nrow(data), 1),
      labels = rownames(data)[seq(1, nrow(data), 1)],
      labels_gp = gpar(fontsize = 5, fontface = 'plain'),
      padding = 0.75),
    width = unit(2.0, 'cm') +
      
      max_text_width(
        rownames(data)[seq(1, nrow(data), 1)],
        gp = gpar(fontsize = 5,  fontface = 'plain')))
  
  
  # plot the heatmap
  hmap=Heatmap(data,
               # split the genes / rows according to the PAM clusters
               #split = NA, # kmeans_clus / hc_clus
               cluster_row_slices = F,
               
               name = 'Methylation\nZ-score',
               
               col = colorRamp2(myBreaks, myCol),
               
               # parameters for the colour-bar that represents gradient of expression
               heatmap_legend_param = list(
                 color_bar = 'continuous',
                 legend_direction = 'vertical',
                 legend_width = unit(12, 'cm'),
                 legend_height = unit(5.0, 'cm'),
                 title_position = 'topleft',
                 title_gp=gpar(fontsize = 10, fontface = 'bold'),
                 labels_gp=gpar(fontsize = 10, fontface = 'bold')),
               
               
               # row (gene) parameters
               cluster_rows = T,
               show_row_dend = TRUE,
               row_title = '',
               row_title_side = 'left',
               row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
               row_title_rot = 90,
               show_row_names = show_row_names,
               row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
               row_names_side = 'right',
               row_dend_width = unit(20,'mm'),
               
               # column (sample) parameters
               cluster_columns = TRUE,
               show_column_dend = TRUE,
               column_title = paste(clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value),
               column_title_side = 'top',
               column_title_gp = gpar(fontsize = 12, fontface = 'plain'),
               column_title_rot = 0,
               show_column_names = show_column_names,
               column_names_gp = gpar(fontsize = 10, fontface = 'plain'),
               column_names_side = "top",
               column_names_max_height = unit(20, 'cm'),
               column_dend_height = unit(20,'mm'),
               
               # cluster methods for rows and columns
               clustering_distance_columns = clustering_distance_columns_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_columns = clustering_method_columns_value, #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               clustering_distance_rows = clustering_distance_rows_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_rows = clustering_method_rows_value,#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               
               # specify top and bottom annotations
               top_annotation = colAnn,
               bottom_annotation = boxplotCol)
  
  return(hmap)
  
  
  
}

heat_map_for_survival_annotation_with_multi_label <- function(raw_data, scale_or_not = T, show_row_names = T, show_column_names = T, 
                                                     clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value){
  #raw_data = common_XGBoost_RF_predictor_target
  #rownames(raw_data) <- raw_data$patient
  #raw_data$surg_pgleasnt = factor(raw_data$surg_pgleasnt, levels = c(0,1), labels = c("Blood Healthy","Blood Tumor"))
  data_methy <- subset(raw_data, select = -c(patient,recurrence, Death, Bad_outcome, PSA_recurrence, risky,surg_lnmets))
  if (scale_or_not == T) {data <- t(apply(data_methy, 2, scale))} else {data = t(data_methy)}
  colnames(data) <- raw_data$patient
  
  
  # set color
  #set main color
  if (scale_or_not == T) { 
    myBreaks <- seq(-5, 5, length.out = 100)
    myCol <- colorRampPalette(c('#712b8f', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  } else {
    myBreaks <- seq(0, 1, length.out = 100)
    myCol <- colorRampPalette(c('#712b8f', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  }
  
  #set annotation color
  # surg_pgleasnt <- raw_data$surg_pgleasnt
  # pick.col <- brewer.pal(9, 'Greens') # 可以设置一个主色（比如绿色），然后从中挑出对应的颜色.allowed maximum for palette Greens is 9
  # col.surg_pgleasnt <- colorRampPalette(pick.col)(length(unique(surg_pgleasnt)))
  
  # create annotation object
  ann <- data.frame(Tumor_Recurrence = raw_data$recurrence, 
                    Lymph_Node_Mets = raw_data$surg_lnmets,
                    PSA_recurrence = raw_data$PSA_recurrence, 
                    All_Death = raw_data$Death, 
                    Tumor_Progression = raw_data$Bad_outcome,
                    Risk = raw_data$risky)
  colors=list(Tumor_Recurrence = c("1" = "black", "0" = "gray80"),
              Lymph_Node_Mets = c("1" = "black", "0" = "gray80"),
              PSA_recurrence = c("1" = "black", "0" = "gray80"),
              All_Death = c("Dead of disease" = "black",  "Dead of other cause" = "gray50", "Alive" = "gray80"),
              Tumor_Progression = c("1" = "black", "0" = "gray80"),
              Risk = c("Low" = "#5b8efd", "Middle" = "goldenrod1", "High"= "#dd217d")
  )

  colAnn <- HeatmapAnnotation(
    df = ann,
    col = colors,
    which = 'col', # set 'col' (samples) or 'row' (gene) annotation
    na_col = 'white', # NA color, default is white
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'),
    # Correct the legend parameter key to match the column names in `ann`
    annotation_legend_param = list(
      Risk = list(
        nrow = 2, # Display legend in 2 rows
        title = 'Risk', # Set the title of the legend
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 10, fontface = 'bold'),
        labels_gp = gpar(fontsize = 10, fontface = 'plain')
      )
    )
  )
  
  # create bottom boxplot
  boxplotCol <- HeatmapAnnotation(
    boxplot = anno_boxplot(
      data,
      border = FALSE,
      gp = gpar(fill = '#CCCCCC'),
      pch = '.',
      size = unit(1, 'mm'),
      axis = TRUE,
      axis_param = list(
        gp = gpar(fontsize = 10),
        side = 'left')),
    annotation_width = unit(c(1.0), 'cm'),
    which = 'col')
  
  # Mark the rows we want
  genelabels <- rowAnnotation(
    Genes = anno_mark(
      at = seq(1, nrow(data), 1),
      labels = rownames(data)[seq(1, nrow(data), 1)],
      labels_gp = gpar(fontsize = 5, fontface = 'plain'),
      padding = 0.75),
    width = unit(2.0, 'cm') +
      
      max_text_width(
        rownames(data)[seq(1, nrow(data), 1)],
        gp = gpar(fontsize = 5,  fontface = 'plain')))
  
  
  # plot the heatmap
  hmap=Heatmap(data,
               # split the genes / rows according to the PAM clusters
               #split = NA, # kmeans_clus / hc_clus
               cluster_row_slices = F,
               
               name = 'Expression\nZ-score',
               
               col = colorRamp2(myBreaks, myCol),
               
               # parameters for the colour-bar that represents gradient of expression
               heatmap_legend_param = list(
                 color_bar = 'continuous',
                 legend_direction = 'vertical',
                 legend_width = unit(12, 'cm'),
                 legend_height = unit(5.0, 'cm'),
                 title_position = 'topleft',
                 title_gp=gpar(fontsize = 10, fontface = 'bold'),
                 labels_gp=gpar(fontsize = 10, fontface = 'bold')),
               
               
               # row (gene) parameters
               cluster_rows = T,
               show_row_dend = TRUE,
               row_title = '',
               row_title_side = 'left',
               row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
               row_title_rot = 90,
               show_row_names = show_row_names,
               row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
               row_names_side = 'right',
               row_dend_width = unit(20,'mm'),
               
               # column (sample) parameters
               cluster_columns = TRUE,
               show_column_dend = TRUE,
               column_title = paste(clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value),
               column_title_side = 'top',
               column_title_gp = gpar(fontsize = 12, fontface = 'plain'),
               column_title_rot = 0,
               show_column_names = show_column_names,
               column_names_gp = gpar(fontsize = 10, fontface = 'plain'),
               column_names_side = "top",
               column_names_max_height = unit(20, 'cm'),
               column_dend_height = unit(2,'cm'),
               
               # cluster methods for rows and columns
               clustering_distance_columns = clustering_distance_columns_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_columns = clustering_method_columns_value, #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               clustering_distance_rows = clustering_distance_rows_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_rows = clustering_method_rows_value,#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               
               # specify top and bottom annotations
               top_annotation = colAnn,
               bottom_annotation = boxplotCol)
  
  return(hmap)
  
}

heat_map_for_clinics_annotation_with_multi_label <- function(raw_data, scale_or_not = T, show_row_names = T, show_column_names = T, 
                                                             clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value){
  rownames(raw_data) <- raw_data$patient
  #raw_data$surg_pgleasnt = factor(raw_data$surg_pgleasnt, levels = c(0,1), labels = c("Blood Healthy","Blood Tumor"))
  data_methy <- subset(raw_data, select = -c(patient, recurrence, Death, Bad_outcome, PSA_recurrence, surg_lnmets,
                                             surg_pgleasnt, surg_pstagt, clinicalt, risky))
  if (scale_or_not == T) {data <- t(apply(data_methy, 2, scale))} else {data = t(data_methy)}
  colnames(data) <- raw_data$patient
  
  
  # set color
  #set main color
  if (scale_or_not == T) { 
    myBreaks <- seq(-5, 5, length.out = 100)
    myCol <- colorRampPalette(c('cornflowerblue', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  } else {
    myBreaks <- seq(0, 1, length.out = 100)
    myCol <- colorRampPalette(c('cornflowerblue', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  }
  
  #set annotation color
  # surg_pgleasnt <- raw_data$surg_pgleasnt
  # pick.col <- brewer.pal(9, 'Greens') # 可以设置一个主色（比如绿色），然后从中挑出对应的颜色.allowed maximum for palette Greens is 9
  # col.surg_pgleasnt <- colorRampPalette(pick.col)(length(unique(surg_pgleasnt)))
  
  # create annotation object
  ann <- data.frame(Tumor_Recurrence = raw_data$recurrence, 
                    Lymph_Node_Mets = raw_data$surg_lnmets,
                    PSA_recurrence = raw_data$PSA_recurrence, 
                    All_Death = raw_data$Death, 
                    Tumor_Progression = raw_data$Bad_outcome,
                    
                    Gleason_score = raw_data$surg_pgleasnt, 
                    Surgical_Stage_T = raw_data$surg_pstagt, 
                    Clinical_Stage_T = raw_data$clinicalt, 
                    Risk = raw_data$risky)
  colors=list(Tumor_Recurrence = c("1" = "black", "0" = "gray80"),
              Lymph_Node_Mets = c("1" = "black", "0" = "gray80"),
              PSA_recurrence = c("1" = "black", "0" = "gray80"),
              All_Death = c("Dead of disease" = "black",  "Dead of other cause" = "gray50", "Alive" = "gray80"),
              Tumor_Progression = c("1" = "black", "0" = "gray80"),
              
              Gleason_score = c("7" = "darkorange", "8" = "darkorange2", "9" = "darkorange3","10"= "darkorange4"),
              Surgical_Stage_T = c("2a" = "deepskyblue", "2b" = "dodgerblue1", "2c" = "dodgerblue2", "3a"="dodgerblue3", "3b"= "dodgerblue4" ),
              Clinical_Stage_T = c("1c" = "darkorchid1", "2a" = "darkorchid2", "2b" = "darkorchid3", "2c"="darkorchid4"),
              Risk = c("Low" = "forestgreen", "Middle" = "goldenrod1", "High"= "firebrick3")
  )
  
  colAnn <- HeatmapAnnotation(
    df = ann,
    col = colors,
    which = 'col', # set 'col' (samples) or 'row' (gene) annotation
    na_col = 'white', # NA颜色，默认白色
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'),
    # 下面这个都是类似的设置：占几行、标题、字体等
    annotation_legend_param = list(
      Gleason_score = list(
        nrow = 2, # 这个legend显示几行
        title = 'Gleason_score',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 10, fontface = 'bold'),
        labels_gp = gpar(fontsize = 10, fontface = 'plain')
      )
    )
  )
  
  # create bottom boxplot
  boxplotCol <- HeatmapAnnotation(
    boxplot = anno_boxplot(
      data,
      border = FALSE,
      gp = gpar(fill = '#CCCCCC'),
      pch = '.',
      size = unit(1, 'mm'),
      axis = TRUE,
      axis_param = list(
        gp = gpar(fontsize = 10),
        side = 'left')),
    annotation_width = unit(c(1.0), 'cm'),
    which = 'col')
  
  # Mark the rows we want
  genelabels <- rowAnnotation(
    Genes = anno_mark(
      at = seq(1, nrow(data), 1),
      labels = rownames(data)[seq(1, nrow(data), 1)],
      labels_gp = gpar(fontsize = 5, fontface = 'plain'),
      padding = 0.75),
    width = unit(2.0, 'cm') +
      
      max_text_width(
        rownames(data)[seq(1, nrow(data), 1)],
        gp = gpar(fontsize = 5,  fontface = 'plain')))
  
  
  # plot the heatmap
  hmap=Heatmap(data,
               # split the genes / rows according to the PAM clusters
               #split = NA, # kmeans_clus / hc_clus
               cluster_row_slices = F,
               
               name = 'Expression\nZ-score',
               
               col = colorRamp2(myBreaks, myCol),
               
               # parameters for the colour-bar that represents gradient of expression
               heatmap_legend_param = list(
                 color_bar = 'continuous',
                 legend_direction = 'vertical',
                 legend_width = unit(12, 'cm'),
                 legend_height = unit(5.0, 'cm'),
                 title_position = 'topleft',
                 title_gp=gpar(fontsize = 10, fontface = 'bold'),
                 labels_gp=gpar(fontsize = 10, fontface = 'bold')),
               
               
               # row (gene) parameters
               cluster_rows = T,
               show_row_dend = TRUE,
               row_title = '',
               row_title_side = 'left',
               row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
               row_title_rot = 90,
               show_row_names = show_row_names,
               row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
               row_names_side = 'right',
               row_dend_width = unit(20,'mm'),
               
               # column (sample) parameters
               cluster_columns = TRUE,
               show_column_dend = TRUE,
               column_title = paste(clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value),
               column_title_side = 'top',
               column_title_gp = gpar(fontsize = 12, fontface = 'plain'),
               column_title_rot = 0,
               show_column_names = show_column_names,
               column_names_gp = gpar(fontsize = 10, fontface = 'plain'),
               column_names_side = "top",
               column_names_max_height = unit(20, 'cm'),
               column_dend_height = unit(2,'cm'),
               
               # cluster methods for rows and columns
               clustering_distance_columns = clustering_distance_columns_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_columns = clustering_method_columns_value, #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               clustering_distance_rows = clustering_distance_rows_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_rows = clustering_method_rows_value,#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               
               # specify top and bottom annotations
               top_annotation = colAnn,
               bottom_annotation = boxplotCol)
  
  return(hmap)
  
}

# 6.1.1 survival infor only----
# setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
setwd("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")

common_XGBoost_RF_predictor_target <- fread("common_XGBoost_RF_predictor_target.csv")
# All_sig_DMR_methy_group_tumor_infor = fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")
All_sig_DMR_methy_group_tumor_infor = fread("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")

common_XGBoost_RF_predictor_target <- subset(All_sig_DMR_methy_group_tumor_infor, select = c("patient", "recurrence", "Death", "Bad_outcome", "PSA_recurrence", "surg_lnmets",
                                                                                             "risky", common_XGBoost_RF_predictor_target$Common_DMR_predictor_target))

common_XGBoost_RF_predictor_target_survival_annotation_heatmap_scaled_certain_order <- heat_map_for_survival_annotation_with_multi_label_and_certain_order(raw_data = common_XGBoost_RF_predictor_target, 
                                                                                                                           scale_or_not = T,show_row_names = T, show_column_names = T,
                                                                                                                           clustering_distance_columns_value = "minkowski",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                                                           clustering_method_columns_value = 'ward.D', #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                                                                                                           clustering_distance_rows_value = "pearson",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                                                           clustering_method_rows_value = 'ward.D'#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
)
pdf("common_XGBoost_RF_predictor_target_survival_annotation_heatmap_scaled_certain_order.pdf", width = 17, height = 7)
common_XGBoost_RF_predictor_target_survival_annotation_heatmap_scaled_certain_order
dev.off()




common_XGBoost_RF_predictor_target_survival_annotation_heatmap_scaled <- heat_map_for_survival_annotation_with_multi_label(raw_data = common_XGBoost_RF_predictor_target, 
                                                                                                         scale_or_not = T,show_row_names = T, show_column_names = T,
                                                                                        clustering_distance_columns_value = "minkowski",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                        clustering_method_columns_value = 'ward.D', #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                                                                        clustering_distance_rows_value = "pearson",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                        clustering_method_rows_value = 'ward.D'#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
)
pdf("common_XGBoost_RF_predictor_target_survival_annotation_heatmap_scaled.pdf", width = 20, height = 7)
common_XGBoost_RF_predictor_target_survival_annotation_heatmap_scaled
dev.off()

common_XGBoost_RF_predictor_target_survival_annotation_heatmap_scaled <- heat_map_for_survival_annotation_with_multi_label(raw_data = common_XGBoost_RF_predictor_target[common_XGBoost_RF_predictor_target$risky != "Middle", ], 
                                                                                                                           scale_or_not = T,show_row_names = T, show_column_names = T,
                                                                                                                           clustering_distance_columns_value = "minkowski",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                                                           clustering_method_columns_value = 'ward.D', #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                                                                                                           clustering_distance_rows_value = "pearson",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                                                           clustering_method_rows_value = 'ward.D'#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
)
pdf("common_XGBoost_RF_predictor_no_Middle_group_target_survival_annotation_heatmap_scaled.pdf", width = 13, height = 7)
common_XGBoost_RF_predictor_target_survival_annotation_heatmap_scaled
dev.off()

# scaled
# Loop through all combinations of the parameters
clustering_distances_columns <- c("euclidean", "manhattan", "canberra", "minkowski", "pearson", "spearman", "kendall")
clustering_methods_columns <- c("ward.D", "ward.D2", "complete", "average", "mcquitty")
clustering_distances_rows <- c("euclidean", "manhattan", "canberra", "minkowski", "pearson", "spearman", "kendall")
clustering_methods_rows <- c("ward.D", "ward.D2", "complete", "average", "mcquitty")

plot_list <- list()
N = 1
for (i in 1:length(clustering_distances_columns)) {
  for (j in 1:length(clustering_methods_columns)) {
    for (k in 1:length(clustering_distances_rows)){
      for (l in 1:length(clustering_methods_rows)){
        plot_list[[N]] <- heat_map_for_survival_annotation_with_multi_label(raw_data = common_XGBoost_RF_predictor_target, 
                                                                            scale_or_not = T,show_row_names = T, show_column_names = T,
                                                  clustering_distance_columns_value = clustering_distances_columns[i],
                                                  clustering_method_columns_value = clustering_methods_columns[j],
                                                  clustering_distance_rows_value = clustering_distances_rows[k],
                                                  clustering_method_rows_value = clustering_methods_rows[l])
        print(paste(N, "of", length(clustering_distances_columns)*length(clustering_methods_columns)*length(clustering_distances_rows)*length(clustering_methods_rows), 
                    clustering_distances_columns[i], clustering_methods_columns[j], clustering_distances_rows[k], clustering_methods_rows[l]))
        N = N + 1
      }
    }
  }
}

# plot all the figures (it will take a while)
pdf("discovery_all_heatmaps_scaled.pdf", height = 5, width = 15)
for (i in 1:length(plot_list)) {
  print(plot_list[[i]])
}
dev.off()

#minkowski ward.D canberra mcquitty
#minkowski ward.D pearson ward.D

# 6.1.2 survival and tumor infor ----
setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
common_XGBoost_RF_predictor_target <- fread("common_XGBoost_RF_predictor_target.csv")
All_sig_DMR_methy_group_tumor_infor = fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")
common_XGBoost_RF_predictor_target <- subset(All_sig_DMR_methy_group_tumor_infor, select = c("patient", "surg_pgleasnt", "surg_pstagt", "clinicalt", "surg_lnmets",
                                                                                             "recurrence", "Death", "Bad_outcome", "PSA_recurrence", "risky", common_XGBoost_RF_predictor_target$Common_DMR_predictor_target))



common_XGBoost_RF_predictor_target_clinics_annotation_heatmap_scaled <- heat_map_for_clinics_annotation_with_multi_label(raw_data = common_XGBoost_RF_predictor_target, 
                                                                                                                           scale_or_not = F,show_row_names = T, show_column_names = T,
                                                                                                                           clustering_distance_columns_value = "euclidean",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                                                           clustering_method_columns_value = 'ward.D2', #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                                                                                                           clustering_distance_rows_value = "euclidean",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                                                           clustering_method_rows_value = 'ward.D2'#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
)
pdf("common_XGBoost_RF_predictor_target_clinics_annotation_heatmap_scaled.pdf", width = 20, height = 8)
common_XGBoost_RF_predictor_target_clinics_annotation_heatmap_scaled
dev.off()



# 6.2 boxplot ----
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
#resesign the function
box_plot_site_annotation <- function(raw_data, pvalue = T, feature_names = "Severity"){
  raw_data = raw_data[raw_data$severity %in% c(1, 0), ]
  raw_data$severity = factor(raw_data$severity)
  data_long <- reshape2::melt(subset(raw_data, select = -c(patient)), id.vars = "severity", variable.name = "DMR", value.name = "Methylation values")
  names(data_long)[1] = "Severity"
  
  # Calculate p-values for each target_site
  calculate_p_value <- function(subset_data) {
    wilcox_test <- wilcox.test(`Methylation values` ~ Severity, data = subset_data)
    return(wilcox_test$p.value)
  }
  
  p_values <- data_long %>% group_by(`DMR`) %>% summarise(p_value = calculate_p_value(cur_data()))
  
  if (pvalue == TRUE) {
    p_values$annotation <- ifelse(p_values$p_value < 0.001, paste0('***\nPvalue=\n', format(p_values$p_value, scientific = TRUE)), 
                                  ifelse(p_values$p_value < 0.01, paste0('**\nPvalue=\n', format(p_values$p_value, scientific = TRUE)), 
                                         ifelse(p_values$p_value < 0.05, paste0('*\nPvalue=\n', format(p_values$p_value, scientific = TRUE)), 
                                                paste0('ns\nPvalue=\n', format(p_values$p_value, scientific = TRUE)))))
  } else {
    p_values$annotation <- ifelse(p_values$p_value < 0.001, paste0('***'), 
                                  ifelse(p_values$p_value < 0.01, paste0('**'), 
                                         ifelse(p_values$p_value < 0.05, paste0('*'), 
                                                paste0('ns'))))
  }
  
  #fwrite(p_values, "Ttest_pvalue_raw.csv")
  
  data_long <- merge(data_long, p_values, by = 'DMR')
  
  # Create the boxplot with shared y-axis and separate comparisons for each target_site
  p <- ggplot(data_long, aes(x = `DMR`, y = `Methylation values`, fill = Severity)) +
    geom_boxplot() +
    facet_wrap(~ `DMR`, scales = "free_x", nrow = 1) + # All facets in one row
    theme_bw()+ #背景变为白色
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), #不显示网格线
      panel.grid.minor = element_blank(),#不显示网格线
      strip.background = element_blank(),
      
      legend.position = "right",
      panel.spacing = unit(0.5, "lines"), # Adjust space between facets if needed
      # strip.text.x = element_text(margin = margin(2,0,2,0)),
      strip.text.x = element_blank(), # Hide facet label text
      panel.background = element_rect(colour = "black", linewidth = 0.5),#添加边框
      
      axis.text.x=element_text(colour="black",family="Arial",face="bold", size=10, angle = 60, hjust = 1), #设置x轴刻度标签的字体属性
      axis.title.x=element_text(family="Arial",size = 15,face="bold"), #设置x轴的标题的字体属性
      axis.text.y=element_text(family="Arial",size=10,face="bold"), #设置x轴刻度标签的字体属性
      axis.title.y=element_text(family="Arial",size = 15,face="bold"), #设置y轴的标题的字体属性
      legend.text=element_text(colour="black",family="Arial",size=10),#图例分类标签设置
      legend.title=element_text(colour="black",family="Arial",size=15),#图例标题设置
      
      plot.title = element_text(family="Arial",size=23,face="bold",hjust = 0.5), #设置总标题的字体属性
    ) +
    labs(title = feature_names) +
    scale_fill_manual(name = feature_names,
      values = c("1" = "black", "0" = "grey80", "2" = "grey50"))+
    coord_cartesian(ylim = c(0, 1), clip = "off") # Extend the y-axis without clipping
  # Note: You would need to adjust the y position and label accordingly
  p <- p + geom_text(aes(label = annotation, y = 1.07), family="Arial", color = 'black', size = 4, vjust = 0, data = unique(data_long[, c('Severity', 'annotation', "DMR")]))
  return(list(p, p_values))
}


# setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
setwd("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
# All_sig_DMR_methy_group_tumor_infor = fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")
All_sig_DMR_methy_group_tumor_infor = fread("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")
common_XGBoost_RF_predictor_target <- fread("common_XGBoost_RF_predictor_target.csv")



# summarize the plots and pvalues
plot_list = list()
p_values_list <- list()
for (feature in c("recurrence", "Death", "Bad_outcome", "PSA_recurrence", "surg_lnmets")) {
  boxplot_data <- subset(All_sig_DMR_methy_group_tumor_infor, select = c("patient", feature, 
                                                                         common_XGBoost_RF_predictor_target$Common_DMR_predictor_target))
  if (feature == "Death") {
    boxplot_data$Death[boxplot_data$Death == "Alive"] = 0
    boxplot_data$Death[boxplot_data$Death == "Dead of disease"] = 1
    boxplot_data$Death[boxplot_data$Death == "Dead of other cause"] = 2
  }
  
  names(boxplot_data)[which(names(boxplot_data) == feature)] = "severity"
  plot_and_pvalue = box_plot_site_annotation(raw_data = boxplot_data, pvalue = F, feature_names = feature)
  plot_list[[feature]] = plot_and_pvalue[[1]]
  p_values_list[[feature]] = plot_and_pvalue[[2]]
}

# output the plots 
combined_plot <- wrap_plots(plot_list, ncol = 1)  # Combine plots in a single column
pdf("boxplot_for_survival_annotation.pdf", width = 10, height = 20)
combined_plot
dev.off()
# output the pvalues 
combined_pvalue <- p_values_list[[1]]$p_value
names(combined_pvalue) = p_values_list[[1]]$DMR
for (i in 2:length(p_values_list)) {
  combined_pvalue <- rbind(combined_pvalue, p_values_list[[i]]$p_value)
}
rownames(combined_pvalue) = names(p_values_list)
fwrite(combined_pvalue, "Pvalue_for_boxplot_for_survival_annotation.csv")


# 6.3 PCA plots ----
library(patchwork)
library(ggplot2)
complexed_PCA_plot <- function(raw_data, title = NULL){
  df_for_PCA <- subset(raw_data, select = names(raw_data)[!names(raw_data) %in% c("patient", "Bad_outcome", "risky")])
  PCA <- prcomp(df_for_PCA)
  PC_plot_df <- data.frame(PCA$x, Bad_outcome = raw_data$Bad_outcome, risky = raw_data$risky)
  summ <- summary(PCA)
  PCA_variance <- PCA$sdev^2
  PCA_variance_per <- round(PCA_variance/sum(PCA_variance)*100,2)
  length(PCA_variance_per)
  xlab <- paste0("PC1(", round(summ$importance[2,1]*100,2), "%)")
  ylab <- paste0("PC2(", round(summ$importance[2,2]*100,2), "%)")
  
  # Combine color and shape into a single legend
  
  COLS <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  plot <- ggplot() +
    stat_ellipse(data = PC_plot_df, aes(x = PC1, y = PC2, fill = Bad_outcome),
                 type = "norm", geom = "polygon", alpha = 0.2, color = NA) +
    scale_fill_manual(name = "Tumor Progression",
                      values = c("1" = "black", "0" = "grey")) + # Adjust the ellipse color
    geom_point(data = PC_plot_df, aes(x = PC1, y = PC2, color = risky), shape = 1, size = 5) +
    scale_colour_manual(name = "Risk",
                        values = c("Low" = "#5b8efd", "Middle" = "goldenrod1", "High"= "#dd217d") 
    ) + 
    #    labs(title = title, x = xlab, y = ylab, color = "CDX2 expression\nZ-score", shape = "Bad_outcome") +
    labs(title = title, x = xlab, y = ylab) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), # 不显示网格线
      panel.grid.minor = element_blank(), # 不显示网格线
      strip.background = element_blank(),
      plot.background = element_blank(), # 无背景
      legend.position = "right",
      panel.background = element_rect(colour = "black", linewidth = 0.5), # 添加边框
      axis.text.x = element_text(colour = "black", face = "bold", size = 10, angle = 0, hjust = 0.5), # 设置x轴刻度标签的字体属性
      axis.title.x = element_text(size = 15, face = "bold"), # 设置x轴的标题的字体属性
      axis.text.y = element_text(size = 10, face = "bold"), # 设置y轴刻度标签的字体属性
      axis.title.y = element_text(size = 15, face = "bold"), # 设置y轴的标题的字体属性
      legend.text = element_text(colour = "black", size = 10), # 图例分类标签设置
      legend.title = element_text(colour = "black", size = 15), # 图例标题设置
      plot.title = element_text(size = 23, face = "bold", hjust = 0.5) # 设置总标题的字体属性
    )
  
  return(plot)
}

PCA_plot <- function(raw_data, feature_names = "Severity"){
  PCA <- prcomp(raw_data[,-c("patient", "severity")])
  PC_plot_df <- data.frame(PCA$x, severity = raw_data$severity)
  summ <- summary(PCA)
  PCA_variance <- PCA$sdev^2
  PCA_variance_per <- round(PCA_variance/sum(PCA_variance)*100,2)
  length(PCA_variance_per)
  xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
  ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
  
  plot = ggplot(data = PC_plot_df,aes(x=PC1,y=PC2,color=severity))+
    stat_ellipse(aes(fill=severity),
                 type = "norm", geom ="polygon",alpha=0.2,color=NA)+
    geom_point()+labs(title = "", x=xlab,y=ylab,color="")+
    scale_fill_manual(name = feature_names, 
                      values = c("#dd217d","goldenrod1", "#5b8efd"))+   
    scale_colour_manual(name = feature_names, 
                        values = c("#dd217d","goldenrod1", "#5b8efd"))+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(), #??????????????????
      panel.grid.minor = element_blank(),#??????????????????
      strip.background = element_blank(),
      plot.background = element_blank(), #?????????
      
      legend.position = "right",
      panel.background = element_rect(colour = "black", linewidth = 0.5),#????????????
      
      axis.text.x=element_text(#family="Times",
        colour="black",face="bold", size=10, angle = 0, hjust = 0.5), #??????x??????????????????????????????
      axis.title.x=element_text(#family="Times",
        size = 15,face="bold"), #??????x???????????????????????????
      axis.text.y=element_text(#family="Times",
        size=10,face="bold"), #??????x??????????????????????????????
      axis.title.y=element_text(#family="Times",
        size = 15,face="bold"), #??????y???????????????????????????
      legend.text=element_text(#family="Times",
        colour="black",size=10),#????????????????????????
      legend.title=element_text(#family="Times",
        colour="black",size=15),#??????????????????
      
      plot.title = element_text(#family="Times",
        size=23,face="bold",hjust = 0.5), #??????????????????????????????
    )
  return(plot)
}

# 6.3.1 plot a complexed PCA ----
setwd("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")

common_XGBoost_RF_predictor_target <- fread("common_XGBoost_RF_predictor_target.csv")
# All_sig_DMR_methy_group_tumor_infor = fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")
All_sig_DMR_methy_group_tumor_infor = fread("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")

common_XGBoost_RF_predictor_target <- subset(All_sig_DMR_methy_group_tumor_infor, select = c("patient", 
                                                                                             "Bad_outcome", "risky", common_XGBoost_RF_predictor_target$Common_DMR_predictor_target))
common_XGBoost_RF_predictor_target$Bad_outcome <- factor(common_XGBoost_RF_predictor_target$Bad_outcome)
str(common_XGBoost_RF_predictor_target)
complext_PCAplot = complexed_PCA_plot(common_XGBoost_RF_predictor_target)

pdf("complext_PCA_plot_survival_annotation.pdf", width = 6, height = 4)
complext_PCAplot
dev.off()

# 6.3.2 plot general PCA plots ----
setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
common_XGBoost_RF_predictor_target <- fread("common_XGBoost_RF_predictor_target.csv")
All_sig_DMR_methy_group_tumor_infor = fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")

plot_list = list()
for (feature in c("recurrence", "Death", "Bad_outcome", "PSA_recurrence", "surg_lnmets")) {
  plot_data <- subset(All_sig_DMR_methy_group_tumor_infor, select = c("patient", feature, 
                                                                         common_XGBoost_RF_predictor_target$Common_DMR_predictor_target))

  names(plot_data)[which(names(plot_data) == feature)] = "severity"
  plot_data$severity <- factor(plot_data$severity)
  PCAplot = PCA_plot(raw_data = plot_data, feature_names = feature)
  plot_list[[feature]] = PCAplot

}
combined_plot <- wrap_plots(plot_list, ncol = 3)  # Combine plots in a single column

pdf("merged_PCA_plot_survival_annotation.pdf", width = 18, height = 10)
combined_plot
dev.off()



# 7.0 Prediction on middle group ----
library(circlize)
library(readxl)
library(readr)
library(caret)
library(pROC)
library(glmmTMB)
library(data.table)
rm(list = ls())
logistic_regression_for_unknown_prediction <- function(train, test, target_var, random_effects_var, other_not_predictor_variables, identification_variable,
                                                       num_iterations = 100, family = "binomial"){
  results <- data.frame()
  DMR_predictor <- data.frame(matrix(ncol = 1, nrow = ncol(train)))
  prediction_result <- data.frame(patient = test[[identification_variable]], test_real = test[[target_var]])
  for (k in 1:num_iterations) {
    set.seed(k)
    # level1 <- subset(data, data[[target_var]] == 0)
    # level2 <- subset(data, data[[target_var]] == 1)
    # index_level1 <- sample(nrow(level1), 0.7*nrow(level1))
    # index_level2 <- sample(nrow(level2), 0.7*nrow(level2))
    # train <- data.frame(rbind(level1[index_level1, ], level2[index_level2, ]))
    # test <- data.frame(rbind(level1[-index_level1, ], level2[-index_level2, ]))
    
    train_x <- subset(train, select = names(train)[-which(names(train) %in% c(target_var, random_effects_var, other_not_predictor_variables))])
    train_y <- as.numeric(train[[target_var]])
    test_x <- subset(test, select = names(test)[-which(names(test) %in% c(target_var, random_effects_var, other_not_predictor_variables))])
    test_y <- as.numeric(test[[target_var]])
    
    scaled_train_x <- scale(train_x)
    scaled_test_x <- scale(test_x)
    
    train <- cbind(subset(train, select = c(target_var, random_effects_var)), scaled_train_x)
    test <- cbind(subset(test, select = c(target_var, random_effects_var)), scaled_test_x)
    str(train)
    
    predictors <- colnames(train)[-which(names(train) %in% c(target_var, random_effects_var))]
    formula_str <- paste0(target_var, "~", paste(predictors, collapse = " + ")) #, "+ (1 | ", random_effects_var, ")"
    model_formula <- as.formula(formula_str)
    
    fit <- glmmTMB(model_formula, data = train, family=family)
    print(summary(fit))
    
    
    # Store results
    coef_values <- fixef(fit)[["cond"]][!is.na(fixef(fit)[["cond"]])]
    coef_names <- names(coef_values)
    
    fit_summary <- summary(fit)
    p_values <- fit_summary$coefficients$cond[, "Pr(>|z|)"]
    odds_ratios <- exp(coef_values)  # Calculate odds ratios by exponentiating the coefficients
    ci <- exp(confint(fit, parm = coef_names, level = 0.95, method = "Wald"))
    # coef_names <- c(coef_names, rep(NA, nrow(DMR_predictor) - length(coef_names)))
    # coef_values <- c(coef_values, rep(NA, nrow(DMR_predictor) - length(coef_values)))
    iteration_coef <- data.frame(coef_names, coef_values, p_values, odds_ratios, OR_lower = ci[, 1], OR_Upper = ci[, 2])
    DMR_predictor <- cbind(DMR_predictor[1:length(coef_values),], iteration_coef)
    
    # Prediction 
    train_pre <- predict(fit, newdata = train, type='response', re.form=NA)
    train_modelroc <- roc(train[[target_var]], as.numeric(train_pre)) #计算ROC绘图数据
    roc_train = train_modelroc$auc
    tab_train <- table(train_y, ifelse(train_pre < 0.5, 0, 1), dnn = c("true", "pre"))
    accuracy_train <- sum(diag(tab_train)) / length(train_y)
    
    test_pre <- predict(fit, newdata = test, type='response', re.form=NA)
    test_pre_01 <- ifelse(test_pre < 0.5, 0, 1)
    prediction_result <- cbind(prediction_result, test_pre_01)
    
    ## ***** For prediction ability summary
    test_modelroc <- roc(test_y, as.numeric(test_pre)) # Calculate ROC for test data
    roc_test <- test_modelroc$auc
    tab_test <- table(test_y, ifelse(test_pre < 0.5, 0, 1), dnn = c("true", "pre"))
    accuracy_test <- sum(diag(tab_test)) / length(test_y)
    
    # pdf("Logistic_ROC.pdf", width = 5, height = 5)
    # plot(test_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
    #      grid.col=c("green", "red"), max.auc.polygon=TRUE,
    #      auc.polygon.col="skyblue", print.thres=TRUE) #绘制ROC曲线图
    # plot(ci(test_modelroc, of = "thresholds", thresholds = "best")) #添加最佳截断值位置
    # dev.off()
    
    ## *****
    
    result <- data.frame(accuracy_train, accuracy_test, roc_train, roc_test)
    print(result)
    print(k)
  }
  return(list(prediction_result, DMR_predictor, test_modelroc, result)) #test_modelroc will be from the last iteration
}


# 7.1 Tumor Progression ----
# setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
setwd("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")

common_XGBoost_RF_predictor_target <- fread("common_XGBoost_RF_predictor_target.csv")
# All_sig_DMR_methy_group_tumor_infor = fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")
All_sig_DMR_methy_group_tumor_infor = fread("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")



common_XGBoost_RF_predictor_target <- subset(All_sig_DMR_methy_group_tumor_infor, select = c("patient", "surg_pgleasnt", "surg_pstagt", "recurrence", "Death", "Bad_outcome", "PSA_recurrence", "surg_lnmets",
                                                                                             "risky", common_XGBoost_RF_predictor_target$Common_DMR_predictor_target))
str(common_XGBoost_RF_predictor_target)
# 7.1.1 By signature ----
# common_XGBoost_RF_predictor_target$Bad_outcome = factor(common_XGBoost_RF_predictor_target$Bad_outcome, levels = c(0, 1), labels = c("Control", "Case"))

train_on_High_Low_and_predic_Middle_by_signature <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                          test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                          target_var = "Bad_outcome",
                                                                                          random_effects_var = "patient", 
                                                                                          other_not_predictor_variables = c("surg_pgleasnt", "surg_pstagt", "subtype", "recurrence", "surg_lnmets", "Death", "PSA_recurrence", "risky"), 
                                                                                          identification_variable = "patient",
                                                                                          num_iterations = 1, family = "binomial")

signature_ROC = train_on_High_Low_and_predic_Middle_by_signature[[3]]
train_on_High_Low_and_predic_Middle_by_signature[[2]]
table(train_on_High_Low_and_predic_Middle_by_signature[[1]]$test_pre_01)

# *** make heatmap with prediction ---- 
train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),]
test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),]
test$subtype = train_on_High_Low_and_predic_Middle_by_signature[[1]]$test_pre_01
train$subtype = NA
full_data = rbind(train,test)

library("pheatmap")
library("RColorBrewer")
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

heat_map_for_survival_annotation_with_multi_label_and_certain_order_and_prediction <- function(raw_data, scale_or_not = T, show_row_names = T, show_column_names = T, 
                                                                                clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value){
  #raw_data = common_XGBoost_RF_predictor_target
  #rownames(raw_data) <- raw_data$patient
  #raw_data$surg_pgleasnt = factor(raw_data$surg_pgleasnt, levels = c(0,1), labels = c("Blood Healthy","Blood Tumor"))
  raw_data$risky <- factor(raw_data$risky, 
                           levels = c("Low", "Middle", "High"))
  raw_data$subtype <- factor(raw_data$subtype, 
                           levels = c(0, 1))
  
  raw_data <- raw_data[order(raw_data$risky, 
                             ifelse(raw_data$risky == "Middle", 
                                    raw_data$subtype, 
                                    NA)), ]
  
  data_methy <- subset(raw_data, select = -c(patient,recurrence, Death, Bad_outcome, PSA_recurrence, risky, surg_lnmets, subtype, surg_pgleasnt, surg_pstagt))
  if (scale_or_not == T) {data <- t(apply(data_methy, 2, scale))} else {data = t(data_methy)}
  colnames(data) <- raw_data$patient
  
  
  # set color
  #set main color
  if (scale_or_not == T) { 
    myBreaks <- seq(-5, 5, length.out = 100)
    myCol <- colorRampPalette(c('#712b8f', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  } else {
    myBreaks <- seq(0, 1, length.out = 100)
    myCol <- colorRampPalette(c('#712b8f', 'white', 'darksalmon'))(100)
    col = colorRamp2(myBreaks, myCol)
  }
  
  #set annotation color
  # surg_pgleasnt <- raw_data$surg_pgleasnt
  # pick.col <- brewer.pal(9, 'Greens') # 可以设置一个主色（比如绿色），然后从中挑出对应的颜色.allowed maximum for palette Greens is 9
  # col.surg_pgleasnt <- colorRampPalette(pick.col)(length(unique(surg_pgleasnt)))
  
  # create annotation object
  ann <- data.frame(Tumor_Recurrence = raw_data$recurrence, 
                    Lymph_Node_Mets = raw_data$surg_lnmets,
                    PSA_recurrence = raw_data$PSA_recurrence, 
                    All_Death = raw_data$Death, 
                    Tumor_Progression = raw_data$Bad_outcome,
                    Risk = raw_data$risky,
                    subtype = raw_data$subtype
                    )
  colors=list(Tumor_Recurrence = c("1" = "black", "0" = "gray80"),
              Lymph_Node_Mets = c("1" = "black", "0" = "gray80"),
              PSA_recurrence = c("1" = "black", "0" = "gray80"),
              All_Death = c("Dead of disease" = "black",  "Dead of other cause" = "gray50", "Alive" = "gray80"),
              Tumor_Progression = c("1" = "black", "0" = "gray80"),
              Risk = c("Low" = "#5b8efd", "Middle" = "goldenrod1", "High"= "#dd217d"),
              subtype = c("1" = "#87CEFA", "0" = "#FFC0CB")
              )
  
  colAnn <- HeatmapAnnotation(
    df = ann,
    col = colors,
    which = 'col', # set 'col' (samples) or 'row' (gene) annotation
    na_col = 'white', # NA color, default is white
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'),
    # Correct the legend parameter key to match the column names in `ann`
    annotation_legend_param = list(
      Risk = list(
        nrow = 2, # Display legend in 2 rows
        title = 'Risk', # Set the title of the legend
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 10, fontface = 'bold'),
        labels_gp = gpar(fontsize = 10, fontface = 'plain')
      )
    )
  )
  
  # create bottom boxplot
  boxplotCol <- HeatmapAnnotation(
    boxplot = anno_boxplot(
      data,
      border = FALSE,
      gp = gpar(fill = '#CCCCCC'),
      pch = '.',
      size = unit(1, 'mm'),
      axis = TRUE,
      axis_param = list(
        gp = gpar(fontsize = 10),
        side = 'left')),
    annotation_width = unit(c(1.0), 'cm'),
    which = 'col')
  
  # Mark the rows we want
  genelabels <- rowAnnotation(
    Genes = anno_mark(
      at = seq(1, nrow(data), 1),
      labels = rownames(data)[seq(1, nrow(data), 1)],
      labels_gp = gpar(fontsize = 5, fontface = 'plain'),
      padding = 0.75),
    width = unit(2.0, 'cm') +
      
      max_text_width(
        rownames(data)[seq(1, nrow(data), 1)],
        gp = gpar(fontsize = 5,  fontface = 'plain')))
  
  
  # plot the heatmap
  hmap=Heatmap(data,
               # split the genes / rows according to the PAM clusters
               #split = NA, # kmeans_clus / hc_clus
               cluster_row_slices = F,
               column_split = raw_data$risky, # 按subtype分组
               
               name = 'Methylation\nZ-score',
               
               col = colorRamp2(myBreaks, myCol),
               
               # parameters for the colour-bar that represents gradient of expression
               heatmap_legend_param = list(
                 color_bar = 'continuous',
                 legend_direction = 'vertical',
                 legend_width = unit(12, 'cm'),
                 legend_height = unit(5.0, 'cm'),
                 title_position = 'topleft',
                 title_gp=gpar(fontsize = 10, fontface = 'bold'),
                 labels_gp=gpar(fontsize = 10, fontface = 'bold')),
               
               
               # row (gene) parameters
               cluster_rows = T,
               show_row_dend = TRUE,
               row_title = '',
               row_title_side = 'left',
               row_title_gp = gpar(fontsize = 12,  fontface = 'plain'),
               row_title_rot = 90,
               show_row_names = show_row_names,
               row_names_gp = gpar(fontsize = 10, fontface = 'plain'),
               row_names_side = 'right',
               row_dend_width = unit(20,'mm'),
               
               # column (sample) parameters
               cluster_columns = T,
               show_column_dend = TRUE,
               column_title = paste(clustering_distance_columns_value, clustering_method_columns_value, clustering_distance_rows_value, clustering_method_rows_value),
               column_title_side = 'top',
               column_title_gp = gpar(fontsize = 12, fontface = 'plain'),
               column_title_rot = 0,
               show_column_names = show_column_names,
               column_names_gp = gpar(fontsize = 10, fontface = 'plain'),
               column_names_side = "top",
               column_names_max_height = unit(20, 'cm'),
               column_dend_height = unit(2,'cm'),
               
               # cluster methods for rows and columns
               clustering_distance_columns = clustering_distance_columns_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_columns = clustering_method_columns_value, #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               clustering_distance_rows = clustering_distance_rows_value,#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
               clustering_method_rows = clustering_method_rows_value,#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
               
               # specify top and bottom annotations
               top_annotation = colAnn,
               bottom_annotation = boxplotCol)
  
  return(hmap)
  
}


Heatmap_with_prediction_and_certain_order <- heat_map_for_survival_annotation_with_multi_label_and_certain_order_and_prediction(raw_data = full_data,
                                                                                                                 scale_or_not = T,show_row_names = T, show_column_names = T,
                                                                                                                 clustering_distance_columns_value = "spearman",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                                                 clustering_method_columns_value = 'ward.D', #ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                                                                                                 clustering_distance_rows_value = "pearson",#euclidean, maximum, manhattan, canberra, binary, minkowski, pearson, spearman, kendall
                                                                                                                 clustering_method_rows_value = 'ward.D'#ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
                                                                                                                 )

pdf("Heatmap_with_prediction_and_certain_order.pdf", width = 19, height = 7)
Heatmap_with_prediction_and_certain_order
dev.off()

# 7.1.1.1 Risk-score----
coef_for_signature <- train_on_High_Low_and_predic_Middle_by_signature[[2]][, -1]
# Define the coefficients and intercept
coef_values <- coef_for_signature$coef_values
names(coef_values) <- coef_for_signature$coef_names
# Separate intercept from other coefficients
intercept <- coef_values["(Intercept)"]
dmr_coefficients <- coef_values[!names(coef_values) %in% "(Intercept)"]
dmr_coefficients <- dmr_coefficients[!is.na(dmr_coefficients)]

# 7.1.1.2 Risk-score calculation ----
str(common_XGBoost_RF_predictor_target)
dmr_matrix <- as.matrix(common_XGBoost_RF_predictor_target[, names(dmr_coefficients), with = FALSE])
coef_vector <- as.numeric(dmr_coefficients)
common_XGBoost_RF_predictor_target$risk_score <- intercept + dmr_matrix %*% coef_vector
fwrite(common_XGBoost_RF_predictor_target, "All_sig_DMR_methy_group_tumor_infor_with_risk_score.csv")

# 7.1.1.3 Risk-Score Distribution ----
# setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
setwd("D:/OneDrive - Washington University in St. Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")

common_XGBoost_RF_predictor_target = fread("All_sig_DMR_methy_group_tumor_infor_with_risk_score.csv")
common_XGBoost_RF_predictor_target = common_XGBoost_RF_predictor_target[order(common_XGBoost_RF_predictor_target$risk_score),]
common_XGBoost_RF_predictor_target$Bad_outcome = factor(common_XGBoost_RF_predictor_target$Bad_outcome)

plot_data <- common_XGBoost_RF_predictor_target %>%
  arrange(risk_score)
str(plot_data)
# Create the plot
plot_data$risky <- factor(plot_data$risky, 
                         levels = c("Low", "Middle", "High"))
plot_data <- plot_data[order(plot_data$risky), ]  # 按新因子顺序排序数据

ggplot(plot_data, 
                  aes(x = seq_along(risk_score), y = risk_score, fill = Bad_outcome)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(name = "Tumor Progression",
                    values = c("blue", "red"),
                    labels = c("No","Yes"))+ # adjust the legend
  # scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
  labs(
    title = "",
    x = "Paitents", y = "Risk Score", fill = "Tumor Progression"
  ) +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),  # Hide x-axis text to avoid clutter
    axis.text.y = element_text(size = 14, color = "black"),  # Increase font size and change color to black
    axis.title = element_text(size = 20, color = "black"),  # Increase axis title size and color
    legend.title = element_text(#family="Times",
      size = 16,face="bold"), # Adjust legend title size here
    legend.text = element_text(#family="Times",
      size = 14,face="plain"), # Adjust legend keys (levels) text size here
    panel.border = element_blank(),
    axis.line.y = element_line(size = 0.5, color = "black"),  # Make axes lines thicker and black
    panel.grid = element_blank(),  # Remove grid lines
    plot.title = element_text(hjust = 0.5, size = 20, color = "black")  # Increase title size and color
  ) +
  scale_y_continuous(limits = c(-60, 40),
                     #y轴刻度，seq(0,2,0.2)函数，分别是最小值0，最大值2，间距0.2
                     breaks = seq(-60, 40, 20))

ggsave("Risk_score_distribution.pdf", width = 14, height = 4)


# 7.1.2 By Risk score ----
train_on_High_Low_and_predic_Middle_by_RS <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                        test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                        target_var = "Bad_outcome",
                                                                                        random_effects_var = "patient", 
                                                                                        other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("risk_score", "patient")], 
                                                                                        identification_variable = "patient",
                                                                                        num_iterations = 1, family = "binomial")

RS_ROC = train_on_High_Low_and_predic_Middle_by_RS[[3]]
train_on_High_Low_and_predic_Middle_by_RS[[1]]


# 7.1.2 By gleason score ----
common_XGBoost_RF_predictor_target$surg_pgleasnt[common_XGBoost_RF_predictor_target$surg_pgleasnt > 7] = 1
common_XGBoost_RF_predictor_target$surg_pgleasnt[common_XGBoost_RF_predictor_target$surg_pgleasnt == 7] = 0

train_on_High_Low_and_predic_Middle_by_GS <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                        test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                        target_var = "Bad_outcome",
                                                                                        random_effects_var = "patient", 
                                                                                        other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("surg_pgleasnt", "patient")], 
                                                                                        identification_variable = "patient",
                                                                                        num_iterations = 1, family = "binomial")

GS_ROC = train_on_High_Low_and_predic_Middle_by_GS[[3]]
train_on_High_Low_and_predic_Middle_by_GS[[1]]

# 7.1.3 By tumor T stage ----
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2a"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2b"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2c"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "3a"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "3b"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "4"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt = as.numeric(common_XGBoost_RF_predictor_target$surg_pstagt)

train_on_High_Low_and_predic_Middle_by_Tstage <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                            test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                            target_var = "Bad_outcome",
                                                                                            random_effects_var = "patient", 
                                                                                            other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("surg_pstagt", "patient")], 
                                                                                            identification_variable = "patient",
                                                                                            num_iterations = 1, family = "binomial")

Tstage_ROC = train_on_High_Low_and_predic_Middle_by_Tstage[[3]]

# 7.1.4 ROC curve with multiple lines ----

pdf("ROC_for_tumor_progression.pdf", width = 5, height = 5)
# Plot the first ROC curve
plot(signature_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "blue", max.auc.polygon = FALSE, print.thres = FALSE, main = "Tumor Progression",
     cex.axis = 1,
     cex.lab = 1.2,
     font.lab = 1, #1 = Plain, 2 = Bold, 3 = Italic, 4 = Bold Italic
     font.main = 2,
)
# Add the second ROC curve
# plot(RS_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
#      col = "red", max.auc.polygon = FALSE, print.thres = FALSE, add = TRUE)
# Add the second ROC curve
plot(GS_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "purple", max.auc.polygon = FALSE, print.thres = FALSE, add = TRUE)
# Add the third ROC curve
plot(Tstage_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "green", max.auc.polygon = FALSE, print.thres = FALSE, add = TRUE)
# Add a legend for the ROC curves
legend("bottomright", title = "AUC:",
       legend = c(paste0("Signature: ", round(auc(signature_ROC), 2)), 
                  #paste0("Risk Score: ", round(auc(RS_ROC), 2)), 
                  paste0("Gleason Score: ", round(auc(GS_ROC), 2)),
                  paste0("T Stage: ", round(auc(Tstage_ROC), 2))),
       col = c("blue", #"red", 
               "purple", "green"), lwd = 2, cex = 0.8)
dev.off()

test_risk_gs <- roc.test(signature_ROC, GS_ROC, method = "delong")
test_risk_tstage <- roc.test(signature_ROC, Tstage_ROC, method = "delong")
test_risk_gs$p.value
test_risk_tstage$p.value


# 7.2 recurrence ----
setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
common_XGBoost_RF_predictor_target <- fread("common_XGBoost_RF_predictor_target.csv")
All_sig_DMR_methy_group_tumor_infor = fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")
common_XGBoost_RF_predictor_target <- subset(All_sig_DMR_methy_group_tumor_infor, select = c("patient", "surg_pgleasnt", "surg_pstagt", "recurrence", "Death", "Bad_outcome", "PSA_recurrence", "surg_lnmets",
                                                                                             "risky", common_XGBoost_RF_predictor_target$Common_DMR_predictor_target))

str(All_sig_DMR_methy_group_tumor_infor)
# 7.2.1 By signature ----

train_on_High_Low_and_predic_Middle_by_risk <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                              test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                              target_var = "recurrence",
                                                                              random_effects_var = "patient", 
                                                                              other_not_predictor_variables = c("surg_pgleasnt", "surg_pstagt", "subtype", "Death", "Bad_outcome", "PSA_recurrence", "surg_lnmets", "risky"), 
                                                                              identification_variable = "patient",
                                                                              num_iterations = 1, family = "binomial")

signature_ROC = train_on_High_Low_and_predic_Middle_by_risk[[3]]
train_on_High_Low_and_predic_Middle_by_risk[[1]]
table(train_on_High_Low_and_predic_Middle_by_risk[[1]]$test_pre_01)

# 7.2.2 By gleason score ----
common_XGBoost_RF_predictor_target$surg_pgleasnt[common_XGBoost_RF_predictor_target$surg_pgleasnt > 7] = 1
common_XGBoost_RF_predictor_target$surg_pgleasnt[common_XGBoost_RF_predictor_target$surg_pgleasnt == 7] = 0


train_on_High_Low_and_predic_Middle_by_GS <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                  test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                  target_var = "recurrence",
                                                                                  random_effects_var = "patient", 
                                                                                  other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("surg_pgleasnt", "patient")], 
                                                                                  identification_variable = "patient",
                                                                                  num_iterations = 1, family = "binomial")

GS_ROC = train_on_High_Low_and_predic_Middle_by_GS[[3]]
train_on_High_Low_and_predic_Middle_by_GS[[2]]

# 7.2.3 By tumor T stage ----
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2a"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2b"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2c"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "3a"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "3b"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "4"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt = as.numeric(common_XGBoost_RF_predictor_target$surg_pstagt)

train_on_High_Low_and_predic_Middle_by_Tstage <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                        test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                        target_var = "recurrence",
                                                                                        random_effects_var = "patient", 
                                                                                        other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("surg_pstagt", "patient")], 
                                                                                        identification_variable = "patient",
                                                                                        num_iterations = 1, family = "binomial")

Tstage_ROC = train_on_High_Low_and_predic_Middle_by_Tstage[[3]]

# 7.2.4 ROC curve with multiple lines ----
roc_data <- bind_rows(
  data.frame(sensitivity = signature_ROC$sensitivities, specificity = signature_ROC$specificities, Model = "Risk", AUC = rep(auc(signature_ROC))),
  data.frame(sensitivity = GS_ROC$sensitivities, specificity = GS_ROC$specificities, Model = "Gleason Score", AUC = rep(auc(GS_ROC))),
  data.frame(sensitivity = Tstage_ROC$sensitivities, specificity = Tstage_ROC$specificities, Model = "T Stage", AUC = rep(auc(Tstage_ROC)))
)


pdf("ROC_for_Recurrence.pdf", width = 5, height = 5)
# Plot the first ROC curve
plot(signature_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "blue", max.auc.polygon = FALSE, print.thres = FALSE, main = "Tumor Recurrence",
     cex.axis = 1,
     cex.lab = 1.2,
     font.lab = 1, #1 = Plain, 2 = Bold, 3 = Italic, 4 = Bold Italic
     font.main = 2,
     )
# Add the second ROC curve
plot(GS_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "purple", max.auc.polygon = FALSE, print.thres = FALSE, add = TRUE)
# Add the third ROC curve
plot(Tstage_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "green", max.auc.polygon = FALSE, print.thres = FALSE, add = TRUE)
# Add a legend for the ROC curves
legend("bottomright", title = "AUC:",
       legend = c(paste0("Signature: ", round(auc(signature_ROC), 2)), 
                  paste0("Gleason Score***: ", round(auc(GS_ROC), 2)),
                  paste0("T Stage***: ", round(auc(Tstage_ROC), 2))),
       col = c("blue", "purple", "green"), lwd = 2, cex = 0.8)
dev.off()

test_risk_gs <- roc.test(signature_ROC, GS_ROC, method = "delong")
test_risk_tstage <- roc.test(signature_ROC, Tstage_ROC, method = "delong")
test_risk_gs$p.value
test_risk_tstage$p.value



# 7.3 Death ----
setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
common_XGBoost_RF_predictor_target <- fread("common_XGBoost_RF_predictor_target.csv")
All_sig_DMR_methy_group_tumor_infor = fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")
common_XGBoost_RF_predictor_target <- subset(All_sig_DMR_methy_group_tumor_infor, select = c("patient", "surg_pgleasnt", "surg_pstagt", "recurrence", "Death", "Bad_outcome", "PSA_recurrence", "surg_lnmets",
                                                                                             "risky", common_XGBoost_RF_predictor_target$Common_DMR_predictor_target))

str(All_sig_DMR_methy_group_tumor_infor)
# 7.3.1 By signature ----
table(common_XGBoost_RF_predictor_target$Death)
common_XGBoost_RF_predictor_target$Death[common_XGBoost_RF_predictor_target$Death == "Alive"] = 0
common_XGBoost_RF_predictor_target$Death[common_XGBoost_RF_predictor_target$Death == "Dead of disease"] = 1
common_XGBoost_RF_predictor_target$Death[common_XGBoost_RF_predictor_target$Death == "Dead of other cause"] = 1 #There is no "Dead of disease" in "middle risk" group, so we have to keep all dead samples
common_XGBoost_RF_predictor_target$Death = as.numeric(common_XGBoost_RF_predictor_target$Death)

train_on_High_Low_and_predic_Middle_by_risk <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                          test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                          target_var = "Death",
                                                                                          random_effects_var = "patient", 
                                                                                          other_not_predictor_variables = c("surg_pgleasnt", "surg_pstagt", "subtype", "recurrence", "Bad_outcome", "PSA_recurrence", "surg_lnmets", "risky"), 
                                                                                          identification_variable = "patient",
                                                                                          num_iterations = 1, family = "binomial")

signature_ROC = train_on_High_Low_and_predic_Middle_by_risk[[3]]
train_on_High_Low_and_predic_Middle_by_risk[[1]]

# 7.3.2 By gleason score ----
common_XGBoost_RF_predictor_target$surg_pgleasnt[common_XGBoost_RF_predictor_target$surg_pgleasnt > 7] = 1
common_XGBoost_RF_predictor_target$surg_pgleasnt[common_XGBoost_RF_predictor_target$surg_pgleasnt == 7] = 0

train_on_High_Low_and_predic_Middle_by_GS <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                        test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                        target_var = "Death",
                                                                                        random_effects_var = "patient", 
                                                                                        other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("surg_pgleasnt", "patient")], 
                                                                                        identification_variable = "patient",
                                                                                        num_iterations = 1, family = "binomial")

GS_ROC = train_on_High_Low_and_predic_Middle_by_GS[[3]]
train_on_High_Low_and_predic_Middle_by_GS[[1]]

# 7.3.3 By tumor T stage ----
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2a"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2b"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2c"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "3a"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "3b"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "4"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt = as.numeric(common_XGBoost_RF_predictor_target$surg_pstagt)

train_on_High_Low_and_predic_Middle_by_Tstage <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                            test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                            target_var = "Death",
                                                                                            random_effects_var = "patient", 
                                                                                            other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("surg_pstagt", "patient")], 
                                                                                            identification_variable = "patient",
                                                                                            num_iterations = 1, family = "binomial")

Tstage_ROC = train_on_High_Low_and_predic_Middle_by_Tstage[[3]]

# 7.3.4 ROC curve with multiple lines ----
roc_data <- bind_rows(
  data.frame(sensitivity = signature_ROC$sensitivities, specificity = signature_ROC$specificities, Model = "Risk", AUC = rep(auc(signature_ROC))),
  data.frame(sensitivity = GS_ROC$sensitivities, specificity = GS_ROC$specificities, Model = "Gleason Score", AUC = rep(auc(GS_ROC))),
  data.frame(sensitivity = Tstage_ROC$sensitivities, specificity = Tstage_ROC$specificities, Model = "T Stage", AUC = rep(auc(Tstage_ROC)))
)


pdf("ROC_for_Death.pdf", width = 5, height = 5)
# Plot the first ROC curve
plot(signature_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "blue", max.auc.polygon = FALSE, print.thres = FALSE, main = "Death",
     cex.axis = 1,
     cex.lab = 1.2,
     font.lab = 1, #1 = Plain, 2 = Bold, 3 = Italic, 4 = Bold Italic
     font.main = 2,
)
# Add the second ROC curve
plot(GS_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "purple", max.auc.polygon = FALSE, print.thres = FALSE, add = TRUE)
# Add the third ROC curve
plot(Tstage_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "green", max.auc.polygon = FALSE, print.thres = FALSE, add = TRUE)
# Add a legend for the ROC curves
legend("bottomright", title = "AUC:",
       legend = c(paste0("Signature: ", round(auc(signature_ROC), 2)), 
                  paste0("Gleason Score: ", round(auc(GS_ROC), 2)),
                  paste0("T Stage: ", round(auc(Tstage_ROC), 2))),
       col = c("blue", "purple", "green"), lwd = 2, cex = 0.8)
dev.off()


test_risk_gs <- roc.test(signature_ROC, GS_ROC, method = "delong")
test_risk_tstage <- roc.test(signature_ROC, Tstage_ROC, method = "delong")
test_risk_gs$p.value
test_risk_tstage$p.value

# 7.4 PSA_recurrence ----
setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
common_XGBoost_RF_predictor_target <- fread("common_XGBoost_RF_predictor_target.csv")
All_sig_DMR_methy_group_tumor_infor = fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")
common_XGBoost_RF_predictor_target <- subset(All_sig_DMR_methy_group_tumor_infor, select = c("patient", "surg_pgleasnt", "surg_pstagt", "recurrence", "Death", "Bad_outcome", "PSA_recurrence", "surg_lnmets",
                                                                                             "risky", common_XGBoost_RF_predictor_target$Common_DMR_predictor_target))

str(All_sig_DMR_methy_group_tumor_infor)
# 7.4.1 By signature ----
table(common_XGBoost_RF_predictor_target$PSA_recurrence)
train_on_High_Low_and_predic_Middle_by_risk <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                          test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                          target_var = "PSA_recurrence",
                                                                                          random_effects_var = "patient", 
                                                                                          other_not_predictor_variables = c("surg_pgleasnt", "surg_pstagt", "subtype", "recurrence", "Bad_outcome", "Death", "surg_lnmets", "risky"), 
                                                                                          identification_variable = "patient",
                                                                                          num_iterations = 1, family = "binomial")

signature_ROC = train_on_High_Low_and_predic_Middle_by_risk[[3]]
train_on_High_Low_and_predic_Middle_by_risk[[1]]

# 7.4.2 By gleason score ----
common_XGBoost_RF_predictor_target$surg_pgleasnt[common_XGBoost_RF_predictor_target$surg_pgleasnt > 7] = 1
common_XGBoost_RF_predictor_target$surg_pgleasnt[common_XGBoost_RF_predictor_target$surg_pgleasnt == 7] = 0

train_on_High_Low_and_predic_Middle_by_GS <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                        test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                        target_var = "PSA_recurrence",
                                                                                        random_effects_var = "patient", 
                                                                                        other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("surg_pgleasnt", "patient")], 
                                                                                        identification_variable = "patient",
                                                                                        num_iterations = 1, family = "binomial")

GS_ROC = train_on_High_Low_and_predic_Middle_by_GS[[3]]
train_on_High_Low_and_predic_Middle_by_GS[[1]]

# 7.4.3 By tumor T stage ----
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2a"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2b"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2c"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "3a"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "3b"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "4"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt = as.numeric(common_XGBoost_RF_predictor_target$surg_pstagt)

train_on_High_Low_and_predic_Middle_by_Tstage <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                            test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                            target_var = "PSA_recurrence",
                                                                                            random_effects_var = "patient", 
                                                                                            other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("surg_pstagt", "patient")], 
                                                                                            identification_variable = "patient",
                                                                                            num_iterations = 1, family = "binomial")

Tstage_ROC = train_on_High_Low_and_predic_Middle_by_Tstage[[3]]

# 7.4.4 ROC curve with multiple lines ----
roc_data <- bind_rows(
  data.frame(sensitivity = signature_ROC$sensitivities, specificity = signature_ROC$specificities, Model = "Risk", AUC = rep(auc(signature_ROC))),
  data.frame(sensitivity = GS_ROC$sensitivities, specificity = GS_ROC$specificities, Model = "Gleason Score", AUC = rep(auc(GS_ROC))),
  data.frame(sensitivity = Tstage_ROC$sensitivities, specificity = Tstage_ROC$specificities, Model = "T Stage", AUC = rep(auc(Tstage_ROC)))
)


pdf("ROC_for_PSA_recurrence.pdf", width = 5, height = 5)
# Plot the first ROC curve
plot(signature_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "blue", max.auc.polygon = FALSE, print.thres = FALSE, main = "PSA Recurrence",
     cex.axis = 1,
     cex.lab = 1.2,
     font.lab = 1, #1 = Plain, 2 = Bold, 3 = Italic, 4 = Bold Italic
     font.main = 2,
)
# Add the second ROC curve
plot(GS_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "purple", max.auc.polygon = FALSE, print.thres = FALSE, add = TRUE)
# Add the third ROC curve
plot(Tstage_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "green", max.auc.polygon = FALSE, print.thres = FALSE, add = TRUE)
# Add a legend for the ROC curves
legend("bottomright", title = "AUC:",
       legend = c(paste0("Signature: ", round(auc(signature_ROC), 2)), 
                  paste0("Gleason Score: ", round(auc(GS_ROC), 2)),
                  paste0("T Stage: ", round(auc(Tstage_ROC), 2))),
       col = c("blue", "purple", "green"), lwd = 2, cex = 0.8)
dev.off()


test_risk_gs <- roc.test(signature_ROC, GS_ROC, method = "delong")
test_risk_tstage <- roc.test(signature_ROC, Tstage_ROC, method = "delong")
test_risk_gs$p.value
test_risk_tstage$p.value


# 7.5 surg_lnmets ----
setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
common_XGBoost_RF_predictor_target <- fread("common_XGBoost_RF_predictor_target.csv")
All_sig_DMR_methy_group_tumor_infor = fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_group_clinics_infor.csv")
common_XGBoost_RF_predictor_target <- subset(All_sig_DMR_methy_group_tumor_infor, select = c("patient", "surg_pgleasnt", "surg_pstagt", "recurrence", "Death", "Bad_outcome", "PSA_recurrence", "surg_lnmets",
                                                                                             "risky", common_XGBoost_RF_predictor_target$Common_DMR_predictor_target))

str(All_sig_DMR_methy_group_tumor_infor)
# 7.5.1 By signature ----
table(common_XGBoost_RF_predictor_target$PSA_recurrence)
train_on_High_Low_and_predic_Middle_by_risk <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                          test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                          target_var = "surg_lnmets",
                                                                                          random_effects_var = "patient", 
                                                                                          other_not_predictor_variables = c("surg_pgleasnt", "surg_pstagt", "subtype", "recurrence", "Bad_outcome", "Death", "PSA_recurrence", "risky"), 
                                                                                          identification_variable = "patient",
                                                                                          num_iterations = 1, family = "binomial")

signature_ROC = train_on_High_Low_and_predic_Middle_by_risk[[3]]
train_on_High_Low_and_predic_Middle_by_risk[[1]]

# 7.5.2 By gleason score ----
common_XGBoost_RF_predictor_target$surg_pgleasnt[common_XGBoost_RF_predictor_target$surg_pgleasnt > 7] = 1
common_XGBoost_RF_predictor_target$surg_pgleasnt[common_XGBoost_RF_predictor_target$surg_pgleasnt == 7] = 0

train_on_High_Low_and_predic_Middle_by_GS <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                        test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                        target_var = "surg_lnmets",
                                                                                        random_effects_var = "patient", 
                                                                                        other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("surg_pgleasnt", "patient")], 
                                                                                        identification_variable = "patient",
                                                                                        num_iterations = 1, family = "binomial")

GS_ROC = train_on_High_Low_and_predic_Middle_by_GS[[3]]
train_on_High_Low_and_predic_Middle_by_GS[[1]]

# 7.5.3 By tumor T stage ----
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2a"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2b"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2c"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "3a"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "3b"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "4"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt = as.numeric(common_XGBoost_RF_predictor_target$surg_pstagt)

train_on_High_Low_and_predic_Middle_by_Tstage <- logistic_regression_for_unknown_prediction(train = common_XGBoost_RF_predictor_target[risky %in% c("High", "Low"),],
                                                                                            test = common_XGBoost_RF_predictor_target[risky %in% c("Middle"),],
                                                                                            target_var = "surg_lnmets",
                                                                                            random_effects_var = "patient", 
                                                                                            other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("surg_pstagt", "patient")], 
                                                                                            identification_variable = "patient",
                                                                                            num_iterations = 1, family = "binomial")

Tstage_ROC = train_on_High_Low_and_predic_Middle_by_Tstage[[3]]

# 7.5.4 ROC curve with multiple lines ----

roc_data <- bind_rows(
  data.frame(sensitivity = signature_ROC$sensitivities, specificity = signature_ROC$specificities, Model = "Risk", AUC = rep(auc(signature_ROC))),
  data.frame(sensitivity = GS_ROC$sensitivities, specificity = GS_ROC$specificities, Model = "Gleason Score", AUC = rep(auc(GS_ROC))),
  data.frame(sensitivity = Tstage_ROC$sensitivities, specificity = Tstage_ROC$specificities, Model = "T Stage", AUC = rep(auc(Tstage_ROC)))
)


pdf("ROC_for_LN_mets.pdf", width = 5, height = 5)
# Plot the first ROC curve
plot(signature_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "blue", max.auc.polygon = FALSE, print.thres = FALSE, main = "Lymph Node Mets",
     cex.axis = 1,
     cex.lab = 1.2,
     font.lab = 1, #1 = Plain, 2 = Bold, 3 = Italic, 4 = Bold Italic
     font.main = 2,
)
# Add the second ROC curve
plot(GS_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "purple", max.auc.polygon = FALSE, print.thres = FALSE, add = TRUE)
# Add the third ROC curve
plot(Tstage_ROC, print.auc = F, auc.polygon = FALSE, legacy.axes = TRUE,
     col = "green", max.auc.polygon = FALSE, print.thres = FALSE, add = TRUE)
# Add a legend for the ROC curves
legend("bottomright", title = "AUC:",
       legend = c(paste0("Signature: ", round(auc(signature_ROC), 2)), 
                  paste0("Gleason Score: ", round(auc(GS_ROC), 2)),
                  paste0("T Stage: ", round(auc(Tstage_ROC), 2))),
       col = c("blue", "purple", "green"), lwd = 2, cex = 0.8)
dev.off()


test_risk_gs <- roc.test(signature_ROC, GS_ROC, method = "delong")
test_risk_tstage <- roc.test(signature_ROC, Tstage_ROC, method = "delong")
test_risk_gs$p.value
test_risk_tstage$p.value


# 8.0 For Odds Ratio ----
library(circlize)
library(readxl)
library(readr)
library(caret)
library(pROC)
library(glmmTMB)

logistic_regression_for_unknown_prediction_without_random_effect <- function(train, test, target_var, random_effects_var, other_not_predictor_variables, identification_variable,
                                                                             num_iterations = 100, family = "binomial"){
  
  
  results <- data.frame()
  DMR_predictor <- data.frame(matrix(ncol = 1, nrow = ncol(train)))
  #DMR_predictor <- data.frame(matrix(ncol = 1, nrow = 29))
  prediction_result <- data.frame(patient = test[[identification_variable]], test_real = test[[target_var]])
  for (k in 1:num_iterations) {
    set.seed(k)
    # level1 <- subset(data, data[[target_var]] == 0)
    # level2 <- subset(data, data[[target_var]] == 1)
    # index_level1 <- sample(nrow(level1), 0.7*nrow(level1))
    # index_level2 <- sample(nrow(level2), 0.7*nrow(level2))
    # train <- data.frame(rbind(level1[index_level1, ], level2[index_level2, ]))
    # test <- data.frame(rbind(level1[-index_level1, ], level2[-index_level2, ]))
    
    train_x <- subset(train, select = names(train)[-which(names(train) %in% c(target_var, random_effects_var, other_not_predictor_variables))])
    train_y <- as.numeric(train[[target_var]])
    test_x <- subset(test, select = names(test)[-which(names(test) %in% c(target_var, random_effects_var, other_not_predictor_variables))])
    test_y <- as.numeric(test[[target_var]])
    
    # we are test the odds ratio of the original values of factor, their values should not be scaled
    # scaled_train_x <- apply(train_x, 2, scale)
    # scaled_test_x <- apply(test_x, 2, scale)
    scaled_train_x <- train_x
    scaled_test_x <- test_x
    
    train <- cbind(subset(train, select = c(target_var)), scaled_train_x)
    test <- cbind(subset(test, select = c(target_var)), scaled_test_x)
    
    predictors <- colnames(train)[-which(names(train) %in% c(target_var))]
    #formula_str <- paste0(target_var, "~", paste(predictors, collapse = " + "), "+ (1 | ", random_effects_var, ")")
    formula_str <- paste0(target_var, "~", paste(predictors, collapse = " + ")) # no random effects
    model_formula <- as.formula(formula_str)
    train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)# Set up cross-validation
    fit <- train(model_formula, data = train, method = "glm", family = family, trControl = train_control, metric = "ROC")# Train the model with cross-validation
    print(summary(fit))
    #summary(fit)


    # Store results
    # 从最终模型统一提取结果
    final_model <- fit$finalModel
    model_summary <- summary(final_model)
    
    # 统一数据获取方式
    coef_values <- coef(final_model)
    coef_names <- names(coef_values)
    p_values <- model_summary$coefficients[, "Pr(>|z|)"]
    p_values_full <- rep(NA, length(coef_names))  # 创建一个与 coef_names 长度相同的向量，填充 NA
    names(p_values_full) <- coef_names
    p_values_full[names(p_values)] <- p_values  # 将提取的 p 值填充到对应位置
    odds_ratios <- exp(coef_values)
    
    # 计算置信区间（使用轮廓似然方法）
    ci <- exp(confint(final_model))
    
    iteration_coef <- data.frame(coef_names, coef_values, p_values = p_values_full, odds_ratios, OR_lower = ci[, 1], OR_Upper = ci[, 2])
    DMR_predictor <- cbind(DMR_predictor[1:length(coef_values),], iteration_coef)
    
    # Prediction
    train_pre <- predict(fit, newdata = train, type = "prob")[,2]
    train_modelroc <- roc(train[[target_var]], as.numeric(train_pre))
    roc_train = train_modelroc$auc
    tab_train <- table(train_y, ifelse(train_pre < 0.5, 0, 1), dnn = c("true", "pre"))
    accuracy_train <- sum(diag(tab_train))/length(train_y)
    
    
    test_pre <- predict(fit, newdata = test, type='prob', re.form=NA)
    test_pre_01 <- ifelse(test_pre[,2] < 0.5, 0, 1)
    prediction_result <- cbind(prediction_result, test_pre_01)
    
    # ***** For prediction ability summary
    test_modelroc <- roc(test_y, as.numeric(test_pre_01)) # Calculate ROC for test data
    roc_test <- test_modelroc$auc
    tab_test <- table(test_y, test_pre_01, dnn = c("true", "pre"))
    accuracy_test <- sum(diag(tab_test)) / length(test_y)
    
    # pdf("Logistic_ROC.pdf", width = 5, height = 5)
    # plot(test_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
    #      grid.col=c("green", "red"), max.auc.polygon=TRUE,
    #      auc.polygon.col="skyblue", print.thres=TRUE) #绘制ROC曲线图
    # plot(ci(test_modelroc, of = "thresholds", thresholds = "best")) #添加最佳截断值位置
    # dev.off()
    
    # result <- data.frame(accuracy_train, accuracy_test, roc_train, roc_test)
    # results <- rbind(results, result)
    ## *****
    
    result <- data.frame(accuracy_train, accuracy_test, roc_train, roc_test)
    print(result)
    print(k)
  }
  return(list(prediction_result, DMR_predictor, test_modelroc, result)) #test_modelroc will be from the last iteration
}


# cor.test(train[,-"Bad_outcome"])
# 
# cor(train[,-"Bad_outcome"])

setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
DMRs <- fread("common_XGBoost_RF_predictor_target.csv")
All_sig_DMR_methy_group_tumor_infor_with_risk_score = fread("All_sig_DMR_methy_group_tumor_infor_with_risk_score.csv")
common_XGBoost_RF_predictor_target <- subset(All_sig_DMR_methy_group_tumor_infor_with_risk_score, select = c("patient", "surg_pgleasnt", "surg_pstagt", "recurrence", "Death", "Bad_outcome", "PSA_recurrence", "surg_lnmets",
                                                                                                             "risky", "risk_score", DMRs$Common_DMR_predictor_target))
str(common_XGBoost_RF_predictor_target)


# 8.1 Univariable logistic regression ----
table(common_XGBoost_RF_predictor_target$risky)

common_XGBoost_RF_predictor_target$surg_pgleasnt[common_XGBoost_RF_predictor_target$surg_pgleasnt > 7] = 1
common_XGBoost_RF_predictor_target$surg_pgleasnt[common_XGBoost_RF_predictor_target$surg_pgleasnt == 7] = 0

common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2a"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2b"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "2c"] = 0
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "3a"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "3b"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt[common_XGBoost_RF_predictor_target$surg_pstagt == "4"] = 1
common_XGBoost_RF_predictor_target$surg_pstagt = as.numeric(common_XGBoost_RF_predictor_target$surg_pstagt)

# common_XGBoost_RF_predictor_target$risk_score <- cut(common_XGBoost_RF_predictor_target$risk_score,
#                                                       breaks = seq(-60, 60, by = 20),
#                                                       labels = seq(0, length(seq(-60, 60, by = 20)) - 2),
#                                                       #labels = seq(0, length(seq(-60, 50, by = 20)) - 2),
#                                                       right = FALSE,  # Exclude the right endpoint (e.g., [-60, -50))
#                                                       include.lowest = TRUE  # Include the lowest point
#                                                      )
# common_XGBoost_RF_predictor_target$risk_score <- as.numeric(common_XGBoost_RF_predictor_target$risk_score)

# common_XGBoost_RF_predictor_target$risk_score[common_XGBoost_RF_predictor_target$risk_score > 0] = 1
# common_XGBoost_RF_predictor_target$risk_score[common_XGBoost_RF_predictor_target$risk_score < 0] = 0

# Select DMR variables by identifying columns that are not target or non-predictor variables
dmr_columns <- names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("Bad_outcome", "patient", "recurrence", "surg_lnmets", "Death", 
                                                                                                           "PSA_recurrence", "risky", "surg_pgleasnt", "surg_pstagt", "risk_score")]
common_XGBoost_RF_predictor_target[, (dmr_columns) := lapply(.SD, function(x) x * 100), .SDcols = dmr_columns]

common_XGBoost_RF_predictor_target$Bad_outcome = factor(common_XGBoost_RF_predictor_target$Bad_outcome, levels = c(0, 1), labels = c("Control", "Case"))

Odds_univariable <- data.frame()
for (DMR in c(DMRs$Common_DMR_predictor_target, "risk_score", "surg_pgleasnt", "surg_pstagt")) { #DMRs$Common_DMR_predictor_target, 
  train_on_High_Low_and_predic_Middle_by_single_var <- logistic_regression_for_unknown_prediction_without_random_effect(train = common_XGBoost_RF_predictor_target,
                                                                                                  test = common_XGBoost_RF_predictor_target,
                                                                                                  target_var = "Bad_outcome",
                                                                                                  random_effects_var = "patient", 
                                                                                                  other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c(DMR, "patient")], 
                                                                                                  identification_variable = "patient",
                                                                                                  num_iterations = 1, 
                                                                                                  family = "binomial")
  
  Odds_univariable <- rbind(Odds_univariable, train_on_High_Low_and_predic_Middle_by_single_var[[2]])
}
fwrite(Odds_univariable, "Odds_univariable.csv")



# 8.2 Multivariable logistic regression ----

# train_on_High_Low_and_predic_Middle_by_multi <- logistic_regression_for_unknown_prediction_without_random_effect(train = common_XGBoost_RF_predictor_target[, -c("chr9_1444688_1444802")],
#                                                                                                test = common_XGBoost_RF_predictor_target[, -c("chr9_1444688_1444802")],
#                                                                                                target_var = "Bad_outcome",
#                                                                                                random_effects_var = "patient", 
#                                                                                                other_not_predictor_variables = c("recurrence", "surg_lnmets", "Death", "PSA_recurrence", "risky"), 
#                                                                                                identification_variable = "patient",
#                                                                                                num_iterations = 1, family = "binomial")
# Odds_multi = train_on_High_Low_and_predic_Middle_by_multi[[2]]
# fwrite(Odds_multi, "Odds_multivariable.csv")

train_on_High_Low_and_predic_Middle_by_multi_DMR <- logistic_regression_for_unknown_prediction_without_random_effect(train = common_XGBoost_RF_predictor_target,
                                                                                               test = common_XGBoost_RF_predictor_target,
                                                                                               target_var = "Bad_outcome",
                                                                                               random_effects_var = "patient",
                                                                                               other_not_predictor_variables = c("surg_pgleasnt", "surg_pstagt", "recurrence", "surg_lnmets", "Death", "PSA_recurrence", "risky"),#, "risk_score"
                                                                                               identification_variable = "patient",
                                                                                               num_iterations = 1, family = "binomial")
Odds_multivariable_DMR = train_on_High_Low_and_predic_Middle_by_multi_DMR[[2]]

# summary(common_XGBoost_RF_predictor_target$risk_score)
train_on_High_Low_and_predic_Middle_by_multi_clinical <- logistic_regression_for_unknown_prediction_without_random_effect(train = common_XGBoost_RF_predictor_target,
                                                                                               test = common_XGBoost_RF_predictor_target,
                                                                                               target_var = "Bad_outcome",
                                                                                               random_effects_var = "patient",
                                                                                               other_not_predictor_variables = names(common_XGBoost_RF_predictor_target)[!names(common_XGBoost_RF_predictor_target) %in% c("surg_pstagt", "risk_score", "surg_pgleasnt", "patient")],#
                                                                                               identification_variable = "patient",
                                                                                               num_iterations = 1, family = "binomial")
Odds_multivariable_clinical = train_on_High_Low_and_predic_Middle_by_multi_clinical[[2]]
Odds_multivariable = rbind(Odds_multivariable_DMR, Odds_multivariable_clinical)
fwrite(Odds_multivariable, "Odds_multivariable.csv")



# 8.3 combined annotation ----
library(ggplot2)
library(scales)
library(patchwork)
library(grid)
library(data.table)
Odds_univariable = fread("Odds_univariable.csv")
Odds_univariable = Odds_univariable[!Odds_univariable$coef_names %in% c("surg_pstagt", "surg_pgleasnt"), ]
Odds_multi = fread("Odds_multivariable.csv")
Odds_multi = Odds_multi[!Odds_multi$coef_names %in% c("surg_pstagt", "surg_pgleasnt"), ]


Odds_univariable = Odds_univariable[order(Odds_univariable$odds_ratios, decreasing = T), ]
Odds_multi = Odds_multi[order(Odds_multi$odds_ratios, decreasing = T), ]

Odds_univariable$Sig=NA
Odds_univariable$Sig[Odds_univariable$p_values < 0.05]="Significant"
Odds_univariable$Sig[Odds_univariable$p_values >= 0.05]="Non-Significant"
Odds_multi$Sig=NA
Odds_multi$Sig[Odds_multi$p_values < 0.05]="Significant"
Odds_multi$Sig[Odds_multi$p_values >= 0.05]="Non-Significant"

Odds_univariable = Odds_univariable[Odds_univariable$coef_names != "(Intercept)", ]
Odds_multi = Odds_multi[Odds_multi$coef_names != "(Intercept)", ]



# custom_log_trans <- trans_new(
#   name = "custom_log",
#   transform = function(x) log(x + 1, base = 1000), # Add 1 to avoid log of zero issues
#   inverse = function(x) 100^x - 1, # Subtract 1 to revert the transformation correctly
#   breaks = function(x) c(0, 1, 2, 5,10, 20, 100), # Custom breaks for desired spacing
#   domain = c(0, Inf)
# )

# Define plot for univariable
plot_univariable <- ggplot(Odds_univariable, aes(x = odds_ratios, y = coef_names, color = Sig)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = `OR_Upper`, xmin = `OR_lower`), size = 1, height = 0.2) +
  geom_point(size = 2, shape = 16) +
  scale_color_manual(name = "",
                     values = c("Significant" = "red3", 
                                "Non-Significant" = "mediumblue")) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(colour = "black", size = 14),
    axis.text.y = element_text(size = 14, face = "plain"),
    axis.title.x = element_text(size = 20, face = "plain"),
    legend.title = element_text(size = 14, face = "plain"),
    legend.text = element_text(size = 14, face = "plain"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    axis.ticks.length.x = unit(0.3, "cm"),  # 加长刻度线以匹配截断标记
    axis.line = element_line(colour = "black")
  ) + 
  scale_x_continuous(
    breaks = c(0, 1, 2),       # 只显示这几个刻度
    limits = c(0, 2)           # 如果想把坐标范围固定在 0~2，设置 limits
  ) +
  # scale_x_continuous(
  #   trans = custom_log_trans,
  #   breaks = c(0, 1, 2, 5, 10),
  #   labels = c("0", "1","2" , "5", "10"),
  #   limits = c(0, 10)
  # ) +
  labs(x = "Odds Ratio (Univariable)", y = "")

# Define plot for multivariable
plot_multivariable <- ggplot(Odds_multi, aes(x = odds_ratios, y = coef_names, color = Sig)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = `OR_Upper`, xmin = `OR_lower`), size = 1, height = 0.2) +
  geom_point(size = 2, shape = 16) +
  scale_color_manual(name = "",
                     values = c("Significant" = "red3", 
                                "Non-Significant" = "mediumblue")) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(colour = "black", size = 14),
    axis.text.y = element_blank(),  # Remove Y-axis text on this plot for alignment
    axis.title.x = element_text(size = 20, face = "plain"),
    legend.title = element_text(size = 14, face = "plain"),
    legend.text = element_text(size = 14, face = "plain"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  ) + 
  scale_x_continuous(
    breaks = c(0, 1, 2),       # 只显示这几个刻度
    limits = c(0, 2)           # 如果想把坐标范围固定在 0~2，设置 limits
  )+
  # scale_x_continuous(
  #   trans = custom_log_trans,
  #   breaks = c(0, 1, 5, 20),
  #   labels = c("0", "1", "5", "20"),
  #   limits = c(0, 20)
  # ) +
  labs(x = "Odds Ratio (Multivariable)", y = "")

# Combine plots with a shared Y-axis
combined_plot <- plot_univariable + plot_multivariable + plot_layout(ncol = 2)
print(combined_plot)

ggsave("Odds_Ratio_summary.pdf", width = 16, height = 8)


for (DMR in DMRs$Common_DMR_predictor_target){
  print(DMR)
  print(summary(common_XGBoost_RF_predictor_target[[DMR]]))
}



# 9.0 DMR signature information summary ----
setwd("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_selection")
All_DMRs <- data.frame(fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/DMR_summary_by_different_strategy/All_merged_DMRs_sig_methy_matrix_noNA.csv"))
All_DMR_gene_annotation = fread("/Users/lmh/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/WashU_Graduate_Study/Maher_Lab/ProstateCancer/output/All_DM_annotation_for_120_dangerous/Overall_methylation_summary/All_DMR_gene_annotation.csv")
DMRs <- fread("common_XGBoost_RF_predictor_target.csv")
Odds_univariable = fread("Odds_univariable.csv")
Odds_multi = fread("Odds_multivariable.csv")

signature_basic_information <- All_DMRs[All_DMRs$DMR %in% DMRs$Common_DMR_predictor_target, ]
signature_basic_information <- signature_basic_information[, c("DMR", "chr", "start", "end", "width", 
                                                               "High_Avg", "Low_Avg", "Differential")]
signature_gene_information <- All_DMR_gene_annotation[All_DMR_gene_annotation$DMR %in% DMRs$Common_DMR_predictor_target, ]
signature_gene_information <- signature_gene_information[signature_gene_information$type == "gene", c("DMR", "type", "gene_id", "gene_name", "gene_biotype")]

signature_odds_information <- left_join(Odds_univariable[Odds_univariable$coef_names != "(Intercept)", c("coef_names", "p_values", "odds_ratios", "OR_lower", "OR_Upper")],
                                        Odds_multi[Odds_multi$coef_names != "(Intercept)", c("coef_names", "p_values", "odds_ratios", "OR_lower", "OR_Upper")], by = "coef_names" )
names(signature_odds_information)[1] <- "DMR"


signature_basic_and_biological_information <- left_join(signature_basic_information, signature_gene_information)
signature_basic_and_biological_and_OR_information <- left_join(signature_basic_and_biological_information, signature_odds_information)
fwrite(signature_basic_and_biological_information, "signature_basic_and_biological_information.csv")
fwrite(signature_basic_and_biological_and_OR_information, "signature_basic_and_biological_and_OR_information.csv")




