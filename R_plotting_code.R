# Plotting script for population genomics analyses
# ----------------------------------------------------
# This script contains ggplot2/ComplexHeatmap-based visualization code.
# The original analysis logic and plotting parameters are retained.


# Plot violin distribution of Fis across continents
plot_fis_by_continent <- ggplot(
  mat_type_239[which(mat_type_239[, 2] != "South America"), ],
  aes(y = Fis, x = continent, fill = continent)
) +
  geom_violin(trim = FALSE, bw = 0.15, width = 0.9, color = "white") +
  # geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = color_value_country) +
  theme(legend.position = "none")


# Plot relationship between missing-site counts and Fis across continents
plot_fis_missing_site <- ggplot(
  mat_type_239[which(mat_type_239[, 2] != "South America"), ],
  aes(y = Fis, x = V22, colour = nuclear)
) +
  geom_point() +
  # geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = c("#EE7671", "#51BC5F", "#719DD1")) +
  xlim(0, 60000) +
  ylim(-1.2, 0.5) +
  scale_colour_manual(values = c("#EE7671", "#51BC5F", "#719DD1")) +
  xlab("Miss-site")


# Plot PCA of all samples colored by continent/nuclear type
plot_pca_by_nuclear <- ggplot(data = df1, aes(x = PC1, y = PC2)) +
  labs(x = xlab, y = ylab, color = "") +
  geom_point(size = 2, alpha = 0.5, aes(colour = mat_type_239[, 16])) +
  theme(axis.text = element_text(size = 20)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values = c("#EE7671", "#51BC5F", "#719DD1")) +
  theme(legend.position = "none") +
  geom_vline(aes(xintercept = 0), colour = "grey60", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "grey60", linetype = "dashed")
# stat_ellipse(lwd = 1, level = 0.8, aes(colour = mat_type_239[, 2]))
# scale_colour_manual(values = c("Inbred pop" = "#FAE620", "Natural pop" = "#36B779", "RNA-seq" = "#32678E"))


# Plot PCA of all samples colored by Fis
plot_pca_by_fis <- ggplot(data = df1, aes(x = PC1, y = PC2)) +
  labs(x = xlab, y = ylab, color = "") +
  geom_point(size = 2, alpha = 0.5, aes(colour = mat_type_239[, 6])) +
  # , shape = mat_type_239[, 2])) +
  theme(axis.text = element_text(size = 20)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_vline(aes(xintercept = 0), colour = "grey60", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "grey60", linetype = "dashed")
# theme(legend.position = "none")


# Plot density distribution of Fis
plot_fis_density <- ggplot() +
  geom_density(data = mat_type_239, aes(x = Fis), fill = "steelblue", color = "steelblue", bw = 0.05) +
  theme_bw() +
  scale_colour_manual(values = "steelblue")


# Plot MAF distribution from VCF
plot_maf_density <- ggplot() +
  geom_density(data = snp_maf, aes(x = V4), fill = "steelblue", color = "steelblue", bw = 0.01) +
  theme_bw() +
  scale_colour_manual(values = "steelblue") +
  xlab("MAF")


# Plot missing-rate distribution of SNPs in VCF
plot_missing_rate_density <- ggplot() +
  geom_density(data = miss_info, aes(x = F_MISS), fill = "steelblue", color = "steelblue", bw = 0.01) +
  theme_bw() +
  scale_colour_manual(values = "steelblue") +
  xlab("MISS")


# Plot IBS among nuclear-type groups
plot_ibs_by_nuclear_pair <- ggplot(
  pst239_ibs,
  aes(y = DST, x = reorder(V17, -DST), fill = reorder(V17, -DST))
) +
  geom_violin(trim = FALSE, width = 0.9, bw = 0.02, aes(color = reorder(V17, -DST))) +
  geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = color_value_country) +
  scale_color_manual(values = color_value_country) +
  theme(legend.position = "none") +
  xlab("nuclear-nuclear") +
  ylab("IBS")


# Plot SNP Fst among nuclear-type groups
plot_fst_by_nuclear_pair <- ggplot(
  all_fst,
  aes(y = WEIR_AND_COCKERHAM_FST, x = reorder(V4, -WEIR_AND_COCKERHAM_FST), fill = reorder(V4, -WEIR_AND_COCKERHAM_FST))
) +
  geom_violin(trim = FALSE, color = "white", bw = 0.15) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = c("#F37E7C", "#FCBF07", "#F7EC18")) +
  scale_color_manual(values = color_value_country) +
  theme(legend.position = "none") +
  xlab("nuclear-nuclear") +
  ylab("Fst")


# Plot mash k-mer density distributions of four nuclei
plot_mash_kmer_density <- ggplot(all_mash3) +
  geom_density_ridges(
    aes(x = shared_kmer, y = nuclear, fill = nuclear, color = nuclear),
    scale = 1.2,
    bandwidth = 15000
  ) +
  scale_fill_manual(values = c("#6C88BE", "#98A7C9", "#B8BED2", "#D4E3A3")) +
  scale_color_manual(values = c("#6C88BE", "#98A7C9", "#B8BED2", "#D4E3A3")) +
  xlim(300000, 530000) +
  theme_bw() +
  theme(legend.position = "none")


genome_snp_inpop[1, 1] <- length(which(pst239_data[, 242] %in% pri_all_snp[which(pri_all_snp[, 2] != "."), 7]))
genome_snp_inpop[2, 1] <- length(which(pst239_data[, 242] %in% pri_all_snp[which(pri_all_snp[, 3] != "."), 7]))
genome_snp_inpop[3, 1] <- length(which(pst239_data[, 242] %in% pri_all_snp[which(pri_all_snp[, 4] != "."), 7]))
genome_snp_inpop[4, 1] <- length(which(pst239_data[, 242] %in% pri_all_snp[which(pri_all_snp[, 2] != pri_all_snp[, 3]), 7]))
genome_snp_inpop[5, 1] <- length(which(pst239_data[, 242] %in% pri_all_snp[which(pri_all_snp[, 2] != pri_all_snp[, 4]), 7]))
genome_snp_inpop[6, 1] <- length(which(pst239_data[, 242] %in% pri_all_snp[which(pri_all_snp[, 3] != pri_all_snp[, 4]), 7]))


# Plot distribution counts of SNPs identified from nuclear genome comparisons in the population
plot_genome_snp_inpop <- ggplot(
  genome_snp_inpop,
  aes(y = snp_inpop, x = reorder(V2, -snp_inpop), fill = reorder(V2, -snp_inpop))
) +
  geom_col(width = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = color_value_country) +
  scale_color_manual(values = color_value_country) +
  theme(legend.position = "none") +
  xlab("nuclear-nuclear") +
  ylab("snp_inpop")


temporary_assignment <- az2a


# Plot mash heatmap of all samples
pst239_mash[, 5:6] <- read.table(
  "/dt3/project/006.pst_revision_ybwang_20250927/004.info/all_mash_pst104E.tab",
  header = FALSE
)

pst239_ha <- HeatmapAnnotation(
  continent = pst239_metadata[which(!is.na(pst239_metadata[, 5])), 2],
  Fis = pst239_metadata[which(!is.na(pst239_metadata[, 5])), 5],
  col = list(
    continent = c(
      "Africa" = "#7D5EAA",
      "Asia" = "#7ED3F6",
      "Europe" = "#9BCC39",
      "North America" = "#F7EC18",
      "Oceania" = "#FCBF07",
      "South America" = "#F37E7C"
    ),
    Fis = colorRamp2(
      c(
        max(pst239_metadata[which(!is.na(pst239_metadata[, 5])), 5]),
        min(pst239_metadata[which(!is.na(pst239_metadata[, 5])), 5])
      ),
      c("#002B7A", "white")
    )
  )
)

Heatmap(
  as.matrix(t(pst239_mash[, 1:6])),
  show_column_names = FALSE,
  top_annotation = pst239_ha,
  # left_annotation = pst239_ha2,
  column_split = as.matrix(pst239_metadata[which(!is.na(pst239_metadata[, 5])), 8, drop = FALSE]),
  col = colorRamp2(breaks = seq(380001, 500000, 20000), colors = viridis(6)),
  column_gap = unit(30, "mm")
)


# Plot density distribution of Fis_cor across site types
plot_fis_cor_density_by_type <- ggplot(data = fis_cor, aes(x = V1^2)) +
  geom_density(color = "white") +
  theme_bw() +
  geom_area(data = subset(dense, x >= 0 & x < 0.25), aes(x, y, fill = "1")) +
  geom_area(data = subset(dense, x >= 0.25 & x < 0.75), aes(x, y, fill = "2")) +
  geom_area(data = subset(dense, x >= 0.75 & x < 1), aes(x, y, fill = "3")) +
  scale_fill_manual(
    values = c("#D5E4A8", "#D0E7ED", "#AFB6D2"),
    breaks = c("1", "2", "3"),
    labels = c("type1", "type2", "type3")
  ) +
  xlab("Fis_cor") +
  ylab("density")


# Plot admixture results
plot_admixture <- plotQ(
  alist[c(1:9)],
  exportplot = FALSE,
  returnplot = TRUE,
  barsize = 1,
  returndata = TRUE,
  ordergrp = TRUE,
  sortind = "all",
  grplabsize = 3.5,
  pointsize = 6,
  linesize = 7,
  linealpha = 0.2,
  pointcol = "white",
  grplabpos = 0.5,
  linepos = 0.5,
  grplabheight = 0.75,
  grplab = onegrpset[, 1, drop = FALSE],
  imgoutput = "join",
  sharedindlab = FALSE,
  clustercol = brewer.pal(10, "Set3")
)


plot_variant_per_bp_by_area <- ggplot(
  all_bed[which(all_bed[, 10] != 4), ],
  aes(y = V14, x = area, fill = area)
) +
  geom_violin(trim = FALSE, bw = 0.0008, width = 0.9, aes(color = area)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.9), outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  scale_color_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  geom_signif(
    comparisons = list(c("area_1", "area_3"), c("area_1", "area_2"), c("area_2", "area_3")),
    y_position = c(0.058, 0.045, 0.051),
    map_signif_level = FALSE,
    tip_length = c(0.005, 0.005, 0.005, 0.005, 0.005, 0.005)
  ) +
  theme(legend.position = "none") +
  xlab("variant per bp") +
  ylab("area")


plot_gene_per_bp_by_area <- ggplot(
  all_bed[which(all_bed[, 10] != 4), ],
  aes(y = V16, x = area, fill = area)
) +
  # geom_violin(trim = FALSE, bw = 0.000001, width = 0.9, aes(color = V15)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.9), outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  # scale_y_continuous(labels = scientific_format()) +
  scale_fill_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  scale_color_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  geom_signif(
    comparisons = list(c("area_1", "area_3"), c("area_1", "area_2"), c("area_2", "area_3")),
    y_position = c(0.00023, 0.00017, 0.00020),
    map_signif_level = FALSE,
    tip_length = c(0.005, 0.005, 0.005, 0.005, 0.005, 0.005)
  ) +
  theme(legend.position = "none") +
  xlab("gene per bp") +
  ylab("area")


plot_het_snp_counts_by_nuclear <- ggplot(snp_distr, aes(y = V1, x = V2, fill = V3)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = c("#AFB6D2", "#D0E7ED", "#D5E4A8")) +
  scale_y_continuous(labels = scientific_format()) +
  xlab("nuclear") +
  ylab("het snp counts") +
  theme(legend.position = "none")


plot_effector_distribution_by_area <- ggplot(
  effector_distr[which(effector_distr[, 2] != "area_4"), ],
  aes(y = V4, x = area, fill = gene_type)
) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("distribution of effectors")


plot_cross_pos_chr1_counts <- ggplot(
  data = cross_pos[which(cross_pos[, 3] == "chr1"), ],
  aes(x = POS / 1000, y = V1)
) +
  geom_col(width = 10) +
  theme_bw() +
  ylab("counts")


plot_gene_count_by_area <- ggplot(
  all_bed[which(all_bed[, 10] != 4), ],
  aes(y = V21, x = area, fill = area)
) +
  geom_violin(trim = FALSE, bw = 3, width = 0.9, aes(color = area)) +
  # geom_boxplot(width = 0.3, position = position_dodge(0.9), outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  scale_color_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  geom_signif(
    comparisons = list(c("area_1", "area_3"), c("area_1", "area_2"), c("area_2", "area_3")),
    y_position = c(23, 17, 20),
    map_signif_level = TRUE,
    tip_length = c(0.005, 0.005, 0.005, 0.005, 0.005, 0.005)
  ) +
  theme(legend.position = "none") +
  xlab("gene per bp") +
  ylab("area") +
  ylim(0, 30)


plot_recombination_breakpoints_by_area <- ggplot(effector_distr[1:3, ], aes(y = V6, x = area, fill = area)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  xlab("area") +
  ylab("recombination breakpoints per 5kb") +
  theme(legend.position = "none")


plot_tpm_by_area <- ggplot(
  data = gene_exp_250_avg_tpm[which(gene_exp_250_avg_tpm[, 2] != "area_4"), ],
  aes(x = area, y = log2(TPM + 1), fill = area)
) +
  geom_violin(aes(color = area)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.9), outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  scale_color_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  geom_signif(
    comparisons = list(c("area_1", "area_3"), c("area_1", "area_2"), c("area_2", "area_3")),
    y_position = c(16, 12, 14),
    map_signif_level = FALSE,
    tip_length = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
  ) +
  theme(legend.position = "none") +
  xlab("area") +
  ylab("log2(TPM+1)") +
  ylim(0, 18)


plot_fpkm_by_area <- ggplot(
  data = gene_exp_250_avg_fpkm[which(gene_exp_250_avg_fpkm[, 2] != "area_4"), ],
  aes(x = area, y = log2(FPKM + 1), fill = area)
) +
  geom_violin(aes(color = area)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.9), outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  scale_color_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  geom_signif(
    comparisons = list(c("area_1", "area_3"), c("area_1", "area_2"), c("area_2", "area_3")),
    y_position = c(19, 14, 16),
    map_signif_level = FALSE,
    tip_length = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
  ) +
  theme(legend.position = "none") +
  xlab("area") +
  ylab("log2(FPKM+1)") +
  ylim(0, 21)


plot_fis_cor_by_area <- ggplot(data = fis_cor, aes(x = V12, y = V10, fill = V12)) +
  geom_violin(aes(color = V12)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  scale_color_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  geom_signif(
    comparisons = list(c("area_1", "area_3"), c("area_1", "area_2"), c("area_2", "area_3")),
    y_position = c(0.8, 0.6, 0.7),
    map_signif_level = FALSE,
    tip_length = c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02)
  ) +
  theme(legend.position = "none") +
  xlab("area") +
  ylab("fis_cor") +
  ylim(0, 0.9)


plot_missing_rate_by_area <- ggplot() +
  geom_density(data = fis_cor, aes(x = V11, color = V12), bw = 0.05, alpha = 0.4, linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = c("#D5E4A8", "#D0E7ED", "#AFB6D2")) +
  xlab("F_MISS")
