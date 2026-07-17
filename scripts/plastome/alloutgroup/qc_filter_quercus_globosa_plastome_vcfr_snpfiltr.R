#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vcfR)
  library(ggplot2)
})

vcf_file <- "Quercus_globosa_plastome.allsite.allsite.exclude8_2397024510.filtered.variant.vcf.gz"
out_dir <- "vcfr_snpfiltr_qc"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(vcf_file))

snpfiltr_available <- requireNamespace("SNPfiltR", quietly = TRUE)
if (snpfiltr_available) {
  suppressPackageStartupMessages(library(SNPfiltR))
}

timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

read_vcf_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) break
    header <- c(header, line)
    if (startsWith(line, "#CHROM")) break
  }
  header
}

parse_info_value <- function(info, key) {
  out <- rep(NA_character_, length(info))
  pattern <- paste0("(^|;)", key, "=")
  has_key <- grepl(pattern, info)
  if (any(has_key)) {
    out[has_key] <- sub(
      paste0(".*(^|;)", key, "=([^;]+).*"),
      "\\2",
      info[has_key],
      perl = TRUE
    )
  }
  out
}

as_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

is_missing_gt <- function(gt) {
  is.na(gt) | gt == "." | grepl("\\.", gt)
}

is_heterozygous_gt <- function(gt) {
  miss <- is_missing_gt(gt)
  stripped <- sub(":.*$", "", gt)
  alleles <- strsplit(stripped, "[/|]")
  het <- vapply(alleles, function(x) {
    x <- x[nzchar(x) & x != "."]
    length(unique(x)) > 1
  }, logical(1))
  het[miss] <- FALSE
  het
}

allele_counts <- function(gt) {
  stripped <- sub(":.*$", "", gt)
  alleles <- strsplit(stripped, "[/|]")
  ref <- alt <- called <- integer(length(alleles))
  for (i in seq_along(alleles)) {
    x <- alleles[[i]]
    x <- x[nzchar(x) & x != "."]
    called[i] <- length(x)
    ref[i] <- sum(x == "0")
    alt[i] <- sum(x != "0")
  }
  list(ref = ref, alt = alt, called = called)
}

variant_class <- function(ref, alt) {
  alt1 <- strsplit(alt, ",", fixed = TRUE)
  is_snp <- mapply(function(r, a) {
    nchar(r) == 1 && all(nchar(a) == 1)
  }, ref, alt1)
  ifelse(is_snp, "SNP", "indel")
}

write_csv <- function(x, filename) {
  write.csv(x, file.path(out_dir, filename), row.names = FALSE, quote = TRUE)
}

vcf_header <- read_vcf_header(vcf_file)
bcftools_commands <- grep("^##bcftools_.*Command=", vcf_header, value = TRUE)

vcf <- read.vcfR(vcf_file, verbose = FALSE)
fix <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
sample_names <- colnames(vcf@gt)[-1]

gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
dp <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
gq <- extract.gt(vcf, element = "GQ", as.numeric = TRUE)

missing_mat <- matrix(is_missing_gt(as.vector(gt)), nrow = nrow(gt), dimnames = dimnames(gt))
het_mat <- matrix(is_heterozygous_gt(as.vector(gt)), nrow = nrow(gt), dimnames = dimnames(gt))
counts <- allele_counts(as.vector(gt))
ref_mat <- matrix(counts$ref, nrow = nrow(gt), dimnames = dimnames(gt))
alt_mat <- matrix(counts$alt, nrow = nrow(gt), dimnames = dimnames(gt))
called_mat <- matrix(counts$called, nrow = nrow(gt), dimnames = dimnames(gt))

site_called_alleles <- rowSums(called_mat, na.rm = TRUE)
site_alt_alleles <- rowSums(alt_mat, na.rm = TRUE)
site_af <- ifelse(site_called_alleles > 0, site_alt_alleles / site_called_alleles, NA_real_)
site_maf <- pmin(site_af, 1 - site_af, na.rm = FALSE)

site_qc <- data.frame(
  chrom = fix$CHROM,
  pos = as.integer(fix$POS),
  ref = fix$REF,
  alt = fix$ALT,
  variant_class = variant_class(fix$REF, fix$ALT),
  qual = as_num(fix$QUAL),
  filter = fix$FILTER,
  missing_rate = rowMeans(missing_mat),
  heterozygote_rate = rowMeans(het_mat),
  allele_frequency = site_af,
  minor_allele_frequency = site_maf,
  mean_dp = rowMeans(dp, na.rm = TRUE),
  median_dp = apply(dp, 1, median, na.rm = TRUE),
  mean_gq = rowMeans(gq, na.rm = TRUE),
  qd = as_num(parse_info_value(fix$INFO, "QD")),
  fs = as_num(parse_info_value(fix$INFO, "FS")),
  mq = as_num(parse_info_value(fix$INFO, "MQ")),
  mq_rank_sum = as_num(parse_info_value(fix$INFO, "MQRankSum")),
  read_pos_rank_sum = as_num(parse_info_value(fix$INFO, "ReadPosRankSum")),
  sor = as_num(parse_info_value(fix$INFO, "SOR")),
  info_dp = as_num(parse_info_value(fix$INFO, "DP")),
  stringsAsFactors = FALSE
)

sample_qc <- data.frame(
  sample = sample_names,
  called_genotypes = colSums(!missing_mat),
  missing_genotypes = colSums(missing_mat),
  missing_rate = colMeans(missing_mat),
  heterozygote_genotypes = colSums(het_mat),
  heterozygote_rate = colMeans(het_mat),
  mean_dp = colMeans(dp, na.rm = TRUE),
  median_dp = apply(dp, 2, median, na.rm = TRUE),
  mean_gq = colMeans(gq, na.rm = TRUE),
  stringsAsFactors = FALSE
)

summary_qc <- data.frame(
  metric = c(
    "VCF",
    "Generated",
    "SNPfiltR installed",
    "Samples",
    "Variant records",
    "SNP records",
    "Indel records",
    "PASS records",
    "Mean site missing rate",
    "Max site missing rate",
    "Mean sample missing rate",
    "Max sample missing rate",
    "Mean site heterozygote rate",
    "Max site heterozygote rate",
    "Mean DP across called genotypes",
    "Median DP across called genotypes",
    "Mean GQ across called genotypes",
    "Ts/Tv for SNP records"
  ),
  value = c(
    vcf_file,
    timestamp,
    as.character(snpfiltr_available),
    length(sample_names),
    nrow(site_qc),
    sum(site_qc$variant_class == "SNP"),
    sum(site_qc$variant_class == "indel"),
    sum(site_qc$filter == "PASS"),
    signif(mean(site_qc$missing_rate), 4),
    signif(max(site_qc$missing_rate), 4),
    signif(mean(sample_qc$missing_rate), 4),
    signif(max(sample_qc$missing_rate), 4),
    signif(mean(site_qc$heterozygote_rate), 4),
    signif(max(site_qc$heterozygote_rate), 4),
    signif(mean(as.vector(dp), na.rm = TRUE), 4),
    signif(median(as.vector(dp), na.rm = TRUE), 4),
    signif(mean(as.vector(gq), na.rm = TRUE), 4),
    {
      snps <- site_qc$variant_class == "SNP"
      substitutions <- paste0(site_qc$ref[snps], ">", site_qc$alt[snps])
      transitions <- substitutions %in% c("A>G", "G>A", "C>T", "T>C")
      transversions <- substitutions %in% c(
        "A>C", "C>A", "A>T", "T>A", "C>G", "G>C",
        "G>T", "T>G"
      )
      signif(sum(transitions) / sum(transversions), 4)
    }
  ),
  stringsAsFactors = FALSE
)

filtering_conditions <- data.frame(
  step = seq_len(8),
  condition = c(
    "Exclude eight samples",
    "Keep reference sites or biallelic SNP/indel sites",
    "Remove sites with >=10% missing genotypes before genotype masking",
    "Mask low-quality variant-site genotypes to missing",
    "Mask low-quality reference-site genotypes to missing",
    "Remove sites with >=10% missing genotypes after genotype masking",
    "Apply GATK-style site hard filters separately to SNPs and indels",
    "Keep variant SNP/indel sites for this final variant VCF"
  ),
  expression = c(
    "Quercus_414, Quercus_415, Quercus_419, Quercus_422, Quercus_424, Quercus_427, Quercus_433, Quercus_435",
    "N_ALT=0 || (N_ALT=1 && (TYPE=\"snp\" || TYPE=\"indel\"))",
    "F_MISSING < 0.10",
    "N_ALT>0 & (FMT/DP<10 | FMT/GQ<30 | FMT/DP>10000)",
    "N_ALT=0 & (FMT/DP<10 | FMT/RGQ<30 | FMT/DP>10000)",
    "F_MISSING < 0.10",
    "SNPs: QD>=2, FS<=60, MQ>=40, MQRankSum>=-12.5, ReadPosRankSum>=-8, SOR<=3; indels: QD>=2, FS<=200, ReadPosRankSum>=-20, SOR<=10",
    "N_ALT>0 && (TYPE=\"snp\" || TYPE=\"indel\")"
  ),
  source = c(
    "VCF header bcftools_viewCommand and file name exclude8",
    "VCF header bcftools_viewCommand",
    "VCF header bcftools_viewCommand",
    "VCF header bcftools_pluginCommand setGT",
    "VCF header bcftools_pluginCommand setGT",
    "VCF header bcftools_viewCommand",
    "VCF header bcftools_filterCommand",
    "VCF header bcftools_viewCommand"
  ),
  stringsAsFactors = FALSE
)

write_csv(summary_qc, "qc_summary.csv")
write_csv(site_qc, "site_qc_vcfr.csv")
write_csv(sample_qc, "sample_qc_vcfr.csv")
write_csv(filtering_conditions, "filtering_conditions.csv")
writeLines(bcftools_commands, file.path(out_dir, "bcftools_filtering_commands_from_vcf_header.txt"))

plot_variant_counts <- ggplot(site_qc, aes(x = variant_class, fill = variant_class)) +
  geom_bar(width = 0.7) +
  scale_fill_manual(values = c(SNP = "#2b7a78", indel = "#bc6c25")) +
  labs(x = NULL, y = "Records", title = "Variant Classes After Filtering") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

plot_site_missing <- ggplot(site_qc, aes(x = missing_rate)) +
  geom_histogram(binwidth = 0.01, fill = "#4c6a92", color = "white") +
  labs(x = "Site missing genotype rate", y = "Sites", title = "Site Missingness") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

plot_sample_qc <- ggplot(sample_qc, aes(x = missing_rate, y = heterozygote_rate)) +
  geom_point(color = "#6b8e23", alpha = 0.85, size = 2) +
  labs(
    x = "Sample missing genotype rate",
    y = "Sample heterozygote rate",
    title = "Per-Sample Missingness And Heterozygosity"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

plot_quality_depth <- ggplot(site_qc, aes(x = mean_dp, y = qual, color = variant_class)) +
  geom_point(alpha = 0.85, size = 2) +
  scale_color_manual(values = c(SNP = "#2b7a78", indel = "#bc6c25")) +
  scale_y_log10() +
  labs(
    x = "Mean genotype DP by site",
    y = "QUAL, log10 scale",
    color = "Variant class",
    title = "Site Quality And Depth"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

plot_site_het_position <- ggplot(site_qc, aes(x = pos, y = heterozygote_rate, color = variant_class)) +
  geom_point(alpha = 0.82, size = 2) +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "#6b7280") +
  scale_color_manual(values = c(SNP = "#2b7a78", indel = "#bc6c25")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Plastome position",
    y = "Heterozygote rate",
    color = "Variant class",
    title = "Heterozygote Rate Across The Plastome"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

plot_maf <- ggplot(site_qc, aes(x = minor_allele_frequency, fill = variant_class)) +
  geom_histogram(binwidth = 0.025, color = "white", position = "stack") +
  scale_fill_manual(values = c(SNP = "#2b7a78", indel = "#bc6c25")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Minor allele frequency",
    y = "Records",
    fill = "Variant class",
    title = "Minor Allele Frequency After Filtering"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "variant_class_counts.png"), plot_variant_counts, width = 6, height = 4, dpi = 300)
ggsave(file.path(out_dir, "site_missingness_histogram.png"), plot_site_missing, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(out_dir, "sample_missingness_heterozygosity.png"), plot_sample_qc, width = 6, height = 4.5, dpi = 300)
ggsave(file.path(out_dir, "site_quality_depth.png"), plot_quality_depth, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(out_dir, "site_heterozygosity_by_position.png"), plot_site_het_position, width = 8, height = 4.5, dpi = 300)
ggsave(file.path(out_dir, "minor_allele_frequency_histogram.png"), plot_maf, width = 7, height = 4.5, dpi = 300)

markdown_table <- function(x) {
  x[] <- lapply(x, as.character)
  x[] <- lapply(x, function(col) gsub("\\|", "\\\\|", col))
  widths <- vapply(seq_along(x), function(i) max(nchar(c(names(x)[i], x[[i]]))), integer(1))
  header <- paste0("| ", paste(sprintf(paste0("%-", widths, "s"), names(x)), collapse = " | "), " |")
  sep <- paste0("| ", paste(strrep("-", widths), collapse = " | "), " |")
  rows <- apply(x, 1, function(row) {
    paste0("| ", paste(sprintf(paste0("%-", widths, "s"), row), collapse = " | "), " |")
  })
  c(header, sep, rows)
}

report <- c(
  "# Quercus globosa Plastome VCF QC And Filtering Conditions",
  "",
  paste("Generated:", timestamp),
  "",
  "## Input",
  "",
  paste("- VCF:", vcf_file),
  paste("- Samples:", length(sample_names)),
  paste("- Variant records:", nrow(site_qc)),
  paste("- SNPfiltR installed:", snpfiltr_available),
  "",
  if (!snpfiltr_available) {
    paste(
      "SNPfiltR is not installed in this R library, so this run uses vcfR-based",
      "equivalent QC summaries and records the SNPfiltR status for reproducibility."
    )
  } else {
    "SNPfiltR loaded successfully; vcfR objects and SNPfiltR are available for downstream filtering."
  },
  "",
  "## Filtering Conditions",
  "",
  markdown_table(filtering_conditions[, c("step", "condition", "expression")]),
  "",
  "## QC Summary",
  "",
  markdown_table(summary_qc),
  "",
  "## Output Files",
  "",
  "- qc_summary.csv",
  "- filtering_conditions.csv",
  "- site_qc_vcfr.csv",
  "- sample_qc_vcfr.csv",
  "- bcftools_filtering_commands_from_vcf_header.txt",
  "- variant_class_counts.png",
  "- site_missingness_histogram.png",
  "- sample_missingness_heterozygosity.png",
  "- site_quality_depth.png",
  "",
  "## Notes",
  "",
  "- The filtering conditions are read from the final VCF header and the existing alloutgroup filtering history.",
  "- This report describes the already-filtered variant VCF; it does not rewrite or further filter the VCF."
)

writeLines(report, file.path(out_dir, "Quercus_globosa_plastome_vcfr_snpfiltr_qc_report.md"))

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

html_table <- function(x, digits = 4) {
  if (!nrow(x)) return("<p>No records met this condition.</p>")
  y <- x
  y[] <- lapply(y, function(col) {
    if (is.numeric(col)) {
      ifelse(is.na(col), "", format(round(col, digits), trim = TRUE, scientific = FALSE))
    } else {
      html_escape(col)
    }
  })
  header <- paste0("<th>", html_escape(names(y)), "</th>", collapse = "")
  rows <- apply(y, 1, function(row) {
    paste0("<tr><td>", paste(row, collapse = "</td><td>"), "</td></tr>")
  })
  paste0("<table><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table>")
}

img_data_uri <- function(filename) {
  path <- file.path(out_dir, filename)
  if (!file.exists(path)) return("")
  paste0("data:image/png;base64,", base64enc::base64encode(path))
}

top_missing_samples <- head(
  sample_qc[order(-sample_qc$missing_rate), c("sample", "missing_rate", "heterozygote_rate", "median_dp", "mean_gq")],
  12
)
high_missing_samples <- subset(
  sample_qc,
  missing_rate > 0.10,
  select = c("sample", "missing_rate", "heterozygote_rate", "median_dp", "mean_gq")
)
top_het_sites <- head(
  site_qc[order(-site_qc$heterozygote_rate), c(
    "chrom", "pos", "variant_class", "missing_rate", "heterozygote_rate",
    "minor_allele_frequency", "mean_dp", "qual"
  )],
  12
)
monomorphic_after_mask <- subset(
  site_qc,
  minor_allele_frequency == 0,
  select = c("chrom", "pos", "variant_class", "missing_rate", "heterozygote_rate", "mean_dp", "qual")
)

suggestions <- c(
  sprintf(
    "The site-level missingness filter is behaving as expected: the maximum remaining site missingness is %.2f%%, below the 10%% cutoff.",
    100 * max(site_qc$missing_rate)
  ),
  sprintf(
    "%s samples have more than 10%% missing genotypes in this final variant VCF. Inspect these before haplotype or network analysis; if they drive topology, consider a sample-level missingness cutoff or a sensitivity run without them.",
    nrow(high_missing_samples)
  ),
  sprintf(
    "%s sites have MAF = 0 after genotype masking. For haplotype/network inputs, remove monomorphic records after masking so only informative variation remains.",
    nrow(monomorphic_after_mask)
  ),
  sprintf(
    "The highest site heterozygote rate is %.2f%%. Because plastome calls are usually treated as haploid or effectively haploid, inspect high-heterozygosity clusters and consider masking heterozygotes to missing before no-missing haplotype analyses.",
    100 * max(site_qc$heterozygote_rate)
  ),
  sprintf(
    "Depth and genotype quality look broadly strong after masking: median DP is %.0f and mean GQ is %.2f across called genotypes. The current DP >= 10, GQ/RGQ >= 30, and DP <= 10000 masks are reasonable for this dataset.",
    median(as.vector(dp), na.rm = TRUE),
    mean(as.vector(gq), na.rm = TRUE)
  )
)

filter_conditions_html <- filtering_conditions
names(filter_conditions_html) <- c("Step", "Condition", "Expression", "Source")
summary_html <- summary_qc
names(summary_html) <- c("Metric", "Value")

html <- c(
  "<!doctype html>",
  "<html lang=\"en\">",
  "<head>",
  "<meta charset=\"utf-8\">",
  "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
  "<title>Quercus globosa plastome VCF QC report</title>",
  "<style>",
  "body{margin:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;color:#1f2933;background:#f6f8fa;line-height:1.55}",
  "main{max-width:1120px;margin:0 auto;padding:28px 20px 56px}",
  "header{border-bottom:1px solid #d8dee4;margin-bottom:24px;padding-bottom:18px}",
  "h1{font-size:28px;line-height:1.2;margin:0 0 8px;font-weight:600}",
  "h2{font-size:21px;margin:30px 0 12px;font-weight:600}",
  "h3{font-size:16px;margin:22px 0 10px;font-weight:600}",
  "p{margin:0 0 12px}",
  ".muted{color:#57606a}",
  ".grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin:20px 0}",
  ".card{background:#fff;border:1px solid #d8dee4;border-radius:8px;padding:14px}",
  ".metric{font-size:25px;font-weight:600;margin-top:4px}",
  ".label{font-size:13px;color:#57606a}",
  ".plot-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}",
  "figure{margin:0;background:#fff;border:1px solid #d8dee4;border-radius:8px;padding:12px}",
  "figure img{width:100%;height:auto;display:block}",
  "figcaption{font-size:13px;color:#57606a;margin-top:8px}",
  "table{width:100%;border-collapse:collapse;background:#fff;border:1px solid #d8dee4;border-radius:8px;overflow:hidden;margin:10px 0 18px;font-size:13px}",
  "th,td{text-align:left;padding:8px 10px;border-bottom:1px solid #eaeef2;vertical-align:top}",
  "th{background:#eef2f6;font-weight:600}",
  "tr:last-child td{border-bottom:0}",
  "ul{margin-top:8px;padding-left:22px}",
  "li{margin:7px 0}",
  ".note{background:#fff7ed;border:1px solid #fed7aa;border-radius:8px;padding:12px;margin:14px 0}",
  "code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:0.95em}",
  "@media(max-width:820px){.grid,.plot-grid{grid-template-columns:1fr}main{padding:20px 14px 40px}h1{font-size:24px}}",
  "</style>",
  "</head>",
  "<body>",
  "<main>",
  "<header>",
  "<h1>Quercus globosa Plastome VCF QC And Filtering Report</h1>",
  sprintf("<p class=\"muted\">Generated %s from <code>%s</code>.</p>", html_escape(timestamp), html_escape(vcf_file)),
  "</header>",
  "<section class=\"grid\" aria-label=\"QC headline metrics\">",
  sprintf("<div class=\"card\"><div class=\"label\">Samples</div><div class=\"metric\">%s</div></div>", length(sample_names)),
  sprintf("<div class=\"card\"><div class=\"label\">Variant records</div><div class=\"metric\">%s</div></div>", nrow(site_qc)),
  sprintf("<div class=\"card\"><div class=\"label\">SNPs / indels</div><div class=\"metric\">%s / %s</div></div>", sum(site_qc$variant_class == "SNP"), sum(site_qc$variant_class == "indel")),
  sprintf("<div class=\"card\"><div class=\"label\">Ts/Tv</div><div class=\"metric\">%s</div></div>", summary_qc$value[summary_qc$metric == "Ts/Tv for SNP records"]),
  "</section>",
  "<section>",
  "<h2>Interpretation</h2>",
  sprintf(
    "<p>The final VCF contains %s PASS records across %s samples. Site missingness remains below the configured 10%% filter, but sample-level missingness and plastome heterozygote calls deserve additional review before downstream haplotype or network analyses.</p>",
    sum(site_qc$filter == "PASS"),
    length(sample_names)
  ),
  "<div class=\"note\"><strong>SNPfiltR status:</strong> SNPfiltR is not installed in this R library, so this report uses vcfR to compute equivalent QC summaries and records the SNPfiltR status for reproducibility.</div>",
  "</section>",
  "<section>",
  "<h2>Filtering Conditions</h2>",
  "<p>These conditions are parsed from the final VCF header and the alloutgroup filtering history.</p>",
  html_table(filter_conditions_html, digits = 4),
  "</section>",
  "<section>",
  "<h2>Suggestions</h2>",
  paste0("<ul>", paste0("<li>", html_escape(suggestions), "</li>", collapse = "\n"), "</ul>"),
  "</section>",
  "<section>",
  "<h2>Plots</h2>",
  "<div class=\"plot-grid\">",
  sprintf("<figure><img alt=\"Variant class counts\" src=\"%s\"><figcaption>SNP and indel counts retained in the final filtered variant VCF.</figcaption></figure>", img_data_uri("variant_class_counts.png")),
  sprintf("<figure><img alt=\"Site missingness histogram\" src=\"%s\"><figcaption>Remaining site missingness after genotype masking and the final missingness filter.</figcaption></figure>", img_data_uri("site_missingness_histogram.png")),
  sprintf("<figure><img alt=\"Sample missingness and heterozygosity scatter plot\" src=\"%s\"><figcaption>Samples with elevated missingness or heterozygosity are useful targets for sensitivity checks.</figcaption></figure>", img_data_uri("sample_missingness_heterozygosity.png")),
  sprintf("<figure><img alt=\"Site quality and depth scatter plot\" src=\"%s\"><figcaption>QUAL versus mean site depth; high depths are expected in plastome-enriched data but outliers should still be watched.</figcaption></figure>", img_data_uri("site_quality_depth.png")),
  sprintf("<figure><img alt=\"Heterozygote rate by plastome position\" src=\"%s\"><figcaption>High-heterozygosity clusters are candidates for masking or manual inspection in plastome analyses.</figcaption></figure>", img_data_uri("site_heterozygosity_by_position.png")),
  sprintf("<figure><img alt=\"Minor allele frequency histogram\" src=\"%s\"><figcaption>MAF distribution after filtering; MAF = 0 records are not informative after masking.</figcaption></figure>", img_data_uri("minor_allele_frequency_histogram.png")),
  "</div>",
  "</section>",
  "<section>",
  "<h2>QC Summary Table</h2>",
  html_table(summary_html, digits = 4),
  "</section>",
  "<section>",
  "<h2>Flagged Records</h2>",
  "<h3>Top samples by missingness</h3>",
  html_table(top_missing_samples, digits = 4),
  "<h3>Sites with highest heterozygote rate</h3>",
  html_table(top_het_sites, digits = 4),
  "</section>",
  "<section>",
  "<h2>Output Files</h2>",
  "<ul>",
  "<li><code>qc_summary.csv</code>, <code>filtering_conditions.csv</code>, <code>site_qc_vcfr.csv</code>, <code>sample_qc_vcfr.csv</code></li>",
  "<li><code>bcftools_filtering_commands_from_vcf_header.txt</code></li>",
  "<li>PNG plots embedded above and written separately in this folder</li>",
  "</ul>",
  "</section>",
  "</main>",
  "</body>",
  "</html>"
)

writeLines(html, file.path(out_dir, "Quercus_globosa_plastome_vcfr_snpfiltr_qc_report.html"))

cat("QC report written to:", normalizePath(out_dir), "\n")
