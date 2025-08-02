#10 bin figure

import gzip

input_file = "RS_AllSites_MethylState.tsv.gz"
output_file = "RS_stacked_values.tsv"

with gzip.open(input_file, 'rt') as fin, open(output_file, 'w') as fout:
    header = fin.readline().strip().split("\t")
    sample_cols = ["RS_008", "RS_016", "RS_017"]
    sample_indices = [header.index(c) for c in sample_cols]
    chr_idx = header.index("chr")

    for line in fin:
        fields = line.strip().split("\t")
        chr_raw = fields[chr_idx]
        if "." not in chr_raw:
            continue
        chrom, pos = chr_raw.split(".", 1)
        for idx in sample_indices:
            try:
                val = float(fields[idx])
                fout.write(f"{chrom}\t{pos}\t{val}\n")
            except ValueError:
                continue


library(OneR)
library(dplyr)
library(ggplot2)

# 读入三列数据：chr, start, value
input_file <- "RS_stacked_values.tsv"
data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(data) <- c("chr", "start", "value")

# 添加唯一 ID
data$ID <- paste0("ID_", seq_len(nrow(data)))

# Step 1: 分 bin（等频分箱为 10）
meth_bins <- OneR::bin(data[, c("ID", "value")], method = "content", nbins = 10)
meth_bins$mean_bin <- as.numeric(meth_bins$value)

# Step 2: 合并 bin 信息回原始数据
data_result <- merge(data, meth_bins[, c("ID", "mean_bin")], by = "ID")

# Step 3: 每个 bin 的均值与 site 数
mean_binmeth <- data_result %>%
  group_by(mean_bin) %>%
  summarise(
    binmean_mean = mean(value, na.rm = TRUE),
    bin_n = n(),
    .groups = "drop"
  )

# Step 4: 合并统计信息到每个 site
data_result <- merge(data_result, mean_binmeth, by = "mean_bin")

# Step 5: 准备 GTF 所需字段
data_result$source <- "none"
data_result$feature <- "transcript"
data_result$score <- "."
data_result$strand <- "+"
data_result$frame <- "."
data_result$seqname <- data_result$chr
data_result$gene_id <- paste0('gene_id "bin', data_result$mean_bin, "_", round(data_result$binmean_mean, 2), "_", data_result$bin_n, '"')
data_result$attribute <- paste(data_result$gene_id, ";")

# Step 6: 构建 GTF 表格
GTF <- data.frame(
  seqname = data_result$seqname,
  source = data_result$source,
  feature = data_result$feature,
  start = data_result$start,
  end = data_result$start,
  score = data_result$score,
  strand = data_result$strand,
  frame = data_result$frame,
  attribute = data_result$attribute
)

# Step 7: 写出 GTF 文件
out_file <- "RS_stacked_meth_value10bin.gtf"
write.table(GTF, file = out_file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

message(paste("✅ GTF written to", out_file))


sed -E 's/gene_id "([^"]+)"[[:space:]]*;/gene_id "\1";/' RS_stacked_meth_value10bin.fixed.gtf > RS_stacked_meth_value10bin.cleaned.gtf


perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/18_GOS.pileup \
        --gtf /work/cyu/meth/Gosling/GOS_stacked_meth_value10bin.cleaned.gtf  \
        --output /work/cyu/meth/Gosling/GOS_10bin.pi.txt \
        --measure pi

18_GOS.pileup



perl Variance-at-position.pl \
  --pool-size 400 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/19_ROB.pileup \
  --gtf /work/cyu/meth/Roberts/RS_stacked_meth_value10bin.cleaned.gtf \
  --output /work/cyu/meth/Roberts/ROB_10bin.pi.txt \
  --measure pi



library(OneR)
library(dplyr)

# 读入原始三列数据（chr, start, value）
input_file <- "RS_stacked_values.tsv"
data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(data) <- c("chr", "start", "value")

# 添加唯一 ID
data$ID <- paste0("ID_", seq_len(nrow(data)))

# Step 1: 使用 OneR 进行频率分箱（50 个 bin）
meth_bins <- OneR::bin(data[, c("ID", "value")], method = "content", nbins = 10)
meth_bins$mean_bin <- as.numeric(meth_bins$value)

# Step 2: 合并 bin 信息回原始数据，确保 value 被保留
data_result <- merge(data, meth_bins[, c("ID", "mean_bin")], by = "ID")

# Step 3: 每个 bin 的均值与 site 数
mean_binmeth <- data_result %>%
  group_by(mean_bin) %>%
  summarise(
    binmean_mean = mean(value, na.rm = TRUE),
    bin_n = n(),
    .groups = "drop"
  )

# Step 4: 合并统计信息到每个 site
data_result <- merge(data_result, mean_binmeth, by = "mean_bin")

# Step 5: 准备 GTF 所需字段
sitedata <- data_result
sitedata$source <- "none"
sitedata$feature <- "transcript"
sitedata$score <- "."
sitedata$strand <- "+"
sitedata$frame <- "."
sitedata$seqname <- sitedata$chr
sitedata$gene_id <- paste0('gene_id "bin', sitedata$mean_bin, "_", round(sitedata$binmean_mean, 2), "_", sitedata$bin_n, '"')
sitedata$attribute <- paste(sitedata$gene_id, ";")

# Step 6: 构建 GTF 表格
GTF <- data.frame(
  seqname = sitedata$seqname,
  source = sitedata$source,
  feature = sitedata$feature,
  start = sitedata$start,
  end = sitedata$start,
  score = sitedata$score,
  strand = sitedata$strand,
  frame = sitedata$frame,
  attribute = sitedata$attribute
)

# Step 7: 写出 GTF 文件
out_file <- "RS_stacked_meth_value10bin.gtf"
write.table(GTF, file = out_file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

message(paste("✅ GTF written to", out_file))




library(ggplot2)

# 读取 Popoolation 输出文件（按 bin）
pi_data <- read.table("/work/cyu/meth/RS_10bin.pi.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(pi_data) <- c("bin_id", "snp_count", "coverage_ratio", "pi")

# 添加一个标记列（后续可以用于分组，比如 Marine / Freshwater）
pi_data$Population <- "RS"

# 拆分 bin_id，例如 bin1_27.35_159925
bin_parts <- strsplit(pi_data$bin_id, "_")
bin_matrix <- do.call(rbind, bin_parts)

# 添加拆分后的列
pi_data$methylation_status <- bin_matrix[, 1]                        # bin编号（如 bin1）
pi_data$mean_methylation <- as.numeric(bin_matrix[, 2])             # bin内平均甲基化
pi_data$n_sites <- as.numeric(bin_matrix[, 3])                      # bin中CpG数量
pi_data$snp_count <- as.numeric(pi_data$snp_count)
pi_data$coverage_ratio <- as.numeric(pi_data$coverage_ratio)
pi_data$pi <- as.numeric(pi_data$pi)

# 生成图形对象
p <- ggplot(pi_data, aes(x = mean_methylation, y = pi)) +
  geom_point(shape = 21, fill = "steelblue4", color = "black", size = 2.5, na.rm = TRUE) +
  labs(
    x = "Rank Mean % Methylation",
    y = expression(paste(pi, " (per bin)")),
    title = "RS Population π per OneR Bin"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "none"
  )

# 保存图像
ggsave("RS__10bin_pi_plot.png", plot = p, width = 6, height = 4, dpi = 300)
