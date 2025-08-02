#Sensitive analysis for finding appropriate cut-off of meth files
#sensitive pi，which cut off should be used?
#rs
for n in {1..20}; do
  col_start=$(( (n-1)*3 + 4 ))
  col1=$col_start
  col2=$((col_start + 1))
  col3=$((col_start + 2))
  
  zcat RS_Meth_Sensitivity.tsv.gz | awk -F'\t' -v a=$col1 -v b=$col2 -v c=$col3 \
    '($a==1 && $b==1 && $c==1) {print $1, $2, $3, $a, $b, $c, $(NF)}' OFS='\t' > meth_"$n"_intersect.tsv
  
  echo "完成 meth_${n}_intersect.tsv"
done

# 定义函数
writeGTF_SiteClass <- function(input_file, output_file, gene_name) {
  sitedata <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  chr_pos <- strsplit(sitedata[,7], "\\.")  
  chr <- sapply(chr_pos, `[`, 1)
  pos <- sapply(chr_pos, `[`, 2)
  
  GTF <- data.frame(
    seqname = chr,
    source = "none",
    feature = "transcript",
    start = pos,
    end = pos,
    score = ".",
    strand = "+",
    frame = ".",
    attribute = paste0('gene_id "', gene_name, '";')
  )
  
  write.table(GTF, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# 批量处理 unmeth_1 到 unmeth_20
for (n in 1:20) {
  input_file <- paste0("meth_", n, "_intersect.tsv")
  output_file <- paste0("meth_", n, ".gtf")
  writeGTF_SiteClass(input_file, output_file, paste0("meth_", n))
  cat("完成 ", output_file, "\n")
}
#修改id带chr
import os

def modify_gtf(input_gtf, output_gtf):
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            chr_name = fields[0]
            # 修改 gene_id
            gene_id = f'gene_id "{chr_name}_{os.path.basename(input_gtf).replace(".gtf", "")}";'
            fields[8] = gene_id
            outfile.write('\t'.join(fields) + '\n')

# 循环处理 unmeth_2.gtf 到 unmeth_20.gtf
for n in range(1, 21):
    input_file = f"meth_{n}.gtf"
    output_file = f"meth_{n}_modified.gtf"
    modify_gtf(input_file, output_file)
    print(f"完成 {output_file}")



for n in {1..20}; do
  perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
    --gtf /work/cyu/meth/meth_${n}_modified.gtf \
    --output /work/cyu/meth/RS_meth_chrlevel_${n}_pi.txt \
    --measure pi
  echo "完成 meth_${n}_modified.gtf"
done






library(ggplot2)
library(dplyr)

# 初始化总数据表
all_data <- data.frame(Pi = numeric(), Group = character())

# 循环读取 RS_meth_chrlevel_1_pi.txt 到 RS_meth_chrlevel_20_pi.txt
for (n in 1:20) {
  file <- paste0("RS_meth_chrlevel_", n, "_pi.txt")
  if (file.exists(file)) {
    dat <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    dat <- dat[!is.na(dat[, 4]), ]  # 过滤 NA
    # 创建一个临时数据框
    temp <- data.frame(Pi = as.numeric(dat[, 4]), Group = paste0("meth_", n))
    all_data <- rbind(all_data, temp)
  } else {
    warning(paste("缺失文件:", file))
  }
}

# 明确设置 Group 的顺序
all_data$Group <- factor(all_data$Group, levels = paste0("meth_", 1:20))

# 绘图
p <- ggplot(all_data, aes(x = Group, y = Pi)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, color = "steelblue") +
  labs(title = "Chr-level π per methylated group (1–20)",
       x = "Methylated group",
       y = "Nucleotide Diversity (π)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 保存图像
ggsave("RS_meth_chrlevel_boxplot_fixed_final.png", p, width = 14, height = 6, dpi = 300)



for n in {1..20}; do
  perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
    --gtf /work/cyu/meth/unmeth_${n}_filtered.gtf \
    --output /work/cyu/meth/RS_unmeth_chrlevelfilter_${n}_pi.txt \
    --measure pi
  echo "完成 unmeth_${n}_filtered.gtf"
done
RS_unmeth_chrlevelfilter_1_pi.txt





library(ggplot2)
library(dplyr)

# 初始化总数据表
all_data <- data.frame(Pi = numeric(), Group = character())

# 循环读取 RS_meth_chrlevel_1_pi.txt 到 RS_meth_chrlevel_20_pi.txt
for (n in 1:20) {
  file <- paste0("RS_unmeth_chrlevelfilter_", n, "_pi.txt")
  if (file.exists(file)) {
    dat <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    dat <- dat[!is.na(dat[, 4]), ]  # 过滤 NA
    # 创建一个临时数据框
    temp <- data.frame(Pi = as.numeric(dat[, 4]), Group = paste0("unmeth_", n))
    all_data <- rbind(all_data, temp)
  } else {
    warning(paste("缺失文件:", file))
  }
}

# 明确设置 Group 的顺序
all_data$Group <- factor(all_data$Group, levels = paste0("unmeth_", 1:20))

# 绘图
p <- ggplot(all_data, aes(x = Group, y = Pi)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, color = "steelblue") +
  labs(title = "Chr-level π per methylated group (1–20)",
       x = "Methylated group",
       y = "Nucleotide Diversity (π)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 保存图像
ggsave("RS_unmeth_chrlevel_boxplot_filtered_final.png", p, width = 14, height = 6, dpi = 300)







#!/bin/bash

for cutoff in $(seq 85 -1 68); do
  # 计算每组的起始列号
  offset=$(( (85 - cutoff) * 3 ))
  col1=$((4 + offset))
  col2=$((5 + offset))
  col3=$((6 + offset))
  
  zcat RS_Meth_Sensitivity_Part2.tsv.gz | awk -F'\t' -v a=$col1 -v b=$col2 -v c=$col3 \
    '($a==1 && $b==1 && $c==1) {print $1, $2, $3, $a, $b, $c, $(NF)}' OFS='\t' > meth_"$cutoff"_intersect.tsv

  echo "✅ 完成 meth_${cutoff}_intersect.tsv"
done
meth_68_intersect.tsv




# 定义函数
writeGTF_SiteClass <- function(input_file, output_file, gene_name) {
  sitedata <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  chr_pos <- strsplit(sitedata[,7], "\\.")  
  chr <- sapply(chr_pos, `[`, 1)
  pos <- sapply(chr_pos, `[`, 2)
  
  GTF <- data.frame(
    seqname = chr,
    source = "none",
    feature = "transcript",
    start = pos,
    end = pos,
    score = ".",
    strand = "+",
    frame = ".",
    attribute = paste0('gene_id "', gene_name, '";')
  )
  
  write.table(GTF, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# 批量处理 unmeth_1 到 unmeth_20
for (n in 68:85) {
  input_file <- paste0("meth_", n, "_intersect.tsv")
  output_file <- paste0("meth_", n, ".gtf")
  writeGTF_SiteClass(input_file, output_file, paste0("meth_", n))
  cat("完成 ", output_file, "\n")
}


import os

def modify_gtf(input_gtf, output_gtf):
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            chr_name = fields[0]
            # 修改 gene_id
            gene_id = f'gene_id "{chr_name}_{os.path.basename(input_gtf).replace(".gtf", "")}";'
            fields[8] = gene_id
            outfile.write('\t'.join(fields) + '\n')

# 循环处理 unmeth_2.gtf 到 unmeth_20.gtf
for n in range(85,86):
    input_file = f"meth_{n}.gtf"
    output_file = f"meth_{n}_modified.gtf"
    modify_gtf(input_file, output_file)
    print(f"完成 {output_file}")



for n in {68..85}; do
  perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
    --gtf /work/cyu/meth/meth_${n}_modified.gtf \
    --output /work/cyu/meth/RS_meth_chrlevel_part2_${n}_pi.txt \
    --measure pi
  echo "完成 meth_${n}_modified.gtf"
done

library(ggplot2)
library(dplyr)

# 初始化数据框
all_data <- data.frame(Pi = numeric(), Group = character())

# 处理 cutoff 85 到 68（倒序）
for (cutoff in 68:85) {
  file <- paste0("RS_meth_chrlevel_part2_", cutoff, "_pi.txt")
  if (file.exists(file)) {
    dat <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    dat <- dat[!is.na(dat[, 4]), ]
    temp <- data.frame(Pi = as.numeric(dat[, 4]), Group = paste0("meth_", cutoff))
    all_data <- rbind(all_data, temp)
  } else {
    warning(paste("❌ 缺失文件:", file))
  }
}

# 设置 Group 的顺序
all_data$Group <- factor(all_data$Group, levels = paste0("meth_", 85:68))

# 绘图
p <- ggplot(all_data, aes(x = Group, y = Pi)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, color = "firebrick") +
  labs(title = "Chr-level π per methylated group (cutoff 85–68)",
       x = "Methylation cutoff group",
       y = "Nucleotide Diversity (π)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 保存图像
ggsave("RS_meth_chrlevel_boxplot_part2_cutoff85to68.png", p, width = 14, height = 6, dpi = 300)





#rob


#!/bin/bash
# ----------------------------------------
# ROB 群体未甲基化位点（unmeth） π 值计算流程
# ----------------------------------------

cd /work/cyu/meth/Roberts

for n in {1..20}; do
  # 每个 bin 占据 7 列，起始列为第 8 列（One 列）
  col1=$(( (n - 1) * 7 + 8 ))
  col2=$((col1 + 1))
  col3=$((col1 + 2))
  col4=$((col1 + 3))
  col5=$((col1 + 4))
  col6=$((col1 + 5))
  col7=$((col1 + 6))

  zcat ROB_Unmeth_Sensitivity.tsv.gz | awk -F'\t' -v a=$col1 -v b=$col2 -v c=$col3 -v d=$col4 -v e=$col5 -v f=$col6 -v g=$col7 \
    '($a==1 && $b==1 && $c==1 && $d==1 && $e==1 && $f==1 && $g==1) {print $a, $b, $c, $d, $e, $f, $g, $(NF)}' OFS='\t' > unmeth_${n}_intersect.tsv

  echo "✅ 完成 unmeth_${n}_intersect.tsv"
done

# Step 2: 生成每个 bin 的 GTF 文件
# R 脚本，建议保存为 generate_gtf.R 并在 shell 中用 Rscript 执行

writeGTF_SiteClass <- function(input_file, output_file, gene_name) {
  sitedata <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

  chr_pos <- strsplit(sitedata[, 8], "\\.")  # 最后一列是 chr.position
  chr <- sapply(chr_pos, `[`, 1)
  pos <- sapply(chr_pos, `[`, 2)

  GTF <- data.frame(
    seqname = chr,
    source = "none",
    feature = "transcript",
    start = pos,
    end = pos,
    score = ".",
    strand = "+",
    frame = ".",
    attribute = paste0('gene_id "', chr, "_", gene_name, '";')
  )

  write.table(GTF, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# 批量生成 GTF
for (n in 1:20) {
  input_file <- paste0("unmeth_", n, "_intersect.tsv")
  output_file <- paste0("unmeth_", n, ".gtf")
  writeGTF_SiteClass(input_file, output_file, paste0("unmeth_", n))
  cat("✅ 写出", output_file, "\n")
}
# Python 脚本
import os

def modify_gtf(input_gtf, output_gtf):
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            chr_name = fields[0]
            bin_name = os.path.basename(input_gtf).replace(".gtf", "")
            fields[8] = f'gene_id "{chr_name}_{bin_name}";'
            outfile.write('\t'.join(fields) + '\n')

# 批量处理
for n in range(1, 21):
    input_file = f"unmeth_{n}.gtf"
    output_file = f"unmeth_{n}_modified.gtf"
    modify_gtf(input_file, output_file)
    print(f"✅ 完成 {output_file}")


for n in {1..20}; do
  perl Variance-at-position.pl \
    --pool-size 400 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/19_ROB.pileup \
    --gtf /work/cyu/meth/Roberts/unmeth_${n}_modified.gtf \
    --output /work/cyu/meth/Roberts/ROB_unmeth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 unmeth_${n} π"
done


#rob meth
#!/bin/bash
cd /work/cyu/meth/Roberts
cd /work/cyu/meth/Roberts

# 生成从 99 到 70 共 30 个 bin 名称的数组
bins=($(seq 99 -1 70))

for i in "${!bins[@]}"; do
  bin="${bins[$i]}"
  n=$((i + 1))  # 索引从 0 开始

  # 每个 bin 对应 7 个样本列，起始列为 8 + (n-1)*7
  col_start=$(( (n - 1) * 7 + 8 ))
  col1=$col_start
  col2=$((col1 + 1))
  col3=$((col1 + 2))
  col4=$((col1 + 3))
  col5=$((col1 + 4))
  col6=$((col1 + 5))
  col7=$((col1 + 6))

  zcat ROB_Meth_Sensitivity.tsv.gz | \
    awk -F'\t' -v a=$col1 -v b=$col2 -v c=$col3 -v d=$col4 -v e=$col5 -v f=$col6 -v g=$col7 \
      '($a==1 && $b==1 && $c==1 && $d==1 && $e==1 && $f==1 && $g==1) {print $(NF)}' \
      > meth_${bin}_intersect.tsv

  echo "✅ 提取完成 meth_${bin}_intersect.tsv"
done

# 保存为 make_meth_gtf.R 并在 R 中运行
writeGTF_SiteClass <- function(input_file, output_file, gene_name) {
  sitedata <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  chr_pos <- strsplit(sitedata[, 1], "\\.")
  chr <- sapply(chr_pos, `[`, 1)
  pos <- sapply(chr_pos, `[`, 2)

  GTF <- data.frame(
    seqname = chr,
    source = "none",
    feature = "transcript",
    start = pos,
    end = pos,
    score = ".",
    strand = "+",
    frame = ".",
    attribute = paste0('gene_id "', gene_name, '";')
  )

  write.table(GTF, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# 生成 70 到 99 bin 的 GTF
bins <- 70:99
for (bin in bins) {
  infile <- sprintf("meth_%d_intersect.tsv", bin)
  outfile <- sprintf("meth_%d.gtf", bin)
  writeGTF_SiteClass(infile, outfile, paste0("meth_", bin))
  cat("✅ 写出", outfile, "\n")
}



import os

def modify_gtf(input_gtf, output_gtf):
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            chr_name = fields[0]
            bin_name = os.path.basename(input_gtf).replace(".gtf", "")
            fields[8] = f'gene_id "{chr_name}_{bin_name}";'
            outfile.write('\t'.join(fields) + '\n')

# 批量处理：meth_99 到 meth_70
for n in range(99, 69, -1):
    input_file = f"meth_{n}.gtf"
    output_file = f"meth_{n}_modified.gtf"
    modify_gtf(input_file, output_file)
    print(f"✅ 完成 {output_file}")




for n in {70..99}; do
  perl Variance-at-position.pl \
    --pool-size 400 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/19_ROB.pileup \
    --gtf /work/cyu/meth/Roberts/meth_${n}_modified.gtf \
    --output /work/cyu/meth/Roberts/ROB_meth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 meth_${n} π"
done


#gos

unmeth
for n in {1..20}; do
  col_start=$(( (n - 1) * 7 + 16 ))  # binary 状态的起始列，从第16列开始
  col1=$col_start
  col2=$((col1 + 1))
  col3=$((col1 + 2))
  col4=$((col1 + 3))
  col5=$((col1 + 4))
  col6=$((col1 + 5))
  col7=$((col1 + 6))

  zcat GOS_Unmeth_Sensitivity.tsv.gz | \
  awk -F'\t' -v a=$col1 -v b=$col2 -v c=$col3 -v d=$col4 -v e=$col5 -v f=$col6 -v g=$col7 \
      '($a==1 && $b==1 && $c==1 && $d==1 && $e==1 && $f==1 && $g==1) {print $a, $b, $c, $d, $e, $f, $g, $15}' \
      OFS='\t' > unmeth_"$n"_intersect.tsv

  echo "✅ 提取完成 unmeth_${n}_intersect.tsv"
done

# 写一个函数将 unmeth_n_intersect.tsv 转为 GTF
writeGTF_SiteClass <- function(input_file, output_file, gene_name) {
  sitedata <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

  chr_pos <- strsplit(sitedata[, 8], "\\.")  # 第 8 列是 chr.position
  chr <- sapply(chr_pos, `[`, 1)
  pos <- sapply(chr_pos, `[`, 2)

  GTF <- data.frame(
    seqname = chr,
    source = "none",
    feature = "transcript",
    start = pos,
    end = pos,
    score = ".",
    strand = "+",
    frame = ".",
    attribute = paste0('gene_id "', chr, "_", gene_name, '";')
  )

  write.table(GTF, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# 批量生成 unmeth_1 到 unmeth_20 的 GTF
for (n in 1:20) {
  input_file <- paste0("unmeth_", n, "_intersect.tsv")
  output_file <- paste0("unmeth_", n, ".gtf")
  writeGTF_SiteClass(input_file, output_file, paste0("unmeth_", n))
  cat("✅ 写出 ", output_file, "\n")
}

import os

def modify_gtf(input_gtf, output_gtf):
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            chr_name = fields[0]  # e.g., chrI
            bin_name = os.path.basename(input_gtf).replace(".gtf", "")
            fields[8] = f'gene_id "{chr_name}_{bin_name}";'
            outfile.write('\t'.join(fields) + '\n')

# 批量处理 unmeth_1.gtf 到 unmeth_20.gtf
for n in range(1, 21):
    input_file = f"unmeth_{n}.gtf"
    output_file = f"unmeth_{n}_modified.gtf"
    modify_gtf(input_file, output_file)
    print(f"✅ 完成 {output_file}")

for n in {1..20}; do
  perl Variance-at-position.pl \
    --pool-size 400 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/18_GOS.pileup \
    --gtf /work/cyu/meth/Gosling/unmeth_${n}_modified.gtf \
    --output /work/cyu/meth/Gosling/GOS_unmeth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 unmeth_${n} π"
done

#gos meth

#!/bin/bash
cd /work/cyu/meth/Gosling

# 生成从 99 到 70 共 30 个 bin 名称的数组
bins=($(seq 99 -1 70))

for i in "${!bins[@]}"; do
  bin="${bins[$i]}"
  n=$((i + 1))  # 索引从 0 开始

  # 每个 bin 对应 7 个样本列，起始列为 8 + (n-1)*7
  col_start=$(( (n - 1) * 7 + 8 ))
  col1=$col_start
  col2=$((col1 + 1))
  col3=$((col1 + 2))
  col4=$((col1 + 3))
  col5=$((col1 + 4))
  col6=$((col1 + 5))
  col7=$((col1 + 6))

  zcat GOS_Meth_Sensitivity.tsv.gz | \
    awk -F'\t' -v a=$col1 -v b=$col2 -v c=$col3 -v d=$col4 -v e=$col5 -v f=$col6 -v g=$col7 \
      '($a==1 && $b==1 && $c==1 && $d==1 && $e==1 && $f==1 && $g==1) {print $(NF)}' \
      > meth_${bin}_intersect.tsv

  echo "✅ 提取完成 meth_${bin}_intersect.tsv"
done

writeGTF_SiteClass <- function(input_file, output_file, gene_name) {
  sitedata <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  chr_pos <- strsplit(sitedata[, 1], "\\.")
  chr <- sapply(chr_pos, `[`, 1)
  pos <- sapply(chr_pos, `[`, 2)

  GTF <- data.frame(
    seqname = chr,
    source = "none",
    feature = "transcript",
    start = pos,
    end = pos,
    score = ".",
    strand = "+",
    frame = ".",
    attribute = paste0('gene_id "', gene_name, '";')
  )

  write.table(GTF, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# 生成 GTF 文件 meth_99 到 meth_70
bins <- 99:70
for (bin in bins) {
  infile <- sprintf("meth_%d_intersect.tsv", bin)
  outfile <- sprintf("meth_%d.gtf", bin)
  writeGTF_SiteClass(infile, outfile, paste0("meth_", bin))
  cat("✅ 写出", outfile, "\n")
}

import os

def modify_gtf(input_gtf, output_gtf):
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            chr_name = fields[0]
            bin_name = os.path.basename(input_gtf).replace(".gtf", "")
            fields[8] = f'gene_id "{chr_name}_{bin_name}";'
            outfile.write('\t'.join(fields) + '\n')

# 处理 meth_70 到 meth_99
for n in range(70, 100):
    input_file = f"meth_{n}.gtf"
    output_file = f"meth_{n}_modified.gtf"
    modify_gtf(input_file, output_file)
    print(f"✅ 完成 {output_file}")

for n in {70..99}; do
  perl Variance-at-position.pl \
    --pool-size 400 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/18_GOS.pileup \
    --gtf /work/cyu/meth/Gosling/meth_${n}_modified.gtf \
    --output /work/cyu/meth/Gosling/GOS_meth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 meth_${n} π"
done




#wik unmeth 

#!/bin/bash

# 进入工作目录
cd /work/cyu/meth/Wik

# 每个 bin 有 8 个 binary 列，从第 18 列开始
for n in {1..20}; do
  # 每组8列，从第18列开始
  col_start=$(( (n - 1) * 8 + 18 ))
  col1=$col_start
  col2=$((col1 + 1))
  col3=$((col1 + 2))
  col4=$((col1 + 3))
  col5=$((col1 + 4))
  col6=$((col1 + 5))
  col7=$((col1 + 6))
  col8=$((col1 + 7))

  # 提取满足所有8列都为1的行，并输出 chr.position（在第17列）
  zcat WK_Unmeth_Sensitivity.tsv.gz | \
    awk -F'\t' -v a=$col1 -v b=$col2 -v c=$col3 -v d=$col4 -v e=$col5 -v f=$col6 -v g=$col7 -v h=$col8 \
      '($a==1 && $b==1 && $c==1 && $d==1 && $e==1 && $f==1 && $g==1 && $h==1) {print $17}' \
      > unmeth_"$n"_intersect.tsv

  echo "✅ 提取完成 unmeth_${n}_intersect.tsv"
done

# 写一个函数将 unmeth_n_intersect.tsv 转为 GTF
writeGTF_SiteClass <- function(input_file, output_file, gene_name) {
  sitedata <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

  # 拆分 chr.position 列（第1列）
  chr_pos <- strsplit(sitedata[, 1], "\\.")
  chr <- sapply(chr_pos, `[`, 1)
  pos <- sapply(chr_pos, `[`, 2)

  # 构建 GTF 表格
  GTF <- data.frame(
    seqname = chr,
    source = "none",
    feature = "transcript",
    start = pos,
    end = pos,
    score = ".",
    strand = "+",
    frame = ".",
    attribute = paste0('gene_id "unmeth_', gene_name, '";')
  )

  # 写出 GTF 文件
  write.table(GTF, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# 批量生成 unmeth_1 到 unmeth_20 的 GTF
for (n in 1:20) {
  input_file <- paste0("unmeth_", n, "_intersect.tsv")
  output_file <- paste0("unmeth_", n, ".gtf")
  writeGTF_SiteClass(input_file, output_file, as.character(n))
  cat("✅ 写出 ", output_file, "\n")
}


import os

def modify_gtf(input_gtf, output_gtf):
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            chr_name = fields[0]  # e.g., chrI
            bin_name = os.path.basename(input_gtf).replace(".gtf", "")
            fields[8] = f'gene_id "{chr_name}_{bin_name}";'
            outfile.write('\t'.join(fields) + '\n')

# 批量处理 unmeth_1.gtf 到 unmeth_20.gtf
for n in range(1, 21):
    input_file = f"unmeth_{n}.gtf"
    output_file = f"unmeth_{n}_modified.gtf"
    modify_gtf(input_file, output_file)
    print(f"✅ 完成 {output_file}")

for n in {1..20}; do
  perl Variance-at-position.pl \
    --pool-size 198 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/8_WK.pileup \
    --gtf /work/cyu/meth/Wik/unmeth_${n}_modified.gtf \
    --output /work/cyu/meth/Wik/WK_unmeth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 unmeth_${n} π"
done







#wik meth
#!/bin/bash
cd /work/cyu/meth/Wik

# 生成从 99 到 70 的 bin
bins=($(seq 99 -1 70))

for i in "${!bins[@]}"; do
  bin="${bins[$i]}"
  n=$((i + 1))  # bin 的编号 1 到 30

  # 每个 bin 8 列，起始列为第9列（第1个bin）
  col_start=$((8 * (n - 1) + 9))
  col1=$col_start
  col2=$((col1 + 1))
  col3=$((col1 + 2))
  col4=$((col1 + 3))
  col5=$((col1 + 4))
  col6=$((col1 + 5))
  col7=$((col1 + 6))
  col8=$((col1 + 7))

  zcat WK_Meth_Sensitivity.tsv.gz | \
    awk -F'\t' -v a=$col1 -v b=$col2 -v c=$col3 -v d=$col4 \
                    -v e=$col5 -v f=$col6 -v g=$col7 -v h=$col8 \
                    '($a==1 && $b==1 && $c==1 && $d==1 && $e==1 && $f==1 && $g==1 && $h==1) {print $(NF)}' \
    > meth_${bin}_intersect.tsv

  echo "✅ 提取完成 meth_${bin}_intersect.tsv"
done


writeGTF_SiteClass <- function(input_file, output_file, gene_name) {
  sitedata <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  chr_pos <- strsplit(sitedata[, 1], "\\.")
  chr <- sapply(chr_pos, `[`, 1)
  pos <- sapply(chr_pos, `[`, 2)

  GTF <- data.frame(
    seqname = chr,
    source = "none",
    feature = "transcript",
    start = pos,
    end = pos,
    score = ".",
    strand = "+",
    frame = ".",
    attribute = paste0('gene_id "', gene_name, '";')
  )

  write.table(GTF, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# 批量生成 meth_70 到 meth_99 的 GTF 文件
for (bin in 70:99) {
  infile <- sprintf("meth_%d_intersect.tsv", bin)
  outfile <- sprintf("meth_%d.gtf", bin)
  writeGTF_SiteClass(infile, outfile, paste0("meth_", bin))
  cat("✅ 写出 ", outfile, "\n")
}




import os

def modify_gtf(input_gtf, output_gtf):
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            chr_name = fields[0]
            bin_name = os.path.basename(input_gtf).replace(".gtf", "")
            fields[8] = f'gene_id "{chr_name}_{bin_name}";'
            outfile.write('\t'.join(fields) + '\n')

# 处理 meth_70 到 meth_99
for n in range(70, 100):
    input_file = f"meth_{n}.gtf"
    output_file = f"meth_{n}_modified.gtf"
    modify_gtf(input_file, output_file)
    print(f"✅ 完成 {output_file}")



for n in {70..99}; do
  perl Variance-at-position.pl \
    --pool-size 198 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/8_WK.pileup \
    --gtf /work/cyu/meth/Wik/meth_${n}_modified.gtf \
    --output /work/cyu/meth/Wik/WK_meth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 meth_${n} π"
done





library(ggplot2)
library(dplyr)

# 初始化总数据表
all_data <- data.frame(Pi = numeric(), Group = character())

# GOS
for (n in 1:20) {
  file <- paste0("WK_unmeth_chrlevel_", n, "_pi.txt")
  if (file.exists(file)) {
    dat <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    dat <- dat[!is.na(dat[, 4]), ]  # 过滤 NA
    # 创建一个临时数据框
    temp <- data.frame(Pi = as.numeric(dat[, 4]), Group = paste0("unmeth_", n))
    all_data <- rbind(all_data, temp)
  } else {
    warning(paste("缺失文件:", file))
  }
}

# 明确设置 Group 的顺序
all_data$Group <- factor(all_data$Group, levels = paste0("unmeth_", 1:20))

# 绘图
p <- ggplot(all_data, aes(x = Group, y = Pi)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, color = "steelblue") +
  labs(title = "Chr-level π per methylated group (1–20)",
       x = "Methylated group",
       y = "Nucleotide Diversity (π)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 保存图像
ggsave("WK_unmeth_chrlevel_boxplot.png", p, width = 14, height = 6, dpi = 300)


library(ggplot2)
library(dplyr)

# 初始化总数据表
all_data <- data.frame(Pi = numeric(), Group = character())

# GOS
for (n in 70:99) {
  file <- paste0("WK_meth_chrlevel_", n, "_pi.txt")
  if (file.exists(file)) {
    dat <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    dat <- dat[!is.na(dat[, 4]), ]  # 过滤 NA
    # 创建一个临时数据框
    temp <- data.frame(Pi = as.numeric(dat[, 4]), Group = paste0("meth_", n))
    all_data <- rbind(all_data, temp)
  } else {
    warning(paste("缺失文件:", file))
  }
}

# 明确设置 Group 的顺序
all_data$Group <- factor(all_data$Group, levels = paste0("meth_", 70:99))

# 绘图
p <- ggplot(all_data, aes(x = Group, y = Pi)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, color = "steelblue") +
  labs(title = "Chr-level π per methylated group (1–20)",
       x = "Methylated group",
       y = "Nucleotide Diversity (π)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 保存图像
ggsave("WK_meth_chrlevel_boxplot.png", p, width = 14, height = 6, dpi = 300)





#reduce rob
for n in {1..20}; do
  # 每个 bin 占据 5 列，起始列为第 6 列
  col1=$(( (n - 1) * 5 + 6 ))
  col2=$((col1 + 1))
  col3=$((col1 + 2))
  col4=$((col1 + 3))
  col5=$((col1 + 4))

  zcat Rob_Unmeth_Sensitivity.tsv.gz | awk -F'\t' -v a=$col1 -v b=$col2 -v c=$col3 -v d=$col4 -v e=$col5 \
    '($a==1 && $b==1 && $c==1 && $d==1 && $e==1) {print $a, $b, $c, $d, $e, $(NF)}' OFS='\t' > unmeth_${n}_intersect.tsv

  echo "✅ 完成 unmeth_${n}_intersect.tsv"
done


writeGTF_SiteClass <- function(input_file, output_file, gene_name) {
  sitedata <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

  # 拆分最后一列 chr.position 为 chr 和 pos
  chr_pos <- strsplit(sitedata[, 6], "\\.")  # 第6列是 chr.position
  chr <- sapply(chr_pos, `[`, 1)
  pos <- sapply(chr_pos, `[`, 2)

  # 创建 GTF 格式
  GTF <- data.frame(
    seqname = chr,
    source = "none",
    feature = "transcript",
    start = pos,
    end = pos,
    score = ".",
    strand = "+",
    frame = ".",
    attribute = paste0('gene_id "', chr, "_", gene_name, '";')
  )

  # 写入文件
  write.table(GTF, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# 循环处理 1–20 个 bin
for (n in 1:20) {
  input_file <- paste0("unmeth_", n, "_intersect.tsv")
  output_file <- paste0("unmeth_", n, ".gtf")
  writeGTF_SiteClass(input_file, output_file, paste0("unmeth_", n))
  cat("✅ 写出", output_file, "\n")
}






# Python 脚本
import os

def modify_gtf(input_gtf, output_gtf):
    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            chr_name = fields[0]
            bin_name = os.path.basename(input_gtf).replace(".gtf", "")
            fields[8] = f'gene_id "{chr_name}_{bin_name}";'
            outfile.write('\t'.join(fields) + '\n')

# 批量处理
for n in range(1, 21):
    input_file = f"unmeth_{n}.gtf"
    output_file = f"unmeth_{n}_modified.gtf"
    modify_gtf(input_file, output_file)
    print(f"✅ 完成 {output_file}")


for n in {1..20}; do
  perl Variance-at-position.pl \
    --pool-size 400 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/19_ROB.pileup \
    --gtf /work/cyu/meth/Roberts/ReducedSamples/unmeth_${n}_modified.gtf \
    --output /work/cyu/meth/Roberts/ReducedSamples/reduced_ROB_unmeth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 unmeth_${n} π"
done








for n in {68..99}; do
  perl Variance-at-position.pl \
    --pool-size 400 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/19_ROB.pileup \
    --gtf /work/cyu/meth/Roberts/ReducedSamples/meth_${n}_modified.gtf \
    --output /work/cyu/meth/Roberts/ReducedSamples/reduced_ROB_meth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 meth_${n} π"
done



for n in {1..20}; do
  perl Variance-at-position.pl \
    --pool-size 400 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/18_GOS.pileup \
    --gtf /work/cyu/meth/Gosling/ReducedSamples/GOS_unmeth_${n}_modified.gtf \
    --output /work/cyu/meth/Gosling/ReducedSamples/GOS_unmeth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 unmeth_${n} π"
done

for n in {68..99}; do
  perl Variance-at-position.pl \
    --pool-size 400 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/18_GOS.pileup \
    --gtf /work/cyu/meth/Gosling/ReducedSamples/GOS_meth_${n}_modified.gtf \
    --output /work/cyu/meth/Gosling/ReducedSamples/GOS_meth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 unmeth_${n} π"
done

#sayward
for n in {1..20}; do
  perl Variance-at-position.pl \
    --pool-size 200 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/17_SAY.pileup \
    --gtf /work/cyu/meth/Sayward/ReducedSamples/SAY_unmeth_${n}_modified.gtf \
    --output /work/cyu/meth/Sayward/ReducedSamples/SAY_unmeth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 unmeth_${n} π"
done

for n in {68..99}; do
  perl Variance-at-position.pl \
    --pool-size 200 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/17_SAY.pileup \
    --gtf /work/cyu/meth/Sayward/ReducedSamples/SAY_meth_${n}_modified.gtf \
    --output /work/cyu/meth/Sayward/ReducedSamples/SAY_meth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 unmeth_${n} π"
done



for n in {1..20}; do
  perl Variance-at-position.pl \
    --pool-size 198 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/8_WK.pileup \
    --gtf /work/cyu/meth/Wik/WGBS/WK_unmeth_${n}_modified.gtf \
    --output /work/cyu/meth/Wik/WGBS/WK_unmeth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 unmeth_${n} π"
done

for n in {68..99}; do
  perl Variance-at-position.pl \
    --pool-size 198 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/8_WK.pileup \
    --gtf /work/cyu/meth/Wik/WGBS/WK_meth_${n}_modified.gtf \
    --output /work/cyu/meth/Wik/WGBS/WK_meth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 unmeth_${n} π"
done


for n in {1..20}; do
  perl Variance-at-position.pl \
    --pool-size 184 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/7_WT.pileup \
    --gtf /work/cyu/meth/Watson/WGBS/WT_unmeth_${n}_modified.gtf \
    --output /work/cyu/meth/Watson/WGBS/WT_unmeth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 unmeth_${n} π"
done

for n in {68..99}; do
  perl Variance-at-position.pl \
    --pool-size 184 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/7_WT.pileup \
    --gtf /work/cyu/meth/Watson/WGBS/WT_meth_${n}_modified.gtf \
    --output /work/cyu/meth/Watson/WGBS/WT_meth_chrlevel_${n}_pi.txt \
    --measure pi

  echo "✅ 计算完 unmeth_${n} π"
done


library(ggplot2)
library(dplyr)

# 初始化总数据表
all_data <- data.frame(Pi = numeric(), Group = character())

# 循环读取 reduced_ROB_unmeth_chrlevel_1_pi.txt 到 _20
for (n in 1:20) {
  file <- sprintf("/work/cyu/meth/Roberts/ReducedSamples/reduced_ROB_unmeth_chrlevel_%d_pi.txt", n)
  if (file.exists(file)) {
    dat <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    dat <- dat[!is.na(dat[, 4]), ]  # 移除 NA
    temp <- data.frame(Pi = as.numeric(dat[, 4]), Group = paste0("unmeth_", n))
    all_data <- rbind(all_data, temp)
  } else {
    warning(paste("❌ 缺失文件:", file))
  }
}

# 设置 Group 因子顺序
all_data$Group <- factor(all_data$Group, levels = paste0("unmeth_", 1:20))

# 绘图
p <- ggplot(all_data, aes(x = Group, y = Pi)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, color = "tomato") +
  labs(title = "Roberts unmethylated chr-level π (bins 1–20)",
       x = "Unmethylated bin",
       y = "Nucleotide Diversity (π)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 保存图像
ggsave("ROB_unmeth_chrlevel_boxplot_filtered_final.png", p, width = 14, height = 6, dpi = 300)

cat("✅ 图像保存成功：ROB_unmeth_chrlevel_boxplot_filtered_final.png\n")

library(ggplot2)
library(dplyr)

# 初始化空数据框
all_data <- data.frame(Pi = numeric(), Group = character())

# 循环读取 reduced_ROB_meth_chrlevel_68_pi.txt 到 99
for (n in 68:99) {
  file <- sprintf("/work/cyu/meth/Roberts/ReducedSamples/reduced_ROB_meth_chrlevel_%d_pi.txt", n)
  if (file.exists(file)) {
    dat <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    dat <- dat[!is.na(dat[, 4]), ]
    temp <- data.frame(Pi = as.numeric(dat[, 4]), Group = paste0("meth_", n))
    all_data <- rbind(all_data, temp)
  } else {
    warning(paste("❌ 缺失文件:", file))
  }
}

# 设置 Group 的顺序
all_data$Group <- factor(all_data$Group, levels = paste0("meth_", 68:99))

# 绘制箱线图
p <- ggplot(all_data, aes(x = Group, y = Pi)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, color = "darkgreen") +
  labs(title = "Roberts methylated chr-level π (bins 68–99)",
       x = "Methylated bin",
       y = "Nucleotide Diversity (π)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 保存图像
ggsave("ROB_meth_chrlevel_boxplot_filtered_final.png", p, width = 14, height = 6, dpi = 300)

cat("✅ 图像保存成功：ROB_meth_chrlevel_boxplot_filtered_final.png\n")


