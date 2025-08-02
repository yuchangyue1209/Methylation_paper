#population level unmeth-inter-meth pi boxplot

#intermethunmeth.gtf for each population 
#RS
import pandas as pd
import gzip
import os

# 输入文件及对应的甲基化类型
methylation_files = {
    "RS_Intersection_Inter.tsv.gz": "inter",
    "RS_Intersection_Meth.tsv.gz": "meth",
    "RS_Intersection_Unmeth.tsv.gz": "unmeth"
}

gtf_rows = []

for fname, label in methylation_files.items():
    if os.path.exists(fname):
        print(f"Reading {fname}...")
        with gzip.open(fname, "rt") as f:
            df = pd.read_csv(f, sep="\t", header=0)  # ✅ 显式设为第一行为 header
            for idx, row in df.iterrows():
                chrom_pos = row["chr"]               # ✅ 用列名访问而不是 .iloc
                if "." not in chrom_pos:
                    print(f"⚠️ Skipping malformed row: {chrom_pos}")
                    continue
                chrom, pos = chrom_pos.split(".")
                start = int(pos)
                end = start

                gtf_row = [
                    chrom,
                    "none",
                    "transcript",
                    start,
                    end,
                    ".",
                    "+",
                    ".",
                    f'gene_id "{chrom}_{label}";'
                ]
                gtf_rows.append(gtf_row)
    else:
        print(f"❌ File not found: {fname}")

# 写出 GTF 文件
if gtf_rows:
    gtf_df = pd.DataFrame(gtf_rows)
    gtf_df.to_csv("RS_intermethunmeth.gtf", sep="\t", index=False, header=False, quoting=3)
    print("✅ GTF 文件已生成：RS_intermethunmeth.gtf")
else:
    print("⚠️ 没有生成任何 GTF 内容，请检查输入文件")



perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/RabbitSlough/RS_intermethunmeth.gtf \
        --output /work/cyu/meth/RabbitSlough/RS_intermethunmeth.pi.txt \
        --measure pi



#r script for boxplot
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggsignif)

############################
# 1. 读取数据并预处理
all_data <- read.table("RS_intermethunmeth.pi.txt", 
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# 添加列名
colnames(all_data) <- c("chr_result", "col2", "col3", "statValue")

# 拆分 chr_result 为 chr 和 result（unmeth/inter/meth）
all_data <- separate(all_data, 
                     col = "chr_result", 
                     into = c("chr", "result"), 
                     sep = "_", 
                     remove = TRUE)

# 排除 chrUn 和 chrM
all_data <- filter(all_data, !(chr %in% c("chrUn", "chrM")))

# 设置 methylation 分类顺序
all_data$result <- factor(all_data$result, levels = c("unmeth", "inter", "meth"))

# 标记性染色体（仅 chrY）
all_data$sex_chr <- ifelse(all_data$chr == "chrY", "sex", "autosome")

# 清理 NaN 并转为数值
all_data$statValue[all_data$statValue == "NaN"] <- NA
all_data$statValue <- as.numeric(all_data$statValue)
all_data <- all_data[!is.na(all_data$statValue), ]

############################
# 2. Wilcoxon 检验（两两比较）
comparisons <- combn(levels(all_data$result), 2, simplify = FALSE)
wt_results <- data.frame()

for (comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  
  x <- all_data %>% filter(result == group1) %>% pull(statValue)
  y <- all_data %>% filter(result == group2) %>% pull(statValue)
  
  wtest <- wilcox.test(x, y, paired = FALSE)
  
  wt_results <- rbind(wt_results, 
                      data.frame(group1 = group1, 
                                 group2 = group2, 
                                 W = wtest$statistic, 
                                 p.value = wtest$p.value))
}

print(wt_results)

############################
# 3. 可视化
bg_pi <- 0.00483154
max_val <- max(all_data$statValue, na.rm = TRUE)

p <- ggplot(all_data, aes(x = result, y = statValue, fill = result)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.8, size = 2,
              aes(shape = sex_chr, color = sex_chr)) +
  geom_hline(yintercept = bg_pi, linetype = "dotted", color = "black", size = 0.6) +
  scale_shape_manual(values = c("autosome" = 21, "sex" = 24)) +
  scale_color_manual(values = c("autosome" = "black", "sex" = "red")) +
  theme_bw() +
  labs(x = "Methylation Category", 
       y = expression(pi),
       title = "Nucleotide Diversity by Methylation State (RS)") +
  theme(legend.position = "right")

# 添加显著性标注
p <- p + geom_signif(comparisons = comparisons,
                     annotations = format.pval(wt_results$p.value, digits = 2, eps = 0.001),
                     y_position = seq(from = max_val * 1.05, length.out = length(comparisons), by = 0.002),
                     tip_length = 0.02, textsize = 4)

# 显示图形
print(p)

# 保存图像
ggsave("RS_methylation_pi_boxplot.png", p, width = 8, height = 6, dpi = 300)







#Rob
import pandas as pd
import gzip
import os

# 输入文件及对应的甲基化类型
methylation_files = {
    "Intersection_ROB_Inter.tsv.gz": "inter",
    "Intersection_ROB_Meth.tsv.gz": "meth",
    "ROB_Intersection_Unmeth.tsv.gz": "unmeth"
}

gtf_rows = []

for fname, label in methylation_files.items():
    if os.path.exists(fname):
        print(f"Reading {fname}...")
        with gzip.open(fname, "rt") as f:
            df = pd.read_csv(f, sep="\t", header=0)  # ✅ 显式设为第一行为 header
            for idx, row in df.iterrows():
                chrom_pos = row["chr"]               # ✅ 用列名访问而不是 .iloc
                if "." not in chrom_pos:
                    print(f"⚠️ Skipping malformed row: {chrom_pos}")
                    continue
                chrom, pos = chrom_pos.split(".")
                start = int(pos)
                end = start

                gtf_row = [
                    chrom,
                    "none",
                    "transcript",
                    start,
                    end,
                    ".",
                    "+",
                    ".",
                    f'gene_id "{chrom}_{label}";'
                ]
                gtf_rows.append(gtf_row)
    else:
        print(f"❌ File not found: {fname}")

# 写出 GTF 文件
if gtf_rows:
    gtf_df = pd.DataFrame(gtf_rows)
    gtf_df.to_csv("Rob_intermethunmeth.gtf", sep="\t", index=False, header=False, quoting=3)
    print("✅ GTF 文件已生成：Rob_intermethunmeth.gtf")
else:
    print("⚠️ 没有生成任何 GTF 内容，请检查输入文件")


perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/19_ROB.pileup \
        --gtf /work/cyu/meth/Roberts/Rob_intermethunmeth.gtf \
        --output /work/cyu/meth/Roberts/Rob_intermethunmeth.pi.txt \
        --measure pi


library(tidyr)
library(ggplot2)
library(dplyr)
library(ggsignif)

############################
# 1. 读取数据并预处理
all_data <- read.table("RS_intermethunmeth.pi.txt", 
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# 添加列名
colnames(all_data) <- c("chr_result", "col2", "col3", "statValue")

# 拆分 chr_result 为 chr 和 result（unmeth/inter/meth）
all_data <- separate(all_data, 
                     col = "chr_result", 
                     into = c("chr", "result"), 
                     sep = "_", 
                     remove = TRUE)

# 排除 chrUn 和 chrM
all_data <- filter(all_data, !(chr %in% c("chrUn", "chrM")))

# 设置 methylation 分类顺序
all_data$result <- factor(all_data$result, levels = c("unmeth", "inter", "meth"))

# 标记性染色体（仅 chrY）
all_data$sex_chr <- ifelse(all_data$chr == "chrY", "sex", "autosome")

# 清理 NaN 并转为数值
all_data$statValue[all_data$statValue == "NaN"] <- NA
all_data$statValue <- as.numeric(all_data$statValue)
all_data <- all_data[!is.na(all_data$statValue), ]

############################
# 2. Wilcoxon 检验（两两比较）
comparisons <- combn(levels(all_data$result), 2, simplify = FALSE)
wt_results <- data.frame()

for (comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  
  x <- all_data %>% filter(result == group1) %>% pull(statValue)
  y <- all_data %>% filter(result == group2) %>% pull(statValue)
  
  wtest <- wilcox.test(x, y, paired = FALSE)
  
  wt_results <- rbind(wt_results, 
                      data.frame(group1 = group1, 
                                 group2 = group2, 
                                 W = wtest$statistic, 
                                 p.value = wtest$p.value))
}

print(wt_results)

############################
# 3. 可视化
bg_pi <- 0.00456532
max_val <- max(all_data$statValue, na.rm = TRUE)

p <- ggplot(all_data, aes(x = result, y = statValue, fill = result)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.8, size = 2,
              aes(shape = sex_chr, color = sex_chr)) +
  geom_hline(yintercept = bg_pi, linetype = "dotted", color = "black", size = 0.6) +
  scale_shape_manual(values = c("autosome" = 21, "sex" = 24)) +
  scale_color_manual(values = c("autosome" = "black", "sex" = "red")) +
  theme_bw() +
  labs(x = "Methylation Category", 
       y = expression(pi),
       title = "Nucleotide Diversity by Methylation State (RS)") +
  theme(legend.position = "right")

# 添加显著性标注
p <- p + geom_signif(comparisons = comparisons,
                     annotations = format.pval(wt_results$p.value, digits = 2, eps = 0.001),
                     y_position = seq(from = max_val * 1.05, length.out = length(comparisons), by = 0.002),
                     tip_length = 0.02, textsize = 4)

# 显示图形
print(p)

# 保存图像
ggsave("RS_methylation_pi_boxplot.png", p, width = 8, height = 6, dpi = 300)

#gos
import pandas as pd
import gzip
import os

# 输入文件及对应的甲基化类型
methylation_files = {
    "Gos_Intersection_Inter.tsv.gz": "inter",
    "Gos_Intersection_Meth.tsv.gz": "meth",
    "Gos_Intersection_Unmeth.tsv.gz": "unmeth"
}

gtf_rows = []

for fname, label in methylation_files.items():
    if os.path.exists(fname):
        print(f"Reading {fname}...")
        with gzip.open(fname, "rt") as f:
            df = pd.read_csv(f, sep="\t", header=0)  # ✅ 显式设为第一行为 header
            for idx, row in df.iterrows():
                chrom_pos = row["chr"]               # ✅ 用列名访问而不是 .iloc
                if "." not in chrom_pos:
                    print(f"⚠️ Skipping malformed row: {chrom_pos}")
                    continue
                chrom, pos = chrom_pos.split(".")
                start = int(pos)
                end = start

                gtf_row = [
                    chrom,
                    "none",
                    "transcript",
                    start,
                    end,
                    ".",
                    "+",
                    ".",
                    f'gene_id "{chrom}_{label}";'
                ]
                gtf_rows.append(gtf_row)
    else:
        print(f"❌ File not found: {fname}")

# 写出 GTF 文件
if gtf_rows:
    gtf_df = pd.DataFrame(gtf_rows)
    gtf_df.to_csv("Gos_intermethunmeth.gtf", sep="\t", index=False, header=False, quoting=3)
    print("✅ GTF 文件已生成：Gos_intermethunmeth.gtf")
else:
    print("⚠️ 没有生成任何 GTF 内容，请检查输入文件")


perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/18_GOS.pileup \
        --gtf /work/cyu/meth/Gosling/Gos_intermethunmeth.gtf \
        --output /work/cyu/meth/Gosling/Gos_intermethunmeth.pi.txt \
        --measure pi











#wik
import pandas as pd
import gzip
import os

# 输入文件及对应的甲基化类型
methylation_files = {
    "WK_Intersection_Inter.tsv.gz": "inter",
    "WK_Intersection_Meth.tsv.gz": "meth",
    "WK_Intersection_Unmeth.tsv.gz": "unmeth"
}

gtf_rows = []

for fname, label in methylation_files.items():
    if os.path.exists(fname):
        print(f"Reading {fname}...")
        with gzip.open(fname, "rt") as f:
            df = pd.read_csv(f, sep="\t", header=0)  # ✅ 显式设为第一行为 header
            for idx, row in df.iterrows():
                chrom_pos = row["chr"]               # ✅ 用列名访问而不是 .iloc
                if "." not in chrom_pos:
                    print(f"⚠️ Skipping malformed row: {chrom_pos}")
                    continue
                chrom, pos = chrom_pos.split(".")
                start = int(pos)
                end = start

                gtf_row = [
                    chrom,
                    "none",
                    "transcript",
                    start,
                    end,
                    ".",
                    "+",
                    ".",
                    f'gene_id "{chrom}_{label}";'
                ]
                gtf_rows.append(gtf_row)
    else:
        print(f"❌ File not found: {fname}")

# 写出 GTF 文件
if gtf_rows:
    gtf_df = pd.DataFrame(gtf_rows)
    gtf_df.to_csv("WK_intermethunmeth.gtf", sep="\t", index=False, header=False, quoting=3)
    print("✅ GTF 文件已生成：WK_intermethunmeth.gtf")
else:
    print("⚠️ 没有生成任何 GTF 内容，请检查输入文件")


perl Variance-at-position.pl --pool-size 198 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/8_WK.pileup \
        --gtf /work/cyu/meth/Wik/WK_intermethunmeth.gtf \
        --output /work/cyu/meth/Wik/WK_intermethunmeth.pi.txt \
        --measure pi


#wik
import pandas as pd
import gzip
import os

# 输入文件及其对应 methylation 类型标签
methylation_files = {
    "WGBS_WK_Intersection_Inter.tsv.gz": "inter",
    "WGBS_WK_Intersection_Meth.tsv.gz": "meth",
    "WGBS_WK_Intersection_Unmeth.tsv.gz": "unmeth"
}

gtf_rows = []

for fname, label in methylation_files.items():
    if os.path.exists(fname):
        print(f"Reading {fname}...")
        with gzip.open(fname, "rt") as f:
            df = pd.read_csv(f, sep="\t", header=0)
            for idx, row in df.iterrows():
                chrom_pos = row["chr"]  # e.g. "chrI.1001661"
                if "." not in chrom_pos:
                    print(f"⚠️ Skipping malformed row: {chrom_pos}")
                    continue
                chrom, pos = chrom_pos.split(".")
                start = int(pos)
                end = start

                gtf_row = [
                    chrom,
                    "none",
                    "transcript",
                    start,
                    end,
                    ".",
                    "+",
                    ".",
                    f'gene_id "{chrom}_{label}";'
                ]
                gtf_rows.append(gtf_row)
    else:
        print(f"❌ File not found: {fname}")

# 写入 GTF 文件
if gtf_rows:
    gtf_df = pd.DataFrame(gtf_rows)
    gtf_df.to_csv("WK_intermethunmeth_WGBS.gtf", sep="\t", index=False, header=False, quoting=3)
    print("✅ GTF 文件已生成：WK_intermethunmeth_WGBS.gtf")
else:
    print("⚠️ 没有生成任何 GTF 内容，请检查输入文件")


perl Variance-at-position.pl --pool-size 198 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/8_WK.pileup \
        --gtf /work/cyu/meth/Wik/WK_intermethunmeth_WGBS.gtf \
        --output /work/cyu/meth/Wik/WK_intermethunmeth.pi.txt \
        --measure pi


ackground pi (WK): 0.00405137



#watson
import pandas as pd
import gzip
import os

# 文件名与甲基化类别标签
methylation_files = {
    "WGBS_WT_Union_Inter.tsv.gz": "inter",
    "WGBS_WT_Union_Meth.tsv.gz": "meth",
    "WGBS_WT_Union_Unmeth.tsv.gz": "unmeth"
}

gtf_rows = []

for fname, label in methylation_files.items():
    if os.path.exists(fname):
        print(f"Reading {fname}...")
        with gzip.open(fname, "rt") as f:
            df = pd.read_csv(f, sep="\t", header=0)
            for idx, row in df.iterrows():
                chrom_pos = row["chr"]  # e.g. "chrI.10000771"
                if "." not in chrom_pos:
                    print(f"⚠️ Skipping malformed row: {chrom_pos}")
                    continue
                chrom, pos = chrom_pos.split(".")
                start = int(pos)
                end = start

                gtf_row = [
                    chrom,
                    "none",
                    "transcript",
                    start,
                    end,
                    ".",
                    "+",
                    ".",
                    f'gene_id "{chrom}_{label}";'
                ]
                gtf_rows.append(gtf_row)
    else:
        print(f"❌ File not found: {fname}")

# 写出 GTF 文件
if gtf_rows:
    gtf_df = pd.DataFrame(gtf_rows)
    gtf_df.to_csv("WT_intermethunmeth_WGBS.gtf", sep="\t", index=False, header=False, quoting=3)
    print("✅ GTF 文件已生成：WT_intermethunmeth_WGBS.gtf")
else:
    print("⚠️ 没有生成任何 GTF 内容，请检查输入文件")

perl Variance-at-position.pl --pool-size 184 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/7_WT.pileup \
        --gtf /work/cyu/meth/Watson/WT_intermethunmeth_WGBS.gtf \
        --output /work/cyu/meth/Watson/WT_intermethunmeth.pi.txt \
        --measure pi

Background pi (WT): 0.00468986