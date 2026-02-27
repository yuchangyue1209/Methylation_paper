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
        --gtf /work/cyu/meth/Roberts/ReducedSamples/Rob_intermethunmeth.gtf \
        --output /work/cyu/meth/Roberts/ReducedSamples/Rob_intermethunmeth.pi.txt \
        --measure pi

Background pi (ROB): 0.00473019


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
        --gtf /work/cyu/meth/Gosling/ReducedSamples/Gos_intermethunmeth.gtf \
        --output /work/cyu/meth/Gosling/ReducedSamples/Gos_intermethunmeth.pi.txt \
        --measure pi

Background pi (GOS): 0.00316497









#wik
import pandas as pd
import gzip
import os

# 输入文件及对应的甲基化类型
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
    gtf_df.to_csv("WGBS_WK_intermethunmeth.gtf", sep="\t", index=False, header=False, quoting=3)
    print("✅ GTF 文件已生成：WGBS_WK_intermethunmeth.gtf")
else:
    print("⚠️ 没有生成任何 GTF 内容，请检查输入文件")


perl Variance-at-position.pl --pool-size 198 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/8_WK.pileup \
        --gtf /work/cyu/meth/Wik/WGBS_WK_intermethunmeth.gtf \
        --output /work/cyu/meth/Wik/WGBS_WK_intermethunmeth.pi.txt \
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
        --gtf /work/cyu/meth/Wik/WK_intermethunmeth.gtf \
        --output /work/cyu/meth/Wik/WGBS_WK_intermethunmeth.pi.txt \
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
        --output /work/cyu/meth/Watson/WGBS_WT_intermethunmeth.pi.txt \
        --measure pi

Background pi (WT): 0.00468986


#sayward
import pandas as pd
import gzip
import os

# 文件名与甲基化类别标签
methylation_files = {
    "Say_Intersection_Inter.tsv.gz": "inter",
    "Say_Intersection_Meth.tsv.gz": "meth",
    "Say_Intersection_Unmeth.tsv.gz": "unmeth"
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
    gtf_df.to_csv("SAY_intermethunmeth.gtf", sep="\t", index=False, header=False, quoting=3)
    print("✅ GTF 文件已生成：SAY_intermethunmeth.gtf")
else:
    print("⚠️ 没有生成任何 GTF 内容，请检查输入文件")

perl Variance-at-position.pl --pool-size 200 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/17_SAY.pileup \
        --gtf /work/cyu/meth/Sayward/ReducedSamples/SAY_intermethunmeth.gtf \
        --output /work/cyu/meth/Sayward/ReducedSamples/SAY_intermethunmeth.pi.txt \
        --measure pi

Background pi (SAY): 0.00514214













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
  theme(
    axis.title = element_text(size = 16),  # 坐标轴标签字体加大
    axis.text = element_text(size = 14),   # 坐标轴刻度字体加大
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.position = "right"
    # 不加 plot.title = element_blank()，保留标题并使用默认字体
  )

# 添加显著性标注
p <- p + geom_signif(comparisons = comparisons,
                     annotations = format.pval(wt_results$p.value, digits = 2, eps = 0.001),
                     y_position = seq(from = max_val * 1.05, length.out = length(comparisons), by = 0.002),
                     tip_length = 0.02, textsize = 4)

# 显示图形
print(p)

# 保存图像
ggsave("RS_bigtext_methylation_pi_boxplot.png", p, width = 8, height = 6, dpi = 300)








#wgbs wik
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
        --gtf /work/cyu/meth/Wik/WK_intermethunmeth.gtf \
        --output /work/cyu/meth/Wik/WGBS_WK_intermethunmeth.pi.txt \
        --measure pi


ackground pi (WK): 0.00405137



#wgbs wt

perl Variance-at-position.pl --pool-size 184 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/7_WT.pileup \
        --gtf /work/cyu/meth/Watson/WT_intermethunmeth_WGBS.gtf \
        --output /work/cyu/meth/Watson/WGBS_WT_intermethunmeth.pi.txt \
        --measure pi

Background pi (WT): 0.00468986


# ---- reusable boxplot function for π by methylation category ----
library(tidyr)
library(ggplot2)
library(dplyr)

# 可选：显著性标注包；没有也能跑
has_signif <- requireNamespace("ggsignif", quietly = TRUE)
if (has_signif) library(ggsignif)

plot_meth_pi <- function(infile, outfile, bg_pi, title_text) {
  # 1) 读取
  dat <- read.table(infile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(dat) <- c("chr_result", "col2", "col3", "statValue")

  # 2) 预处理
  dat <- separate(dat, col = "chr_result", into = c("chr","result"),
                  sep = "_", remove = TRUE)
  dat <- dplyr::filter(dat, !(chr %in% c("chrUn","chrM")))
  dat$result  <- factor(dat$result, levels = c("unmeth","inter","meth"))
  dat$sex_chr <- ifelse(dat$chr == "chrY", "sex", "autosome")
  dat$statValue[dat$statValue == "NaN"] <- NA
  dat$statValue <- as.numeric(dat$statValue)
  dat <- dat[!is.na(dat$statValue), ]

  # 3) Wilcoxon 两两比较
  comps <- combn(levels(dat$result), 2, simplify = FALSE)
  wtres <- lapply(comps, function(comp){
    x <- dplyr::filter(dat, result==comp[1]) |> dplyr::pull(statValue)
    y <- dplyr::filter(dat, result==comp[2]) |> dplyr::pull(statValue)
    w <- wilcox.test(x, y, paired = FALSE, exact = FALSE)
    data.frame(group1=comp[1], group2=comp[2], W=unname(w$statistic), p.value=w$p.value)
  }) |> dplyr::bind_rows()
  print(wtres)

  # 4) 作图
  max_val <- max(dat$statValue, na.rm = TRUE)

  p <- ggplot(dat, aes(x = result, y = statValue, fill = result)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.8, size = 2,
                aes(shape = sex_chr, color = sex_chr)) +
    geom_hline(yintercept = bg_pi, linetype = "dotted", color = "black", size = 0.6) +
    scale_shape_manual(values = c("autosome" = 21, "sex" = 24)) +
    scale_color_manual(values = c("autosome" = "black", "sex" = "red")) +
    theme_bw() +
    labs(x = "Methylation Category", y = expression(pi),
         title = title_text) +
    theme(legend.position = "right")

  if (has_signif) {
    p <- p + geom_signif(
      comparisons = comps,
      annotations = format.pval(wtres$p.value, digits = 2, eps = 0.001),
      y_position = seq(from = max_val * 1.05,
                       length.out = length(comps),
                       by = max_val * 0.03),
      tip_length = 0.02, textsize = 4
    )
  } else {
    # 无 ggsignif：手动画横线与 p 值
    levs  <- levels(dat$result)
    pos_y <- seq(from = max_val * 1.05,
                 length.out = nrow(wtres),
                 by = max_val * 0.03)
    for (i in seq_len(nrow(wtres))) {
      g1 <- wtres$group1[i]; g2 <- wtres$group2[i]
      x1 <- which(levs == g1); x2 <- which(levs == g2)
      y  <- pos_y[i]
      lab <- format.pval(wtres$p.value[i], digits = 2, eps = 0.001)
      p <- p +
        geom_segment(aes(x = x1, xend = x2, y = y, yend = y)) +
        geom_segment(aes(x = x1, xend = x1, y = y, yend = y - max_val*0.01)) +
        geom_segment(aes(x = x2, xend = x2, y = y, yend = y - max_val*0.01)) +
        annotate("text", x = (x1 + x2)/2, y = y + max_val*0.01, label = lab, size = 4)
    }
  }

  print(p)
  ggsave(filename = outfile, plot = p, width = 8, height = 6, dpi = 300)
  message("Saved: ", outfile)
}

# =========================
# 使用示例（与你WT脚本完全等价）
plot_meth_pi("/work/cyu/meth/Wik/WGBS_WK_intermethunmeth.pi.txt",
              "/work/cyu/meth/Wik/WK_methylation_pi_boxplot.png",
             0.00405137,
            "Nucleotide Diversity by Methylation State (WK)")
# plot_meth_pi("/work/cyu/meth/Sayward/ReducedSamples/SAY_intermethunmeth.pi.txt",
#              "/work/cyu/meth/Sayward/ReducedSamples/SAY_methylation_pi_boxplot.png",
#              0.00514214,
#              "Nucleotide Diversity by Methylation State (SAY)")
# plot_meth_pi("RS_intermethunmeth.pi.txt",
#              "RS_methylation_pi_boxplot.png",
#              0.00483154,
#              "Nucleotide Diversity by Methylation State (RS)")

