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
        --output  RS_intermethunmeth.pi.txt \
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









#say for publishing
suppressPackageStartupMessages({
  library(tidyr)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

# =========================
# 1. read pi data for boxplot
# =========================
all_data <- read.table(
  "/work/cyu/meth/Sayward/ReducedSamples/SAY_intermethunmeth.pi.txt",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)

colnames(all_data) <- c("chr_result", "col2", "col3", "statValue")

all_data <- separate(
  all_data,
  col = "chr_result",
  into = c("chr", "result"),
  sep = "_",
  remove = TRUE
)

all_data <- all_data %>%
  filter(!(chr %in% c("chrUn", "chrM")))

all_data$result <- factor(all_data$result, levels = c("unmeth", "inter", "meth"))
all_data$sex_chr <- ifelse(all_data$chr == "chrY", "sex", "autosome")

all_data$statValue[all_data$statValue == "NaN"] <- NA
all_data$statValue <- as.numeric(all_data$statValue)
all_data <- all_data %>% filter(!is.na(statValue))

# =========================
# 2. count CpG sites used in each category
#    (print only, not shown on plot)
# =========================
count_valid_cpg <- function(file, label) {
  dt <- fread(cmd = paste("zcat", shQuote(file)))
  dt[, chr_only := sub("\\..*$", "", chr)]
  dt <- dt[!chr_only %in% c("chrUn", "chrM")]
  
  data.frame(
    result = label,
    n_cpg = nrow(dt)
  )
}

count_df <- bind_rows(
  count_valid_cpg("/work/cyu/meth/Sayward/ReducedSamples/Say_Intersection_Unmeth.tsv.gz", "unmeth"),
  count_valid_cpg("/work/cyu/meth/Sayward/ReducedSamples/Say_Intersection_Inter.tsv.gz",  "inter"),
  count_valid_cpg("/work/cyu/meth/Sayward/ReducedSamples/Say_Intersection_Meth.tsv.gz",   "meth")
)

count_df$result <- factor(count_df$result, levels = c("unmeth", "inter", "meth"))

cat("\nCpG site counts used for each methylation category:\n")
print(count_df)

# =========================
# 3. pairwise Wilcoxon tests
# =========================
comparisons <- list(
  c("unmeth", "inter"),
  c("unmeth", "meth"),
  c("inter", "meth")
)

wt_results <- lapply(comparisons, function(comp) {
  x <- all_data %>% filter(result == comp[1]) %>% pull(statValue)
  y <- all_data %>% filter(result == comp[2]) %>% pull(statValue)
  wtest <- wilcox.test(x, y, paired = FALSE, exact = FALSE)
  
  data.frame(
    group1 = comp[1],
    group2 = comp[2],
    W = unname(wtest$statistic),
    p.value = wtest$p.value
  )
}) %>% bind_rows()

wt_results$label <- ifelse(
  wt_results$p.value < 0.001,
  "p < 0.001",
  paste0("p = ", signif(wt_results$p.value, 3))
)

cat("\nWilcoxon test results:\n")
print(wt_results)

# =========================
# 4. y positions
# =========================
y_min <- min(all_data$statValue, na.rm = TRUE)
y_max <- max(all_data$statValue, na.rm = TRUE)
y_rng <- y_max - y_min
if (y_rng == 0) y_rng <- 0.001

sig_df <- data.frame(
  xmin  = c(1, 1, 2),
  xmax  = c(2, 3, 3),
  y     = c(
    y_max + 0.2 * y_rng,
    y_max + 0.3 * y_rng,
    y_max + 0.4 * y_rng
  ),
  label = wt_results$label
)

upper_limit <- y_max + 0.50 * y_rng

# =========================
# 5. plot
# =========================
bg_pi <- 0.00514214

p <- ggplot(all_data, aes(x = result, y = statValue, fill = result)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.75, linewidth = 0.9) +
  geom_jitter(
    width = 0.15, alpha = 0.85, size = 3,
    aes(shape = sex_chr, color = sex_chr)
  ) +
  geom_hline(
    yintercept = bg_pi,
    linetype = "dotted",
    color = "black",
    linewidth = 0.8
  ) +
  scale_shape_manual(values = c("autosome" = 21, "sex" = 24)) +
  scale_color_manual(values = c("autosome" = "black", "sex" = "red")) +
  scale_fill_manual(values = c(
    "unmeth" = "#F8766D",
    "inter"  = "#00BA38",
    "meth"   = "#619CFF"
  )) +
  coord_cartesian(ylim = c(y_min, upper_limit), clip = "off") +
  theme_bw() +
  labs(
    title = "Sayward-Ancestral Marine",
    x = "Methylation Category",
    y = expression(pi),
    fill = "Result",
    shape = "Chr",
    color = "Chr"
  ) +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 17),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position = "right",
    plot.margin = margin(15, 40, 15, 15)
  )

# =========================
# 6. add significance brackets
# =========================
for (i in 1:nrow(sig_df)) {
  x1  <- sig_df$xmin[i]
  x2  <- sig_df$xmax[i]
  y   <- sig_df$y[i]
  lab <- sig_df$label[i]
  
  p <- p +
    annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.7) +
    annotate("segment", x = x1, xend = x1, y = y, yend = y - 0.015 * y_rng, linewidth = 0.7) +
    annotate("segment", x = x2, xend = x2, y = y, yend = y - 0.015 * y_rng, linewidth = 0.7) +
    annotate("text", x = (x1 + x2) / 2, y = y + 0.04 * y_rng, label = lab, size = 6)
}

print(p)

ggsave(
  "/work/cyu/meth/Sayward/ReducedSamples/SAY_methylation_pi_boxplot_clean.png",
  p, width = 10, height = 8, dpi = 300
)
cat(
  "\nCpG sites used for π calculation:\n",
  "unmethylated =", count_df$n_cpg[count_df$result=="unmeth"], "\n",
  "intermediate =", count_df$n_cpg[count_df$result=="inter"], "\n",
  "methylated   =", count_df$n_cpg[count_df$result=="meth"], "\n"
)










suppressPackageStartupMessages({
  library(tidyr)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

# =========================================================
# Global settings for publication-quality consistency
# =========================================================
Y_LIM <- c(0, 0.013)   # unified y-axis across all populations
FIG_W <- 8
FIG_H <- 8
FIG_DPI <- 300

COL_FILL <- c(
  "unmeth" = "#F8766D",
  "inter"  = "#00BA38",
  "meth"   = "#619CFF"
)

COL_SHAPE <- c("autosome" = 21, "sex" = 24)
COL_OUTLINE <- c("autosome" = "black", "sex" = "red")

# =========================================================
# Population metadata
# =========================================================
pop_configs <- list(
  list(
    pop = "SAY",
    title = "Sayward-Ancestral Marine",
    pi_file = "/work/cyu/meth/Sayward/ReducedSamples/SAY_intermethunmeth.pi.txt",
    unmeth_file = "/work/cyu/meth/Sayward/ReducedSamples/Say_Intersection_Unmeth.tsv.gz",
    inter_file  = "/work/cyu/meth/Sayward/ReducedSamples/Say_Intersection_Inter.tsv.gz",
    meth_file   = "/work/cyu/meth/Sayward/ReducedSamples/Say_Intersection_Meth.tsv.gz",
    bg_pi = 0.00514214,
    out_png = "/work/cyu/meth/Sayward/ReducedSamples/SAY_methylation_pi_boxplot_pub.png"
  ),
  list(
    pop = "RS",
    title = "Rabbit Slough-Ancestral Marine",
    pi_file = "/work/cyu/meth/RabbitSlough/RS_intermethunmeth.pi.txt",
    unmeth_file = "/work/cyu/meth/RabbitSlough/RS_Intersection_Unmeth.tsv.gz",
    inter_file  = "/work/cyu/meth/RabbitSlough/RS_Intersection_Inter.tsv.gz",
    meth_file   = "/work/cyu/meth/RabbitSlough/RS_Intersection_Meth.tsv.gz",
    bg_pi = 0.00483154,
    out_png = "/work/cyu/meth/RabbitSlough/RS_methylation_pi_boxplot_pub.png"
  ),
  list(
    pop = "ROB",
    title = "Roberts-Freshwater",
    pi_file = "/work/cyu/meth/Roberts/ReducedSamples/Rob_intermethunmeth.pi.txt",
    unmeth_file = "/work/cyu/meth/Roberts/ReducedSamples/ROB_Intersection_Unmeth.tsv.gz",
    inter_file  = "/work/cyu/meth/Roberts/ReducedSamples/Intersection_ROB_Inter.tsv.gz",
    meth_file   = "/work/cyu/meth/Roberts/ReducedSamples/Intersection_ROB_Meth.tsv.gz",
    bg_pi = 0.00473019,
    out_png = "/work/cyu/meth/Roberts/ReducedSamples/ROB_methylation_pi_boxplot_pub.png"
  ),
  list(
    pop = "GOS",
    title = "Gosling-Freshwater",
    pi_file = "/work/cyu/meth/Gosling/ReducedSamples/Gos_intermethunmeth.pi.txt",
    unmeth_file = "/work/cyu/meth/Gosling/ReducedSamples/Gos_Intersection_Unmeth.tsv.gz",
    inter_file  = "/work/cyu/meth/Gosling/ReducedSamples/Gos_Intersection_Inter.tsv.gz",
    meth_file   = "/work/cyu/meth/Gosling/ReducedSamples/Gos_Intersection_Meth.tsv.gz",
    bg_pi = 0.00316497,
    out_png = "/work/cyu/meth/Gosling/ReducedSamples/GOS_methylation_pi_boxplot_pub.png"
  ),
  list(
    pop = "WK",
    title = "Wik-Freshwater",
    pi_file = "/work/cyu/meth/Wik/WGBS_WK_intermethunmeth.pi.txt",
    unmeth_file = "/work/cyu/meth/Wik/WGBS_WK_Intersection_Unmeth.tsv.gz",
    inter_file  = "/work/cyu/meth/Wik/WGBS_WK_Intersection_Inter.tsv.gz",
    meth_file   = "/work/cyu/meth/Wik/WGBS_WK_Intersection_Meth.tsv.gz",
    bg_pi = 0.00405137,
    out_png = "/work/cyu/meth/Wik/WK_methylation_pi_boxplot_pub.png"
  ),
  list(
    pop = "WT",
    title = "Watson-Freshwater",
    pi_file = "/work/cyu/meth/Watson/WGBS_WT_intermethunmeth.pi.txt",
    unmeth_file = "/work/cyu/meth/Watson/WGBS_WT_Union_Unmeth.tsv.gz",
    inter_file  = "/work/cyu/meth/Watson/WGBS_WT_Union_Inter.tsv.gz",
    meth_file   = "/work/cyu/meth/Watson/WGBS_WT_Union_Meth.tsv.gz",
    bg_pi = 0.00468986,
    out_png = "/work/cyu/meth/Watson/WT_methylation_pi_boxplot_pub.png"
  )
)

# =========================================================
# Helper functions
# =========================================================
count_valid_cpg <- function(file, label) {
  dt <- fread(cmd = paste("zcat", shQuote(file)))
  dt[, chr_only := sub("\\..*$", "", chr)]
  dt <- dt[!chr_only %in% c("chrUn", "chrM")]
  
  data.frame(
    result = label,
    n_cpg = nrow(dt)
  )
}

format_p_label <- function(p) {
  ifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 3)))
}

# =========================================================
# Main plotting function
# =========================================================
plot_one_population <- function(cfg) {
  cat("\n==============================\n")
  cat("Processing:", cfg$pop, "\n")
  cat("==============================\n")
  
  # -------------------------
  # 1. Read pi data
  # -------------------------
  dat <- read.table(
    cfg$pi_file,
    sep = "\t", header = FALSE, stringsAsFactors = FALSE
  )
  
  colnames(dat) <- c("chr_result", "col2", "col3", "statValue")
  
  dat <- separate(
    dat,
    col = "chr_result",
    into = c("chr", "result"),
    sep = "_",
    remove = TRUE
  )
  
  dat <- dat %>%
    filter(!(chr %in% c("chrUn", "chrM")))
  
  dat$result <- factor(dat$result, levels = c("unmeth", "inter", "meth"))
  dat$sex_chr <- ifelse(dat$chr == "chrY", "sex", "autosome")
  
  dat$statValue[dat$statValue == "NaN"] <- NA
  dat$statValue <- as.numeric(dat$statValue)
  dat <- dat %>% filter(!is.na(statValue))
  
  # -------------------------
  # 2. Count CpG sites
  # -------------------------
  count_df <- bind_rows(
    count_valid_cpg(cfg$unmeth_file, "unmeth"),
    count_valid_cpg(cfg$inter_file,  "inter"),
    count_valid_cpg(cfg$meth_file,   "meth")
  )
  
  count_df$result <- factor(count_df$result, levels = c("unmeth", "inter", "meth"))
  
  cat("\nCpG sites used for π calculation:\n")
  print(count_df)
  
  cat(
    sprintf(
      "\nCpG sites used for π calculation: unmeth=%s, inter=%s, meth=%s\n",
      format(count_df$n_cpg[count_df$result=="unmeth"], big.mark=","),
      format(count_df$n_cpg[count_df$result=="inter"], big.mark=","),
      format(count_df$n_cpg[count_df$result=="meth"], big.mark=",")
    )
  )
  
  # -------------------------
  # 3. Pairwise Wilcoxon tests
  # -------------------------
  comparisons <- list(
    c("unmeth", "inter"),
    c("unmeth", "meth"),
    c("inter", "meth")
  )
  
  wt_results <- lapply(comparisons, function(comp) {
    x <- dat %>% filter(result == comp[1]) %>% pull(statValue)
    y <- dat %>% filter(result == comp[2]) %>% pull(statValue)
    wtest <- wilcox.test(x, y, paired = FALSE, exact = FALSE)
    
    data.frame(
      pop = cfg$pop,
      group1 = comp[1],
      group2 = comp[2],
      W = unname(wtest$statistic),
      p.value = wtest$p.value,
      label = format_p_label(wtest$p.value)
    )
  }) %>% bind_rows()
  
  cat("\nWilcoxon test results:\n")
  print(wt_results)
  
  # -------------------------
  # 4. Fixed-from-top bracket positions
  # -------------------------
  y_top <- Y_LIM[2]
  y_rng <- diff(Y_LIM)
  
  sig_df <- data.frame(
    xmin = c(1, 1, 2),
    xmax = c(2, 3, 3),
    y = c(
      y_top - 0.15 * y_rng,
      y_top - 0.08 * y_rng,
      y_top - 0.01 * y_rng
    ),
    label = wt_results$label
  )
  
  # -------------------------
  # 5. Plot
  # -------------------------
  p <- ggplot(dat, aes(x = result, y = statValue, fill = result)) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.65,
      alpha = 0.80,
      linewidth = 1.0
    ) +
    geom_jitter(
      width = 0.12,
      alpha = 0.90,
      size = 3.2,
      aes(shape = sex_chr, color = sex_chr)
    ) +
    geom_hline(
      yintercept = cfg$bg_pi,
      linetype = "dotted",
      color = "black",
      linewidth = 0.9
    ) +
    scale_shape_manual(values = COL_SHAPE) +
    scale_color_manual(values = COL_OUTLINE) +
    scale_fill_manual(values = COL_FILL) +
    coord_cartesian(ylim = Y_LIM, clip = "off") +
    theme_bw() +
    labs(
      title = cfg$title,
      x = "Methylation Category",
      y = expression(pi),
      fill = "Result",
      shape = "Chr",
      color = "Chr"
    ) +
    theme(
      plot.title = element_text(size = 24, face = "bold"),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 17),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.position = "right",
      panel.grid.major = element_line(linewidth = 0.35),
      panel.grid.minor = element_line(linewidth = 0.20),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  # -------------------------
  # 6. Add fixed-from-top brackets
  # -------------------------
  for (i in 1:nrow(sig_df)) {
    x1  <- sig_df$xmin[i]
    x2  <- sig_df$xmax[i]
    y   <- sig_df$y[i]
    lab <- sig_df$label[i]
    
    p <- p +
      annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.7) +
      annotate("segment", x = x1, xend = x1, y = y, yend = y - 0.012 * y_rng, linewidth = 0.7) +
      annotate("segment", x = x2, xend = x2, y = y, yend = y - 0.012 * y_rng, linewidth = 0.7) +
      annotate("text", x = (x1 + x2) / 2, y = y + 0.03 * y_rng, label = lab, size = 4.5)
  }
  
  print(p)
  
  ggsave(
    filename = cfg$out_png,
    plot = p,
    width = FIG_W,
    height = FIG_H,
    dpi = FIG_DPI
  )
  
  cat("\nSaved plot to:\n", cfg$out_png, "\n")
  
  list(
    counts = count_df %>% mutate(pop = cfg$pop, .before = 1),
    wilcox = wt_results,
    plot = p
  )
}

# =========================================================
# Run all populations
# =========================================================
all_counts <- list()
all_wilcox <- list()

for (cfg in pop_configs) {
  res <- plot_one_population(cfg)
  all_counts[[cfg$pop]] <- res$counts
  all_wilcox[[cfg$pop]] <- res$wilcox
}

all_counts_df <- bind_rows(all_counts)
all_wilcox_df <- bind_rows(all_wilcox)

# =========================================================
# Save summary tables
# =========================================================
fwrite(
  all_counts_df,
  "/work/cyu/meth/all_populations_methylation_pi_CpG_counts.tsv",
  sep = "\t"
)

fwrite(
  all_wilcox_df,
  "/work/cyu/meth/all_populations_methylation_pi_wilcox.tsv",
  sep = "\t"
)

cat("\n========================================\n")
cat("Finished all populations.\n")
cat("Saved summary tables:\n")
cat("/work/cyu/meth/all_populations_methylation_pi_CpG_counts.tsv\n")
cat("/work/cyu/meth/all_populations_methylation_pi_wilcox.tsv\n")
cat("========================================\n")





#may 14 redo fig1

mkdir -p /work/shared/cyu/meth/NoSNPs_methyl/gtf
cd /work/shared/cyu/meth/NoSNPs_methyl

cat > make_all_gtf.py <<'PY'
import pandas as pd
import gzip
import os

BASE_DIR = "/work/shared/cyu/meth/NoSNPs_methyl"
OUT_DIR  = f"{BASE_DIR}/gtf"

os.makedirs(OUT_DIR, exist_ok=True)

pop_configs = {
    "RS": {
        "inter": "RS_Intersection_Inter.tsv.gz",
        "meth":  "RS_Intersection_Meth.tsv.gz",
        "unmeth":"RS_Intersection_Unmeth.tsv.gz"
    },

    "ROB": {
        "inter": "Intersection_ROB_Inter.tsv.gz",
        "meth":  "Intersection_ROB_Meth.tsv.gz",
        "unmeth":"ROB_Intersection_Unmeth.tsv.gz"
    },

    "GOS": {
        "inter": "Gos_Intersection_Inter.tsv.gz",
        "meth":  "Gos_Intersection_Meth.tsv.gz",
        "unmeth":"Gos_Intersection_Unmeth.tsv.gz"
    },

    "SAY": {
        "inter": "Say_Intersection_Inter.tsv.gz",
        "meth":  "Say_Intersection_Meth.tsv.gz",
        "unmeth":"Say_Intersection_Unmeth.tsv.gz"
    },

    "WK": {
        "inter": "WGBS_WK_Intersection_Inter.tsv.gz",
        "meth":  "WGBS_WK_Intersection_Meth.tsv.gz",
        "unmeth":"WGBS_WK_Intersection_Unmeth.tsv.gz"
    },

    "WT": {
        "inter": "WGBS_WT_Intersection_Inter.tsv.gz",
        "meth":  "WGBS_WT_Intersection_Meth.tsv.gz",
        "unmeth":"WGBS_WT_Intersection_Unmeth.tsv.gz"
    }
}

for pop, files in pop_configs.items():

    print(f"\n======================")
    print(f"Processing {pop}")
    print(f"======================")

    gtf_rows = []

    for label, fname in files.items():

        fullpath = f"{BASE_DIR}/{fname}"

        print(f"Reading {fullpath}")

        if not os.path.exists(fullpath):
            print(f"ERROR: missing file {fullpath}")
            continue

        with gzip.open(fullpath, "rt") as f:
            df = pd.read_csv(f, sep="\t", header=0)

        if "chr" not in df.columns:
            raise ValueError(f"'chr' column missing in {fname}")

        n_before = len(df)

        for _, row in df.iterrows():

            chrom_pos = str(row["chr"])

            if "." not in chrom_pos:
                continue

            chrom, pos = chrom_pos.split(".", 1)

            try:
                pos = int(pos)
            except:
                continue

            gtf_rows.append([
                chrom,
                "none",
                "transcript",
                pos,
                pos,
                ".",
                "+",
                ".",
                f'gene_id "{chrom}_{label}";'
            ])

        print(f"Loaded {n_before:,} rows")

    out_gtf = f"{OUT_DIR}/{pop}_intermethunmeth.gtf"

    gtf_df = pd.DataFrame(gtf_rows)

    gtf_df.to_csv(
        out_gtf,
        sep="\t",
        index=False,
        header=False,
        quoting=3
    )

    print(f"\nWrote: {out_gtf}")
    print(f"Total GTF rows: {len(gtf_rows):,}")

print("\nFinished all populations.")
PY

python3 make_all_gtf.py



