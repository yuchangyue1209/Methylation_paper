#run bwa-meth
#!/bin/bash

for sample_id in /work/ekerns/trimmed_reads/Leder/*_val_1.fq.gz
 do
        base_id=$(basename $sample_id "_val_1.fq.gz")
        /home/cyu/bwa-meth/bwameth.py --reference /work/ekerns/LederData/stickleback_v5_assembly.fa ${base_id}_val_1.fq.gz ${base_id}_val_2.fq.gz | samtools view -b | samtools sort --output-fmt BAM > ${base_id}.bam
 done


 #!/bin/bash

for fq1 in /work/ekerns/trimmed_reads/Leder/*_val_1.fq.gz
do
    # base_id
    base_id=$(basename "$fq1" "_1_val_1.fq.gz")
    
    # fq2 path
    fq2="/work/ekerns/trimmed_reads/Leder/${base_id}_2_val_2.fq.gz"

    # check fq2 
    if [[ -f "$fq1" && -f "$fq2" ]]; then
        echo "Processing: $base_id"
        /home/cyu/bwa-meth/bwameth.py --reference /work/ekerns/LederData/stickleback_v5_assembly.fa \
        "$fq1" "$fq2" | samtools view -b | samtools sort --output-fmt BAM > "${base_id}.bam"
    else
        echo "Skipping $base_id: Paired file not found!"
    fi
done







#Poolseq data treatment: from trimming to pileup

1. #Trim raw_data with bbduk 38.9

#!/bin/bash

# Set input and output directories and parameters
INPUT_DIR="/work/cyu/poolseq/raw_data"
OUTPUT_DIR="/work/cyu/poolseq/PPalign_output/trimmed"
ADAPTER_REF="/home/cyu/adapters.fa"
THREADS=48

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all paired R1 and R2 files
for R1_FILE in "$INPUT_DIR"/*_R1_001.fastq; do
    # Get the corresponding R2 file
    R2_FILE="${R1_FILE/_R1_/_R2_}"
    
    # Extract the file prefix (e.g., 16_AMO_S16)
    PREFIX=$(basename "$R1_FILE" | sed 's/_R1_001.fastq//')
    
    # Set output file paths
    OUT_R1="$OUTPUT_DIR/trimmed_R1_${PREFIX}.fastq"
    OUT_R2="$OUTPUT_DIR/trimmed_R2_${PREFIX}.fastq"
    LOG_FILE="$OUTPUT_DIR/${PREFIX}_bbduk_log.txt"
    
    # Run bbduk.sh
    bbduk.sh \
        in1="$R1_FILE" \
        in2="$R2_FILE" \
        out1="$OUT_R1" \
        out2="$OUT_R2" \
        ref="$ADAPTER_REF" \
        ktrim=rl \
        trimq=20 \
        minlength=25 \
        ftl=10 \
        tossbrokenreads=t \
        threads="$THREADS" > "$LOG_FILE" 2>&1
done

2. Mapping with bowtie2
#!/bin/bash
# Path
TRIMMED_DIR="/work/cyu/poolseq/PPalign_output/trimmed"
MAPPED_DIR="/work/cyu/meth/pooldata/mapped"
REFERENCE_GZ="/work/cyu/meth/stickleback_v5_assembly.fa.gz"
REFERENCE="/work/cyu/meth/stickleback_v5_assembly.fa"

# output
mkdir -p "$MAPPED_DIR"

# unzip ref
if [[ ! -f "$REFERENCE" ]]; then
    echo "Decompressing reference genome..."
    gunzip -c "$REFERENCE_GZ" > "$REFERENCE"
fi

# index ref
if [[ ! -f "${REFERENCE}.1.bt2" ]]; then
    echo "Reference index not found, creating index..."
    bowtie2-build "$REFERENCE" "$REFERENCE"
fi

# samples input
declare -A SAMPLES=(
    ["25_RS"]="trimmed_R1_25_RS.fastq trimmed_R2_25_RS.fastq"
)

# map
for SAMPLE in "${!SAMPLES[@]}"; do
    read -r R1_FILE R2_FILE <<< "${SAMPLES[$SAMPLE]}"
    
    R1_PATH="$TRIMMED_DIR/$R1_FILE"
    R2_PATH="$TRIMMED_DIR/$R2_FILE"
    
    # check 
    if [[ -f "$R1_PATH" && -f "$R2_PATH" ]]; then
        # output path
        SAM_FILE="$MAPPED_DIR/${SAMPLE}.sam"
        LOG_FILE="$MAPPED_DIR/${SAMPLE}_bowtie2.log"

        # bowtie2 default 
        echo "Mapping $SAMPLE..."
        bowtie2 -x "$REFERENCE" -1 "$R1_PATH" -2 "$R2_PATH" -p 48 -S "$SAM_FILE" > "$LOG_FILE" 2>&1
    else
        echo "Warning: R1 or R2 file for $SAMPLE not found, skipping." | tee -a "$MAPPED_DIR/mapping_warnings.log"
    fi
done

echo "Mapping completed."


3. Convert to bam and sort
#!/bin/bash

# path
SAM_DIR="/work/cyu/meth/pooldata/mapped"
BAM_DIR="/work/cyu/meth/pooldata/bam"

# output
mkdir -p "$BAM_DIR"

for SAM_FILE in "$SAM_DIR"/*.sam; do
    # name
    BASENAME=$(basename "$SAM_FILE" .sam)
    BAM_FILE="$BAM_DIR/${BASENAME}.bam"
    SORTED_BAM_FILE="$BAM_DIR/${BASENAME}_sorted.bam"

    # check
    if [[ -f "$SAM_FILE" ]]; then
        echo "Converting $SAM_FILE to BAM with MAPQ ≥ 20..."
        samtools view -b -q 20 "$SAM_FILE" > "$BAM_FILE"

        echo "Sorting $BAM_FILE..."
        samtools sort -o "$SORTED_BAM_FILE" "$BAM_FILE"
        echo "Finished processing $BASENAME."
    else
        echo "Warning: No SAM files found in $SAM_DIR."
    fi
done

echo "Conversion and sorting completed."





4. dedup sambamba 1.0.1
#!/bin/bash

# Increase the open file limit
ulimit -n 100000

# Set BAM file path
BAM_DIR="/work/cyu/meth/pooldata/bam"

# Iterate through sorted BAM files
for BAM in "$BAM_DIR"/*SL*_sorted.bam; do
    SAMPLE=$(basename "$BAM" _sorted.bam)
    DEDUP_BAM="${BAM_DIR}/${SAMPLE}_sorted_dedup.bam"

    echo "Removing duplicates for ${SAMPLE}..."
    sambamba markdup -r -t 15 "$BAM" "$DEDUP_BAM"

    echo "Finished processing $SAMPLE."
done

echo "Duplicate removal completed."



5. stat alignment numbers 
#!/bin/bash

# Set the directory containing BAM files
BAM_DIR="/work/cyu/meth/pooldata/bam"
OUTPUT_DIR="/work/cyu/meth/pooldata/alignment_stats"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through all *_sorted.bam files
for BAM in "$BAM_DIR"/*_sorted.bam; do
    SAMPLE=$(basename "$BAM" _sorted.bam)
    
    echo "Processing $SAMPLE..."
    
    # Calculate alignment statistics using samtools
    samtools stats "$BAM" > "$OUTPUT_DIR/${SAMPLE}_stats.txt"
    
    echo "Saved alignment stats for $SAMPLE."
done

echo "All alignment stats saved in $OUTPUT_DIR."



6. mileup dedup files 
#!/bin/bash

BAM_DIR="/work/cyu/meth/pooldata/bam"
PILEUP_DIR="/work/cyu/meth/pooldata/pileup"
REFERENCE="/work/cyu/meth/stickleback_v5_assembly.fa"

# Ensure the pileup directory exists
mkdir -p "$PILEUP_DIR"

# Loop through all _sorted_dedup.bam files
for BAM in "$BAM_DIR"/*SL*_sorted_dedup.bam; do
    SAMPLE=$(basename "$BAM" _sorted_dedup.bam)

    echo "Indexing $BAM..."
    samtools index "$BAM"

    echo "Generating pileup for $SAMPLE..."
    samtools mpileup -B -q 20 -Q 20 -f "$REFERENCE" "$BAM" > "$PILEUP_DIR/${SAMPLE}.pileup"

    echo "Finished processing $SAMPLE."
done




#create non hyper hypo fig2
#trans bed files to gtf，prep for popoolation input
import pandas as pd
import csv

# Read the BED files (assuming the file contains a header row)
dmc_sites = pd.read_csv("RSRob_DMCs_minDepth5.bed", sep=" ", header=0, quotechar='"')
all_cpg_sites = pd.read_csv("RSRob_CpGs_minDepth5.bed", sep=" ", header=0, quotechar='"')

# Split the 'chr' column to extract chromosome name and starting position 
# (for entries formatted like "chr1.1000")
dmc_sites[['chr', 'chr_start']] = dmc_sites['chr'].str.split('.', expand=True)
all_cpg_sites[['chr', 'chr_start']] = all_cpg_sites['chr'].str.split('.', expand=True)

# Convert the split starting positions to integers
dmc_sites['start'] = dmc_sites['chr_start'].astype(int)
all_cpg_sites['start'] = all_cpg_sites['chr_start'].astype(int)

# Retain necessary columns and remove duplicate entries
# Note: 'result_GS' is used as the column name because the actual column in the file is named result_GS
dmc_sites_clean = dmc_sites[['chr', 'start', 'result_GS']].drop_duplicates(['chr', 'start'])
all_cpg_sites_clean = all_cpg_sites[['chr', 'start']].drop_duplicates(['chr', 'start'])

# Merge the data and fill in missing values (fill with 'non' if there is no corresponding result_GS)
merged_sites = all_cpg_sites_clean.merge(dmc_sites_clean, on=['chr', 'start'], how='left')
merged_sites['result_GS'] = merged_sites['result_GS'].fillna('non')

# Construct a DataFrame in GTF format. The attribute column is formatted as: gene_id "chr_result_GS";
final_gtf_df = pd.DataFrame({
    'seqname': merged_sites['chr'],
    'source': 'none',
    'feature': 'transcript',
    'start': merged_sites['start'],
    'end': merged_sites['start'],  # Adjust the end value if needed
    'score': '.',
    'strand': '+',
    'frame': '.',
    'attribute': 'gene_id "' + merged_sites['chr'] + '_' + merged_sites['result_GS'] + '";'
})

# Export the DataFrame as a GTF file, ensuring that no quotes are automatically added
final_output_file = "RSRob_DMCs_final_output.gtf"
final_gtf_df.to_csv(final_output_file, sep='\t', header=False, index=False,
                    quoting=csv.QUOTE_NONE, escapechar="\\")

# Print the counts of sites by type to verify the distribution of hypo/hyper/non
print("Counts of sites by type:")
print(merged_sites['result_GS'].value_counts())


RSGos_CpGs_minDepth5
non      127302
hyper      2058
hypo        479

non      112877
hyper      1565
hypo        430


#popoolation for meth

for k in pi theta D; do
    perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/19_ROB.pileup \
        --gtf /work/cyu/meth/RSvsRob/RSRob_DMCs_final_output.gtf \
        --output /work/cyu/meth/RSvsRob/ROB_RSvsROB_${k}_final.txt \
        --measure ${k}
done

for k in pi theta D; do
    perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/RSvsRob/RSRob_DMCs_final_output.gtf \
        --output /work/cyu/meth/RSvsRob/RS_RSvsROB_${k}_final.txt \
        --measure ${k}
done
perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/19_ROB.pileup \
        --gtf /work/cyu/meth/RSvsRob/RSRob_DMCs_final_output.gtf \
        --output /work/cyu/meth/RSvsRob/ROB_RSvsROB_pi_final.txt\
        --measure pi

perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/RSvsRob/RSRob_DMCs_final_output.gtf \
        --output /work/cyu/meth/RSvsRob/RS_RSvsROB_pi_final.txt\
        --measure pi

perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/RSvsGos/RSGos_DMCs_final_output.gtf \
        --output /work/cyu/meth/RSvsGos/RS_RSvsGOS_pi_final.txt \
        --measure pi



perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/18_GOS.pileup \
        --gtf /work/cyu/meth/RSvsGos/RSGos_DMCs_final_output.gtf \
        --output /work/cyu/meth/RSvsGos/GOS_RSvsGOS_pi_final.txt \
        --measure pi




#barplot hypohypernon-RSGOS RSROB

library(tidyr)
library(ggplot2)
library(dplyr)

############################
# 1. Read and merge data
marine <- read.table("RS_RSvsGOS_pi_final.txt", 
                     sep = "\t", header = FALSE, stringsAsFactors = FALSE)
freshwater <- read.table("GOS_RSvsGOS_pi_final.txt", 
                         sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# Add a column to distinguish the populations
marine$Population <- "Marine"
freshwater$Population <- "Freshwater"

# Combine the data
all_data <- rbind(marine, freshwater)

# Rename the first four columns
colnames(all_data)[1:4] <- c("chr_result", "col2", "col3", "statValue")

# Remove rows where col2 and statValue are both 0
all_data <- all_data[!(all_data$col2 == 0 & all_data$statValue == 0), ]

# Split the chr_result column into chr and result
all_data <- separate(all_data, 
                     col = "chr_result", 
                     into = c("chr", "result"), 
                     sep = "_", 
                     remove = TRUE)

# Set the result column as a factor with specified order
all_data$result <- factor(all_data$result, levels = c("non", "hypo", "hyper"))

# Convert statValue to numeric
all_data$statValue <- as.numeric(all_data$statValue)

# Exclude rows with NA in statValue
all_data <- all_data[!is.na(all_data$statValue), ]

############################
# 2. Run the Wilcoxon test
wt_results <- data.frame()  # Data frame to store results
for(cat in c("non", "hypo", "hyper")){
  # Subset Marine group
  x <- subset(all_data, result == cat & Population == "Marine")$statValue
  # Subset Freshwater group
  y <- subset(all_data, result == cat & Population == "Freshwater")$statValue
  
  # Perform unpaired Wilcoxon test
  wtest <- wilcox.test(x, y, paired = FALSE, alternative = "two.sided")
  
  # Save the results
  wt_results <- rbind(wt_results, 
                      data.frame(result = cat,
                                 W = wtest$statistic,
                                 p.value = wtest$p.value))
}

print(wt_results)

############################
# 3. Plot boxplot + scatter plot

p <- ggplot(all_data, aes(x = result, y = statValue, fill = Population)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.7, pch = 21) +
  theme_bw() +
  labs(x = "", 
       y = expression(pi~"(chromosome-level)"),
       title = "RSvsGOS") +
  theme(legend.position = "bottom")

# Save Wilcoxon test p-values
wt_data <- wt_results %>% select(result, p.value)

# Extract p-values
p_non   <- wt_data$p.value[wt_data$result == "non"]
p_hypo  <- wt_data$p.value[wt_data$result == "hypo"]
p_hyper <- wt_data$p.value[wt_data$result == "hyper"]

# Max value for placing p-value labels
max_val <- max(all_data$statValue, na.rm = TRUE)

# Function to format p-values
format_p <- function(p) {
  ifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 3)))
}

# Add annotations
p <- p +
  annotate("text", x = 1, y = max_val * 1.05, label = format_p(p_non), size = 4) +
  annotate("text", x = 2, y = max_val * 1.05, label = format_p(p_hypo), size = 4) +
  annotate("text", x = 3, y = max_val * 1.05, label = format_p(p_hyper), size = 4)

# Display plot
print(p)

# Save plot to file
ggsave("pi_RSvsGOS_boxplot.png", p, width = 8, height = 6, dpi = 300)















#figre3a density
# 加载必要的包
library(dplyr)
library(ggplot2)

# 1. 读取数据文件
# 读取 DMCs 数据
dmc_data <- read.table(gzfile("RSvsGos_PercMeth_DMCs.tsv.gz"),
                       header = TRUE, sep = "\t", check.names = FALSE)

# 读取 AllCpGs 数据
allcpg_data <- read.table(gzfile("RSvsGos_PercMeth_AllCpGs.tsv.gz"),
                          header = TRUE, sep = "\t", check.names = FALSE)

# 2. 直接提取 AllCpGs 数据中的前12列，并手动指定列名
allcpg_data_subset <- allcpg_data[, 1:12]
colnames(allcpg_data_subset) <- c("chr", "start", "end", "strand", "pvalue", "qvalue",
                                  "meth.diff", "result_GS", "mean_Gos", "SD_Gos",
                                  "mean_RS", "SD_RS")

# 3. 转换 join 所需的列类型为字符型
allcpg_data_subset$chr <- as.character(allcpg_data_subset$chr)
dmc_data$chr <- as.character(dmc_data$chr)

# 如果 start、end 也可能存在类型差异，可考虑统一转换为字符或数值，
# 例如这里统一转换为数值：
allcpg_data_subset$start <- as.numeric(allcpg_data_subset$start)
allcpg_data_subset$end   <- as.numeric(allcpg_data_subset$end)
dmc_data$start <- as.numeric(dmc_data$start)
dmc_data$end   <- as.numeric(dmc_data$end)

# 如果 strand 列也有差异，统一转换为字符：
allcpg_data_subset$strand <- as.character(allcpg_data_subset$strand)
dmc_data$strand <- as.character(dmc_data$strand)

# 4. 利用 anti_join 得到 non 位点（即 AllCpGs 中有而 DMCs 中没有的位点）
non_sites <- anti_join(allcpg_data_subset, dmc_data, by = c("chr", "start", "end", "strand"))
non_sites$result_GS <- "non"  # 为 non 位点添加 result_GS 标签

# 5. 针对 DMCs 数据，分别提取 Marine（RS）和淡水（Gos）数据
# Marine（RS）数据
levels_M_dmc <- dmc_data[, c("chr", "start", "strand", "mean_RS", "SD_RS", "result_GS")]
colnames(levels_M_dmc)[4:5] <- c("mean", "SD")
levels_M_dmc$Population <- "Marine"

# 淡水（Gos）数据
levels_F_dmc <- dmc_data[, c("chr", "start", "strand", "mean_Gos", "SD_Gos", "result_GS")]
colnames(levels_F_dmc)[4:5] <- c("mean", "SD")
levels_F_dmc$Population <- "Freshwater"

# 6. 针对 non 数据，从 AllCpGs 数据中得到，分别提取 Marine（RS）和淡水（Gos）数据
levels_M_non <- non_sites[, c("chr", "start", "strand", "mean_RS", "SD_RS", "result_GS")]
colnames(levels_M_non)[4:5] <- c("mean", "SD")
levels_M_non$Population <- "Marine"

levels_F_non <- non_sites[, c("chr", "start", "strand", "mean_Gos", "SD_Gos", "result_GS")]
colnames(levels_F_non)[4:5] <- c("mean", "SD")
levels_F_non$Population <- "Freshwater"

# 7. 合并所有数据
combined_data <- rbind(levels_M_dmc, levels_F_dmc, levels_M_non, levels_F_non)

# 设置因子水平，保证绘图时分面顺序正确
combined_data$result_GS <- factor(combined_data$result_GS, levels = c("non", "hypo", "hyper"))
combined_data$Population <- factor(combined_data$Population, levels = c("Marine", "Freshwater"))

# 设置分面标签（可根据需要自定义）
res.labs <- c("Non-DMC", "FW-Hypo", "FW-Hyper")
names(res.labs) <- c("non", "hypo", "hyper")

# 8. 绘制甲基化百分比密度图
methlevels_hist <- ggplot(combined_data, aes(x = mean, fill = Population, alpha = Population)) +
  theme_classic() +
  geom_density() +
  scale_alpha_manual(values = c(1, 0.5), name = "Population: ") +
  labs(x = "Mean % methylation", y = "Site density", title = "A") +
  facet_wrap(result_GS ~ ., scales = "free_y", drop = FALSE, labeller = labeller(result_GS = res.labs)) +
  theme(strip.background = element_rect(color = "black"),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = "white")) +
  scale_fill_manual(values = c("steelblue4", "palegreen"), name = "Population: ")

# 显示图形
print(methlevels_hist)

# 9. 保存图形为 PNG 文件
ggsave("figure3_methlevels_hist.png", plot = methlevels_hist, width = 8, height = 6, dpi = 300)











# -------------------------------
# 1. 数据读取与预处理
# -------------------------------
library(dplyr)
library(ggplot2)

# 读取 DMCs 数据
dmc_data <- read.table(gzfile("RSvsGos_PercMeth_DMCs.tsv.gz"),
                       header = TRUE, sep = "\t", check.names = FALSE)

# 读取 AllCpGs 数据
allcpg_data <- read.table(gzfile("RSvsGos_PercMeth_AllCpGs.tsv.gz"),
                          header = TRUE, sep = "\t", check.names = FALSE)

# 提取 AllCpGs 数据中的前12列，并重新命名
allcpg_data_subset <- allcpg_data[, 1:12]
colnames(allcpg_data_subset) <- c("chr", "start", "end", "strand", "pvalue", "qvalue",
                                  "meth.diff", "result_GS", "mean_Gos", "SD_Gos",
                                  "mean_RS", "SD_RS")

# 转换关键字段类型（确保两边一致）
allcpg_data_subset$chr <- as.character(allcpg_data_subset$chr)
dmc_data$chr <- as.character(dmc_data$chr)

allcpg_data_subset$start <- as.numeric(allcpg_data_subset$start)
allcpg_data_subset$end   <- as.numeric(allcpg_data_subset$end)
dmc_data$start <- as.numeric(dmc_data$start)
dmc_data$end   <- as.numeric(dmc_data$end)

allcpg_data_subset$strand <- as.character(allcpg_data_subset$strand)
dmc_data$strand <- as.character(dmc_data$strand)

# 得到 non 位点（AllCpGs 中有而 DMCs 中没有的位点）
non_sites <- anti_join(allcpg_data_subset, dmc_data, by = c("chr", "start", "end", "strand"))
non_sites$result_GS <- "non"  # 为 non 位点设置标签

# 分别提取 Marine（RS）和 Freshwater（Gos）数据
# Marine 数据（DMCs）
levels_M_dmc <- dmc_data[, c("chr", "start", "strand", "mean_RS", "SD_RS", "result_GS")]
colnames(levels_M_dmc)[4:5] <- c("mean", "SD")
levels_M_dmc$Population <- "Marine"

# Freshwater 数据（DMCs）
levels_F_dmc <- dmc_data[, c("chr", "start", "strand", "mean_Gos", "SD_Gos", "result_GS")]
colnames(levels_F_dmc)[4:5] <- c("mean", "SD")
levels_F_dmc$Population <- "Freshwater"

# Marine 数据（non）
levels_M_non <- non_sites[, c("chr", "start", "strand", "mean_RS", "SD_RS", "result_GS")]
colnames(levels_M_non)[4:5] <- c("mean", "SD")
levels_M_non$Population <- "Marine"

# Freshwater 数据（non）
levels_F_non <- non_sites[, c("chr", "start", "strand", "mean_Gos", "SD_Gos", "result_GS")]
colnames(levels_F_non)[4:5] <- c("mean", "SD")
levels_F_non$Population <- "Freshwater"

# 合并所有数据
combined_data <- rbind(levels_M_dmc, levels_F_dmc, levels_M_non, levels_F_non)
combined_data$result_GS <- factor(combined_data$result_GS, levels = c("non", "hypo", "hyper"))
combined_data$Population <- factor(combined_data$Population, levels = c("Marine", "Freshwater"))

# 为后续 merge 添加唯一 ID（若原数据已有唯一标识，此步可删除）
combined_data$ID <- paste0("site", 1:nrow(combined_data))


# -------------------------------
# 2. 绘制密度图
# -------------------------------
res.labs <- c("Non-DMC", "FW-Hypo", "FW-Hyper")
names(res.labs) <- c("non", "hypo", "hyper")

methlevels_hist <- ggplot(combined_data, aes(x = mean, fill = Population, alpha = Population)) +
  theme_classic() +
  geom_density() +
  scale_alpha_manual(values = c(1, 0.5), name = "Population: ") +
  labs(x = "Mean % methylation", y = "Site density", title = "A") +
  facet_wrap(result_GS ~ ., scales = "free_y", drop = FALSE,
             labeller = labeller(result_GS = res.labs)) +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = "white")) +
  scale_fill_manual(values = c("steelblue4", "palegreen"), name = "Population: ")

print(methlevels_hist)
ggsave("figure3_methlevels_hist.png", plot = methlevels_hist, width = 8, height = 6, dpi = 300)







#bin 50

#combine tsv
library(dplyr)
library(ggplot2)

###############################
### 1. 生成 combined_data  ###
###############################

# 读取 DMCs 数据
dmc_data <- read.table(gzfile("RSvsGos_PercMeth_DMCs.tsv.gz"),
                       header = TRUE, sep = "\t", check.names = FALSE)

# 读取 AllCpGs 数据
allcpg_data <- read.table(gzfile("RSvsGos_PercMeth_AllCpGs.tsv.gz"),
                          header = TRUE, sep = "\t", check.names = FALSE)

# 提取 AllCpGs 数据中的前12列并重命名
allcpg_data_subset <- allcpg_data[, 1:12]
colnames(allcpg_data_subset) <- c("chr", "start", "end", "strand", "pvalue", "qvalue",
                                  "meth.diff", "result_GS", "mean_Gos", "SD_Gos",
                                  "mean_RS", "SD_RS")

# 转换关键字段类型
allcpg_data_subset$chr <- as.character(allcpg_data_subset$chr)
dmc_data$chr <- as.character(dmc_data$chr)
allcpg_data_subset$start <- as.numeric(allcpg_data_subset$start)
allcpg_data_subset$end   <- as.numeric(allcpg_data_subset$end)
dmc_data$start <- as.numeric(dmc_data$start)
dmc_data$end   <- as.numeric(dmc_data$end)
allcpg_data_subset$strand <- as.character(allcpg_data_subset$strand)
dmc_data$strand <- as.character(dmc_data$strand)

# 得到 non 位点（AllCpGs 有而 DMCs 没有）
non_sites <- anti_join(allcpg_data_subset, dmc_data, by = c("chr", "start", "end", "strand"))
non_sites$result_GS <- "non"

# 分别提取 Marine（RS）和 Freshwater（Gos）数据
# Marine 数据（DMCs）使用 mean_RS, SD_RS
levels_M_dmc <- dmc_data[, c("chr", "start", "end", "strand", "mean_RS", "SD_RS", "result_GS")]
colnames(levels_M_dmc)[5:6] <- c("mean", "SD")
levels_M_dmc$Population <- "Marine"

# Freshwater 数据（DMCs）使用 mean_Gos, SD_Gos
levels_F_dmc <- dmc_data[, c("chr", "start", "end", "strand", "mean_Gos", "SD_Gos", "result_GS")]
colnames(levels_F_dmc)[5:6] <- c("mean", "SD")
levels_F_dmc$Population <- "Freshwater"

# Marine 数据（non）使用 mean_RS, SD_RS
levels_M_non <- non_sites[, c("chr", "start", "end", "strand", "mean_RS", "SD_RS", "result_GS")]
colnames(levels_M_non)[5:6] <- c("mean", "SD")
levels_M_non$Population <- "Marine"

# Freshwater 数据（non）使用 mean_Gos, SD_Gos
levels_F_non <- non_sites[, c("chr", "start", "end", "strand", "mean_Gos", "SD_Gos", "result_GS")]
colnames(levels_F_non)[5:6] <- c("mean", "SD")
levels_F_non$Population <- "Freshwater"

# 合并所有数据
combined_data <- rbind(levels_M_dmc, levels_F_dmc, levels_M_non, levels_F_non)
combined_data$result_GS <- factor(combined_data$result_GS, levels = c("non", "hypo", "hyper"))
combined_data$Population <- factor(combined_data$Population, levels = c("Marine", "Freshwater"))

# 添加唯一 ID（如果原数据中没有唯一标识）
combined_data$ID <- paste0("site", 1:nrow(combined_data))

# 检查生成的数据（选取部分关键字段）
head(combined_data[, c("ID", "chr", "start", "end", "strand", "Population", "result_GS", "mean", "SD")])

#####################################
### 2. 定义辅助函数             ###
#####################################

# bin() 函数：对指定列进行分箱
bin <- function(df, col = "mean", nbins = 50) {
  if (!(col %in% names(df))) {
    stop("Data frame中未找到列：", col)
  }
  bin_col <- paste0("bin_", col)
  df[[bin_col]] <- as.numeric(cut(df[[col]], breaks = nbins))
  return(df)
}

# 写出 GTF 的函数（此处写为制表符分隔文本文件；实际GTF可根据需要调整）
writeGTF_SiteClass <- function(df, output) {
  write.table(df, file = output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# 定义转换函数：将 chr 字段仅保留“.”前的部分，例如 "chrVIII.8665226" -> "chrVIII"
convert_chr <- function(x) {
  sub("\\..*", "", x)
}

##################################################
### 3. 对均值（mean）进行分箱并生成 GTF 文件 ###
##################################################
mean_bins_list <- list()  # 用于存储每个群体的均值分箱结果

for (pop in c("Marine", "Freshwater")) {
  pop_data <- subset(combined_data, Population == pop)
  mean_bins_pop <- NULL
  
  # 针对每个结果类别
  for (res in c("hyper", "hypo", "non")) {
    sub_data <- subset(pop_data, result_GS == res)
    if(nrow(sub_data) == 0) next
    
    # 对 mean 列进行分箱
    mean_data <- bin(sub_data, col = "mean", nbins = 50)
    mean_data <- rename(mean_data, mean_bin = bin_mean)
    
    # 合并分箱信息（保留 ID 用于合并所有其他字段）
    mean_bins <- merge(mean_data[, c("ID", "mean_bin")], sub_data, by = "ID")
    
    # 统计每个分箱内均值的平均值和站点数
    mean_stats <- mean_bins %>%
      group_by(mean_bin) %>%
      summarise(binmean_mean = mean(as.numeric(as.character(mean)), na.rm = TRUE),
                bin_n = n())
    
    mean_bins <- merge(mean_bins, mean_stats, by = "mean_bin")
    # 构造分箱标签：结果类别_箱内均值_箱内站点数
    mean_bins$Result_mean_n <- paste(mean_bins$result_GS,
                                     mean_bins$binmean_mean,
                                     mean_bins$bin_n, sep = "_")
    mean_bins_pop <- rbind(mean_bins_pop, mean_bins)
  }
  # 选择输出字段，同时保留位置信息
  out_mean <- mean_bins_pop[, c("ID", "chr", "start", "end", "strand", "Population", "result_GS", "mean", "SD", "Result_mean_n")]
  
  # 转换 chr 字段，仅保留点前的部分
  out_mean$chr <- sapply(out_mean$chr, convert_chr)
  
  # 将结果转换为标准 GTF 格式（9列）
  gtf_mean <- out_mean %>%
    mutate(
      source = "popbin",
      feature = "site",
      score = ".",
      frame = ".",
      attribute = paste0('ID "', ID, '"; ',
                         'Population "', Population, '"; ',
                         'result_GS "', result_GS, '"; ',
                         'mean "', mean, '"; ',
                         'SD "', SD, '"; ',
                         'Bin_Label "', Result_mean_n, '";')
    ) %>%
    select(seqname = chr, source, feature, start, end, score, strand, frame, attribute)
  
  output_filename <- paste0(pop, "_mean_bins.gtf")
  writeGTF_SiteClass(gtf_mean, output = output_filename)
  mean_bins_list[[pop]] <- gtf_mean
}

##################################################
### 4. 对标准差（SD）进行分箱并生成 GTF 文件 ###
##################################################
sd_bins_list <- list()  # 用于存储每个群体的 SD 分箱结果

for (pop in c("Marine", "Freshwater")) {
  pop_data <- subset(combined_data, Population == pop)
  sd_bins_pop <- NULL
  
  for (res in c("hyper", "hypo", "non")) {
    sub_data <- subset(pop_data, result_GS == res)
    if(nrow(sub_data) == 0) next
    
    # 对 SD 列进行分箱
    sd_data <- bin(sub_data, col = "SD", nbins = 50)
    sd_data <- rename(sd_data, SD_bin = bin_SD)
    
    sd_bins <- merge(sd_data[, c("ID", "SD_bin")], sub_data, by = "ID")
    
    sd_stats <- sd_bins %>%
      group_by(SD_bin) %>%
      summarise(binmean_SD = mean(as.numeric(as.character(SD)), na.rm = TRUE),
                bin_n = n())
    
    sd_bins <- merge(sd_bins, sd_stats, by = "SD_bin")
    sd_bins$Result_SD_n <- paste(sd_bins$result_GS,
                                 sd_bins$binmean_SD,
                                 sd_bins$bin_n, sep = "_")
    sd_bins_pop <- rbind(sd_bins_pop, sd_bins)
  }
  # 选择输出字段：保留位置信息
  out_sd <- sd_bins_pop[, c("ID", "chr", "start", "end", "strand", "Population", "result_GS", "mean", "SD", "Result_SD_n")]
  
  # 转换 chr 字段，仅保留点前的部分
  out_sd$chr <- sapply(out_sd$chr, convert_chr)
  
  # 转换为标准 GTF 格式
  gtf_sd <- out_sd %>%
    mutate(
      source = "popbin",
      feature = "site",
      score = ".",
      frame = ".",
      attribute = paste0('ID "', ID, '"; ',
                         'Population "', Population, '"; ',
                         'result_GS "', result_GS, '"; ',
                         'mean "', mean, '"; ',
                         'SD "', SD, '"; ',
                         'Bin_Label "', Result_SD_n, '";')
    ) %>%
    select(seqname = chr, source, feature, start, end, score, strand, frame, attribute)
  
  output_filename <- paste0(pop, "_SD_bins.gtf")
  writeGTF_SiteClass(gtf_sd, output = output_filename)
  sd_bins_list[[pop]] <- gtf_sd
}

for k in pi theta D; do
    perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/18_GOS.pileup \
        --gtf /work/cyu/meth/RSvsRob/RSRob_DMCs_final_output.gtf \
        --output /work/cyu/meth/RSvsRob/ROB_RSvsROB_${k}_final.txt \
        --measure ${k}
done

for k in pi theta D; do
    perl Variance-at-position.pl --pool-size 200 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/RSvsRob/RSRob_DMCs_final_output.gtf \
        --output /work/cyu/meth/RSvsRob/RS_RSvsROB_${k}_final.txt \
        --measure ${k}
done

    perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/RS_intermethunmeth.gtf \
        --output /work/cyu/meth/RSvsRob/RS_inter_${k}.txt \
        --measure pi



#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=your.email@example.com
#SBATCH --job-name="popoolation_single.sh"
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G

# 设置参数
pop="Freshwater"
bin_type="mean"
pileup_file="/work/cyu/meth/pooldata/pileup/18_GOS.pileup"
pool_size=400
gtf_file="/work/cyu/meth/RSvsGos/binning_outputs/pop_mean_bins_GOS.gtf"


echo "Processing $pop, bin type $bin_type"

# 针对三个多样性指标：pi, theta, D
for measure in pi theta D; do
  perl Variance-at-position.pl \
    --pool-size ${pool_size} \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup ${pileup_file} \
    --gtf ${gtf_file} \
    --output ${pop}_temp.${measure} \
    --measure ${measure}
  
  # 在每一行末尾添加群体标识，便于后续整合
  sed -i "s/$/\t${pop}/" ${pop}_temp.${measure}
done



#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=your.email@example.com
#SBATCH --job-name="popoolation_single.sh"
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G

# 设置参数
pop="Freshwater"
bin_type="mean"
pileup_file="/work/cyu/meth/pooldata/pileup/18_GOS.pileup"
pool_size=400
gtf_file="/work/cyu/meth/RSvsGos/binning_outputs/pop_SD_bins_GOS.gtf"


echo "Processing $pop, bin type $bin_type"

# 针对三个多样性指标：pi, theta, D
for measure in pi theta D; do
  perl Variance-at-position.pl \
    --pool-size ${pool_size} \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup ${pileup_file} \
    --gtf ${gtf_file} \
    --output ${pop}_temp.${measure} \
    --measure ${measure}
  
  # 在每一行末尾添加群体标识，便于后续整合
  sed -i "s/$/\t${pop}/" ${pop}_temp.${measure}
done



# 加载必要库
library(dplyr)
library(readr)

# 读取输入文件
all_cpgs <- read_tsv("RSvsGos_PercMeth_AllCpGs.tsv.gz", show_col_types = FALSE)
dmcs <- read_tsv("RSvsGos_PercMeth_DMCs.tsv.gz", show_col_types = FALSE)

# 确保关键列存在
stopifnot(all(c("chr", "start", "mean_RS", "SD_RS", "mean_Gos", "SD_Gos") %in% colnames(all_cpgs)))
stopifnot(all(c("chr", "start", "result_GS") %in% colnames(dmcs)))

# 对DMC位点建立标记（chr + start）
dmcs_key <- dmcs %>%
  mutate(result_GS = tolower(result_GS)) %>%
  select(chr, start, result_GS)

# 合并 allCpGs 和 dmcs，若没有匹配则 result_GS = "non"
combined <- all_cpgs %>%
  select(chr, start, mean_RS, SD_RS, mean_Gos, SD_Gos) %>%
  left_join(dmcs_key, by = c("chr", "start")) %>%
  mutate(result_GS = ifelse(is.na(result_GS), "non", result_GS))

# 构造 Marine 表格
marine <- combined %>%
  transmute(chr = gsub("\\..*", "", chr),
            start = start,
            end = start,
            strand = "+",
            mean = mean_RS,
            SD = SD_RS,
            result_GS = result_GS,
            Population = "Marine")

# 构造 Freshwater 表格
fresh <- combined %>%
  transmute(chr = gsub("\\..*", "", chr),
            start = start,
            end = start,
            strand = "+",
            mean = mean_Gos,
            SD = SD_Gos,
            result_GS = result_GS,
            Population = "Freshwater")

# 合并 Marine 和 Freshwater
final_combined <- bind_rows(marine, fresh) %>%
  mutate(ID = paste0("site", row_number()))

# 写入文件
write_tsv(final_combined, "combined_non_and_diff.tsv")

message("✅ 合并完成，文件保存为 combined_non_and_diff.tsv")







library(dplyr)
library(readr)
library(OneR)

# ==== 读取数据 ====
meth_long <- read_tsv("combined_non_and_diff.tsv", show_col_types = FALSE)
colnames(meth_long)[which(names(meth_long) == "result_GS")] <- "result_pop"

# ==== 创建输出文件夹 ====
dir.create("binning_outputs", showWarnings = FALSE)

# ==== 分箱函数 ====
bin_content <- function(df, nbins = 50, var = "mean") {
  OneR::bin(df[[var]], nbins = nbins, method = "content")
}

# ==== 写标准 GTF 的函数 ====
write_standard_gtf <- function(df, label_col, output_path) {
  gtf_ready <- df %>%
    select(chr, start, end, !!label_col) %>%
    mutate(
      source = "none",
      feature = "transcript",
      score = ".",
      strand = "+",
      frame = ".",
      attribute = paste0('gene_id "', .[[label_col]], '";')
    ) %>%
    select(chr, source, feature, start, end, score, strand, frame, attribute)
  
  write.table(gtf_ready, file = output_path, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ==== 主循环：每个群体分别做 mean 和 SD 分箱 ====
for (i in c("Marine", "Freshwater")) {
  meth_pop <- subset(meth_long, Population == i)
  pop_acc <- ifelse(i == "Marine", "RS", "GOS")
  
  # ==== MEAN BINS ====
  meth_mean_bins_pop <- NULL
  for (j in c("hyper", "hypo", "non")) {
    meth_result <- subset(meth_pop, result_pop == j)
    if (nrow(meth_result) == 0) next
    
    meth_bins <- bin_content(meth_result, nbins = 50, var = "mean")
    meth_result$mean_bin <- as.character(meth_bins)  # ⬅️ 关键修复
    
    bin_summary <- meth_result %>%
      group_by(mean_bin) %>%
      summarise(binmean_mean = mean(mean), bin_n = n(), .groups = "drop")
    
    meth_result <- left_join(meth_result, bin_summary, by = "mean_bin")
    meth_result$Result_mean_n <- paste(meth_result$result_pop,
                                       round(meth_result$binmean_mean, 2),
                                       meth_result$bin_n, sep = "_")
    
    meth_mean_bins_pop <- rbind(meth_mean_bins_pop, meth_result)
  }
  
  # 写入标准 GTF（mean）
  write_standard_gtf(meth_mean_bins_pop, "Result_mean_n",
                     paste0("binning_outputs/pop_mean_bins_", pop_acc, ".gtf"))
  
  # ==== SD BINS ====
  meth_SD_bins_pop <- NULL
  for (j in c("hyper", "hypo", "non")) {
    meth_result <- subset(meth_pop, result_pop == j)
    if (nrow(meth_result) == 0) next
    
    meth_bins <- bin_content(meth_result, nbins = 50, var = "SD")
    meth_result$SD_bin <- as.character(meth_bins)  # ⬅️ 关键修复
    
    bin_summary <- meth_result %>%
      group_by(SD_bin) %>%
      summarise(binmean_SD = mean(SD), bin_n = n(), .groups = "drop")
    
    meth_result <- left_join(meth_result, bin_summary, by = "SD_bin")
    meth_result$Result_SD_n <- paste(meth_result$result_pop,
                                     round(meth_result$binmean_SD, 2),
                                     meth_result$bin_n, sep = "_")
    
    meth_SD_bins_pop <- rbind(meth_SD_bins_pop, meth_result)
  }
  
  # 写入标准 GTF（SD）
  write_standard_gtf(meth_SD_bins_pop, "Result_SD_n",
                     paste0("binning_outputs/pop_SD_bins_", pop_acc, ".gtf"))
}

message("✅ 所有标准 GTF 文件已成功生成！")



#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=your.email@example.com
#SBATCH --job-name="popoolation_single.sh"
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G

# 设置参数
pop="Marine"
bin_type="mean"
pileup_file="/work/cyu/meth/pooldata/pileup/25_RS.pileup"
pool_size=400
gtf_file="/work/cyu/meth/RSvsGos/binning_outputs/pop_mean_bins_RS.gtf"

echo "Processing $pop, bin type $bin_type"

# 输出前缀
out_prefix="/work/cyu/meth/RSvsGos/binning_outputs/RS_RSvsGOS_mean_bin_temp"

# 针对三个多样性指标：pi, theta, D
for measure in pi theta D; do
  perl Variance-at-position.pl \
    --pool-size ${pool_size} \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup ${pileup_file} \
    --gtf ${gtf_file} \
    --output ${out_prefix}.${measure} \
    --measure ${measure}

  # 在每一行末尾添加群体标识，便于后续整合
  sed -i "s/$/\t${pop}/" ${out_prefix}.${measure}
done









library(ggplot2)

# 读取数据
freshwater_data <- read.table("GOS_RSvsGOS_mean_bin_temp.pi", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(freshwater_data) <- c("bin_id", "snp_count", "coverage_ratio", "pi", "Population")

marine_data <- read.table("RS_RSvsGOS_mean_bin_temp.pi", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(marine_data) <- c("bin_id", "snp_count", "coverage_ratio", "pi", "Population")

combined_data <- rbind(freshwater_data, marine_data)

# 拆分 bin_id
bin_parts <- strsplit(combined_data$bin_id, "_")
bin_matrix <- do.call(rbind, bin_parts)
combined_data$methylation_status <- bin_matrix[, 1]
combined_data$mean_methylation <- as.numeric(bin_matrix[, 2])
combined_data$n_sites <- as.numeric(bin_matrix[, 3])
combined_data$snp_count <- as.numeric(combined_data$snp_count)
combined_data$coverage_ratio <- as.numeric(combined_data$coverage_ratio)
combined_data$pi <- as.numeric(combined_data$pi)

# 设置因子
combined_data$methylation_status <- factor(combined_data$methylation_status, levels = c("non", "hypo", "hyper"))
combined_data$Population <- factor(combined_data$Population, levels = c("Marine", "Freshwater"))

# 标签映射
facet_labels <- c("non" = "Non-DMC", "hypo" = "FW-hypo", "hyper" = "FW-hyper")

# 生成图形对象
p <- ggplot(combined_data, aes(x = mean_methylation, y = pi, shape = Population, fill = Population)) +
  geom_point(color = "black", na.rm = TRUE, size = 2.5) +
  scale_shape_manual(values = c(21, 24), name = "Population") +
  scale_fill_manual(values = c("steelblue4", "palegreen"), name = "Population") +
  facet_wrap(~ methylation_status, labeller = as_labeller(facet_labels), scales = "free_y") +
  labs(
    x = "Rank Mean % Methylation",
    y = expression(paste(pi, " (rank-level)")),
    title = "B"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.background = element_rect(color = "black")
  )

# 保存图像（使用 ggsave）
ggsave("GOS_RSvsGOS_pi_plot.png", plot = p, width = 8, height = 4, dpi = 300)






library(dplyr)
library(readr)

# 读取 all_cpgs 和 dmcs 数据
all_cpgs <- read_tsv("RSvsGos_PercMeth_AllCpGs.tsv.gz", show_col_types = FALSE)
dmcs <- read_tsv("RSvsGos_PercMeth_DMCs.tsv.gz", show_col_types = FALSE)

# 对DMC位点建立标记（chr + start）
dmcs_key <- dmcs %>%
  mutate(result_GS = tolower(result_GS)) %>%
  select(chr, start, result_GS)

# 查找在 all_cpgs 中但不在 dmcs 中的位点，即 'non' 位点
non_sites <- all_cpgs %>%
  anti_join(dmcs_key, by = c("chr", "start")) %>%
  mutate(result_GS = "non")

# 提取 hyper 和 hypo 位点
dmcs_filtered <- dmcs %>%
  filter(result_GS %in% c("hyper", "hypo"))

# 合并 non 位点和 dmcs 中的 hyper 和 hypo 位点
combined <- bind_rows(non_sites, dmcs_filtered)

# 提取 Marine 和 Freshwater 的 mean 和 SD
marine_data <- all_cpgs %>%
  select(chr, start, mean_RS, SD_RS) %>%
  rename(mean = mean_RS, SD = SD_RS)

freshwater_data <- all_cpgs %>%
  select(chr, start, mean_Gos, SD_Gos) %>%
  rename(mean = mean_Gos, SD = SD_Gos)

# 合并 Marine 和 Freshwater 数据
combined_final <- combined %>%
  left_join(marine_data, by = c("chr", "start")) %>%
  left_join(freshwater_data, by = c("chr", "start"))

# 保存最终合并的数据
write_tsv(combined_final, "combined_non_and_diff_final.tsv")

message("✅ 合并完成，文件保存为 combined_non_and_diff_final.tsv")

# --- 检查 combined_non_and_diff_final.tsv 中的 'non'、'hypo' 和 'hyper' 位点数量 ---
result_count <- combined_final %>%
  count(result_GS)

# 输出每类的数量
cat("✅ combined_non_and_diff_final.tsv 中 'non'、'hyper' 和 'hypo' 位点数量：\n")
print(result_count)

# --- 检查 all_cpgs 中 'non' 位点数量 ---
all_non_sites <- all_cpgs %>%
  anti_join(dmcs_key, by = c("chr", "start")) %>%
  mutate(result_GS = "non")

cat("✅ all_cpgs 中 'non' 位点数量：", nrow(all_non_sites), "\n")

# --- 检查 dmcs 中 'hyper' 和 'hypo' 位点数量 ---
dmcs_hypo_count <- dmcs %>%
  filter(result_GS == "hypo") %>%
  nrow()

dmcs_hyper_count <- dmcs %>%
  filter(result_GS == "hyper") %>%
  nrow()

cat("✅ dmcs 中 'hypo' 位点数量：", dmcs_hypo_count, "\n")
cat("✅ dmcs 中 'hyper' 位点数量：", dmcs_hyper_count, "\n")









# 加载必要的库
library(dplyr)
library(readr)

# 读取 all_cpgs 和 dmcs 数据
all_cpgs <- read_tsv("RSvsGos_PercMeth_AllCpGs.tsv.gz", show_col_types = FALSE)
dmcs <- read_tsv("RSvsGos_PercMeth_DMCs.tsv.gz", show_col_types = FALSE)

# 对DMC位点建立标记（chr）
dmcs_key <- dmcs %>%
  mutate(result_GS = tolower(result_GS)) %>%
  select(chr, result_GS)

# 查找在 all_cpgs 中但不在 dmcs 中的位点，即 'non' 位点
non_sites <- all_cpgs %>%
  anti_join(dmcs_key, by = "chr") %>%
  mutate(result_GS = "non")

# 提取 hyper 和 hypo 位点
dmcs_filtered <- dmcs %>%
  filter(result_GS %in% c("hyper", "hypo"))

# 合并 non 位点和 dmcs 中的 hyper 和 hypo 位点
combined <- bind_rows(non_sites, dmcs_filtered)

# 提取 Marine 和 Freshwater 的 mean 和 SD
marine_data <- all_cpgs %>%
  select(chr, mean_RS, SD_RS) %>%
  rename(mean = mean_RS, SD = SD_RS)

freshwater_data <- all_cpgs %>%
  select(chr, mean_Gos, SD_Gos) %>%
  rename(mean = mean_Gos, SD = SD_Gos)

# 合并 Marine 和 Freshwater 数据
combined_final <- combined %>%
  left_join(marine_data, by = "chr") %>%
  left_join(freshwater_data, by = "chr")

# 保存最终合并的数据
write_tsv(combined_final, "combined_non_and_diff_final.tsv")

message("✅ 合并完成，文件保存为 combined_non_and_diff_final.tsv")

# --- 检查 combined_non_and_diff_final.tsv 中的 'non'、'hypo' 和 'hyper' 位点数量 ---
result_count <- combined_final %>%
  count(result_GS)

# 输出每类的数量
cat("✅ combined_non_and_diff_final.tsv 中 'non'、'hyper' 和 'hypo' 位点数量：\n")
print(result_count)

# --- 检查 all_cpgs 中 'non' 位点数量 ---
all_non_sites <- all_cpgs %>%
  anti_join(dmcs_key, by = "chr") %>%
  mutate(result_GS = "non")

cat("✅ all_cpgs 中 'non' 位点数量：", nrow(all_non_sites), "\n")

# --- 检查 dmcs 中 'hyper' 和 'hypo' 位点数量 ---
dmcs_hypo_count <- dmcs %>%
  filter(result_GS == "hypo") %>%
  nrow()

dmcs_hyper_count <- dmcs %>%
  filter(result_GS == "hyper") %>%
  nrow()

cat("✅ dmcs 中 'hypo' 位点数量：", dmcs_hypo_count, "\n")
cat("✅ dmcs 中 'hyper' 位点数量：", dmcs_hyper_count, "\n")




# 加载必要的库
library(dplyr)
library(readr)

# 读取 combined_non_and_diff_final.tsv 数据
combined_data <- read_tsv("combined_non_and_diff_final.tsv", show_col_types = FALSE)

# 提取 chr 列中的后缀部分作为 start 列的值
combined_data <- combined_data %>%
  mutate(start = gsub(".*\\.", "", chr))  # 使用正则表达式提取 chr 列中 '.' 后面的部分作为 start 列

# 查看修改后的数据
head(combined_data)

# 保存修改后的数据为新的文件
write_tsv(combined_data, "modified_combined_non_and_diff_final.tsv")

message("✅ 修改完成，文件保存为 modified_combined_non_and_diff_final.tsv")



library(dplyr)
library(readr)

# 读取 modified_combined_non_and_diff_final.tsv 文件
modified_combined_data <- read_tsv("modified_combined_non_and_diff_final.tsv")

# 查看列名，确保文件加载正确
print(colnames(modified_combined_data))

# 提取 freshwater 和 marine 数据
# 假设您想提取包含 mean_Gos 和 mean_RS 等列的数据
freshwater_data <- modified_combined_data %>%
  select(chr, start, mean_Gos, SD_Gos, result_GS)

marine_data <- modified_combined_data %>%
  select(chr, start, mean_RS, SD_RS, result_GS)

# 查看前几行数据
head(freshwater_data)
head(marine_data)

# 保存提取的数据
write_tsv(freshwater_data, "freshwater_data.tsv")
write_tsv(marine_data, "marine_data.tsv")





# 加载必要的库
library(dplyr)
library(tidyr)
library(readr)
library(OneR)

# 读取数据（请替换成你的文件路径）
freshwater_data <- read_tsv("freshwater_data.tsv")
marine_data <- read_tsv("marine_data.tsv")

# 添加唯一 ID（用于合并）
freshwater_data$ID <- paste0("FW_", seq_len(nrow(freshwater_data)))
marine_data$ID <- paste0("MA_", seq_len(nrow(marine_data)))

# 清理 chr 列后缀
clean_chr_format <- function(chr) {
  gsub("\\..*", "", chr)
}

# 写入 GTF 文件函数（通用）
writeGTF_SiteClass <- function(input, output) {
  sitedata <- input
  sitedata$source <- "none"
  sitedata$feature <- "transcript"
  sitedata$score <- "."
  sitedata$strand <- "+"
  sitedata$frame <- "."
  sitedata$seqname <- clean_chr_format(sitedata$chr)
  sitedata$gene_id <- paste0('gene_id "', sitedata$result_GS, "_", round(sitedata$binmean_sd, 2), "_", sitedata$bin_n, '"')
  sitedata$attribute <- paste(sitedata$gene_id, ";", sep = "")
  
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
  
  write.table(GTF, file = output, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  message(paste("✅ GTF written to", output))
}

# SD 分箱并写入 GTF（使用等频分箱）
bin_and_write_GTF_SD <- function(data, pop_type) {
  all_bins <- data.frame()
  column_to_use <- ifelse(pop_type == "Marine", "SD_RS", "SD_Gos")
  
  for (result_type in c("hyper", "hypo", "non")) {
    message(paste("🔄 Processing", result_type, "SD for", pop_type, "population"))
    
    data_result <- filter(data, result_GS == result_type)
    
    if (nrow(data_result) > 0) {
      # 添加唯一 ID（防冲突）
      data_result$ID <- paste0(pop_type, "_", result_type, "_", seq_len(nrow(data_result)))
      
      # 等频分箱
      meth_bins <- OneR::bin(data_result[, c("ID", column_to_use)], method = "content", nbins = 50)
      meth_bins$sd_bin <- as.numeric(meth_bins[[column_to_use]])
      
      # 合并回数据
      data_result <- merge(data_result, meth_bins[, c("ID", "sd_bin")], by = "ID")
      
      # 统计 bin 均值 + bin 内位点数
      sd_binmeth <- data_result %>%
        group_by(sd_bin) %>%
        summarise(binmean_sd = mean(.data[[column_to_use]], na.rm = TRUE),
                  bin_n = n())
      
      data_result <- merge(data_result, sd_binmeth, by = "sd_bin")
      
      # 打印信息
      message(paste("✅", result_type, "in", pop_type, "(SD) →", n_distinct(data_result$sd_bin), "bins with", nrow(data_result), "sites"))
      
      # 写 GTF 文件
      sitedata <- data_result
      sitedata$source <- "none"
      sitedata$feature <- "transcript"
      sitedata$score <- "."
      sitedata$strand <- "+"
      sitedata$frame <- "."
      sitedata$seqname <- clean_chr_format(sitedata$chr)
      sitedata$gene_id <- paste0('gene_id "', sitedata$result_GS, "_", round(sitedata$binmean_sd, 2), "_", sitedata$bin_n, '"')
      sitedata$attribute <- paste(sitedata$gene_id, ";")
      
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
      
      out_file <- paste0(pop_type, "_", result_type, "_sd_bins.gtf")
      write.table(GTF, file = out_file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
      message(paste("✅ GTF written to", out_file))
      
      all_bins <- rbind(all_bins, data_result)
    } else {
      message(paste("⚠️ No SD data found for", result_type, "in", pop_type))
    }
  }
  
  return(all_bins)
}

# 执行 SD 分箱并输出合并 GTF
marine_sd_bins <- bin_and_write_GTF_SD(marine_data, "Marine")
freshwater_sd_bins <- bin_and_write_GTF_SD(freshwater_data, "Freshwater")

writeGTF_SiteClass(marine_sd_bins, "combined_marine_SD_bins.gtf")
writeGTF_SiteClass(freshwater_sd_bins, "combined_freshwater_SD_bins.gtf")

message("✅ All SD binning GTF files written successfully!")




# 加载必要的库
library(dplyr)
library(tidyr)
library(readr)
library(OneR)

# 读取数据（请替换为你的实际路径）
freshwater_data <- read_tsv("freshwater_data.tsv")
marine_data <- read_tsv("marine_data.tsv")

# 添加唯一 ID（用于 merge）
freshwater_data$ID <- paste0("FW_", seq_len(nrow(freshwater_data)))
marine_data$ID <- paste0("MA_", seq_len(nrow(marine_data)))

# 清理 chr 列后缀
clean_chr_format <- function(chr) {
  gsub("\\..*", "", chr)
}

# 写入 GTF 文件的通用函数
writeGTF_SiteClass <- function(input, output) {
  sitedata <- input
  sitedata$source <- "none"
  sitedata$feature <- "transcript"
  sitedata$score <- "."
  sitedata$strand <- "+"
  sitedata$frame <- "."
  sitedata$seqname <- clean_chr_format(sitedata$chr)
  sitedata$gene_id <- paste0('gene_id "', sitedata$result_GS, "_", round(sitedata$binmean_mean, 2), "_", sitedata$bin_n, '"')
  sitedata$attribute <- paste(sitedata$gene_id, ";", sep = "")
  
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
  
  write.table(GTF, file = output, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  message(paste("✅ GTF written to", output))
}

# ✅ Mean 分箱函数（等频分箱）
bin_and_write_GTF_Mean <- function(data, pop_type) {
  all_bins <- data.frame()
  column_to_use <- ifelse(pop_type == "Marine", "mean_RS", "mean_Gos")
  
  for (result_type in c("hyper", "hypo", "non")) {
    message(paste("🔄 Processing", result_type, "Mean for", pop_type, "population"))
    
    data_result <- filter(data, result_GS == result_type)
    
    if (nrow(data_result) > 0) {
      data_result$ID <- paste0(pop_type, "_", result_type, "_", seq_len(nrow(data_result)))
      
      # 使用 OneR 等频分箱
      meth_bins <- OneR::bin(data_result[, c("ID", column_to_use)], method = "content", nbins = 50)
      meth_bins$mean_bin <- as.numeric(meth_bins[[column_to_use]])
      
      data_result <- merge(data_result, meth_bins[, c("ID", "mean_bin")], by = "ID")
      
      # 统计每个 bin 的均值与大小
      mean_binmeth <- data_result %>%
        group_by(mean_bin) %>%
        summarise(binmean_mean = mean(.data[[column_to_use]], na.rm = TRUE),
                  bin_n = n())
      
      data_result <- merge(data_result, mean_binmeth, by = "mean_bin")
      
      # 打印统计信息
      message(paste("✅", result_type, "in", pop_type, "(Mean) →", n_distinct(data_result$mean_bin), "bins with", nrow(data_result), "sites"))
      
      # 写出 GTF 文件
      sitedata <- data_result
      sitedata$source <- "none"
      sitedata$feature <- "transcript"
      sitedata$score <- "."
      sitedata$strand <- "+"
      sitedata$frame <- "."
      sitedata$seqname <- clean_chr_format(sitedata$chr)
      sitedata$gene_id <- paste0('gene_id "', sitedata$result_GS, "_", round(sitedata$binmean_mean, 2), "_", sitedata$bin_n, '"')
      sitedata$attribute <- paste(sitedata$gene_id, ";")
      
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
      
      out_file <- paste0(pop_type, "_", result_type, "_mean_bins.gtf")
      write.table(GTF, file = out_file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
      message(paste("✅ GTF written to", out_file))
      
      all_bins <- rbind(all_bins, data_result)
    } else {
      message(paste("⚠️ No mean data found for", result_type, "in", pop_type))
    }
  }
  
  return(all_bins)
}

# ✅ 运行 mean 分箱并写入最终合并 GTF
marine_mean_bins <- bin_and_write_GTF_Mean(marine_data, "Marine")
freshwater_mean_bins <- bin_and_write_GTF_Mean(freshwater_data, "Freshwater")

writeGTF_SiteClass(marine_mean_bins, "combined_marine_mean_bins.gtf")
writeGTF_SiteClass(freshwater_mean_bins, "combined_freshwater_mean_bins.gtf")

message("✅ All Mean binning GTF files written successfully!")

🔄 Processing hyper SD for Marine population
✅ hyper in Marine (SD) → 41 bins with 2058 sites
✅ GTF written to Marine_hyper_sd_bins.gtf
🔄 Processing hypo SD for Marine population
✅ hypo in Marine (SD) → 50 bins with 479 sites
✅ GTF written to Marine_hypo_sd_bins.gtf
🔄 Processing non SD for Marine population
✅ non in Marine (SD) → 36 bins with 127302 sites
✅ GTF written to Marine_non_sd_bins.gtf
🔄 Processing hyper SD for Freshwater population
✅ hyper in Freshwater (SD) → 40 bins with 2058 sites
✅ GTF written to Freshwater_hyper_sd_bins.gtf
🔄 Processing hypo SD for Freshwater population
✅ hypo in Freshwater (SD) → 31 bins with 479 sites
✅ GTF written to Freshwater_hypo_sd_bins.gtf
🔄 Processing non SD for Freshwater population
✅ non in Freshwater (SD) → 50 bins with 127302 sites
✅ GTF written to Freshwater_non_sd_bins.gtf
✅ GTF written to combined_marine_sd_bins.gtf
✅ GTF written to combined_freshwater_sd_bins.gtf
✅ All SD binning GTF files written successfully!

🔄 Processing hyper Mean for Marine population
✅ hyper in Marine (Mean) → 41 bins with 2058 sites
✅ GTF written to Marine_hyper_mean_bins.gtf
🔄 Processing hypo Mean for Marine population
✅ hypo in Marine (Mean) → 49 bins with 479 sites
✅ GTF written to Marine_hypo_mean_bins.gtf
🔄 Processing non Mean for Marine population
✅ non in Marine (Mean) → 37 bins with 127302 sites
✅ GTF written to Marine_non_mean_bins.gtf
🔄 Processing hyper Mean for Freshwater population
✅ hyper in Freshwater (Mean) → 39 bins with 2058 sites
✅ GTF written to Freshwater_hyper_mean_bins.gtf
🔄 Processing hypo Mean for Freshwater population
✅ hypo in Freshwater (Mean) → 31 bins with 479 sites
✅ GTF written to Freshwater_hypo_mean_bins.gtf
🔄 Processing non Mean for Freshwater population
✅ non in Freshwater (Mean) → 50 bins with 127302 sites
✅ GTF written to Freshwater_non_mean_bins.gtf





#!/bin/bash

# 设置种群信息
pop="Freshwater"
pileup_file="/work/cyu/meth/pooldata/pileup/18_GOS.pileup"

# 分析类型列表：SD 和 mean
for bin_type in SD mean; do
    gtf_file="/work/cyu/meth/RSvsGos/combined_freshwater_${bin_type}_bins.gtf"
    out_prefix="/work/cyu/meth/RSvsGos/binning_outputs/freshwater_${bin_type}_bin_temp"

    echo "🔄 Processing ${pop} population with ${bin_type} GTF: ${gtf_file}"

    # 针对三个多样性指标：pi, theta, D
    for measure in pi theta D; do
        echo "⚙️ Running Popoolation for ${measure}..."

        perl Variance-at-position.pl \
            --pool-size 400 \
            --min-qual 20 \
            --min-coverage 3 \
            --min-count 2 \
            --fastq-type sanger \
            --pileup ${pileup_file} \
            --gtf ${gtf_file} \
            --output ${out_prefix}.${measure} \
            --measure ${measure}

        # 添加群体标识，方便后续合并
        sed -i "s/$/\t${pop}/" ${out_prefix}.${measure}
    done

    echo "✅ Popoolation analysis completed for ${pop} ${bin_type} bins!"
done



#!/bin/bash

# 设置参数
pop="Marine"
pileup_file="/work/cyu/meth/pooldata/pileup/25_RS.pileup"

# 分析类型：mean 和 SD
for bin_type in "mean" "SD"; do
    gtf_file="/work/cyu/meth/RSvsGos/combined_marine_${bin_type}_bins.gtf"
    out_prefix="/work/cyu/meth/RSvsGos/binning_outputs/marine_${bin_type}_bin_temp"

    echo "🔄 Processing ${pop} population - ${bin_type} bins"
    echo "📄 GTF: ${gtf_file}"
    echo "📂 Pileup: ${pileup_file}"

    # 针对三个指标：pi, theta, D
    for measure in pi theta D; do
        echo "⚙️ Running Popoolation for ${measure}..."

        perl Variance-at-position.pl \
            --pool-size 400 \
            --min-qual 20 \
            --min-coverage 3 \
            --min-count 2 \
            --fastq-type sanger \
            --pileup ${pileup_file} \
            --gtf ${gtf_file} \
            --output ${out_prefix}.${measure} \
            --measure ${measure}

        # 添加群体标签
        sed -i "s/$/\t${pop}/" ${out_prefix}.${measure}
    done

    echo "✅ Popoolation analysis completed for ${pop} - ${bin_type} bins"
    echo "---------------------------------------------"
done

echo "🎉 All Marine Popoolation bin-based analyses completed!"



perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/RS_intermethunmeth.gtf \
        --output /work/cyu/meth/RSvsRob/RS_inter_${k}.txt \
        --measure pi



#!/bin/bash

# 设置参数
pop="Marine"
pileup_file="/work/cyu/meth/pooldata/pileup/25_RS.pileup"

# 分析类型：mean 和 SD
for bin_type in "mean" "SD"; do
    gtf_file="/work/cyu/meth/RSvsGos/combined_marine_${bin_type}_bins.gtf"
    out_prefix="/work/cyu/meth/RSvsGos/binning_outputs/marine_${bin_type}_bin_temp"

    echo "🔄 Processing ${pop} population - ${bin_type} bins"
    echo "📄 GTF: ${gtf_file}"
    echo "📂 Pileup: ${pileup_file}"

    # 只计算 pi
    measure="pi"
    echo "⚙️ Running Popoolation for ${measure}..."

    perl Variance-at-position.pl \
        --pool-size 400 \
        --min-qual 20 \
        --min-coverage 3 \
        --min-count 2 \
        --fastq-type sanger \
        --pileup ${pileup_file} \
        --gtf ${gtf_file} \
        --output ${out_prefix}.${measure} \
        --measure ${measure}

    # 添加群体标签
    sed -i "s/$/\t${pop}/" ${out_prefix}.${measure}

    echo "✅ Popoolation π analysis completed for ${pop} - ${bin_type} bins"
    echo "---------------------------------------------"
done

echo "🎉 All Marine π-only analyses completed!"



library(ggplot2)

# 读取数据
freshwater_data <- read.table("GOS_RSvsGOS_mean_bin_temp.pi", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(freshwater_data) <- c("bin_id", "snp_count", "coverage_ratio", "pi", "Population")

marine_data <- read.table("RS_RSvsGOS_mean_bin_temp.pi", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(marine_data) <- c("bin_id", "snp_count", "coverage_ratio", "pi", "Population")

combined_data <- rbind(freshwater_data, marine_data)

# 拆分 bin_id
bin_parts <- strsplit(combined_data$bin_id, "_")
bin_matrix <- do.call(rbind, bin_parts)
combined_data$methylation_status <- bin_matrix[, 1]
combined_data$mean_methylation <- as.numeric(bin_matrix[, 2])
combined_data$n_sites <- as.numeric(bin_matrix[, 3])
combined_data$snp_count <- as.numeric(combined_data$snp_count)
combined_data$coverage_ratio <- as.numeric(combined_data$coverage_ratio)
combined_data$pi <- as.numeric(combined_data$pi)

# 设置因子
combined_data$methylation_status <- factor(combined_data$methylation_status, levels = c("non", "hypo", "hyper"))
combined_data$Population <- factor(combined_data$Population, levels = c("Marine", "Freshwater"))

# 标签映射
facet_labels <- c("non" = "Non-DMC", "hypo" = "FW-hypo", "hyper" = "FW-hyper")

# 生成图形对象
p <- ggplot(combined_data, aes(x = mean_methylation, y = pi, shape = Population, fill = Population)) +
  geom_point(color = "black", na.rm = TRUE, size = 2.5) +
  scale_shape_manual(values = c(21, 24), name = "Population") +
  scale_fill_manual(values = c("steelblue4", "palegreen"), name = "Population") +
  facet_wrap(~ methylation_status, labeller = as_labeller(facet_labels), scales = "free_y") +
  labs(
    x = "Rank Mean % Methylation",
    y = expression(paste(pi, " (rank-level)")),
    title = "B"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.background = element_rect(color = "black")
  )

# 保存图像（使用 ggsave）
ggsave("GOS_RSvsGOS_pi_plot.png", plot = p, width = 8, height = 4, dpi = 300)








# 加载绘图包
library(ggplot2)

# 读取数据文件
freshwater_data <- read.table("freshwater_mean_bin_temp.pi", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
marine_data <- read.table("marine_mean_bin_temp.pi", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# 添加列名
colnames(freshwater_data) <- colnames(marine_data) <- c("bin_id", "snp_count", "coverage_ratio", "pi", "Population")

# 合并数据
combined_data <- rbind(freshwater_data, marine_data)

# 拆分 bin_id 字段（格式：result_binmean_mean_n）
bin_parts <- strsplit(combined_data$bin_id, "_")
bin_matrix <- do.call(rbind, bin_parts)
combined_data$methylation_status <- bin_matrix[, 1]  # non/hypo/hyper
combined_data$mean_methylation <- as.numeric(bin_matrix[, 2])
combined_data$n_sites <- as.numeric(bin_matrix[, 3])

# 转换其他列为数值型
combined_data$snp_count <- as.numeric(combined_data$snp_count)
combined_data$coverage_ratio <- as.numeric(combined_data$coverage_ratio)
combined_data$pi <- as.numeric(combined_data$pi)

# 设置因子顺序
combined_data$methylation_status <- factor(combined_data$methylation_status, levels = c("non", "hypo", "hyper"))
combined_data$Population <- factor(combined_data$Population, levels = c("Marine", "Freshwater"))

# Facet 标签
facet_labels <- c("non" = "Non-DMC", "hypo" = "FW-hypo", "hyper" = "FW-hyper")

# 生成图形对象
p <- ggplot(combined_data, aes(x = mean_methylation, y = pi, shape = Population, fill = Population)) +
  geom_point(color = "black", na.rm = TRUE, size = 2.5) +
  scale_shape_manual(values = c(21, 24), name = "Population") +
  scale_fill_manual(values = c("steelblue4", "palegreen"), name = "Population") +
  facet_wrap(~ methylation_status, labeller = as_labeller(facet_labels), scales = "free_y") +
  labs(
    x = "Rank Mean % Methylation",
    y = expression(paste(pi, " (rank-level)")),
    title = "B"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.background = element_rect(color = "black")
  )

# 保存图像（PNG 文件）
ggsave("figure3b_pi_plot.png", plot = p, width = 8, height = 4, dpi = 300)









#intermethunmeth.gtf for each population 
import pandas as pd
import os

# 定义文件名和 gene_id 标签
gtf_files = {
    "RS_inter_sites.gtf": "inter",
    "RS_Meth_sites.gtf": "meth",
    "RS_Unmeth_sites.gtf": "unmeth"
}

combined_gtf = []

for file, label in gtf_files.items():
    if os.path.exists(file):
        print(f"Reading {file}...")
        df = pd.read_csv(file, sep='\t', header=None,
                         names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
        
        # 修改 gene_id
        df["attribute"] = 'gene_id "' + df["seqname"] + '_' + label + '";'
        combined_gtf.append(df)
    else:
        print(f"❌ File not found: {file}")

# 合并并输出
if combined_gtf:
    final_gtf = pd.concat(combined_gtf)
    final_gtf.to_csv("RS_intermethunmeth.gtf", sep="\t", index=False, header=False, quoting=3)
    print("✅ 合并完成并保存为 RS_intermethunmeth.gtf")
else:
    print("⚠️ 没有读取到任何 GTF 文件，请检查文件路径")




perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/RS_intermethunmeth.gtf \
        --output /work/cyu/meth/RS_intermethunmeth.pi.txt \
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

# 设置 result 的因子顺序
all_data$result <- factor(all_data$result, levels = c("unmeth", "inter", "meth"))

# 清理 NaN 字符，并转换为数值型
all_data$statValue[all_data$statValue == "NaN"] <- NA
all_data$statValue <- as.numeric(all_data$statValue)

# 去除 statValue 中的 NA 行
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
# 3. 可视化 boxplot + jitter + p 值
# 设置标签位置
max_val <- max(all_data$statValue, na.rm = TRUE)

p <- ggplot(all_data, aes(x = result, y = statValue, fill = result)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2, shape = 21) +
  theme_bw() +
  labs(x = "Methylation Category", 
       y = expression(pi),
       title = "Nucleotide Diversity by Methylation State") +
  theme(legend.position = "none")

# 添加显著性标注
p <- p + geom_signif(comparisons = comparisons,
                     annotations = format.pval(wt_results$p.value, digits = 2, eps = 0.001),
                     y_position = seq(from = max_val*1.05, length.out = length(comparisons), by = 0.002),
                     tip_length = 0.02, textsize = 4)

# 显示图形
print(p)

# 保存图形
ggsave("RS_methylation_pi_boxplot.png", p, width = 8, height = 6, dpi = 300)












#!/bin/bash

# 输入文件
input="RSvsGos_PercMeth_AllCpGs.tsv.gz"
output="RS_selected_chr_split.tsv"

# 提取标题行，找到列号
header=$(zcat "$input" | head -1)

# 获取列号（chr 和3个样本）
chr_col=$(echo "$header" | tr '\t' '\n' | grep -n "^chr$" | cut -d: -f1)
col1=$(echo "$header" | tr '\t' '\n' | grep -n "^RS_008$" | cut -d: -f1)
col2=$(echo "$header" | tr '\t' '\n' | grep -n "^RS_016$" | cut -d: -f1)
col3=$(echo "$header" | tr '\t' '\n' | grep -n "^RS_017$" | cut -d: -f1)

# 检查是否找到列号
if [[ -z $chr_col || -z $col1 || -z $col2 || -z $col3 ]]; then
    echo "❌ 找不到某些列，请检查列名是否正确。"
    exit 1
fi

# 提取列并格式化
zcat "$input" | awk -v chr="$chr_col" -v c1="$col1" -v c2="$col2" -v c3="$col3" 'BEGIN{OFS="\t"} 
NR==1 { print "chr", "start", "RS_008", "RS_016", "RS_017"; next }
{
    split($chr, a, ".")
    print a[1], a[2], $c1, $c2, $c3
}' > "$output"

echo "✅ 已保存到 $output"


#!/bin/bash

input="RS_selected_chr_split.tsv"
output="RS_selected_long.tsv"

# 跳过表头，处理数据
awk 'NR==1 { next } { 
  print $1, $2, $3
  print $1, $2, $4
  print $1, $2, $5
}' OFS="\t" "$input" > "$output"

echo "✅ 已保存为 $output"





# 加载必要包
library(OneR)
library(dplyr)

# 读入原始三列数据（chr, start, value）
input_file <- "RS_selected_long.tsv"
data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(data) <- c("chr", "start", "value")

# 添加唯一 ID
data$ID <- paste0("ID_", seq_len(nrow(data)))

# Step 1: 使用 OneR 进行频率分箱（50 个 bin）
meth_bins <- OneR::bin(data[, c("ID", "value")], method = "content", nbins = 50)
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
out_file <- "RS_long_meth_bins.gtf"
write.table(GTF, file = out_file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

message(paste("✅ GTF written to", out_file))




perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/RSvsGos/RS_stacked_meth_value10bin.fixed.gtf \
        --output /work/cyu/meth/RS_allsite_methper_bin.pi.txt \
        --measure pi




library(ggplot2)

# 读取 RS 数据（来自 popoolation 输出）
marine_data <- read.table("RS_allsite_methper_bin.pi.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(marine_data) <- c("bin_id", "snp_count", "coverage_ratio", "pi")

# 添加一个标记列，作为后续 plot 的分面或颜色因子
marine_data$Population <- "Marine"

# 拆分 bin_id（如 bin1_27.35_159925）
bin_parts <- strsplit(marine_data$bin_id, "_")
bin_matrix <- do.call(rbind, bin_parts)

# 添加拆分后的列
marine_data$methylation_status <- bin_matrix[, 1]  # e.g. bin1
marine_data$mean_methylation <- as.numeric(bin_matrix[, 2])
marine_data$n_sites <- as.numeric(bin_matrix[, 3])
marine_data$snp_count <- as.numeric(marine_data$snp_count)
marine_data$coverage_ratio <- as.numeric(marine_data$coverage_ratio)
marine_data$pi <- as.numeric(marine_data$pi)

# 生成图形对象（Marine only，颜色保持 steelblue4）
p <- ggplot(marine_data, aes(x = mean_methylation, y = pi)) +
  geom_point(shape = 21, fill = "steelblue4", color = "black", size = 2.5, na.rm = TRUE) +
  labs(
    x = "Rank Mean % Methylation",
    y = expression(paste(pi, " (rank-level)")),
    title = "RS Population π per Bin"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "none"
  )

# 保存图像
ggsave("RS_meth_bin_pi_plot.png", plot = p, width = 6, height = 4, dpi = 300)







perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/RSvsGos/RS_long_meth_bins.gtf  \
        --output /work/cyu/meth/RS_allsite_methper_bin.pi.txt \
        --measure pi

perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 2 --min-count 1 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/RSvsGos/RS_long_meth_bins.gtf  \
        --output /work/cyu/meth/RS_allsite_methper_bin_mincover2count1.pi.txt \
        --measure pi


perl Variance-at-position.pl \
  --pileup 25_RS.pileup \
  --gtf test_bin.gtf \
  --output test_bin.pi.txt \
  --measure pi \
  --min-count 2 --min-coverage 3 --min-qual 20 --pool-size 400 --fastq-type sanger
perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 2 --min-count 1 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/test_bin.gtf  \
        --output /work/cyu/meth/test.pi.txt \
        --measure pi

perl Variance-at-position.pl \
  --pileup /work/cyu/meth/RSvsGos/matched_hyper.pileup \
  --gtf /work/cyu/meth/RSvsGos/test_chrXXI_hyper.gtf \
  --output /work/cyu/meth//RSvsGos/test_chrXXI_hyper.pi.txt \
  --measure pi \
  --pool-size 400 --min-count 2 --min-qual 20 --min-coverage 3 --fastq-type sanger



perl variance-sliding.pl --input matched_hyper.pileup --output test_chrXXI_hyper.pi.txt ...（参数）















#计算bg pi exclude chrY and chrUn

grep -v -E '^#|chrUn|chrY' /work/cyu/stickleback_v5_ensembl_genes.gtf > /work/cyu/stickleback_noUnY.gtf


perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
    --gtf /work/cyu/stickleback_noUnY_popoolation.gtf \
    --output /work/cyu/meth/RS_genomewide_pi.txt \
    --measure pi
  echo "完成 unmeth_${n}"

stickleback_wholegenome_cleaned.gtf

perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
    --gtf /work/cyu/stickleback_wholegenome_cleaned.gtf \
    --output /work/cyu/meth/RS_genomewide_pi.txt \
    --measure pi
  echo "完成 unmeth_${n}"



#10kb non-overlap step

grep -v -E 'chrY|chrUn' /work/cyu/meth/pooldata/pileup/25_RS.pileup > /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup

tmux bg
perl Variance-sliding.pl \
    --input /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
    --window-size 10000 \
    --step-size 10000 \
    --min-count 2 \
    --min-coverage 4 \
    --max-coverage 5000 \
    --min-qual 20 \
    --measure pi \
    --output /work/cyu/meth/RS_10kb_nostep_pi.txt \
    --pool-size 400 \
    --fastq-type sanger

awk '{sum += $5; n++} END {if (n > 0) print "Background pi:", sum / n; else print "No data"}' /work/cyu/meth/RS_10kb_nostep_pi.txt
Background pi: 0.00483154



#gos
# 1. 排除 chrY 和 chrUn
grep -v -E 'chrY|chrUn' /work/cyu/meth/pooldata/pileup/18_GOS.pileup > /work/cyu/meth/pooldata/pileup/18_GOS_noYUn.pileup

# 2. 计算窗口 π
perl Variance-sliding.pl \
    --input /work/cyu/meth/pooldata/pileup/18_GOS_noYUn.pileup \
    --window-size 10000 \
    --step-size 10000 \
    --min-count 2 \
    --min-coverage 4 \
    --max-coverage 5000 \
    --min-qual 20 \
    --measure pi \
    --output /work/cyu/meth/GOS_10kb_nostep_pi.txt \
    --pool-size 400 \
    --fastq-type sanger

# 3. 计算 genome-wide 平均 π
awk '{sum += $5; n++} END {if (n > 0) print "Background pi (GOS):", sum / n; else print "No data"}' /work/cyu/meth/GOS_10kb_nostep_pi.txt
Background pi (GOS): 0.00316497

#rob
# 1. 排除 chrY 和 chrUn
grep -v -E 'chrY|chrUn' /work/cyu/meth/pooldata/pileup/19_ROB.pileup > /work/cyu/meth/pooldata/pileup/19_ROB_noYUn.pileup

# 2. 计算窗口 π
perl Variance-sliding.pl \
    --input /work/cyu/meth/pooldata/pileup/19_ROB_noYUn.pileup \
    --window-size 10000 \
    --step-size 10000 \
    --min-count 2 \
    --min-coverage 4 \
    --max-coverage 5000 \
    --min-qual 20 \
    --measure pi \
    --output /work/cyu/meth/ROB_10kb_nostep_pi.txt \
    --pool-size 400 \
    --fastq-type sanger

# 3. 计算 genome-wide 平均 π
awk '{sum += $5; n++} END {if (n > 0) print "Background pi (ROB):", sum / n; else print "No data"}' /work/cyu/meth/ROB_10kb_nostep_pi.txt
Background pi (ROB): 0.00473019


#wik
# 1. 排除 chrY 和 chrUn
grep -v -E 'chrY|chrUn' /work/cyu/meth/pooldata/pileup/8_WK.pileup > /work/cyu/meth/pooldata/pileup/8_WK_noYUn.pileup

# 2. 计算窗口 π
perl Variance-sliding.pl \
    --input /work/cyu/meth/pooldata/pileup/8_WK_noYUn.pileup \
    --window-size 10000 \
    --step-size 10000 \
    --min-count 2 \
    --min-coverage 4 \
    --max-coverage 5000 \
    --min-qual 20 \
    --measure pi \
    --output /work/cyu/meth/WK_10kb_nostep_pi.txt \
    --pool-size 400 \
    --fastq-type sanger

# 3. 计算 genome-wide 平均 π
awk '{sum += $5; n++} END {if (n > 0) print "Background pi (WK):", sum / n; else print "No data"}' /work/cyu/meth/WK_10kb_nostep_pi.txt
Background pi (WK): 0.00405137

#watson
# 1. 排除 chrY 和 chrUn
grep -v -E 'chrY|chrUn' /work/cyu/meth/pooldata/pileup/7_WT.pileup > /work/cyu/meth/pooldata/pileup/7_WT_noYUn.pileup

# 2. 计算窗口 π
perl Variance-sliding.pl \
    --input /work/cyu/meth/pooldata/pileup/7_WT_noYUn.pileup \
    --window-size 10000 \
    --step-size 10000 \
    --min-count 2 \
    --min-coverage 4 \
    --max-coverage 5000 \
    --min-qual 20 \
    --measure pi \
    --output /work/cyu/meth/WT_10kb_nostep_pi.txt \
    --pool-size 400 \
    --fastq-type sanger

# 3. 计算 genome-wide 平均 π
awk '{sum += $5; n++} END {if (n > 0) print "Background pi (WT):", sum / n; else print "No data"}' /work/cyu/meth/WT_10kb_nostep_pi.txt




#sl
# 1. 排除 chrY 和 chrUn
grep -v -E 'chrY|chrUn' /work/cyu/meth/pooldata/pileup/4_SL.pileup > /work/cyu/meth/pooldata/pileup/4_SL_noYUn.pileup

# 2. 计算窗口 π
perl Variance-sliding.pl \
    --input /work/cyu/meth/pooldata/pileup/4_SL_noYUn.pileup \
    --window-size 10000 \
    --step-size 10000 \
    --min-count 2 \
    --min-coverage 4 \
    --max-coverage 5000 \
    --min-qual 20 \
    --measure pi \
    --output /work/cyu/meth/SL_10kb_nostep_pi.txt \
    --pool-size 148 \
    --fastq-type sanger

# 3. 计算 genome-wide 平均 π
awk '{sum += $5; n++} END {if (n > 0) print "Background pi (SL):", sum / n; else print "No data"}' /work/cyu/meth/SL_10kb_nostep_pi.txt
Background pi (SL): 0.00448732








#persite

import pandas as pd
import gzip

# 输入文件路径
input_file = "RS_AllSites_Fig.tsv.gz"
output_file = "RS_AllSites_Fig.gtf"

gtf_rows = []

# 读取压缩文件
with gzip.open(input_file, "rt") as f:
    df = pd.read_csv(f, sep="\t", header=0)

    for idx, row in df.iterrows():
        chrom = row["chr"]
        start = int(row["start"])
        end = int(row["end"])
        pos_id = f"{chrom}.{start}"  # 如 chrI.10002181

        gtf_row = [
            chrom,
            "none",
            "transcript",
            start,
            end,
            ".",
            "+",
            ".",
            f'gene_id "{pos_id}";'
        ]
        gtf_rows.append(gtf_row)

# 写出 GTF 文件
if gtf_rows:
    gtf_df = pd.DataFrame(gtf_rows)
    gtf_df.to_csv(output_file, sep="\t", index=False, header=False, quoting=3)
    print(f"✅ GTF 文件已生成：{output_file}")
else:
    print("⚠️ 没有生成任何 GTF 内容，请检查输入数据格式")

perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS.pileup \
        --gtf /work/cyu/meth/RabbitSlough/RS_AllSites_Fig.gtf \
        --output /work/cyu/meth/RabbitSlough/RS_AllSites_Fig.pi.txt \
        --measure pi