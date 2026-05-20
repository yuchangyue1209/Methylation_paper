#violin hyper hypo non


#RSGos_CpGs_minDepth5.bed.gz  3 island GTF
#!/usr/bin/env python
import pandas as pd
import csv

# ===== 配置 =====
infile = "RSGos_CpGs_minDepth5.bed.gz"   # 你的输入
island_len = 1000                        # 一个 island 最大跨度（从起点算）
cats = ["non", "hypo", "hyper"]          # 三种类别

# ===== 读入 BED =====
# 注意：原始文件用空格分隔并带引号
df = pd.read_csv(
    infile,
    sep=" ",
    header=0,
    quotechar='"',
    compression="gzip"
)

# 拆分 "chrI.10005485" → chrI 和 10005485
df[["chr_clean", "chr_pos"]] = df["chr"].str.split(".", expand=True)
df["chr"] = df["chr_clean"]
df["start"] = df["chr_pos"].astype(int)

# DMC 列就是 non / hypo / hyper
df["DMC"] = df["DMC"].astype(str)

# 去重：同一 chr/start 不要重复
df = df[["chr", "start", "DMC"]].drop_duplicates()

print("总 CpG 数：")
print(df["DMC"].value_counts())

# ===== 对每个类别分别做 island 并输出 GTF =====
for cat in cats:
    sub = df[df["DMC"] == cat].copy()
    if sub.empty:
        print(f"[{cat}] 没有位点，跳过")
        continue

    # 按 chr, start 排序
    sub = sub.sort_values(["chr", "start"])

    # 按染色体分组
    gtf_rows = []

    for chrom, g in sub.groupby("chr"):
        g = g.sort_values("start")
        starts = g["start"].tolist()

        island_idx = 1
        island_start = starts[0]
        current_id = f"island{island_idx}_{cat}"

        # 为每个 CpG 分配 island gene_id
        gene_ids = []
        for pos in starts:
            if pos - island_start >= island_len:
                island_idx += 1
                island_start = pos
                current_id = f"island{island_idx}_{cat}"
            gene_ids.append(current_id)

        # 构建这个 chr 上的 GTF 行
        tmp = pd.DataFrame({
            "seqname": chrom,
            "source": "none",
            "feature": "transcript",
            "start": starts,
            "end": starts,
            "score": ".",
            "strand": "+",
            "frame": ".",
            "attribute": [f'gene_id "{gid}";' for gid in gene_ids]
        })
        gtf_rows.append(tmp)

    gtf_df = pd.concat(gtf_rows, ignore_index=True)

    outname = f"RSGos_{cat}_islands_1kb.gtf"
    gtf_df.to_csv(
        outname,
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
        escapechar="\\"
    )

    print(f"✅ 写出 {cat} GTF: {outname}, 行数 = {len(gtf_df)}")


RSGos_non_islands_1kb.gtf

RSGos_hypo_islands_1kb.gtf

RSGos_hyper_islands_1kb.gtf


for k in pi; do
perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/18_GOS_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsGos/RSGos_non_islands_1kb.gtf \
    --output /work/cyu/meth/RSvsGos/GOS_RSvsGOS_pi_non_islands_1k.txt \
    --measure pi
done

for k in pi; do
perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/18_GOS_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsGos/RSGos_hypo_islands_1kb.gtf \
    --output /work/cyu/meth/RSvsGos/GOS_RSvsGOS_pi_hypo_islands_1k.txt \
    --measure pi
done

for k in pi; do
perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/18_GOS_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsGos/RSGos_hyper_islands_1kb.gtf \
    --output /work/cyu/meth/RSvsGos/GOS_RSvsGOS_pi_hyper_islands_1k.txt \
    --measure pi
done



for k in pi; do
perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsGos/RSGos_non_islands_1kb.gtf \
    --output /work/cyu/meth/RSvsGos/RS_RSvsGOS_pi_non_islands_1k.txt \
    --measure pi
done



for k in pi; do
perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsGos/RSGos_hypo_islands_1kb.gtf \
    --output /work/cyu/meth/RSvsGos/RS_RSvsGOS_pi_hypo_islands_1k.txt \
    --measure pi
done


for k in pi; do
perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsGos/RSGos_hyper_islands_1kb.gtf \
    --output /work/cyu/meth/RSvsGos/RS_RSvsGOS_pi_hyper_islands_1k.txt \
    --measure pi
done



#rsrob
#!/usr/bin/env python
import pandas as pd
import csv

# ===== 配置 =====
infile = "RSRob_CpGs_minDepth5.bed.gz"   # 输入文件
island_len = 1000                        # 一个 island 最大跨度（从起点算）
cats = ["non", "hypo", "hyper"]          # 三种类别

# ===== 读入 BED =====
# 格式类似：
# "chr" "start" "end" "strand" "pvalue" "qvalue" "meth.diff" "DMC"
df = pd.read_csv(
    infile,
    sep=" ",
    header=0,
    quotechar='"',
    compression="gzip"
)

# 拆分 "chrI.10005485" → chrI 和 10005485
df[["chr_clean", "chr_pos"]] = df["chr"].str.split(".", expand=True)
df["chr"] = df["chr_clean"]
df["start"] = df["chr_pos"].astype(int)

# DMC 列就是 non / hypo / hyper
df["DMC"] = df["DMC"].astype(str)

# 去重：同一 chr/start 不要重复
df = df[["chr", "start", "DMC"]].drop_duplicates()

print("RSRob 总 CpG 数：")
print(df["DMC"].value_counts())

# ===== 对每个类别分别做 island 并输出 GTF =====
for cat in cats:
    sub = df[df["DMC"] == cat].copy()
    if sub.empty:
        print(f"[{cat}] 没有位点，跳过")
        continue

    # 按 chr, start 排序
    sub = sub.sort_values(["chr", "start"])

    # 按染色体分组
    gtf_rows = []

    for chrom, g in sub.groupby("chr"):
        g = g.sort_values("start")
        starts = g["start"].tolist()

        island_idx = 1
        island_start = starts[0]
        current_id = f"island{island_idx}_{cat}"

        # 为每个 CpG 分配 island gene_id
        gene_ids = []
        for pos in starts:
            if pos - island_start >= island_len:
                island_idx += 1
                island_start = pos
                current_id = f"island{island_idx}_{cat}"
            gene_ids.append(current_id)

        # 构建这个 chr 上的 GTF 行
        tmp = pd.DataFrame({
            "seqname": chrom,
            "source": "none",
            "feature": "transcript",
            "start": starts,
            "end": starts,
            "score": ".",
            "strand": "+",
            "frame": ".",
            "attribute": [f'gene_id "{gid}";' for gid in gene_ids]
        })
        gtf_rows.append(tmp)

    gtf_df = pd.concat(gtf_rows, ignore_index=True)

    outname = f"RSRob_{cat}_islands_1kb.gtf"
    gtf_df.to_csv(
        outname,
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
        escapechar="\\"
    )

    print(f"✅ 写出 {cat} GTF: {outname}, 行数 = {len(gtf_df)}")


# non
for k in pi; do
    perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/19_ROB_noYUn.pileup \
        --gtf /work/cyu/meth/RSvsRob/RSRob_non_islands_1kb.gtf \
        --output /work/cyu/meth/RSvsRob/ROB_RSvsROB_pi_non_islands_1k.txt \
        --measure ${k}
done

# hypo
for k in pi; do
    perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/19_ROB_noYUn.pileup \
        --gtf /work/cyu/meth/RSvsRob/RSRob_hypo_islands_1kb.gtf \
        --output /work/cyu/meth/RSvsRob/ROB_RSvsROB_pi_hypo_islands_1k.txt \
        --measure ${k}
done

# hyper
for k in pi; do
    perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/19_ROB_noYUn.pileup \
        --gtf /work/cyu/meth/RSvsRob/RSRob_hyper_islands_1kb.gtf \
        --output /work/cyu/meth/RSvsRob/ROB_RSvsROB_pi_hyper_islands_1k.txt \
        --measure ${k}
done

# non
for k in pi; do
    perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
        --gtf /work/cyu/meth/RSvsRob/RSRob_non_islands_1kb.gtf \
        --output /work/cyu/meth/RSvsRob/RS_RSvsROB_pi_non_islands_1k.txt \
        --measure ${k}
done

# hypo
for k in pi; do
    perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
        --gtf /work/cyu/meth/RSvsRob/RSRob_hypo_islands_1kb.gtf \
        --output /work/cyu/meth/RSvsRob/RS_RSvsROB_pi_hypo_islands_1k.txt \
        --measure ${k}
done

# hyper
for k in pi; do
    perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
        --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
        --gtf /work/cyu/meth/RSvsRob/RSRob_hyper_islands_1kb.gtf \
        --output /work/cyu/meth/RSvsRob/RS_RSvsROB_pi_hyper_islands_1k.txt \
        --measure ${k}
done


#rswt
#!/usr/bin/env python
import pandas as pd
import csv

# ===== 配置 =====
infile = "RSWGBSWT_CpGs_minDepth5.bed.gz"   # 输入文件名
island_len = 1000                           # 一个 island 最大跨度（从起点算）
cats = ["non", "hypo", "hyper"]             # 三种类别

# ===== 读入 BED =====
# 你的文件头是：
# "chr" "start" "end" "strand" "pvalue" "qvalue" "meth.diff" "result_GS"
df = pd.read_csv(
    infile,
    sep=" ",
    header=0,
    quotechar='"',
    compression="gzip"
)

print("列名:", list(df.columns))

# 使用 result_GS 作为 non/hypo/hyper 分类列
if "result_GS" not in df.columns:
    raise ValueError("找不到列 'result_GS'，请确认文件格式/列名")

df["DMC"] = df["result_GS"].astype(str)

# 拆分 chr 列： "chrI.1000020" → "chrI" + "1000020"
# 这里只取前半部分作为真正的染色体名
df[["chr_clean", "dummy_id"]] = df["chr"].astype(str).str.replace('"', "", regex=False).str.split(".", expand=True)
df["chr"] = df["chr_clean"]

# start 列已经是数字坐标（第二列）
df["start"] = df["start"].astype(int)

# 去重：同一 chr/start 不要重复
df = df[["chr", "start", "DMC"]].drop_duplicates()

print("RSWGBSWT 总 CpG 数（按类别）:")
print(df["DMC"].value_counts())

# ===== 对每个类别分别做 island 并输出 GTF =====
for cat in cats:
    sub = df[df["DMC"] == cat].copy()
    if sub.empty:
        print(f"[{cat}] 没有位点，跳过")
        continue

    # 按 chr, start 排序
    sub = sub.sort_values(["chr", "start"])

    # 按染色体分组
    gtf_rows = []

    for chrom, g in sub.groupby("chr"):
        g = g.sort_values("start")
        starts = g["start"].tolist()

        island_idx = 1
        island_start = starts[0]
        current_id = f"island{island_idx}_{cat}"

        # 为每个 CpG 分配 island gene_id
        gene_ids = []
        for pos in starts:
            if pos - island_start >= island_len:
                island_idx += 1
                island_start = pos
                current_id = f"island{island_idx}_{cat}"
            gene_ids.append(current_id)

        # 构建这个 chr 上的 GTF 行
        tmp = pd.DataFrame({
            "seqname": chrom,
            "source": "none",
            "feature": "transcript",
            "start": starts,
            "end": starts,
            "score": ".",
            "strand": "+",
            "frame": ".",
            "attribute": [f'gene_id "{gid}";' for gid in gene_ids]
        })
        gtf_rows.append(tmp)

    gtf_df = pd.concat(gtf_rows, ignore_index=True)

    outname = f"RSWGBSWT_{cat}_islands_1kb.gtf"
    gtf_df.to_csv(
        outname,
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
        escapechar="\\"
    )

    print(f"✅ 写出 {cat} GTF: {outname}, 行数 = {len(gtf_df)}")

cd /work/cyu/meth/RSvsWat

########## WT (WGBS) – non islands ##########
for k in pi; do
  perl Variance-at-position.pl \
    --pool-size 184 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/7_WT_noYUn.pileup \
    --gtf    /work/cyu/meth/RSvsWat/RSWGBSWT_non_islands_1kb.gtf \
    --output /work/cyu/meth/RSvsWat/WGBSWT_RSWGBSWT_${k}_non_islands_1k.txt \
    --measure ${k}
done

########## WT (WGBS) – hypo islands ##########
for k in pi; do
  perl Variance-at-position.pl \
    --pool-size 184 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/7_WT_noYUn.pileup \
    --gtf    /work/cyu/meth/RSvsWat/RSWGBSWT_hypo_islands_1kb.gtf \
    --output /work/cyu/meth/RSvsWat/WGBSWT_RSWGBSWT_${k}_hypo_islands_1k.txt \
    --measure ${k}
done

########## WT (WGBS) – hyper islands ##########
for k in pi; do
  perl Variance-at-position.pl \
    --pool-size 184 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/7_WT_noYUn.pileup \
    --gtf    /work/cyu/meth/RSvsWat/RSWGBSWT_hyper_islands_1kb.gtf \
    --output /work/cyu/meth/RSvsWat/WGBSWT_RSWGBSWT_${k}_hyper_islands_1k.txt \
    --measure ${k}
done


########## RS – non islands ##########
for k in pi; do
  perl Variance-at-position.pl \
    --pool-size 400 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
    --gtf    /work/cyu/meth/RSvsWat/RSWGBSWT_non_islands_1kb.gtf \
    --output /work/cyu/meth/RSvsWat/RS_RSWGBSWT_${k}_non_islands_1k.txt \
    --measure ${k}
done

########## RS – hypo islands ##########
for k in pi; do
  perl Variance-at-position.pl \
    --pool-size 400 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
    --gtf    /work/cyu/meth/RSvsWat/RSWGBSWT_hypo_islands_1kb.gtf \
    --output /work/cyu/meth/RSvsWat/RS_RSWGBSWT_${k}_hypo_islands_1k.txt \
    --measure ${k}
done

########## RS – hyper islands ##########
for k in pi; do
  perl Variance-at-position.pl \
    --pool-size 400 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
    --gtf    /work/cyu/meth/RSvsWat/RSWGBSWT_hyper_islands_1kb.gtf \
    --output /work/cyu/meth/RSvsWat/RS_RSWGBSWT_${k}_hyper_islands_1k.txt \
    --measure ${k}
done




#rswik

cd /work/cyu/meth/RSvsWik

########## WK – non islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 198 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/8_WK_noYUn.pileup \
  --gtf    /work/cyu/meth/RSvsWik/RSWGBSWK_non_islands_1kb.gtf \
  --output /work/cyu/meth/RSvsWik/WK_RSWGBSWK_pi_non_islands_1k.txt \
  --measure pi


########## WK – hypo islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 198 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/8_WK_noYUn.pileup \
  --gtf    /work/cyu/meth/RSvsWik/RSWGBSWK_hypo_islands_1kb.gtf \
  --output /work/cyu/meth/RSvsWik/WK_RSWGBSWK_pi_hypo_islands_1k.txt \
  --measure pi


########## WK – hyper islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 198 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/8_WK_noYUn.pileup \
  --gtf    /work/cyu/meth/RSvsWik/RSWGBSWK_hyper_islands_1kb.gtf \
  --output /work/cyu/meth/RSvsWik/WK_RSWGBSWK_pi_hyper_islands_1k.txt \
  --measure pi



########## RS – non islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 400 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
  --gtf    /work/cyu/meth/RSvsWik/RSWGBSWK_non_islands_1kb.gtf \
  --output /work/cyu/meth/RSvsWik/RS_RSWGBSWK_pi_non_islands_1k.txt \
  --measure pi


########## RS – hypo islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 400 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
  --gtf    /work/cyu/meth/RSvsWik/RSWGBSWK_hypo_islands_1kb.gtf \
  --output /work/cyu/meth/RSvsWik/RS_RSWGBSWK_pi_hypo_islands_1k.txt \
  --measure pi


########## RS – hyper islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 400 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
  --gtf    /work/cyu/meth/RSvsWik/RSWGBSWK_hyper_islands_1kb.gtf \
  --output /work/cyu/meth/RSvsWik/RS_RSWGBSWK_pi_hyper_islands_1k.txt \
  --measure pi




#saygos


cd /work/cyu/meth/SayvsGos

########## SAY – non islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 200 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/17_SAT_noYUn.pileup \
  --gtf    /work/cyu/meth/SayvsGos/SayGos_non_islands_1kb.gtf \
  --output /work/cyu/meth/SayvsGos/SAY_SayGos_pi_non_islands_1k.txt \
  --measure pi

########## SAY – hypo islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 200 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/17_SAT_noYUn.pileup \
  --gtf    /work/cyu/meth/SayvsGos/SayGos_hypo_islands_1kb.gtf \
  --output /work/cyu/meth/SayvsGos/SAY_SayGos_pi_hypo_islands_1k.txt \
  --measure pi

########## SAY – hyper islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 200 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/17_SAT_noYUn.pileup \
  --gtf    /work/cyu/meth/SayvsGos/SayGos_hyper_islands_1kb.gtf \
  --output /work/cyu/meth/SayvsGos/SAY_SayGos_pi_hyper_islands_1k.txt \
  --measure pi



########## GOS – non islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 400 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/18_GOS_noYUn.pileup \
  --gtf    /work/cyu/meth/SayvsGos/SayGos_non_islands_1kb.gtf \
  --output /work/cyu/meth/SayvsGos/GOS_SayGos_pi_non_islands_1k.txt \
  --measure pi

########## GOS – hypo islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 400 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/18_GOS_noYUn.pileup \
  --gtf    /work/cyu/meth/SayvsGos/SayGos_hypo_islands_1kb.gtf \
  --output /work/cyu/meth/SayvsGos/GOS_SayGos_pi_hypo_islands_1k.txt \
  --measure pi

########## GOS – hyper islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 400 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/18_GOS_noYUn.pileup \
  --gtf    /work/cyu/meth/SayvsGos/SayGos_hyper_islands_1kb.gtf \
  --output /work/cyu/meth/SayvsGos/GOS_SayGos_pi_hyper_islands_1k.txt \
  --measure pi



#sayrob

cd /work/cyu/meth/SayvsRob

########## SAY – non islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 200 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/17_SAT_noYUn.pileup \
  --gtf    /work/cyu/meth/SayvsRob/SayRob_non_islands_1kb.gtf \
  --output /work/cyu/meth/SayvsRob/SAY_SayRob_pi_non_islands_1k.txt \
  --measure pi

########## SAY – hypo islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 200 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/17_SAT_noYUn.pileup \
  --gtf    /work/cyu/meth/SayvsRob/SayRob_hypo_islands_1kb.gtf \
  --output /work/cyu/meth/SayvsRob/SAY_SayRob_pi_hypo_islands_1k.txt \
  --measure pi

########## SAY – hyper islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 200 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/17_SAT_noYUn.pileup \
  --gtf    /work/cyu/meth/SayvsRob/SayRob_hyper_islands_1kb.gtf \
  --output /work/cyu/meth/SayvsRob/SAY_SayRob_pi_hyper_islands_1k.txt \
  --measure pi



########## ROB – non islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 400 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/19_ROB_noYUn.pileup \
  --gtf    /work/cyu/meth/SayvsRob/SayRob_non_islands_1kb.gtf \
  --output /work/cyu/meth/SayvsRob/ROB_SayRob_pi_non_islands_1k.txt \
  --measure pi

########## ROB – hypo islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 400 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/19_ROB_noYUn.pileup \
  --gtf    /work/cyu/meth/SayvsRob/SayRob_hypo_islands_1kb.gtf \
  --output /work/cyu/meth/SayvsRob/ROB_SayRob_pi_hypo_islands_1k.txt \
  --measure pi

########## ROB – hyper islands (π) ##########
perl Variance-at-position.pl \
  --pool-size 400 \
  --min-qual 20 \
  --min-coverage 3 \
  --min-count 2 \
  --fastq-type sanger \
  --pileup /work/cyu/meth/pooldata/pileup/19_ROB_noYUn.pileup \
  --gtf    /work/cyu/meth/SayvsRob/SayRob_hyper_islands_1kb.gtf \
  --output /work/cyu/meth/SayvsRob/ROB_SayRob_pi_hyper_islands_1k.txt \
  --measure pi






#5k
# #RSGOS
#!/usr/bin/env python
import pandas as pd
import csv

# ===== 配置 =====
infile = "RSGos_CpGs_minDepth5.bed.gz"
island_len = 5000
cats = ["non", "hypo", "hyper"]

df = pd.read_csv(
    infile,
    sep=" ",
    header=0,
    quotechar='"',
    compression="gzip"
)

# 拆分 "chrI.10005485" → chrI 和 10005485
df[["chr_clean", "chr_pos"]] = df["chr"].astype(str).str.replace('"', "", regex=False).str.split(".", expand=True)
df["chr"] = df["chr_clean"]
df["start"] = df["chr_pos"].astype(int)

# 分类列
df["DMC"] = df["DMC"].astype(str)

# 去重
df = df[["chr", "start", "DMC"]].drop_duplicates()

print("RSGos 总 CpG 数：")
print(df["DMC"].value_counts())

for cat in cats:
    sub = df[df["DMC"] == cat].copy()
    if sub.empty:
        print(f"[{cat}] 没有位点，跳过")
        continue

    sub = sub.sort_values(["chr", "start"])
    gtf_rows = []

    for chrom, g in sub.groupby("chr"):
        starts = g["start"].tolist()

        island_idx = 1
        island_start = starts[0]
        current_id = f"island{island_idx}_{cat}"

        gene_ids = []
        for pos in starts:
            if pos - island_start >= island_len:
                island_idx += 1
                island_start = pos
                current_id = f"island{island_idx}_{cat}"
            gene_ids.append(current_id)

        tmp = pd.DataFrame({
            "seqname": chrom,
            "source": "none",
            "feature": "transcript",
            "start": starts,
            "end": starts,
            "score": ".",
            "strand": "+",
            "frame": ".",
            "attribute": [f'gene_id "{gid}";' for gid in gene_ids]
        })
        gtf_rows.append(tmp)

    gtf_df = pd.concat(gtf_rows, ignore_index=True)

    outname = f"RSGos_{cat}_islands_5kb.gtf"
    gtf_df.to_csv(
        outname,
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
        escapechar="\\"
    )
    print(f"✅ 写出 {cat} GTF: {outname}, 行数 = {len(gtf_df)}")
cd /work/cyu/meth/RSvsGos

# ========= GOS (400) =========
for cat in non hypo hyper; do
  perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/18_GOS_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsGos/RSGos_${cat}_islands_5kb.gtf \
    --output /work/cyu/meth/RSvsGos/GOS_RSvsGOS_pi_${cat}_islands_5k.txt \
    --measure pi
done

# ========= RS (400) =========
for cat in non hypo hyper; do
  perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsGos/RSGos_${cat}_islands_5kb.gtf \
    --output /work/cyu/meth/RSvsGos/RS_RSvsGOS_pi_${cat}_islands_5k.txt \
    --measure pi
done


#RSROB
#!/usr/bin/env python
import pandas as pd
import csv

# ===== 配置 =====
infile = "RSRob_CpGs_minDepth5.bed.gz"
island_len = 5000
cats = ["non", "hypo", "hyper"]

df = pd.read_csv(
    infile,
    sep=" ",
    header=0,
    quotechar='"',
    compression="gzip"
)

# 拆分 "chrI.10005485" → chrI 和 10005485
df[["chr_clean", "chr_pos"]] = df["chr"].astype(str).str.replace('"', "", regex=False).str.split(".", expand=True)
df["chr"] = df["chr_clean"]
df["start"] = df["chr_pos"].astype(int)

# 分类列
df["DMC"] = df["DMC"].astype(str)

# 去重
df = df[["chr", "start", "DMC"]].drop_duplicates()

print("RSRob 总 CpG 数：")
print(df["DMC"].value_counts())

for cat in cats:
    sub = df[df["DMC"] == cat].copy()
    if sub.empty:
        print(f"[{cat}] 没有位点，跳过")
        continue

    sub = sub.sort_values(["chr", "start"])
    gtf_rows = []

    for chrom, g in sub.groupby("chr"):
        starts = g["start"].tolist()

        island_idx = 1
        island_start = starts[0]
        current_id = f"island{island_idx}_{cat}"

        gene_ids = []
        for pos in starts:
            if pos - island_start >= island_len:
                island_idx += 1
                island_start = pos
                current_id = f"island{island_idx}_{cat}"
            gene_ids.append(current_id)

        tmp = pd.DataFrame({
            "seqname": chrom,
            "source": "none",
            "feature": "transcript",
            "start": starts,
            "end": starts,
            "score": ".",
            "strand": "+",
            "frame": ".",
            "attribute": [f'gene_id "{gid}";' for gid in gene_ids]
        })
        gtf_rows.append(tmp)

    gtf_df = pd.concat(gtf_rows, ignore_index=True)

    outname = f"RSRob_{cat}_islands_5kb.gtf"
    gtf_df.to_csv(
        outname,
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
        escapechar="\\"
    )
    print(f"✅ 写出 {cat} GTF: {outname}, 行数 = {len(gtf_df)}")

cd /work/cyu/meth/RSvsRob

# ========= ROB (400) =========
for cat in non hypo hyper; do
  perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/19_ROB_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsRob/RSRob_${cat}_islands_5kb.gtf \
    --output /work/cyu/meth/RSvsRob/ROB_RSvsROB_pi_${cat}_islands_5k.txt \
    --measure pi
done

# ========= RS (400) =========
for cat in non hypo hyper; do
  perl Variance-at-position.pl --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsRob/RSRob_${cat}_islands_5kb.gtf \
    --output /work/cyu/meth/RSvsRob/RS_RSvsROB_pi_${cat}_islands_5k.txt \
    --measure pi
done

#rswgbswt

#!/usr/bin/env python
import pandas as pd
import csv

# ===== 配置 =====
infile = "RSWGBSWT_CpGs_minDepth5.bed.gz"
island_len = 5000
cats = ["non", "hypo", "hyper"]

# ===== 读入 BED =====
df = pd.read_csv(
    infile,
    sep=" ",
    header=0,
    quotechar='"',
    compression="gzip"
)

# 分类列
if "result_GS" not in df.columns:
    raise ValueError("找不到列 result_GS")
df["DMC"] = df["result_GS"].astype(str)

# chrI.1000020 → chrI
df[["chr_clean", "dummy"]] = (
    df["chr"].astype(str)
      .str.replace('"', "", regex=False)
      .str.split(".", expand=True)
)
df["chr"] = df["chr_clean"]
df["start"] = df["start"].astype(int)

# 去重
df = df[["chr", "start", "DMC"]].drop_duplicates()

print("RSWGBSWT CpG counts:")
print(df["DMC"].value_counts())

# ===== 构建 islands =====
for cat in cats:
    sub = df[df["DMC"] == cat].copy()
    if sub.empty:
        print(f"[{cat}] 无位点，跳过")
        continue

    sub = sub.sort_values(["chr", "start"])
    gtf_rows = []

    for chrom, g in sub.groupby("chr"):
        starts = g["start"].tolist()

        island_idx = 1
        island_start = starts[0]
        current_id = f"island{island_idx}_{cat}"

        gene_ids = []
        for pos in starts:
            if pos - island_start >= island_len:
                island_idx += 1
                island_start = pos
                current_id = f"island{island_idx}_{cat}"
            gene_ids.append(current_id)

        tmp = pd.DataFrame({
            "seqname": chrom,
            "source": "none",
            "feature": "transcript",
            "start": starts,
            "end": starts,
            "score": ".",
            "strand": "+",
            "frame": ".",
            "attribute": [f'gene_id "{gid}";' for gid in gene_ids]
        })
        gtf_rows.append(tmp)

    gtf_df = pd.concat(gtf_rows, ignore_index=True)
    outname = f"RSWGBSWT_{cat}_islands_5kb.gtf"

    gtf_df.to_csv(
        outname,
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
        escapechar="\\"
    )

    print(f"✅ 写出 {outname} ({len(gtf_df)} rows)")

cd /work/cyu/meth/RSvsWat

for cat in non hypo hyper; do
  perl Variance-at-position.pl \
    --pool-size 184 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/7_WT_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsWat/RSWGBSWT_${cat}_islands_5kb.gtf \
    --output /work/cyu/meth/RSvsWat/WGBSWT_RSWGBSWT_pi_${cat}_islands_5k.txt \
    --measure pi
done

for cat in non hypo hyper; do
  perl Variance-at-position.pl \
    --pool-size 400 \
    --min-qual 20 \
    --min-coverage 3 \
    --min-count 2 \
    --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsWat/RSWGBSWT_${cat}_islands_5kb.gtf \
    --output /work/cyu/meth/RSvsWat/RS_RSWGBSWT_pi_${cat}_islands_5k.txt \
    --measure pi
done


#rswik
#!/usr/bin/env python
import pandas as pd
import csv

infile = "RSWGBSWK_CpGs_minDepth5.bed.gz"
island_len = 5000
cats = ["non", "hypo", "hyper"]

df = pd.read_csv(infile, sep=" ", header=0, quotechar='"', compression="gzip")

# 这类文件通常是 DMC 或 result_GS（二选一）
if "DMC" in df.columns:
    df["DMC2"] = df["DMC"].astype(str)
elif "result_GS" in df.columns:
    df["DMC2"] = df["result_GS"].astype(str)
else:
    raise ValueError(f"找不到分类列：DMC 或 result_GS。当前列：{list(df.columns)}")

# chrI.1000020 → chrI
df[["chr_clean", "dummy"]] = (
    df["chr"].astype(str).str.replace('"', "", regex=False).str.split(".", expand=True)
)
df["chr"] = df["chr_clean"]
df["start"] = df["start"].astype(int)

df = df[["chr", "start", "DMC2"]].drop_duplicates()
print("RSWGBSWK counts:")
print(df["DMC2"].value_counts())

for cat in cats:
    sub = df[df["DMC2"] == cat].copy()
    if sub.empty:
        print(f"[{cat}] 无位点，跳过")
        continue

    sub = sub.sort_values(["chr", "start"])
    gtf_rows = []

    for chrom, g in sub.groupby("chr"):
        starts = g["start"].tolist()
        island_idx = 1
        island_start = starts[0]
        current_id = f"island{island_idx}_{cat}"

        gene_ids = []
        for pos in starts:
            if pos - island_start >= island_len:
                island_idx += 1
                island_start = pos
                current_id = f"island{island_idx}_{cat}"
            gene_ids.append(current_id)

        tmp = pd.DataFrame({
            "seqname": chrom,
            "source": "none",
            "feature": "transcript",
            "start": starts,
            "end": starts,
            "score": ".",
            "strand": "+",
            "frame": ".",
            "attribute": [f'gene_id "{gid}";' for gid in gene_ids]
        })
        gtf_rows.append(tmp)

    gtf_df = pd.concat(gtf_rows, ignore_index=True)
    outname = f"RSWGBSWK_{cat}_islands_5kb.gtf"
    gtf_df.to_csv(outname, sep="\t", header=False, index=False,
                  quoting=csv.QUOTE_NONE, escapechar="\\")
    print(f"✅ wrote {outname}: {len(gtf_df)} rows")
cd /work/cyu/meth/RSvsWik

# WK (pool-size 198)
for cat in non hypo hyper; do
  perl Variance-at-position.pl \
    --pool-size 198 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/8_WK_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsWik/RSWGBSWK_${cat}_islands_5kb.gtf \
    --output /work/cyu/meth/RSvsWik/WK_RSWGBSWK_pi_${cat}_islands_5k.txt \
    --measure pi
done

# RS (pool-size 400)
for cat in non hypo hyper; do
  perl Variance-at-position.pl \
    --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/25_RS_noYUn.pileup \
    --gtf /work/cyu/meth/RSvsWik/RSWGBSWK_${cat}_islands_5kb.gtf \
    --output /work/cyu/meth/RSvsWik/RS_RSWGBSWK_pi_${cat}_islands_5k.txt \
    --measure pi
done


#saygos
#!/usr/bin/env python
import pandas as pd
import csv

infile = "SayGos_CpGs_minDepth5.bed.gz"
island_len = 5000
cats = ["non", "hypo", "hyper"]

df = pd.read_csv(infile, sep=" ", header=0, quotechar='"', compression="gzip")

if "DMC" in df.columns:
    df["DMC2"] = df["DMC"].astype(str)
elif "result_GS" in df.columns:
    df["DMC2"] = df["result_GS"].astype(str)
else:
    raise ValueError(f"找不到分类列：DMC 或 result_GS。当前列：{list(df.columns)}")

df[["chr_clean", "dummy"]] = (
    df["chr"].astype(str).str.replace('"', "", regex=False).str.split(".", expand=True)
)
df["chr"] = df["chr_clean"]
df["start"] = df["start"].astype(int)

df = df[["chr", "start", "DMC2"]].drop_duplicates()
print("SayGos counts:")
print(df["DMC2"].value_counts())

for cat in cats:
    sub = df[df["DMC2"] == cat].copy()
    if sub.empty:
        print(f"[{cat}] 无位点，跳过")
        continue

    sub = sub.sort_values(["chr", "start"])
    gtf_rows = []

    for chrom, g in sub.groupby("chr"):
        starts = g["start"].tolist()
        island_idx = 1
        island_start = starts[0]
        current_id = f"island{island_idx}_{cat}"

        gene_ids = []
        for pos in starts:
            if pos - island_start >= island_len:
                island_idx += 1
                island_start = pos
                current_id = f"island{island_idx}_{cat}"
            gene_ids.append(current_id)

        tmp = pd.DataFrame({
            "seqname": chrom,
            "source": "none",
            "feature": "transcript",
            "start": starts,
            "end": starts,
            "score": ".",
            "strand": "+",
            "frame": ".",
            "attribute": [f'gene_id "{gid}";' for gid in gene_ids]
        })
        gtf_rows.append(tmp)

    gtf_df = pd.concat(gtf_rows, ignore_index=True)
    outname = f"SayGos_{cat}_islands_5kb.gtf"
    gtf_df.to_csv(outname, sep="\t", header=False, index=False,
                  quoting=csv.QUOTE_NONE, escapechar="\\")
    print(f"✅ wrote {outname}: {len(gtf_df)} rows")
cd /work/cyu/meth/SayvsGos

# SAY (pool-size 200; pileup=17_SAT_noYUn)
for cat in non hypo hyper; do
  perl Variance-at-position.pl \
    --pool-size 200 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/17_SAT_noYUn.pileup \
    --gtf /work/cyu/meth/SayvsGos/SayGos_${cat}_islands_5kb.gtf \
    --output /work/cyu/meth/SayvsGos/SAY_SayGos_pi_${cat}_islands_5k.txt \
    --measure pi
done

# GOS (pool-size 400)
for cat in non hypo hyper; do
  perl Variance-at-position.pl \
    --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/18_GOS_noYUn.pileup \
    --gtf /work/cyu/meth/SayvsGos/SayGos_${cat}_islands_5kb.gtf \
    --output /work/cyu/meth/SayvsGos/GOS_SayGos_pi_${cat}_islands_5k.txt \
    --measure pi
done

#sayrob
#!/usr/bin/env python
import pandas as pd
import csv

infile = "SayRob_CpGs_minDepth5.bed.gz"
island_len = 5000
cats = ["non", "hypo", "hyper"]

df = pd.read_csv(infile, sep=" ", header=0, quotechar='"', compression="gzip")

if "DMC" in df.columns:
    df["DMC2"] = df["DMC"].astype(str)
elif "result_GS" in df.columns:
    df["DMC2"] = df["result_GS"].astype(str)
else:
    raise ValueError(f"找不到分类列：DMC 或 result_GS。当前列：{list(df.columns)}")

df[["chr_clean", "dummy"]] = (
    df["chr"].astype(str).str.replace('"', "", regex=False).str.split(".", expand=True)
)
df["chr"] = df["chr_clean"]
df["start"] = df["start"].astype(int)

df = df[["chr", "start", "DMC2"]].drop_duplicates()
print("SayRob counts:")
print(df["DMC2"].value_counts())

for cat in cats:
    sub = df[df["DMC2"] == cat].copy()
    if sub.empty:
        print(f"[{cat}] 无位点，跳过")
        continue

    sub = sub.sort_values(["chr", "start"])
    gtf_rows = []

    for chrom, g in sub.groupby("chr"):
        starts = g["start"].tolist()
        island_idx = 1
        island_start = starts[0]
        current_id = f"island{island_idx}_{cat}"

        gene_ids = []
        for pos in starts:
            if pos - island_start >= island_len:
                island_idx += 1
                island_start = pos
                current_id = f"island{island_idx}_{cat}"
            gene_ids.append(current_id)

        tmp = pd.DataFrame({
            "seqname": chrom,
            "source": "none",
            "feature": "transcript",
            "start": starts,
            "end": starts,
            "score": ".",
            "strand": "+",
            "frame": ".",
            "attribute": [f'gene_id "{gid}";' for gid in gene_ids]
        })
        gtf_rows.append(tmp)

    gtf_df = pd.concat(gtf_rows, ignore_index=True)
    outname = f"SayRob_{cat}_islands_5kb.gtf"
    gtf_df.to_csv(outname, sep="\t", header=False, index=False,
                  quoting=csv.QUOTE_NONE, escapechar="\\")
    print(f"✅ wrote {outname}: {len(gtf_df)} rows")
cd /work/cyu/meth/SayvsRob

# SAY (pool-size 200)
for cat in non hypo hyper; do
  perl Variance-at-position.pl \
    --pool-size 200 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/17_SAT_noYUn.pileup \
    --gtf /work/cyu/meth/SayvsRob/SayRob_${cat}_islands_5kb.gtf \
    --output /work/cyu/meth/SayvsRob/SAY_SayRob_pi_${cat}_islands_5k.txt \
    --measure pi
done

# ROB (pool-size 400)
for cat in non hypo hyper; do
  perl Variance-at-position.pl \
    --pool-size 400 --min-qual 20 --min-coverage 3 --min-count 2 --fastq-type sanger \
    --pileup /work/cyu/meth/pooldata/pileup/19_ROB_noYUn.pileup \
    --gtf /work/cyu/meth/SayvsRob/SayRob_${cat}_islands_5kb.gtf \
    --output /work/cyu/meth/SayvsRob/ROB_SayRob_pi_${cat}_islands_5k.txt \
    --measure pi
done


#SAY ROB
  pop    habitat   cat n_islands
1 ROB Freshwater   non      1854
2 SAY     Marine   non      2091
3 ROB Freshwater  hypo         9
4 SAY     Marine  hypo        11
5 ROB Freshwater hyper         6
6 SAY     Marine hyper        18



#!/usr/bin/env Rscript

library(ggplot2)

# =========================
# 通用：读入 + 过滤
# =========================
read_pi <- function(path, pop, habitat, cat) {
  df <- read.table(path, header=FALSE, stringsAsFactors=FALSE)
  colnames(df) <- c("island", "n_valid", "p3", "pi")

  df$n_valid <- suppressWarnings(as.integer(df$n_valid))
  df$pi      <- suppressWarnings(as.numeric(df$pi))

  # 关键过滤：>=2 valid sites
  df <- df[df$n_valid >= 2 & !is.na(df$pi), ]

  df$pop     <- pop
  df$habitat <- habitat
  df$cat     <- cat
  df
}

# =========================
# 通用：画一组对比
# =========================
plot_one_group <- function(base_dir, files, out_prefix, title_suffix="") {

  all_df <- data.frame()

  for (x in files) {
    fpath <- file.path(base_dir, x$file)
    if (!file.exists(fpath)) stop("Missing file: ", fpath)
    tmp <- read_pi(fpath, x$pop, x$habitat, x$cat)
    all_df <- rbind(all_df, tmp)
  }

  # 因子顺序
  all_df$cat     <- factor(all_df$cat, levels=c("non","hypo","hyper"))
  all_df$habitat <- factor(all_df$habitat, levels=c("Freshwater","Marine"))

  # ---- 关键参数：避免 overlap + 适度瘦身（你上次想“瘦点但别过”）----
  dodge_w <- 0.78
  vio_w   <- 0.72
  box_w   <- 0.12

  # 标题（用 unicode π，避免出 pdf 字体问题的话下面 ggsave 已加 useDingbats=FALSE）
  plot_title <- paste0("Per-island pi (5 kb; valid sites \u2265 2)", title_suffix)

  p <- ggplot(all_df, aes(x=cat, y=pi, fill=habitat)) +
    geom_violin(
      position = position_dodge(width = dodge_w),
      width    = vio_w,
      trim     = TRUE,
      scale    = "width",
      alpha    = 0.9,
      linewidth= 0.3
    ) +
    geom_boxplot(
      width = box_w,
      outlier.shape = NA,
      position = position_dodge(width = dodge_w),
      linewidth = 0.3
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Per-island pi",
      fill = "Population"
    ) +
    theme_classic(base_size=14) +
    theme(
      plot.title = element_text(hjust=0.5),
      legend.position="top"
    )

  # 保存到各自 base_dir（注意 path）
  pdf_out <- file.path(base_dir, paste0(out_prefix, ".pdf"))
  png_out <- file.path(base_dir, paste0(out_prefix, ".png"))

  ggsave(pdf_out, p, width=10, height=5, useDingbats=FALSE)
  ggsave(png_out, p, width=10, height=5, dpi=300)

  # QC：每组剩多少 islands
  qc <- aggregate(island ~ pop + habitat + cat, data=all_df, FUN=length)
  colnames(qc)[4] <- "n_islands"
  message("\n[QC] ", out_prefix)
  print(qc)

  invisible(list(plot=p, qc=qc, data=all_df))
}

# ============================================================
# 1) RS vs GOS  (RS = Marine, GOS = Freshwater)
# ============================================================
base_dir <- "/work/cyu/meth/RSvsGos"
files <- list(
  list(pop="RS",  habitat="Marine",     cat="non",   file="RS_RSvsGOS_pi_non_islands_5k.txt"),
  list(pop="RS",  habitat="Marine",     cat="hypo",  file="RS_RSvsGOS_pi_hypo_islands_5k.txt"),
  list(pop="RS",  habitat="Marine",     cat="hyper", file="RS_RSvsGOS_pi_hyper_islands_5k.txt"),
  list(pop="GOS", habitat="Freshwater", cat="non",   file="GOS_RSvsGOS_pi_non_islands_5k.txt"),
  list(pop="GOS", habitat="Freshwater", cat="hypo",  file="GOS_RSvsGOS_pi_hypo_islands_5k.txt"),
  list(pop="GOS", habitat="Freshwater", cat="hyper", file="GOS_RSvsGOS_pi_hyper_islands_5k.txt")
)
plot_one_group(base_dir, files, out_prefix="violin_RSvsGOS_5kb_pi_validNge2")

# ============================================================
# 2) RS vs ROB  (RS = Marine, ROB = Freshwater)
# ============================================================
base_dir <- "/work/cyu/meth/RSvsRob"
files <- list(
  list(pop="RS",  habitat="Marine",     cat="non",   file="RS_RSvsROB_pi_non_islands_5k.txt"),
  list(pop="RS",  habitat="Marine",     cat="hypo",  file="RS_RSvsROB_pi_hypo_islands_5k.txt"),
  list(pop="RS",  habitat="Marine",     cat="hyper", file="RS_RSvsROB_pi_hyper_islands_5k.txt"),
  list(pop="ROB", habitat="Freshwater", cat="non",   file="ROB_RSvsROB_pi_non_islands_5k.txt"),
  list(pop="ROB", habitat="Freshwater", cat="hypo",  file="ROB_RSvsROB_pi_hypo_islands_5k.txt"),
  list(pop="ROB", habitat="Freshwater", cat="hyper", file="ROB_RSvsROB_pi_hyper_islands_5k.txt")
)
plot_one_group(base_dir, files, out_prefix="violin_RSvsROB_5kb_pi_validNge2")

# ============================================================
# 3) RS vs WT (WGBS)  (RS = Marine, WT = Freshwater)
# 你当前输出名：WGBSWT_RSWGBSWT_pi_*_islands_5k.txt  &  RS_RSWGBSWT_pi_*_islands_5k.txt
# ============================================================
base_dir <- "/work/cyu/meth/RSvsWat"
files <- list(
  list(pop="RS", habitat="Marine",      cat="non",   file="RS_RSWGBSWT_pi_non_islands_5k.txt"),
  list(pop="RS", habitat="Marine",      cat="hypo",  file="RS_RSWGBSWT_pi_hypo_islands_5k.txt"),
  list(pop="RS", habitat="Marine",      cat="hyper", file="RS_RSWGBSWT_pi_hyper_islands_5k.txt"),
  list(pop="WT", habitat="Freshwater",  cat="non",   file="WGBSWT_RSWGBSWT_pi_non_islands_5k.txt"),
  list(pop="WT", habitat="Freshwater",  cat="hypo",  file="WGBSWT_RSWGBSWT_pi_hypo_islands_5k.txt"),
  list(pop="WT", habitat="Freshwater",  cat="hyper", file="WGBSWT_RSWGBSWT_pi_hyper_islands_5k.txt")
)
plot_one_group(base_dir, files, out_prefix="violin_RSvsWT_5kb_pi_validNge2")

# ============================================================
# 4) RS vs WK (WGBS)  (RS = Marine, WK = Freshwater)
# 你当前输出名：WK_RSWGBSWK_pi_*_islands_5k.txt  &  RS_RSWGBSWK_pi_*_islands_5k.txt
# ============================================================
base_dir <- "/work/cyu/meth/RSvsWik"
files <- list(
  list(pop="RS", habitat="Marine",      cat="non",   file="RS_RSWGBSWK_pi_non_islands_5k.txt"),
  list(pop="RS", habitat="Marine",      cat="hypo",  file="RS_RSWGBSWK_pi_hypo_islands_5k.txt"),
  list(pop="RS", habitat="Marine",      cat="hyper", file="RS_RSWGBSWK_pi_hyper_islands_5k.txt"),
  list(pop="WK", habitat="Freshwater",  cat="non",   file="WK_RSWGBSWK_pi_non_islands_5k.txt"),
  list(pop="WK", habitat="Freshwater",  cat="hypo",  file="WK_RSWGBSWK_pi_hypo_islands_5k.txt"),
  list(pop="WK", habitat="Freshwater",  cat="hyper", file="WK_RSWGBSWK_pi_hyper_islands_5k.txt")
)
plot_one_group(base_dir, files, out_prefix="violin_RSvsWK_5kb_pi_validNge2")

# ============================================================
# 5) SAY vs GOS  (SAY = Marine, GOS = Freshwater)
# 你当前输出名：SAY_SayGos_pi_*_islands_5k.txt  &  GOS_SayGos_pi_*_islands_5k.txt
# ============================================================
base_dir <- "/work/cyu/meth/SayvsGos"
files <- list(
  list(pop="SAY", habitat="Marine",     cat="non",   file="SAY_SayGos_pi_non_islands_5k.txt"),
  list(pop="SAY", habitat="Marine",     cat="hypo",  file="SAY_SayGos_pi_hypo_islands_5k.txt"),
  list(pop="SAY", habitat="Marine",     cat="hyper", file="SAY_SayGos_pi_hyper_islands_5k.txt"),
  list(pop="GOS", habitat="Freshwater", cat="non",   file="GOS_SayGos_pi_non_islands_5k.txt"),
  list(pop="GOS", habitat="Freshwater", cat="hypo",  file="GOS_SayGos_pi_hypo_islands_5k.txt"),
  list(pop="GOS", habitat="Freshwater", cat="hyper", file="GOS_SayGos_pi_hyper_islands_5k.txt")
)
plot_one_group(base_dir, files, out_prefix="violin_SAYvsGOS_5kb_pi_validNge2")

cat("\nALL DONE\n")






#!/usr/bin/env Rscript




library(ggplot2)

# =========================
# 基本信息
# =========================
base_dir <- "/work/cyu/meth/SayvsRob"

files <- list(
  list(pop="SAY", habitat="Marine",     cat="non",   file="SAY_SayRob_pi_non_islands_5k.txt"),
  list(pop="SAY", habitat="Marine",     cat="hypo",  file="SAY_SayRob_pi_hypo_islands_5k.txt"),
  list(pop="SAY", habitat="Marine",     cat="hyper", file="SAY_SayRob_pi_hyper_islands_5k.txt"),
  list(pop="ROB", habitat="Freshwater", cat="non",   file="ROB_SayRob_pi_non_islands_5k.txt"),
  list(pop="ROB", habitat="Freshwater", cat="hypo",  file="ROB_SayRob_pi_hypo_islands_5k.txt"),
  list(pop="ROB", habitat="Freshwater", cat="hyper", file="ROB_SayRob_pi_hyper_islands_5k.txt")
)

plot_title <- "Per-island pi (5 kb; valid sites \u2265 2)"

# =========================
# 读入 + 过滤
# =========================
read_pi <- function(path, pop, habitat, cat) {
  df <- read.table(path, header=FALSE, stringsAsFactors=FALSE)
  colnames(df) <- c("island", "n_valid", "p3", "pi")

  df$n_valid <- as.integer(df$n_valid)
  df$pi      <- as.numeric(df$pi)

  # 关键过滤：≥2 valid sites
  df <- df[df$n_valid >= 2 & !is.na(df$pi), ]

  df$pop     <- pop
  df$habitat <- habitat
  df$cat     <- cat
  return(df)
}

# =========================
# 合并所有数据
# =========================
all_df <- data.frame()

for (x in files) {
  fpath <- file.path(base_dir, x$file)
  if (!file.exists(fpath)) {
    stop("Missing file: ", fpath)
  }
  tmp <- read_pi(fpath, x$pop, x$habitat, x$cat)
  all_df <- rbind(all_df, tmp)
}

# 因子顺序（非常重要）
all_df$cat <- factor(all_df$cat, levels=c("non","hypo","hyper"))
all_df$habitat <- factor(all_df$habitat, levels=c("Freshwater","Marine"))

# =========================
# 作图参数（关键微调区）
# =========================
dodge_w <- 0.78   # 两组间距：避免贴太近
vio_w   <- 0.72   # violin 宽度：瘦但不成线
box_w   <- 0.12   # box 稍微宽一点增强可读性

p <- ggplot(all_df, aes(x=cat, y=pi, fill=habitat)) +
  geom_violin(
    position = position_dodge(width = dodge_w),
    width    = vio_w,
    trim     = TRUE,
    scale    = "width",
    alpha    = 0.9,
    linewidth= 0.3
  ) +
  geom_boxplot(
    width = box_w,
    outlier.shape = NA,
    position = position_dodge(width = dodge_w),
    linewidth = 0.3
  ) +
  labs(
    title = plot_title,
    x = NULL,
    y = "Per-island pi",
    fill = "Population"
  ) +
  theme_classic(base_size=14) +
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position="top"
  )

# =========================
# 保存
# =========================
ggsave("violin_SayRob_5kb_pi_validNge2.pdf",
       p, width=10, height=5, useDingbats=FALSE)

ggsave("violin_SayRob_5kb_pi_validNge2.png",
       p, width=10, height=5, dpi=300)

# =========================
# QC：每组 island 数量
# =========================
qc <- aggregate(island ~ pop + habitat + cat,
                data = all_df,
                FUN = length)
colnames(qc)[4] <- "n_islands"
print(qc)

cat("DONE\n")

#



R

R version 4.3.3 (2024-02-29) -- "Angel Food Cake"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

[Previously saved workspace restored]

> #!/usr/bin/env Rscript

library(ggplot2)

# =========================
# 通用：读入 + 过滤
# =========================
read_pi <- function(path, pop, habitat, cat) {
  df <- read.table(path, header=FALSE, stringsAsFactors=FALSE)
  colnames(df) <- c("island", "n_valid", "p3", "pi")

  df$n_valid <- suppressWarnings(as.integer(df$n_valid))
  df$pi      <- suppressWarnings(as.numeric(df$pi))

  # 关键过滤：>=2 valid sites
  df <- df[df$n_valid >= 2 & !is.na(df$pi), ]

  df$pop     <- pop
  df$habitat <- habitat
  df$cat     <- cat
  df
}

# =========================
# 通用：画一组对比
cat("\nALL DONE\n")_dir, files, out_prefix="violin_SAYvsGOS_5kb_pi_validNge2")_islands_5k.txt"),,t")

[QC] violin_RSvsGOS_5kb_pi_validNge2
  pop    habitat   cat n_islands
1 GOS Freshwater   non      1026
2  RS     Marine   non      1969
3 GOS Freshwater  hypo        28
4  RS     Marine  hypo        51
5 GOS Freshwater hyper        20
6  RS     Marine hyper       113
Warning messages:
1: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
2: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
3: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
4: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
9: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
+ 
+ 
  Error: unexpected end of input
> #!/usr/bin/env Rscript

library(ggplot2)

# =========================
# 通用：读入 + 过滤
# =========================
read_pi <- function(path, pop, habitat, cat) {
  df <- read.table(path, header=FALSE, stringsAsFactors=FALSE)
  colnames(df) <- c("island", "n_valid", "p3", "pi")

  df$n_valid <- suppressWarnings(as.integer(df$n_valid))
  df$pi      <- suppressWarnings(as.numeric(df$pi))

  # 关键过滤：>=2 valid sites
  df <- df[df$n_valid >= 2 & !is.na(df$pi), ]

  df$pop     <- pop
  df$habitat <- habitat
  df$cat     <- cat
  df
}

# =========================
# 通用：画一组对比
cat("\nALL DONE\n")_dir, files, out_prefix="violin_SAYvsGOS_5kb_pi_validNge2")_islands_5k.txt"),,t")




[QC] violin_RSvsGOS_5kb_pi_validNge2
  pop    habitat   cat n_islands
1 GOS Freshwater   non      1026
2  RS     Marine   non      1969
3 GOS Freshwater  hypo        28
4  RS     Marine  hypo        51
5 GOS Freshwater hyper        20
6  RS     Marine hyper       113
Warning messages:
1: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
2: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
3: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
4: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
9: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
+ 
  Error: unexpected end of input
> 

Save workspace image? [y/n/c]: N

(base) cyu@stickleback:~$ 
(base) cyu@stickleback:~$ cd
(base) cyu@stickleback:~$ cd /work/cyu/meth
(base) cyu@stickleback:/work/cyu/meth$ nano plot_all_5groups_5kb_valid2.R
(base) cyu@stickleback:/work/cyu/meth$ 
(base) cyu@stickleback:/work/cyu/meth$ Rscript plot_all_5groups_5kb_valid2.R

[QC] violin_RSvsGOS_5kb_pi_validNge2
  pop    habitat   cat n_islands
1 GOS Freshwater   non      1026
2  RS     Marine   non      1969
3 GOS Freshwater  hypo        28
4  RS     Marine  hypo        51
5 GOS Freshwater hyper        20
6  RS     Marine hyper       113
Warning messages:
1: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
2: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
3: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
4: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
9: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>

[QC] violin_RSvsROB_5kb_pi_validNge2
  pop    habitat   cat n_islands
1 ROB Freshwater   non      1778
2  RS     Marine   non      2249
3 ROB Freshwater  hypo        19
4  RS     Marine  hypo        23
5 ROB Freshwater hyper        32
6  RS     Marine hyper        84
Warning messages:
1: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
2: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
3: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
4: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
9: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>

[QC] violin_RSvsWT_5kb_pi_validNge2
  pop    habitat   cat n_islands
1  WT Freshwater   non      6321
2  RS     Marine   non      6426
3  WT Freshwater  hypo      1103
4  RS     Marine  hypo      1003
5  WT Freshwater hyper      2650
6  RS     Marine hyper      3042
Warning messages:
1: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
2: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
3: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
4: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
9: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>

[QC] violin_RSvsWK_5kb_pi_validNge2
  pop    habitat   cat n_islands
1  WK Freshwater   non       854
2  RS     Marine   non      1448
3  WK Freshwater  hypo        18
4  RS     Marine  hypo        16
5  WK Freshwater hyper        28
6  RS     Marine hyper        55
Warning messages:
1: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
2: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
3: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
4: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
9: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>

[QC] violin_SAYvsGOS_5kb_pi_validNge2
  pop    habitat   cat n_islands
1 GOS Freshwater   non      1255
2 SAY     Marine   non      1874
3 GOS Freshwater  hypo        11
4 SAY     Marine  hypo        15
5 GOS Freshwater hyper         9
6 SAY     Marine hyper        23
Warning messages:
1: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
2: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
3: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
4: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>
7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <e2>
8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <89>
9: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'Per-island pi (5 kb; valid sites ≥ 2)' in 'mbcsToSbcs': dot substituted for <a5>

ALL DONE
(base) cyu@stickleback:/work/cyu/meth$ 



#violin 发表版 不加ab clean
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# =========================================================
# helper: read and filter one pi file
# =========================================================
read_pi <- function(path, pop, cat) {
  df <- read.table(path, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("island", "n_valid", "p3", "pi")
  
  df$n_valid <- suppressWarnings(as.integer(df$n_valid))
  df$pi      <- suppressWarnings(as.numeric(df$pi))
  
  # keep only islands used in analysis
  df <- df[df$n_valid >= 2 & !is.na(df$pi), ]
  
  df$pop <- pop
  df$cat <- cat
  df
}

# =========================================================
# helper: make one comparison plot + QC + Tukey
# =========================================================
plot_one_group <- function(base_dir, files, out_prefix, pop_order, pop_colors) {
  
  # -------------------------
  # read all files
  # -------------------------
  all_df <- data.frame()
  
  for (x in files) {
    fpath <- file.path(base_dir, x$file)
    if (!file.exists(fpath)) stop("Missing file: ", fpath)
    tmp <- read_pi(fpath, x$pop, x$cat)
    all_df <- rbind(all_df, tmp)
  }
  
  # factor order
  all_df$cat <- factor(all_df$cat, levels = c("non", "hypo", "hyper"))
  all_df$pop <- factor(all_df$pop, levels = pop_order)
  
  # -------------------------
  # QC table
  # -------------------------
  qc <- all_df %>%
    group_by(pop, cat) %>%
    summarise(
      n_islands = n(),
      sum_valid_sites = sum(n_valid, na.rm = TRUE),
      mean_valid_sites_per_island = mean(n_valid, na.rm = TRUE),
      median_pi = median(pi, na.rm = TRUE),
      mean_pi = mean(pi, na.rm = TRUE),
      .groups = "drop"
    )
  
  message("\n[QC] ", out_prefix)
  print(qc)
  
  # -------------------------
  # overall Tukey across pop × category groups
  # -------------------------
  all_df$group <- interaction(all_df$pop, all_df$cat, sep = ".")
  
  fit_all <- aov(pi ~ group, data = all_df)
  tuk_all <- TukeyHSD(fit_all)
  
  tuk_all_df <- as.data.frame(tuk_all$group)
  tuk_all_df$comparison <- rownames(tuk_all_df)
  rownames(tuk_all_df) <- NULL
  tuk_all_df$out_prefix <- out_prefix
  
  message("\n[Tukey overall] ", out_prefix)
  print(tuk_all_df)
  
  # -------------------------
  # within-category Tukey: compare the two populations inside each cat
  # -------------------------
  tuk_cat_list <- list()
  
  for (cc in levels(all_df$cat)) {
    sub <- all_df[all_df$cat == cc, ]
    
    # only run if both pops exist
    if (length(unique(sub$pop)) >= 2) {
      fit_cat <- aov(pi ~ pop, data = sub)
      tuk_cat <- TukeyHSD(fit_cat)
      
      tmp <- as.data.frame(tuk_cat$pop)
      tmp$comparison <- rownames(tmp)
      rownames(tmp) <- NULL
      tmp$cat <- cc
      tmp$out_prefix <- out_prefix
      tuk_cat_list[[cc]] <- tmp
      
      message("\n[Tukey within category: ", cc, "] ", out_prefix)
      print(tmp)
    }
  }
  
  tuk_cat_df <- bind_rows(tuk_cat_list)
  
  # -------------------------
  # plot
  # -------------------------
  dodge_w <- 0.78
  vio_w   <- 0.72
  box_w   <- 0.12
  
  p <- ggplot(all_df, aes(x = cat, y = pi, fill = pop)) +
    geom_violin(
      position  = position_dodge(width = dodge_w),
      width     = vio_w,
      trim      = TRUE,
      scale     = "width",
      alpha     = 0.90,
      linewidth = 0.3
    ) +
    geom_boxplot(
      width = box_w,
      outlier.shape = NA,
      position = position_dodge(width = dodge_w),
      linewidth = 0.3
    ) +
    scale_fill_manual(values = pop_colors, breaks = pop_order) +
    labs(
      title = NULL,
      x = NULL,
      y = "Per-island pi",
      fill = "Population"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_blank()
    )
  
  pdf_out <- file.path(base_dir, paste0(out_prefix, ".pdf"))
  png_out <- file.path(base_dir, paste0(out_prefix, ".png"))
  qc_out  <- file.path(base_dir, paste0(out_prefix, "_QC.tsv"))
  tuk_all_out <- file.path(base_dir, paste0(out_prefix, "_Tukey_all.tsv"))
  tuk_cat_out <- file.path(base_dir, paste0(out_prefix, "_Tukey_withinCategory.tsv"))
  
  # show plot in interactive session
  print(p)
  
  # save plot files
  ggsave(filename = pdf_out, plot = p, width = 10, height = 5, useDingbats = FALSE)
  ggsave(filename = png_out, plot = p, width = 10, height = 5, dpi = 300)
  
  message("Saved plot: ", pdf_out)
  message("Saved plot: ", png_out)
  
  # save tables
  write.table(qc, qc_out, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(tuk_all_df, tuk_all_out, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(tuk_cat_df, tuk_cat_out, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Saved QC table: ", qc_out)
  message("Saved Tukey table: ", tuk_all_out)
  message("Saved Tukey within-category table: ", tuk_cat_out)
  
  invisible(list(
    plot = p,
    qc = qc,
    tukey_all = tuk_all_df,
    tukey_within_cat = tuk_cat_df,
    data = all_df
  ))
}

# =========================================================
# color palette for named populations
# =========================================================
pop_palette <- c(
  "Rabbit Slough" = "#00BFC4",
  "Sayward"       = "#00BFC4",
  "Gosling"       = "#F8766D",
  "Roberts"       = "#F8766D",
  "Watson"        = "#F8766D",
  "Wik"           = "#F8766D"
)

# =========================================================
# 1) RS vs GOS
# =========================================================
base_dir <- "/work/cyu/meth/RSvsGos"
files <- list(
  list(pop = "Rabbit Slough", cat = "non",   file = "RS_RSvsGOS_pi_non_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hypo",  file = "RS_RSvsGOS_pi_hypo_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hyper", file = "RS_RSvsGOS_pi_hyper_islands_5k.txt"),
  list(pop = "Gosling",       cat = "non",   file = "GOS_RSvsGOS_pi_non_islands_5k.txt"),
  list(pop = "Gosling",       cat = "hypo",  file = "GOS_RSvsGOS_pi_hypo_islands_5k.txt"),
  list(pop = "Gosling",       cat = "hyper", file = "GOS_RSvsGOS_pi_hyper_islands_5k.txt")
)
plot_one_group(
  base_dir = base_dir,
  files = files,
  out_prefix = "violin_RSvsGOS_5kb_pi_validNge2_namedPop",
  pop_order = c("Gosling", "Rabbit Slough"),
  pop_colors = pop_palette[c("Gosling", "Rabbit Slough")]
)

# =========================================================
# 2) RS vs ROB
# =========================================================
base_dir <- "/work/cyu/meth/RSvsRob"
files <- list(
  list(pop = "Rabbit Slough", cat = "non",   file = "RS_RSvsROB_pi_non_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hypo",  file = "RS_RSvsROB_pi_hypo_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hyper", file = "RS_RSvsROB_pi_hyper_islands_5k.txt"),
  list(pop = "Roberts",       cat = "non",   file = "ROB_RSvsROB_pi_non_islands_5k.txt"),
  list(pop = "Roberts",       cat = "hypo",  file = "ROB_RSvsROB_pi_hypo_islands_5k.txt"),
  list(pop = "Roberts",       cat = "hyper", file = "ROB_RSvsROB_pi_hyper_islands_5k.txt")
)
plot_one_group(
  base_dir = base_dir,
  files = files,
  out_prefix = "violin_RSvsROB_5kb_pi_validNge2_namedPop",
  pop_order = c("Roberts", "Rabbit Slough"),
  pop_colors = pop_palette[c("Roberts", "Rabbit Slough")]
)

# =========================================================
# 3) RS vs WT
# =========================================================
base_dir <- "/work/cyu/meth/RSvsWat"
files <- list(
  list(pop = "Rabbit Slough", cat = "non",   file = "RS_RSWGBSWT_pi_non_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hypo",  file = "RS_RSWGBSWT_pi_hypo_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hyper", file = "RS_RSWGBSWT_pi_hyper_islands_5k.txt"),
  list(pop = "Watson",        cat = "non",   file = "WGBSWT_RSWGBSWT_pi_non_islands_5k.txt"),
  list(pop = "Watson",        cat = "hypo",  file = "WGBSWT_RSWGBSWT_pi_hypo_islands_5k.txt"),
  list(pop = "Watson",        cat = "hyper", file = "WGBSWT_RSWGBSWT_pi_hyper_islands_5k.txt")
)
plot_one_group(
  base_dir = base_dir,
  files = files,
  out_prefix = "violin_RSvsWT_5kb_pi_validNge2_namedPop",
  pop_order = c("Watson", "Rabbit Slough"),
  pop_colors = pop_palette[c("Watson", "Rabbit Slough")]
)

# =========================================================
# 4) RS vs WK
# =========================================================
base_dir <- "/work/cyu/meth/RSvsWik"
files <- list(
  list(pop = "Rabbit Slough", cat = "non",   file = "RS_RSWGBSWK_pi_non_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hypo",  file = "RS_RSWGBSWK_pi_hypo_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hyper", file = "RS_RSWGBSWK_pi_hyper_islands_5k.txt"),
  list(pop = "Wik",           cat = "non",   file = "WK_RSWGBSWK_pi_non_islands_5k.txt"),
  list(pop = "Wik",           cat = "hypo",  file = "WK_RSWGBSWK_pi_hypo_islands_5k.txt"),
  list(pop = "Wik",           cat = "hyper", file = "WK_RSWGBSWK_pi_hyper_islands_5k.txt")
)
plot_one_group(
  base_dir = base_dir,
  files = files,
  out_prefix = "violin_RSvsWK_5kb_pi_validNge2_namedPop",
  pop_order = c("Wik", "Rabbit Slough"),
  pop_colors = pop_palette[c("Wik", "Rabbit Slough")]
)

# =========================================================
# 5) SAY vs GOS
# =========================================================
base_dir <- "/work/cyu/meth/SayvsGos"
files <- list(
  list(pop = "Sayward", cat = "non",   file = "SAY_SayGos_pi_non_islands_5k.txt"),
  list(pop = "Sayward", cat = "hypo",  file = "SAY_SayGos_pi_hypo_islands_5k.txt"),
  list(pop = "Sayward", cat = "hyper", file = "SAY_SayGos_pi_hyper_islands_5k.txt"),
  list(pop = "Gosling", cat = "non",   file = "GOS_SayGos_pi_non_islands_5k.txt"),
  list(pop = "Gosling", cat = "hypo",  file = "GOS_SayGos_pi_hypo_islands_5k.txt"),
  list(pop = "Gosling", cat = "hyper", file = "GOS_SayGos_pi_hyper_islands_5k.txt")
)
plot_one_group(
  base_dir = base_dir,
  files = files,
  out_prefix = "violin_SAYvsGOS_5kb_pi_validNge2_namedPop",
  pop_order = c("Gosling", "Sayward"),
  pop_colors = pop_palette[c("Gosling", "Sayward")]
)

# =========================================================
# 6) SAY vs ROB
# =========================================================
base_dir <- "/work/cyu/meth/SayvsRob"
files <- list(
  list(pop = "Sayward", cat = "non",   file = "SAY_SayRob_pi_non_islands_5k.txt"),
  list(pop = "Sayward", cat = "hypo",  file = "SAY_SayRob_pi_hypo_islands_5k.txt"),
  list(pop = "Sayward", cat = "hyper", file = "SAY_SayRob_pi_hyper_islands_5k.txt"),
  list(pop = "Roberts", cat = "non",   file = "ROB_SayRob_pi_non_islands_5k.txt"),
  list(pop = "Roberts", cat = "hypo",  file = "ROB_SayRob_pi_hypo_islands_5k.txt"),
  list(pop = "Roberts", cat = "hyper", file = "ROB_SayRob_pi_hyper_islands_5k.txt")
)
plot_one_group(
  base_dir = base_dir,
  files = files,
  out_prefix = "violin_SAYvsROB_5kb_pi_validNge2_namedPop",
  pop_order = c("Roberts", "Sayward"),
  pop_colors = pop_palette[c("Roberts", "Sayward")]
)

cat("\nALL DONE\n")





#violin final
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# =========================================================
# helper: read and filter one pi file
# =========================================================
read_pi <- function(path, pop, cat) {
  df <- read.table(path, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("island", "n_valid", "p3", "pi")
  
  df$n_valid <- suppressWarnings(as.integer(df$n_valid))
  df$pi      <- suppressWarnings(as.numeric(df$pi))
  
  df <- df[df$n_valid >= 2 & !is.na(df$pi), ]
  
  df$pop <- pop
  df$cat <- cat
  df
}

# =========================================================
# 2-group letters
# =========================================================
letters_for_two_groups <- function(group1, group2, pval, alpha = 0.05) {
  if (is.na(pval)) {
    out <- data.frame(group = c(group1, group2), letter = c("a", "a"))
  } else if (pval < alpha) {
    out <- data.frame(group = c(group1, group2), letter = c("a", "b"))
  } else {
    out <- data.frame(group = c(group1, group2), letter = c("a", "a"))
  }
  out
}

# =========================================================
# 3-category letters: non, hypo, hyper
# =========================================================
letters_for_three_groups <- function(p_non_hypo, p_non_hyper, p_hypo_hyper,
                                     alpha = 0.05,
                                     labs = c("non", "hypo", "hyper")) {
  sig12 <- !is.na(p_non_hypo)   && p_non_hypo   < alpha
  sig13 <- !is.na(p_non_hyper)  && p_non_hyper  < alpha
  sig23 <- !is.na(p_hypo_hyper) && p_hypo_hyper < alpha
  
  letters <- c("a", "a", "a")
  names(letters) <- labs
  
  if (!sig12 && !sig13 && !sig23) {
    return(data.frame(group = labs, letter = letters, stringsAsFactors = FALSE))
  }
  
  if (sig12 && sig13 && sig23) {
    letters <- c("a", "b", "c")
    names(letters) <- labs
    return(data.frame(group = labs, letter = letters, stringsAsFactors = FALSE))
  }
  
  if (!sig12 && sig13 && sig23) {
    letters[c("non", "hypo", "hyper")] <- c("a", "a", "b")
    return(data.frame(group = labs, letter = letters, stringsAsFactors = FALSE))
  }
  if (!sig13 && sig12 && sig23) {
    letters[c("non", "hypo", "hyper")] <- c("a", "b", "a")
    return(data.frame(group = labs, letter = letters, stringsAsFactors = FALSE))
  }
  if (!sig23 && sig12 && sig13) {
    letters[c("non", "hypo", "hyper")] <- c("a", "b", "b")
    return(data.frame(group = labs, letter = letters, stringsAsFactors = FALSE))
  }
  
  if (sig12 && !sig13 && !sig23) {
    letters[c("non", "hypo", "hyper")] <- c("a", "b", "ab")
    return(data.frame(group = labs, letter = letters, stringsAsFactors = FALSE))
  }
  if (sig13 && !sig12 && !sig23) {
    letters[c("non", "hypo", "hyper")] <- c("a", "ab", "b")
    return(data.frame(group = labs, letter = letters, stringsAsFactors = FALSE))
  }
  if (sig23 && !sig12 && !sig13) {
    letters[c("non", "hypo", "hyper")] <- c("ab", "a", "b")
    return(data.frame(group = labs, letter = letters, stringsAsFactors = FALSE))
  }
  
  data.frame(group = labs, letter = c("a", "b", "c"), stringsAsFactors = FALSE)
}

# =========================================================
# helper: extract Tukey adjusted p
# =========================================================
get_tukey_p <- function(df, comparison_name) {
  idx <- which(df$comparison == comparison_name)
  if (length(idx) == 0) return(NA_real_)
  df$`p adj`[idx[1]]
}

# =========================================================
# main function
# =========================================================
plot_one_group <- function(base_dir, files, out_prefix, pop_order, pop_colors, label_nudge = NULL) {
  
  # -------------------------
  # read
  # -------------------------
  all_df <- data.frame()
  
  for (x in files) {
    fpath <- file.path(base_dir, x$file)
    if (!file.exists(fpath)) stop("Missing file: ", fpath)
    tmp <- read_pi(fpath, x$pop, x$cat)
    all_df <- rbind(all_df, tmp)
  }
  
  all_df$cat <- factor(all_df$cat, levels = c("non", "hypo", "hyper"))
  all_df$pop <- factor(all_df$pop, levels = pop_order)
  
  pop1_name <- pop_order[1]
  pop2_name <- pop_order[2]
  
  # -------------------------
  # QC
  # -------------------------
  qc <- all_df %>%
    group_by(pop, cat) %>%
    summarise(
      n_islands = n(),
      sum_valid_sites = sum(n_valid, na.rm = TRUE),
      mean_valid_sites_per_island = mean(n_valid, na.rm = TRUE),
      median_pi = median(pi, na.rm = TRUE),
      mean_pi = mean(pi, na.rm = TRUE),
      max_pi = max(pi, na.rm = TRUE),
      .groups = "drop"
    )
  
  message("\n[QC] ", out_prefix)
  print(qc)
  
  # -------------------------
  # Tukey within pop1
  # -------------------------
  sub_pop1 <- all_df[all_df$pop == pop1_name, ]
  fit_pop1 <- aov(pi ~ cat, data = sub_pop1)
  tuk_pop1 <- TukeyHSD(fit_pop1)
  
  tuk_pop1_df <- as.data.frame(tuk_pop1$cat)
  tuk_pop1_df$comparison <- rownames(tuk_pop1_df)
  rownames(tuk_pop1_df) <- NULL
  tuk_pop1_df$pop <- pop1_name
  tuk_pop1_df$out_prefix <- out_prefix
  
  message("\n[Tukey within ", pop1_name, "] ", out_prefix)
  print(tuk_pop1_df)
  
  # -------------------------
  # Tukey within pop2
  # -------------------------
  sub_pop2 <- all_df[all_df$pop == pop2_name, ]
  fit_pop2 <- aov(pi ~ cat, data = sub_pop2)
  tuk_pop2 <- TukeyHSD(fit_pop2)
  
  tuk_pop2_df <- as.data.frame(tuk_pop2$cat)
  tuk_pop2_df$comparison <- rownames(tuk_pop2_df)
  rownames(tuk_pop2_df) <- NULL
  tuk_pop2_df$pop <- pop2_name
  tuk_pop2_df$out_prefix <- out_prefix
  
  message("\n[Tukey within ", pop2_name, "] ", out_prefix)
  print(tuk_pop2_df)
  
  # -------------------------
  # Tukey between pops within category
  # -------------------------
  tuk_cat_list <- list()
  
  for (cc in levels(all_df$cat)) {
    sub <- all_df[all_df$cat == cc, ]
    if (length(unique(sub$pop)) >= 2) {
      fit_cat <- aov(pi ~ pop, data = sub)
      tuk_cat <- TukeyHSD(fit_cat)
      
      tmp <- as.data.frame(tuk_cat$pop)
      tmp$comparison <- rownames(tmp)
      rownames(tmp) <- NULL
      tmp$cat <- cc
      tmp$out_prefix <- out_prefix
      tuk_cat_list[[cc]] <- tmp
      
      message("\n[Tukey between pops within category: ", cc, "] ", out_prefix)
      print(tmp)
    }
  }
  
  tuk_cat_df <- bind_rows(tuk_cat_list)
  
  # -------------------------
  # letters within pop1
  # -------------------------
  p12_1 <- get_tukey_p(tuk_pop1_df, "hypo-non")
  p13_1 <- get_tukey_p(tuk_pop1_df, "hyper-non")
  p23_1 <- get_tukey_p(tuk_pop1_df, "hyper-hypo")
  
  letters_pop1 <- letters_for_three_groups(
    p_non_hypo = p12_1,
    p_non_hyper = p13_1,
    p_hypo_hyper = p23_1,
    labs = c("non", "hypo", "hyper")
  )
  letters_pop1$pop <- pop1_name
  
  # -------------------------
  # letters within pop2
  # -------------------------
  p12_2 <- get_tukey_p(tuk_pop2_df, "hypo-non")
  p13_2 <- get_tukey_p(tuk_pop2_df, "hyper-non")
  p23_2 <- get_tukey_p(tuk_pop2_df, "hyper-hypo")
  
  letters_pop2 <- letters_for_three_groups(
    p_non_hypo = p12_2,
    p_non_hyper = p13_2,
    p_hypo_hyper = p23_2,
    labs = c("non", "hypo", "hyper")
  )
  letters_pop2$pop <- pop2_name
  
  letters_within_pop <- bind_rows(letters_pop1, letters_pop2)
  colnames(letters_within_pop)[1] <- "cat"
  
  # -------------------------
  # letters between populations within category
  # -------------------------
  letters_between_list <- list()
  
  for (cc in levels(all_df$cat)) {
    tmp <- tuk_cat_df[tuk_cat_df$cat == cc, ]
    if (nrow(tmp) > 0) {
      pval <- tmp$`p adj`[1]
      letters2 <- letters_for_two_groups(pop1_name, pop2_name, pval)
      colnames(letters2)[1] <- "pop"
      letters2$cat <- cc
      letters_between_list[[cc]] <- letters2
    }
  }
  
  letters_between_pop <- bind_rows(letters_between_list)
  
  # -------------------------
  # label table
  # -------------------------
  pos_df <- all_df %>%
    group_by(pop, cat) %>%
    summarise(y = max(pi, na.rm = TRUE), .groups = "drop")
  
  label_df <- pos_df %>%
    left_join(letters_within_pop, by = c("pop", "cat")) %>%
    rename(letter_within_pop = letter) %>%
    left_join(letters_between_pop, by = c("pop", "cat")) %>%
    rename(letter_between_pop = letter)
  
  label_df$label <- paste0(
    ifelse(is.na(label_df$letter_within_pop), "", label_df$letter_within_pop),
    "/",
    ifelse(is.na(label_df$letter_between_pop), "", label_df$letter_between_pop)
  )
  
  global_range <- diff(range(all_df$pi, na.rm = TRUE))
  if (global_range == 0) global_range <- 0.01
  
  # default label position: above each violin max
  label_df$y_lab <- label_df$y + pmax(0.12 * label_df$y, 0.06 * global_range)
  
  # manual label-only nudges
  if (!is.null(label_nudge)) {
    label_df <- label_df %>%
      left_join(label_nudge, by = c("pop", "cat")) %>%
      mutate(
        nudge = ifelse(is.na(nudge), 0, nudge),
        y_lab = y_lab + nudge
      ) %>%
      select(-nudge)
  }
  
  # -------------------------
  # plot
  # -------------------------
  dodge_w <- 0.78
  vio_w   <- 0.72
  box_w   <- 0.12
  
  data_max <- max(all_df$pi, na.rm = TRUE)
  
  p <- ggplot(all_df, aes(x = cat, y = pi, fill = pop)) +
    geom_violin(
      position  = position_dodge(width = dodge_w),
      width     = vio_w,
      trim      = TRUE,
      scale     = "width",
      alpha     = 0.90,
      linewidth = 0.3
    ) +
    geom_boxplot(
      width = box_w,
      outlier.shape = NA,
      position = position_dodge(width = dodge_w),
      linewidth = 0.3
    ) +
    geom_text(
      data = label_df,
      aes(x = cat, y = y_lab, label = label, group = pop),
      position = position_dodge(width = dodge_w),
      inherit.aes = FALSE,
      size = 4.5,
      vjust = 0
    ) +
    scale_fill_manual(values = pop_colors, breaks = pop_order) +
    coord_cartesian(
      ylim = c(min(all_df$pi, na.rm = TRUE), data_max * 1.10),
      clip = "off"
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = "Per-island pi",
      fill = "Population"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  pdf_out <- file.path(base_dir, paste0(out_prefix, ".pdf"))
  png_out <- file.path(base_dir, paste0(out_prefix, ".png"))
  qc_out  <- file.path(base_dir, paste0(out_prefix, "_QC.tsv"))
  tuk_pop1_out <- file.path(base_dir, paste0(out_prefix, "_Tukey_", gsub(" ", "_", pop1_name), "_within.tsv"))
  tuk_pop2_out <- file.path(base_dir, paste0(out_prefix, "_Tukey_", gsub(" ", "_", pop2_name), "_within.tsv"))
  tuk_cat_out  <- file.path(base_dir, paste0(out_prefix, "_Tukey_sameCategory_betweenPops.tsv"))
  letters_out  <- file.path(base_dir, paste0(out_prefix, "_letters.tsv"))
  
  print(p)
  
  ggsave(filename = pdf_out, plot = p, width = 10, height = 5, useDingbats = FALSE)
  ggsave(filename = png_out, plot = p, width = 10, height = 5, dpi = 300)
  
  message("Saved plot: ", pdf_out)
  message("Saved plot: ", png_out)
  
  write.table(qc, qc_out, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(tuk_pop1_df, tuk_pop1_out, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(tuk_pop2_df, tuk_pop2_out, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(tuk_cat_df, tuk_cat_out, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(label_df, letters_out, sep = "\t", quote = FALSE, row.names = FALSE)
  
  invisible(list(
    plot = p,
    qc = qc,
    tukey_pop1_within = tuk_pop1_df,
    tukey_pop2_within = tuk_pop2_df,
    tukey_between_pop_same_cat = tuk_cat_df,
    letters = label_df,
    data = all_df
  ))
}

# =========================================================
# colors
# =========================================================
pop_palette <- c(
  "Rabbit Slough" = "#00BFC4",
  "Sayward"       = "#00BFC4",
  "Gosling"       = "#F8766D",
  "Roberts"       = "#F8766D",
  "Watson"        = "#F8766D",
  "Wik"           = "#F8766D"
)

# =========================================================
# 1) RS vs GOS
# =========================================================
base_dir <- "/work/cyu/meth/RSvsGos"
files <- list(
  list(pop = "Rabbit Slough", cat = "non",   file = "RS_RSvsGOS_pi_non_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hypo",  file = "RS_RSvsGOS_pi_hypo_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hyper", file = "RS_RSvsGOS_pi_hyper_islands_5k.txt"),
  list(pop = "Gosling",       cat = "non",   file = "GOS_RSvsGOS_pi_non_islands_5k.txt"),
  list(pop = "Gosling",       cat = "hypo",  file = "GOS_RSvsGOS_pi_hypo_islands_5k.txt"),
  list(pop = "Gosling",       cat = "hyper", file = "GOS_RSvsGOS_pi_hyper_islands_5k.txt")
)
plot_one_group(
  base_dir = base_dir,
  files = files,
  out_prefix = "violin_RSvsGOS_5kb_pi_validNge2_splitLetters",
  pop_order = c("Gosling", "Rabbit Slough"),
  pop_colors = pop_palette[c("Gosling", "Rabbit Slough")]
)

# =========================================================
# 2) RS vs ROB
# =========================================================
base_dir <- "/work/cyu/meth/RSvsRob"
files <- list(
  list(pop = "Rabbit Slough", cat = "non",   file = "RS_RSvsROB_pi_non_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hypo",  file = "RS_RSvsROB_pi_hypo_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hyper", file = "RS_RSvsROB_pi_hyper_islands_5k.txt"),
  list(pop = "Roberts",       cat = "non",   file = "ROB_RSvsROB_pi_non_islands_5k.txt"),
  list(pop = "Roberts",       cat = "hypo",  file = "ROB_RSvsROB_pi_hypo_islands_5k.txt"),
  list(pop = "Roberts",       cat = "hyper", file = "ROB_RSvsROB_pi_hyper_islands_5k.txt")
)
plot_one_group(
  base_dir = base_dir,
  files = files,
  out_prefix = "violin_RSvsROB_5kb_pi_validNge2_splitLetters",
  pop_order = c("Roberts", "Rabbit Slough"),
  pop_colors = pop_palette[c("Roberts", "Rabbit Slough")],
  label_nudge = data.frame(
    pop   = c("Roberts", "Rabbit Slough", "Roberts", "Rabbit Slough"),
    cat   = c("hypo",    "hypo",          "hyper",   "hyper"),
    nudge = c(0.005,     0.02,           0.05,     0.001)
  )
)

# =========================================================
# 3) RS vs WT
# =========================================================
base_dir <- "/work/cyu/meth/RSvsWat"
files <- list(
  list(pop = "Rabbit Slough", cat = "non",   file = "RS_RSWGBSWT_pi_non_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hypo",  file = "RS_RSWGBSWT_pi_hypo_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hyper", file = "RS_RSWGBSWT_pi_hyper_islands_5k.txt"),
  list(pop = "Watson",        cat = "non",   file = "WGBSWT_RSWGBSWT_pi_non_islands_5k.txt"),
  list(pop = "Watson",        cat = "hypo",  file = "WGBSWT_RSWGBSWT_pi_hypo_islands_5k.txt"),
  list(pop = "Watson",        cat = "hyper", file = "WGBSWT_RSWGBSWT_pi_hyper_islands_5k.txt")
)
plot_one_group(
  base_dir = base_dir,
  files = files,
  out_prefix = "violin_RSvsWT_5kb_pi_validNge2_splitLetters",
  pop_order = c("Watson", "Rabbit Slough"),
  pop_colors = pop_palette[c("Watson", "Rabbit Slough")],
  label_nudge = data.frame(
    pop   = c("Watson", "Rabbit Slough", "Watson", "Rabbit Slough"),
    cat   = c("hypo",   "hypo",          "hyper",  "hyper"),
    nudge = c(-0.1,     0.05,            0.04,     0.02)
  )
)

# =========================================================
# 4) RS vs WK
# =========================================================
base_dir <- "/work/cyu/meth/RSvsWik"
files <- list(
  list(pop = "Rabbit Slough", cat = "non",   file = "RS_RSWGBSWK_pi_non_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hypo",  file = "RS_RSWGBSWK_pi_hypo_islands_5k.txt"),
  list(pop = "Rabbit Slough", cat = "hyper", file = "RS_RSWGBSWK_pi_hyper_islands_5k.txt"),
  list(pop = "Wik",           cat = "non",   file = "WK_RSWGBSWK_pi_non_islands_5k.txt"),
  list(pop = "Wik",           cat = "hypo",  file = "WK_RSWGBSWK_pi_hypo_islands_5k.txt"),
  list(pop = "Wik",           cat = "hyper", file = "WK_RSWGBSWK_pi_hyper_islands_5k.txt")
)
plot_one_group(
  base_dir = base_dir,
  files = files,
  out_prefix = "violin_RSvsWK_5kb_pi_validNge2_splitLetters",
  pop_order = c("Wik", "Rabbit Slough"),
  pop_colors = pop_palette[c("Wik", "Rabbit Slough")],
  label_nudge = data.frame(
    pop   = c("Wik", "Rabbit Slough", "Wik", "Rabbit Slough"),
    cat   = c("hypo",   "hypo",          "hyper",  "hyper"),
    nudge = c(-0.2, 0.2, -0.1, 0.12)
  )
)

# =========================================================
# 5) SAY vs GOS
# =========================================================
base_dir <- "/work/cyu/meth/SayvsGos"
files <- list(
  list(pop = "Sayward", cat = "non",   file = "SAY_SayGos_pi_non_islands_5k.txt"),
  list(pop = "Sayward", cat = "hypo",  file = "SAY_SayGos_pi_hypo_islands_5k.txt"),
  list(pop = "Sayward", cat = "hyper", file = "SAY_SayGos_pi_hyper_islands_5k.txt"),
  list(pop = "Gosling", cat = "non",   file = "GOS_SayGos_pi_non_islands_5k.txt"),
  list(pop = "Gosling", cat = "hypo",  file = "GOS_SayGos_pi_hypo_islands_5k.txt"),
  list(pop = "Gosling", cat = "hyper", file = "GOS_SayGos_pi_hyper_islands_5k.txt")
)
plot_one_group(
  base_dir = base_dir,
  files = files,
  out_prefix = "violin_SAYvsGOS_5kb_pi_validNge2_splitLetters",
  pop_order = c("Gosling", "Sayward"),
  pop_colors = pop_palette[c("Gosling", "Sayward")]
)

# =========================================================
# 6) SAY vs ROB
# =========================================================
base_dir <- "/work/cyu/meth/SayvsRob"
files <- list(
  list(pop = "Sayward", cat = "non",   file = "SAY_SayRob_pi_non_islands_5k.txt"),
  list(pop = "Sayward", cat = "hypo",  file = "SAY_SayRob_pi_hypo_islands_5k.txt"),
  list(pop = "Sayward", cat = "hyper", file = "SAY_SayRob_pi_hyper_islands_5k.txt"),
  list(pop = "Roberts", cat = "non",   file = "ROB_SayRob_pi_non_islands_5k.txt"),
  list(pop = "Roberts", cat = "hypo",  file = "ROB_SayRob_pi_hypo_islands_5k.txt"),
  list(pop = "Roberts", cat = "hyper", file = "ROB_SayRob_pi_hyper_islands_5k.txt")
)
plot_one_group(
  base_dir = base_dir,
  files = files,
  out_prefix = "violin_SAYvsROB_5kb_pi_validNge2_splitLetters",
  pop_order = c("Roberts", "Sayward"),
  pop_colors = pop_palette[c("Roberts", "Sayward")]
)

cat("\nALL DONE\n")






#nano /work/shared/cyu/meth/NoSNPs_methyl/run_5kb_pi_all.sh
#!/bin/bash
set -euo pipefail

BASE="/work/shared/cyu/meth/NoSNPs_methyl"
GTF_DIR="${BASE}/gtf"
PI_DIR="${BASE}/pi_out"
PILEUP_DIR="/work/cyu/meth/pooldata/pileup"

mkdir -p "$GTF_DIR" "$PI_DIR"

cd "$BASE"

echo "Step 1: make 5kb GTF from updated BED files"

python make_all_gtf_5kb.py

echo "Step 2: run Variance-at-position.pl for all 5kb islands"

run_pi () {
  POP=$1
  POOLSIZE=$2
  PILEUP=$3
  PREFIX=$4
  GTF_PREFIX=$5

  for CAT in non hypo hyper; do
    echo "Running ${POP} ${PREFIX} ${CAT}"

    perl Variance-at-position.pl \
      --pool-size ${POOLSIZE} \
      --min-qual 20 \
      --min-coverage 3 \
      --min-count 2 \
      --fastq-type sanger \
      --pileup ${PILEUP_DIR}/${PILEUP} \
      --gtf ${GTF_DIR}/${GTF_PREFIX}_${CAT}_islands_5kb.gtf \
      --output ${PI_DIR}/${POP}_${PREFIX}_pi_${CAT}_islands_5k.txt \
      --measure pi
  done
}

# RS vs GOS
run_pi GOS 400 18_GOS_noYUn.pileup RSvsGOS RSGos
run_pi RS  400 25_RS_noYUn.pileup  RSvsGOS RSGos

# RS vs ROB
run_pi ROB 400 19_ROB_noYUn.pileup RSvsROB RSRob
run_pi RS  400 25_RS_noYUn.pileup  RSvsROB RSRob

# RS vs WT
run_pi WGBSWT 184 7_WT_noYUn.pileup RSWGBSWT RSWGBSWT
run_pi RS     400 25_RS_noYUn.pileup RSWGBSWT RSWGBSWT

# RS vs WK
run_pi WK 198 8_WK_noYUn.pileup RSWGBSWK RSWGBSWK
run_pi RS 400 25_RS_noYUn.pileup RSWGBSWK RSWGBSWK

# SAY vs GOS
run_pi SAY 200 17_SAT_noYUn.pileup SayGos SayGos
run_pi GOS 400 18_GOS_noYUn.pileup SayGos SayGos

# SAY vs ROB
run_pi SAY 200 17_SAT_noYUn.pileup SayRob SayRob
run_pi ROB 400 19_ROB_noYUn.pileup SayRob SayRob

echo "Done. Check:"
echo "find ${PI_DIR} -name '*5k.txt' | wc -l"



nano /work/shared/cyu/meth/NoSNPs_methyl/make_all_gtf_5kb.py
#!/usr/bin/env python3

import pandas as pd
import csv
import os

BASE = "/work/shared/cyu/meth/NoSNPs_methyl"
GTF_DIR = os.path.join(BASE, "gtf")
os.makedirs(GTF_DIR, exist_ok=True)

ISLAND_LEN = 5000
CATS = ["non", "hypo", "hyper"]

BED_FILES = {
    "RSGos": "RSGos_CpGs_minDepth5.bed.gz",
    "RSRob": "RSRob_CpGs_minDepth5.bed.gz",
    "RSWGBSWT": "RSWGBSWT_CpGs_minDepth5.bed.gz",
    "RSWGBSWK": "RSWGBSWK_CpGs_minDepth5.bed.gz",
    "SayGos": "SayGos_CpGs_minDepth5.bed.gz",
    "SayRob": "SayRob_CpGs_minDepth5.bed.gz",
}

def make_gtf(prefix, infile):
    path = os.path.join(BASE, infile)
    if not os.path.exists(path):
        print(f"[SKIP] missing {path}")
        return

    df = pd.read_csv(path, sep=r"\s+", header=0, quotechar='"', compression="gzip")

    if "DMC" in df.columns:
        df["DMC2"] = df["DMC"].astype(str)
    elif "result_GS" in df.columns:
        df["DMC2"] = df["result_GS"].astype(str)
    else:
        raise ValueError(f"{infile}: cannot find DMC or result_GS column. Columns: {list(df.columns)}")

    # handle chr like chrI.10005485
    if df["chr"].astype(str).str.contains(r"\.").any():
        split_chr = df["chr"].astype(str).str.replace('"', "", regex=False).str.split(".", expand=True)
        df["chr"] = split_chr[0]
        if "start" not in df.columns or df["start"].isna().all():
            df["start"] = split_chr[1].astype(int)
        else:
            df["start"] = df["start"].astype(int)
    else:
        df["chr"] = df["chr"].astype(str).str.replace('"', "", regex=False)
        df["start"] = df["start"].astype(int)

    df = df[["chr", "start", "DMC2"]].drop_duplicates()

    print(f"\n{prefix} counts:")
    print(df["DMC2"].value_counts())

    for cat in CATS:
        sub = df[df["DMC2"] == cat].copy()
        if sub.empty:
            print(f"[{prefix} {cat}] no sites, skip")
            continue

        sub = sub.sort_values(["chr", "start"])
        rows = []

        for chrom, g in sub.groupby("chr"):
            starts = g["start"].tolist()
            island_idx = 1
            island_start = starts[0]
            current_id = f"island{island_idx}_{cat}"

            for pos in starts:
                if pos - island_start >= ISLAND_LEN:
                    island_idx += 1
                    island_start = pos
                    current_id = f"island{island_idx}_{cat}"

                rows.append([
                    chrom, "none", "transcript", pos, pos,
                    ".", "+", ".", f'gene_id "{current_id}";'
                ])

        out = os.path.join(GTF_DIR, f"{prefix}_{cat}_islands_5kb.gtf")
        pd.DataFrame(rows).to_csv(
            out, sep="\t", header=False, index=False,
            quoting=csv.QUOTE_NONE, escapechar="\\"
        )
        print(f"wrote {out}: {len(rows)} rows")

for prefix, infile in BED_FILES.items():
    make_gtf(prefix, infile)
    