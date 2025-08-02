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
    --pool-size 184 \
    --fastq-type sanger

# 3. 计算 genome-wide 平均 π
awk '{sum += $5; n++} END {if (n > 0) print "Background pi (WT):", sum / n; else print "No data"}' /work/cyu/meth/WT_10kb_nostep_pi.txt
Background pi (WT): 0.00468986



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





