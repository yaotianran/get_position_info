# get_position_info，一个检查bam文件特定位点的相关信息的脚本

当你拿到一个小panel富集测序结果，或者全外测序结果，你可能只对某些特定的位点或者区域感兴趣。你需要收集特定文件的一些信息，例如测序覆盖度，碱基组成，是否和已知位点信息相符等等，你可以使用这个脚本来检测bam文件特定位点的信息并且输出成tsv文本格式。

输出示例：
[![d7NN54n.png](https://iili.io/d7NN54n.png)](https://freeimage.host/i/d7NN54n)

---
### 1，安装

get_position_info是一个python脚本，所以你必须安装python。除此之外必须安装pysam。
脚本本身无需安装，直接运行

---
### 2，基础用法

`get_position_info.py [bam_file] [locus_file]`

其中 bam_file 是一个经过排序和index的bam文件。locus_file是一个包含位置信息的文件，支持三种格式（BED，VCF和POS）。其中POS格式是一个双列的文本文件，第一列为染色体名称，第二列为位置。而BED格式参考[这里](https://asia.ensembl.org/info/website/upload/bed.html)，VCF格式参考[这里](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

无论哪种格式，都会忽略空行和以'#'开头的行

---
### 3，比较测序结果中位点与标准位点基因型的差异

当你测定一些标准样本，例如GIAB的HG001 ～ HG007的样本，它们提供了一系列已知基因型的标准位点，或者你已经基于已有信息的构建了你自己的标准位点，你可以使用这些信息，get_position_info将自动计算bam中相关位点是否与标准位点的基因型相符。例如：

`get_position_info.py -v <vcf_file> [bam_file] [locus_file]`

输出的文件会附加几列，表示有多少碱基和InDel与标准位点的基因型相符

---
### 4，获取位点在基因组上下游的序列

当你发现某些位点的测序结果与标准基因型不符时，你可能想知道它们的基因组上下游序列，以便确定它们是否在重复序列区域或者存在以下特定的上下游pattern。这时你可以使用-r选项和-c选项，get_position_info会输出位点在参考基因组上下游的序列。例如：

`get_position_info.py -r <reference.fasta> -c 5 [bam_file] [locus_file]`

这里reference.fasta必须使用samtools faidx命令index，-c 5 表示 分别输出上下游5bp，加上位点本身，共11bp序列

---
### 5，使用 -f / --format选项输出位点的额外信息

有时你需要位点的额外信息，例如位点的cycle，测序质量或者read的MAPQ值等等，你可以使用-f或--format选项获取额外的信息。当然有些信息只有同时使用-v或者-r选项时才有意义。例如：

`get_position_info.py -f 'A_seq_quality, T_seq_quality, C_seq_quality, G_seq_quality' [bam_file] [locus_file]`

或者

`get_position_info.py -v <vcf_file> -f 'matched_snp_seq_quality, unmatched_snp_seq_quality' [bam_file] [locus_file]`

目前支持的信息有：
```
# 位点所在位置各种碱基的数量，测序质量，所在read的MAPQ值，cycle循环数
# X为A T C G N miss中的一种
X_count: List[int, int, int, int] 
X_seq_quality: List[int, ...]
X_MAPQ: List[int, ...]
X_cycle: List[int, ...] 

# 位点所在位置与标准位点基因型相同或不相同的base类型的数量，测序质量，所在read的MAPQ值，cycle循环数
matched_snp_count: List[int, int, int, int]
unmatched_snp_count: List[int, int, int, int]
matched_snp_seq_quality: List[int, ...]
unmatched_snp_seq_quality: List[int, ...]
matched_snp_MAPQ: List[int, ...]
unmatched_snp_MAPQ: List[int, ...]
matched_snp_cycle: List[int, ...]
unmatched_snp_cycle: List[int, ...]

# 位点所在位置的insert的数量，测序质量，所在read的MAPQ值，cycle循环数
ins_count: List[int, int, int, int]
ins_seq_quality: List[int, ...]
ins_MAPQ: List[int, ...]
ins_cycle: List[int, ...]

# 位点所在位置的delete的数量，所在read的MAPQ值，cycle循环数
del_count: List[int, int, int, int]
del_MAPQ: List[int, ...]
del_cycle: List[int, ...]

# 位点所在位置与标准位点基因型相同或不相同的InDel类型的数量，测序质量，所在read的MAPQ值，cycle循环数
matched_indel_count: List[int, int, int, int]
unmatched_indel_count: List[int, int, int, int]
matched_ins_seq_quality: List[int, ...]
unmatched_ins_seq_quality: List[int, ...]
matched_indel_MAPQ: List[int, ...]
unmatched_indel_MAPQ: List[int, ...]
matched_indel_cycle: List[int, ...]
unmatched_indel_cycle: List[int, ...]

# 以上各项指标的均值，其中X为A T C G N miss中的一种
X_mean_seq_quality
ins_mean_seq_quality
matched_snp_mean_seq_quality
unmatched_snp_mean_seq_quality
matched_ins_mean_seq_quality
unmatched_ins_mean_seq_quality
X_mean_MAPQ 
del_mean_MAPQ
ins_mean_MAPQ
matched_snp_mean_MAPQ
unmatched_snp_mean_MAPQ
matched_indel_mean_MAPQ
unmatched_indel_mean_MAPQ
X_mean_cycle 
del_mean_cycle
ins_mean_cycle
matched_snp_mean_cycle
unmatched_snp_mean_cycle
matched_indel_mean_cycle
unmatched_indel_mean_cycle

# 位点所在位置的InDel的长度数量统计，正数为insert，负数为delete，0为无插入缺失
indel_length_counter
```

---
### 6，FAQs

- Q：为什么在X_count列不是一个整数，而是四个整数？<br/>
  A：X_count列的的格式为四个以逗号分割的整数，它们依次表示forward 1st read, forward 2nd read, reverse 1st read, reverse 2nd read。如果是单端测序，则forward 2nd read和reverse 2nd read都为0。将不同方向的reads数单独列出，可以帮助识别由一些PCR或者上下游序列造成的测序错误。
<br/>

- Q：miss和delete都是什么意思？<br/>
  A：miss指在查询位点所在位置有read覆盖，但是位于read的delete内部，而delete和insert是指在VCF格式中，发生在查询位点所在位置后1bp的InDel。
<br/>

- Q：insert和delete所在的cycle是指哪个碱基？<br/>
  A：insert和delete所在的cycle，如果在正向read中，为pysam.PileupRead.query_position_or_next 属性 + 1 。如果在反向read中，为read长度 - pysam.PileupRead.query_position_or_next
<br/>

- Q：怎样确定输出结果是正确的？<br/>
  A：可以使用[IGV](https://igv.org/)手动查看位点信息。IGV的坐标是1-based，所以建议使用1-based的POS或VCF格式输入位点位置。

  
