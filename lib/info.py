# 点位信息获取和输出的相关函数
import pysam
from dataclasses import dataclass, field
import statistics
import collections
import sys
from . import utils


BASES = ['A', 'T', 'C', 'G']


@dataclass
class PositionInfo:

   # =====================count=====================
   A_count: list = field(default_factory=list)  # 4 length integer list, 分别表示此核苷酸种类在F1， F2, R1, R2 flag reads 上的数量
   T_count: list = field(default_factory=list)
   C_count: list = field(default_factory=list)
   G_count: list = field(default_factory=list)
   N_count: list = field(default_factory=list)
   miss_count: list = field(default_factory=list) # in-situ deletion
   del_count: list = field(default_factory=list) # next-site deletion
   ins_count: list = field(default_factory=list)

   matched_snp_count : list = field(default_factory=list)
   unmatched_snp_count : list = field(default_factory=list)
   matched_indel_count : list = field(default_factory=list)
   unmatched_indel_count : list = field(default_factory=list)

   background_count : list = field(default_factory=list)

   # =====================seq_quality=====================
   A_seq_quality: list = field(default_factory=list) # A 的测序质量，长度等于A的数量
   T_seq_quality: list = field(default_factory=list)
   C_seq_quality: list = field(default_factory=list)
   G_seq_quality: list = field(default_factory=list)
   N_seq_quality: list = field(default_factory=list)
   ins_seq_quality: list = field(default_factory=list) # 记录全部插入序列碱基的测序质量

   matched_snp_seq_quality : list = field(default_factory=list)
   unmatched_snp_seq_quality : list = field(default_factory=list)
   matched_ins_seq_quality : list = field(default_factory=list)
   unmatched_ins_seq_quality : list = field(default_factory=list)

   # =====================MAPQ=====================
   A_MAPQ: list = field(default_factory=list) # A所在read的MAPQ质量，长度等于A的数量
   T_MAPQ: list = field(default_factory=list)
   C_MAPQ: list = field(default_factory=list)
   G_MAPQ: list = field(default_factory=list)
   N_MAPQ: list = field(default_factory=list)
   miss_MAPQ: list = field(default_factory=list)
   del_MAPQ: list = field(default_factory=list)
   ins_MAPQ: list = field(default_factory=list)

   matched_snp_MAPQ : list = field(default_factory=list)
   unmatched_snp_MAPQ : list = field(default_factory=list)
   matched_indel_MAPQ : list = field(default_factory=list)
   unmatched_indel_MAPQ : list = field(default_factory=list)


   # =====================Cycle=====================
   A_cycle: list = field(default_factory=list)  # A所在read的cycle, 长度等于A的数量
   T_cycle: list = field(default_factory=list)
   C_cycle: list = field(default_factory=list)
   G_cycle: list = field(default_factory=list)
   N_cycle: list = field(default_factory=list)
   miss_cycle: list = field(default_factory=list)
   del_cycle: list = field(default_factory=list) # 只记录发生InDel的起始位置的碱基的cycle
   ins_cycle: list = field(default_factory=list)

   matched_snp_cycle: list = field(default_factory=list)
   unmatched_snp_cycle: list = field(default_factory=list)
   matched_indel_cycle: list = field(default_factory=list)
   unmatched_indel_cycle: list = field(default_factory=list)

   # =====================其他=====================
   chrom: str = None
   pos: int = None
   coverage: int = None  # 覆盖度，deletion计算在内

   indel_length: list = field(default_factory=list) # read在该位点后面的indel长度, 计数结果。来之[pileup_read.indel for pileup_read in pileupcolumn.pileups]
   query_snp: list = field(default_factory=list) # 该位点所在位置的序列的种类和数量， 包括ATCG*,其中*表示缺失
   query_indel: list = field(default_factory=list) # 该位点所在位置后接InDel的序列的种类和数量，以samtools的格式表示，例如‘+2AC’, '-3NNN'等

   # List中的每个元素是一个allele
   real_allele_snp: list = field(default_factory=list)   # 原位allele，例如ATCG或者*
   real_allele_indel: list = field(default_factory=list)   # Indel allele，例如‘+2AC’, '-3NNN'等等。注意*和InDel不能出现在同一位点

   reference: str = None   # 该位置的ref碱基
   context: str = None # 该位置上下游的base
   other: str = None  # 其他信息



# 输入一个PositionInfo对象和它的一个属性，首先将这个属性的字符串转化为人类易懂的格式，然后输出
def __output_attr(pos_info:PositionInfo, attribute: str) -> str:

   '''
   提取一个PositionInfo对象的一个属性，转换为适用于打印的字符串

   Parameters:
   **pos_info**: PositionInfo类
      一个PositionInfo类对象

   **attribute**: str
      属性字符串

   Returns:
       **output**: string
       适用于打印的字符串

   '''
   try:
      value = getattr(pos_info, attribute)
   except AttributeError:
      message = f'__output_attr： PositionInfo 对象不存在 {attribute} 属性'
      print(message)
      sys.exit(message)

   if isinstance(value, list):
      if value == []:
         return ''

      if isinstance(value[0], int): # [12,4,2,9]
         if value == [0, 0, 0, 0]:
            return ''
         else:
            return str(value)[1:-1]
      else:   #  ['A', 'A', '+2AC']
         value = list(set(value))
         return str(value)[1:-1].replace("'", "")

   if isinstance(value, int):
      return str(value)

   if isinstance(value, float):
      return str(round(value, 1))

   if isinstance(value, collections.Counter):
      return str(value)[9:-2].replace("'", '')

   if value is None:
      return ''

   return str(value)


# 输入一个PositionInfo对象和它的一个属性列表，然后输出
def output_attributes_pos_info(pos_info:PositionInfo, attributes: list[str, ...]) -> str:
   '''
   输入一个PositionInfo对象和它的一个属性列表，然后输出

   Parameters:
      **pos_info**: PositionInfo
         一个PositionInfo对象

      **attributes**: list[str, ...]
         列表，每一个item是一个属性字符串

   Returns:
       **output_str**: str
         输出的字符串
   '''
   output_str = ''
   for attr_str in attributes:
      if not isinstance(attr_str, str):
         message = f'output_attributes_pos_info: attribute {attr} ({type(attr)}) is not a string.'
         sys.exit(message)

      output_str += '\t' + __output_attr(pos_info, attr_str)

   return output_str.strip()

# 输入一个PositionInfo对象，利用副作用，设置PositionInfo对象内部的一些其他属性
# 设置了以下属性：
# X_mean_seq_quality    # X in 'A T C G N miss
# ins_mean_seq_quality
# matched_snp_mean_seq_quality
# unmatched_snp_mean_seq_quality
# matched_ins_mean_seq_quality
# unmatched_ins_mean_seq_quality

# X_mean_MAPQ    # X in 'A T C G N miss
# del_mean_MAPQ
# ins_mean_MAPQ
# matched_snp_mean_MAPQ
# unmatched_snp_mean_MAPQ
# matched_indel_mean_MAPQ
# unmatched_indel_mean_MAPQ

# X_mean_cycle    # X in 'A T C G N miss
# del_mean_cycle
# ins_mean_cycle
# matched_snp_mean_cycle
# unmatched_snp_mean_cycle
# matched_indel_mean_cycle
# unmatched_indel_mean_cycle

# indel_length_counter
# query_snp_counter
# query_indel_counter
def add_attributes_pos_info(pos_info:PositionInfo) -> int:
   '''
   输入一个PositionInfo对象，利用副作用，设置PositionInfo对象内部的一些其他属性
   '''

   for base in ['A', 'T', 'C', 'G', 'N', 'miss', 'matched_snp', 'unmatched_snp', 'matched_ins', 'unmatched_ins', 'matched_indel', 'unmatched_indel']:
      for attr in ['seq_quality', 'MAPQ', 'cycle']:
         try:
            origial_attr_lst = getattr(pos_info, f'{base}_{attr}')
         except AttributeError:
            continue
         except Exception as ex:
            print('add_attributes_pos_info:', ex, pos_info)
            return(pos_info)

         if origial_attr_lst == []:
            setattr(pos_info, f'{base}_mean_{attr}', None)
         else:
            try:
               setattr(pos_info, f'{base}_mean_{attr}', statistics.mean(origial_attr_lst))
            except Exception as ex:
               setattr(pos_info, f'{base}_mean_{attr}', None)
               print('add_attributes_pos_info:', ex, f'{base}_mean_{attr}', origial_attr_lst)

   pos_info.indel_length_counter = collections.Counter(pos_info.indel_length) if pos_info.indel_length != [] else None
   pos_info.query_snp_counter = collections.Counter(pos_info.query_snp) if pos_info.query_snp != [] else None
   pos_info.query_indel_counter = collections.Counter(pos_info.query_indel) if pos_info.query_indel != [] else None

   return 0


# 提取给定的位置bam文件的比对信息
# allele为一个List，其中的每个元素是一个allele，tuple格式，包含alt 和对应的genotype。例如（‘A’，‘0|1’） （'*', '1|1'）或者 ('+2AC', '1|1')
# 返回一个PositionInfo对象
# 设置了以下属性：
# coverage: int
# background: List[int, int, int, int]
# indel_length: List[int, ...]
# query_snp: List[str, ...]
# query_indel: List[str, ...]

# X_count: List[int, int, int, int]   # X in 'A T C G N miss'
# X_seq_quality: List[int, ...]   # X in 'A T C G N'
# X_MAPQ: List[int, ...]   # X in 'A T C G N miss'
# X_cycle: List[int, ...]   # X in 'A T C G N miss'

# matched_snp_count: List[int, int, int, int]
# unmatched_snp_count: List[int, int, int, int]

# matched_snp_seq_quality: List[int, ...]
# unmatched_snp_seq_quality: List[int, ...]

# matched_snp_MAPQ: List[int, ...]
# unmatched_snp_MAPQ: List[int, ...]

# matched_snp_cycle: List[int, ...]
# unmatched_snp_cycle: List[int, ...]

# ins_count: List[int, int, int, int]
# ins_seq_quality: List[int, ...]
# ins_MAPQ: List[int, ...]
# ins_cycle: List[int, ...]

# del_count: List[int, int, int, int]
# del_MAPQ: List[int, ...]
# del_cycle: List[int, ...]

# matched_indel_count: List[int, int, int, int]
# unmatched_indel_count: List[int, int, int, int]

# matched_ins_seq_quality: List[int, ...]
# unmatched_ins_seq_quality: List[int, ...]

# matched_indel_MAPQ: List[int, ...]
# unmatched_indel_MAPQ: List[int, ...]

# matched_indel_cycle: List[int, ...]
# unmatched_indel_cycle: List[int, ...]
def get_pos_info(bam_af: pysam.AlignmentFile, chrom: str, pos: int, real_allele_snp: tuple[str, ...] = None, real_allele_indel: tuple[str, ...] = None) -> PositionInfo:
   '''
   提取位点信息，包括位点深度，四种碱基read数（百分比），四种碱基平均测序质量，四种碱基的正反向数量，四种碱基orientation数量等

   Parameters:
      **bam_af**: pysam.AlignmentFile
         一个pysam.AlignmentFile对象，bam_af = pysam.AlignmentFile(bamfile, 'rb')

      **chrom**: str
         染色体

      **pos**: int
         位置（1-based）,与IGV兼容

      **real_allele_snp, real_allele_indel**: tuple
         可选, 用于计算各项matched 和unmatched指标
         其中的每个元素是一个allele，string格式，通常tuple的长度为该物种的染色体倍性（Ploidy），如人类为2。
         real_allele_snp表示该位点所在位置的序列的种类和数量，碱基用ATCG， miss用*表示。
         real_allele_indel表示该位点所在位置后接InDel的序列的种类和数量，以samtools的格式表示，例如‘+2AC’, '-3NNN'等等。

   Returns:
       **PositionInfo**: class
         PositionInfo类
   '''

   result_pos = PositionInfo()


   result_pos.A_count = [0, 0, 0, 0]
   result_pos.T_count = [0, 0, 0, 0]
   result_pos.C_count = [0, 0, 0, 0]
   result_pos.G_count = [0, 0, 0, 0]
   result_pos.N_count = [0, 0, 0, 0]
   result_pos.miss_count = [0, 0, 0, 0]
   result_pos.del_count = [0, 0, 0, 0]
   result_pos.ins_count = [0, 0, 0, 0]

   result_pos.matched_snp_count = [0, 0, 0, 0]
   result_pos.unmatched_snp_count = [0, 0, 0, 0]
   result_pos.matched_indel_count = [0, 0, 0, 0]
   result_pos.unmatched_indel_count = [0, 0, 0, 0]
   result_pos.background_count = [0, 0, 0, 0]


   result_pos.chrom = chrom
   result_pos.pos = pos
   result_pos.real_allele_indel = real_allele_indel
   result_pos.real_allele_snp = real_allele_snp

   for pileupcolumn in bam_af.pileup(contig = chrom, start = pos - 1, stop = pos, truncate = True, max_depth = 999999999, stepper = 'nofilter', ignore_overlaps=False, ignore_orphans = False, min_base_quality=0):

      result_pos.coverage = pileupcolumn.get_num_aligned()
      query_indel_lst = [x.upper()[1:] for x in pileupcolumn.get_query_sequences(add_indels = True)] # 用来在后面获取InDel read的query string ['', '', '', '', '+2AC', '-4NNN']
      # query_snp_lst = [x.upper()[0] for x in pileupcolumn.get_query_sequences(add_indels = True)]

      #temp = []
      #for x in query_indel_lst:
         #if x != '':
            #temp.append(x)
      #print(collections.Counter(query_snp_lst))
      #print(collections.Counter(temp))


      for i, pileup_read in enumerate(pileupcolumn.pileups):
         segment = pileup_read.alignment
         if pileup_read.is_refskip:
            message = 'get_pos_info：read {} is_refskip 为真(flag {})，忽略此read（is_forward:{}, is_reverse:{}, is_read1:{}, is_read2:{}'.format(segment.query_name, segment.flag, segment.is_forward, segment.is_reverse, segment.is_read1, segment.is_read2)
            print(message)
            continue

         flag_index_int = utils.get_index(segment)
         if flag_index_int is None:
            message = 'get_pos_info：无法判断read {} 方向(flag {})，忽略此read（is_forward:{}, is_reverse:{}, is_read1:{}, is_read2:{}'.format(segment.query_name, segment.flag, segment.is_forward, segment.is_reverse, segment.is_read1, segment.is_read2)
            print(message)
            continue


         result_pos.background_count[flag_index_int] += 1
         result_pos.indel_length.append(pileup_read.indel) #  indel length (- 0 +)for the position following the current pileup site.

         mapq_int = segment.mapping_quality

         if segment.is_forward:
            cycle_int = pileup_read.query_position_or_next + 1 # 当前位置在read上面的cycle数，如果是miss，是下一个碱基的cycle。+1 为了将0-based转换成1-based
         else:
            cycle_int = segment.infer_read_length() - pileup_read.query_position_or_next

         # = = = = = = = = = = = = = = 获取并设置当前位置（snp）的信息 = = = = = = = = = = = = = =
         if pileup_read.query_position is not None: #  当前位置不是miss
            base = segment.query_sequence[pileup_read.query_position]
            result_pos.query_snp.append(base)
            seq_quality_int = segment.get_forward_qualities()[pileup_read.query_position]
         else:   #  当前位置为miss
            base = 'miss'
            result_pos.query_snp.append('*')
            # seq_quality_int = -1

         # set count
         attr_str = base + '_count'
         temp_lst = getattr(result_pos, attr_str)
         temp_lst[flag_index_int] += 1
         setattr(result_pos, attr_str, temp_lst)

         # set seq quality
         attr_str = base + '_seq_quality'
         if base != 'miss':  #  只设置 非miss reads的测序质量
            temp_lst = getattr(result_pos, attr_str)
            temp_lst.append(seq_quality_int)
            setattr(result_pos, attr_str, temp_lst)

         # set map quality
         attr_str = base + '_MAPQ'
         temp_lst = getattr(result_pos, attr_str)
         temp_lst.append(mapq_int)
         setattr(result_pos, attr_str, temp_lst)

         # set cycle
         attr_str = base + '_cycle'
         temp_lst = getattr(result_pos, attr_str)
         temp_lst.append(cycle_int)
         setattr(result_pos, attr_str, temp_lst)

         # set matched and unmatched
         if real_allele_snp is not None and base in real_allele_snp: #  matched
            result_pos.matched_snp_count[flag_index_int] += 1
            result_pos.matched_snp_MAPQ.append(mapq_int)
            result_pos.matched_snp_cycle.append(cycle_int)
            if base != 'miss':  # 只设置非miss的matched_snp_seq_quality
               result_pos.matched_snp_seq_quality.append(seq_quality_int)
         elif real_allele_snp is not None and base not in real_allele_snp: #  matched:
            result_pos.unmatched_snp_count[flag_index_int] += 1
            result_pos.unmatched_snp_MAPQ.append(mapq_int)
            result_pos.unmatched_snp_cycle.append(cycle_int)
            if base != 'miss':  # 只设置非miss的unmatched_snp_seq_quality
               result_pos.unmatched_snp_seq_quality.append(seq_quality_int)
         else:
            pass


         # = = = = = = = = = = = = = = 获取并设置当前位置后接InDel的信息 = = = = = = = = = = = = = =
         if pileup_read.indel != 0:  # 如果后接位置存在InDel
            indel_alt_str = query_indel_lst[i]
            result_pos.query_indel.append(indel_alt_str)

            if pileup_read.indel > 0: # 后方有插入
               result_pos.ins_count[flag_index_int] += 1 #  count
               seq_quality_lst = list(segment.query_qualities)[pileup_read.query_position + 1:pileup_read.query_position + pileup_read.indel + 1]
               result_pos.ins_seq_quality.extend(seq_quality_lst)
               result_pos.ins_MAPQ.append(mapq_int)
               result_pos.ins_cycle.append(cycle_int)

               # indel_alt_str = segment.query_sequence[pileup_read.query_position + 1: pileup_read.query_position + pileup_read.indel + 1]
               # indel_alt_str = '+' + str(len(indel_alt_str)) + indel_alt_str

            if pileup_read.indel < 0: # 后方有缺失
               result_pos.del_count[flag_index_int] += 1
               seq_quality_lst = []
               result_pos.del_MAPQ.append(mapq_int)
               result_pos.del_cycle.append(cycle_int)

               # indel_alt_str = str(pileup_read.indel) + 'N' * abs(pileup_read.indel)

            # set matched and unmatched
            if real_allele_indel is not None and indel_alt_str in real_allele_indel:
               result_pos.matched_indel_count[flag_index_int] += 1
               result_pos.matched_ins_seq_quality.extend(seq_quality_lst)
               result_pos.matched_indel_MAPQ.append(mapq_int)
               result_pos.matched_indel_cycle.append(cycle_int)

            elif real_allele_indel is not None and indel_alt_str not in real_allele_indel:
               result_pos.unmatched_indel_count[flag_index_int] += 1
               result_pos.unmatched_ins_seq_quality.extend(seq_quality_lst)
               result_pos.unmatched_indel_MAPQ.append(mapq_int)
               result_pos.unmatched_indel_cycle.append(cycle_int)

            else:
               pass

   return result_pos


if __name__ == '__main__':

   import os.path as path
   import copy

   bam_file = '~/server/result/sequencer/for_partner/zgbio/snvindelLOD_20240717/qua25_unqua10_len50/bam/sample_03/sample_03.filter.bam'
   bam_file = '~/server/result/sequencer/salus/giab/hg003_na24149_father/wgs/pro63_modelopt_20240705/qua25_unqua10_len50/rmdup/hg003_na24149_modelopt/hg003_na24149_modelopt.rmdup.sorted.bam'
   # bam_file = '~/server/result/sequencer/for_partner/sino_us/read2_lowquality_20240724/abnormal_data_SL200/AC340.merge.bam'
   bam_file = path.realpath(path.expanduser(bam_file))
   bam_af = pysam.AlignmentFile(bam_file, 'rb')

   reference = '~/server/data/reference/human/hg38/hg38.fasta'
   reference = '~/server/data/reference/human/hg19/hg19.fasta'
   reference = path.realpath(path.expanduser(reference))
   fasta = pysam.FastaFile(reference)

   chrom = 'chr14'
   pos = 23183147

   chrom = 'chr17'
   pos = 39724814
   # pos = 39724821
   # pos = 39724822

   chrom = '1'
   # pos = 790696  # HG003  1    790696  rs76089329      C       CAT     50      PASS    platforms=4;platformnames=PacBio,Illumina,10X,CG;datase 纯合
   pos = 790758  # HG003    1       790758  rs72558500      GTA     G       50      PASS    platforms=3;platformnames=PacBio,Illumina,10X;datasets=5;datasetnames=CCS15kb_20kb,HiSeqPE300x,10XChromiumL> 纯合
   pos = 31676174 # HG003   1  31676174 .     GAAGGG                                        T        50 PASS  GT:PS:DP:ADALL:AD:GQ 1/1:.:151:0,0:0,0:320
   pos = 31676172
   pos = 240708775 #  HG003 2     1 240708775 .     TAATTGTAA                                     C        50 PASS  GT:PS:DP:ADALL:AD:GQ 0/1:.:45:0,0:0,0:174


   result = get_pos_info(bam_af, chrom, pos)
   result_copy = copy.copy(result)
   r = add_attributes_pos_info(result_copy)

   for attr in list(result_copy.__dict__.keys()):
      print(attr, result_copy.__dict__[attr])

   i = 1
   # print(result)
   # print(result.coverage)
   # print(collections.Counter(result.query_snp))
   # print(collections.Counter(result.query_indel))

