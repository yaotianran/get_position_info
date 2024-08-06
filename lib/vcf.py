#!/usr/bin/env python3
# 处理位点vcf文件的函数
# 也可以直接运行: vcf.py XXXX.vcf.gz 生成realsite文件 XXXX.vcf.realsite
import os
import sys
import os.path as path
import collections
import gzip

BASES = ['A', 'T', 'C', 'G']

# 输入vcf文件的REF和ALT字段, 形成一个genotype list，形如[REF, ALT1, ALT2, ...]
# 然后在每一个ALT后加入空格，将每个字段的长度增加至ref的长度（pad = True）
# 例如['A', 'ACGGGGG', 'TCGGGGG']  或者
# ['TACACACACACACAC', 'TACACACACAC    ', 'T              ']
# 否则不加入空格
# ['TACACACACACACAC', 'TACACACACAC', 'T']
def __pad_alt(ref: str, alt: str, pad = True) -> list[str, ...]:

   genotype_lst = [ref]
   genotype_lst.extend(alt.split(','))
   ref_length_int = len(genotype_lst[0])  # ref的长度

   result_lst = []
   for s in genotype_lst:
      if pad:
         result_lst.append(s.ljust(ref_length_int))
      else:
         result_lst.append(s)


   return result_lst

# 输入vcf文件的POS, REF和ALT字段, 输出一个修正POS和samtools style的字符例如'+2AC'或者'-3NNN'
# 例如 34513115，‘T’, 'TAC' 会返回 34513115 ， +2AC
# 199570773，‘TG’, 'AGTAAATTAT' 会返回 199570774 ， +8TAAATTAT
# 2540473 TACACACACACACAC TACACACACAC 会返回 2540483 ， -4NNNN
# 注意ALT不能有‘,’
# 有错误发生时返回[None, None]
def __indel_ref_and_alt(pos: int, ref: str, alt: str) -> tuple[int, str]:

   if len(ref) == len(alt):
      return [None, None]

   if len(ref) < len(alt):  # 插入
      pos_correct_int = pos + len(ref) - 1
      insert_str =  alt[len(ref):]
      samtools_style_str = '+{}{}'.format(len(insert_str), insert_str)

   if len(ref) > len(alt):  #  缺失
      pos_correct_int = pos + len(alt) - 1
      del_str =  ref[len(alt):]
      samtools_style_str = '-{}{}'.format(len(del_str), 'N'*len(del_str))

   return pos_correct_int, samtools_style_str

# 读入realsite文件
def get_real_variants_from_realsite(real_site_file: str) -> dict:

   real_site_dict = collections.defaultdict(list)
   i = 0
   with open(path.realpath(path.expanduser(real_site_file))) as in_f:
      for line_str in in_f:
         if line_str.strip() == '':
            continue

         line_lst = line_str.strip().split()
         real_site_dict[(line_lst[0], int(line_lst[1]))].extend(line_lst[2:])
         i += 1

         if i % 100000 == 0:
            print(i, line_lst, '                ', end = '\r')
   message = f'read {i} sites. Done.                         '
   print(message)
   return real_site_dict

# 将金标准位点VCF文件，转换为真实位点文件，并读入
# 输出为一个字典 {(chrom, pos):[variant_1, variant_2, ....]}
# 真实位点文件为一个tsv文本文件，chrom pos variant_1 variant_2 ....
# variant的格式为samtools风格， 例如'A', 'T', '*', ‘+2AC’， ‘-3NNN’等
def get_real_variants_from_vcf(vcf_file: str, pass_only = True, qual = 0) -> dict:
   '''
   读取vcf文件中特定位置的pos，genotype ref和alt信息
   返回一个字典，格式为 {(chrom, pos):[variant_1, variant_2, ....]}。用来作为对比pysam输出的结果是否match
   例如：
            ('6', 166650651): ['A', '-1N'],
            ('6', 166650652): ['*'],
            ('1', 150614867): ['A', '+9GTAAATTAT'],
            ('1', 199570773): ['C', '+2CA'],
            ('2', 18391469): ['G', '+31GAAAGTATTAGGGAAAGTATTAGGGAAAGTA'],
            ('2', 218537764): ['G', '+1G'],
            ('2', 220901634): ['A', 'C', '+1T'],

   Parameter:
      **vcf_file**: string
         vcf文件，可以zip压缩

      **pass_only**: bool
         只读入PASS的位点

      **qual**: int
         只读入QUAL字段大于等于该值的位点

   Return:
      **real_site_dict**: dict
         {(chrom, pos):[..., ...], }

   '''

   vcf_file_str = path.realpath(path.expanduser(vcf_file))
   if vcf_file_str.endswith('.gz'):
      real_site_str = vcf_file_str[:-3]+ '.realsite'
   else:
      real_site_str = vcf_file_str + '.realsite'

   # real_site_dict = {('chr1', 1314235):['A', 'T', '-3NNN']}
   real_site_dict = collections.defaultdict(list)

   # read the real site file
   if os.access(real_site_str, os.R_OK):
      message = 'read {}'.format(real_site_str)
      print(message)
      real_site_dict = get_real_variants_from_realsite(real_site_str)
      return real_site_dict

   n = 1
   with open(vcf_file_str) if not vcf_file_str.endswith('.gz') else gzip.open(vcf_file_str) as in_f, open(real_site_str, 'wt') as out_f:
      for line in in_f:

         if isinstance(line, bytes):
            line_str = line.decode().strip()
         else:
            line_str = line.strip()

         if line_str.startswith('#'):
            continue

         try:
            line_lst = line_str.split()
            chrom = line_lst[0]
            pos = int(line_lst[1])
            ref = line_lst[3]
            alt = line_lst[4]
            qual_int = int(line_lst[5])
            filter_str = line_lst[6]
            format_str = line_lst[8]
            sample_str = line_lst[9]

            if pass_only and filter_str != 'PASS':
               raise ValueError('FILTER is not PASS.')

            if qual_int < qual:
               raise ValueError('QUAL is less than {}'.format(qual))

            if 'GT' not in format_str:
               raise TypeError('FORMAT string does not have GT tag.')

         except Exception as ex:
            message = '跳过该条目'
            print(ex, message, line_str.strip())
            continue

         if n % 100000 == 0:
            print(n, chrom, pos, ref, alt, '                            ', end = '\r')

         gt_index = format_str.split(':').index('GT')
         gt_str = sample_str.split(':')[gt_index]
         if '/' in gt_str:
            gt_lst = gt_str.split('/')
         elif '|' in gt_str:
            gt_lst = gt_str.split('|')
         else:
            continue

         gt_lst = list(set(gt_lst))
         gt_lst.sort()
         gt_lst = [int(x) for x in gt_lst] #  在双倍体中，形如 [0, 1] （杂合）或者 [1, 2]（双alt杂合），或者 [1]（纯合）


         pos_in_this_loop_lst = []  # 记录每一行vcf文件衍生的真实位点，用于写入realsite文件
         # 添加点突变 (直接将可能的碱基或者'*'(无论是ref还是alt)添加至相应位置)
         genotype_lst = __pad_alt(ref, alt, pad = True)  #  genotype_lst 形如 ['TACACAC', 'TACACACACAC', 'T      '] 第一个元素为ref
         ref_length_int = len(genotype_lst[0])
         for i in range(ref_length_int):
            for j in gt_lst:
               variant_str = genotype_lst[j][i] if genotype_lst[j][i] != ' ' else '*'
               real_site_dict[(chrom, pos + i)].append(variant_str)
               pos_in_this_loop_lst.append(pos + i)

         # 添加InDel （将alt与ref的比较值添加至相应位置）
         genotype_lst = __pad_alt(ref, alt, pad = False)  #  genotype_lst 形如 ['TACACACACACACAC', 'TACACACACAC', 'T'] 第一个元素为ref
         for j in gt_lst:
            if j == 0:
               continue

            pos_shift, indel_str = __indel_ref_and_alt(pos, ref = genotype_lst[0], alt = genotype_lst[j])

            if pos_shift is not None:
               real_site_dict[(chrom, pos_shift)].append(indel_str)
               pos_in_this_loop_lst.append(pos_shift)


         # 写入realsite文件
         for pos_to_write in set(pos_in_this_loop_lst):
            temp_str =  '\t'.join(real_site_dict[(chrom, pos_to_write)])
            out_f.writelines(f'{chrom}\t{pos_to_write}\t{temp_str}\n')


         # 累加计数器
         n += 1

   print('read', n - 1, 'sites.                 ')
   return real_site_dict

if __name__ == '__main__':
   real_site_dict = get_real_variants_from_vcf(sys.argv[1])
   i = 1



