#!/usr/bin/env python3
# 点突变信息提取
# v0.4
# What's New: 重构所有代码，可以准确和快速地读取金标准vcf文件中的位点信息，不需要再按照SNP和InDel读取，更简练，更结构化

# TODO:
# 1, 输出进度时，显示总数
# 2, 按照MAPQ过滤reads

import os
import sys
import os.path as path
import argparse
import pysam
import re
import multiprocessing as mp

import lib.utils as utils
import lib.info as info
import lib.vcf as vcf


ARGUMENTS_DICT = {}
BASES = ['A', 'T', 'C', 'G']
GC_COMPILE = re.compile(r'C|G')

def get_arguments() -> None:
   '''
   读取命令函参数

   Returns:
       **arguments_dict**: dict
           参数值
   '''
   global ARGUMENTS_DICT

   argvList = sys.argv
   parser_ar = argparse.ArgumentParser(prog = 'PROG',
                                       usage = '{0} [OPTION] <bam file> <locus file>'.format(argvList[0]),
                                       description ='提取指定位置的比对信息。',
                                       epilog='本脚本提取一个bam文件指定位置的比对信息，例如碱基数量，方向，错误率，上下文等。',
                                       formatter_class=argparse.RawTextHelpFormatter)

   parser_ar.add_argument('BAM_FILE', help = 'FILE. 比对文件', metavar = 'bam file')
   parser_ar.add_argument('LOCUS_FILE', help = 'FILE. 位置文件 （VCF, BED, POS）', metavar='locus file')

   parser_ar.add_argument('-l', '--locus-format', default = 'VCF', help = 'STR. 位置文件的格式 （VCF, BED, POS）', dest='LOCUS_FORMAT')
   parser_ar.add_argument('-o', '--output', default='', help= '输出文件', metavar = '', dest='OUTPUT')
   parser_ar.add_argument('-r', '--reference', default='', help= 'FILE. faidx indexed参考基因组文件（.fasta）', metavar = '', dest='REFERENCE_FILE')
   parser_ar.add_argument('-v', '--vcf', default='', help= 'FILE. 标准位点VCF文件', metavar = '', dest='VCF_FILE')
   parser_ar.add_argument('-c', '--context', default=5, type=int, help= 'INT. 提取上下游的各n个碱基写入结果文件，默认值为5', metavar = '', dest='CONTEXT_FLANK')

   parser_ar.add_argument('-f', '--format', default='', help= 'STR. 需要额外输出的位点信息，用,分割，例如matched_snp_cycle,unmatched_snp_cycle', metavar = '', dest='FORAMT_STRING')
   parser_ar.add_argument('-n', '--no-header', action='store_true', default=False, help= '输出文件不需要header', dest='IS_NO_HEADER')
   parser_ar.add_argument('-u', '--locus-as-standard', action='store_true', default=False, help= '如果locus为VCF文件，则直接使用它作为标准位点', dest='LOCUS_AS_STANDARD')
   parser_ar.add_argument('-t', '--threads', default=10, type=int, help= 'INT. 进程数，默认值为10', metavar = '', dest='PROCESS')



   paramters = parser_ar.parse_args()

   ARGUMENTS_DICT['BAM'] = paramters.BAM_FILE
   ARGUMENTS_DICT['LOCUS'] = paramters.LOCUS_FILE

   ARGUMENTS_DICT['LOCUS_FORMAT'] = paramters.LOCUS_FORMAT
   ARGUMENTS_DICT['OUTPUT'] = paramters.OUTPUT
   ARGUMENTS_DICT['REFERENCE'] = paramters.REFERENCE_FILE
   ARGUMENTS_DICT['VCF_FILE'] = paramters.VCF_FILE
   ARGUMENTS_DICT['CONTEXT_FLANK'] = int(paramters.CONTEXT_FLANK)

   ARGUMENTS_DICT['FORAMT_STRING'] = paramters.FORAMT_STRING
   ARGUMENTS_DICT['IS_NO_HEADER'] = paramters.IS_NO_HEADER
   ARGUMENTS_DICT['LOCUS_AS_STANDARD'] = paramters.LOCUS_AS_STANDARD
   ARGUMENTS_DICT['PROCESS'] = paramters.PROCESS

   return None

# 写入结果文件, 输入'#done#'结束
def write_file(q: mp.Queue, output_file: str):
   """
   持续监听queue，写入文件，在queue中输入'#done#'结束
   continue to listen for messages on the queue and writes to file when receive one
   if it receives a '#done#' message it will exit
   """
   with open(output_file, 'a', buffering = 1000) as out_f:
      while True:
         m = q.get()
         if m == '#done#':
            break
         out_f.write(m)
         out_f.flush()

   return None

def multiple_process_helper(bam_file: str, loci_lst: list, format_list: list, q: mp.Queue, reference_file: str = '', real_site_dict: dict = None, flank: int = 5, counter: mp.Value = None, counter_lock: mp.Lock = None) -> int:
   '''
   多线程运行的helper，负责打开bam_file, 返回句柄，收集位点信息，写入StringIO

   Parameters:
      **bam_file**: string
         bam file

      **loci_lst**: list[tuple[str, int, str], ...]
         chrom, pos, other

      **reference_file**: 参考基因组
         indexed fasta file


   Returns:
       **value**: type
           0 on success, other on failure
   '''


   bam_af = pysam.AlignmentFile(path.realpath(path.expanduser(bam_file)))

   if reference_file != '':
      genome_reference_file_handle, index_dict = utils.read_reference(reference_file)
   else:
      genome_reference_file_handle = None
      index_dict = None

   # ==================================================
   i = 0
   for chrom, pos, other_str in loci_lst:
      with counter_lock:
         counter.value += 1
         if counter.value % 1000 == 0:
            print(' '*50, end = '\r')
            print('{}\t{}\t{}'.format(counter.value, chrom, pos), end = '\r')

      if genome_reference_file_handle is not None and index_dict is not None:
         ref_base = utils.get_base_fast(genome_reference_file_handle, index_dict, chrom, pos)
         context = utils.get_base_fast(genome_reference_file_handle, index_dict, chrom, pos - flank, end = pos + flank)
      else:
         ref_base = ''
         context = ''

      if real_site_dict is not None:
         real_allele_snp = []
         real_allele_indel = []
         for real_allele_str in real_site_dict[(chrom, pos)]:
            if '+' in real_allele_str or '-' in real_allele_str:
               real_allele_indel.append(real_allele_str)
            else:
               real_allele_snp.append(real_allele_str)
      else:
         real_allele_snp = None
         real_allele_indel = None

      try:
         # real_allele_snp 和 real_allele_indel 都是 [str, ...]
         pos_PositionInfo = info.get_pos_info(bam_af, chrom, pos, real_allele_snp, real_allele_indel)
         pos_PositionInfo.reference = ref_base
         pos_PositionInfo.context = context
         pos_PositionInfo.real_allele_snp = real_allele_snp
         pos_PositionInfo.real_allele_indel = real_allele_indel
         pos_PositionInfo.other = other_str
      except Exception as ex:
         message = 'multiple_process_helper：位置文件 Line {}: {} {} {}'.format(i, chrom, pos, ex)
         print(message)
         continue

      _ = info.add_attributes_pos_info(pos_PositionInfo)
      line_str = info.output_attributes_pos_info(pos_PositionInfo, format_list)
      q.put(line_str + '\n')

   try:
      bam_af.close()
   except:
      pass

   try:
      genome_reference_file_handle.close()
   except:
      pass

   return 0

def main(argvList = sys.argv, argv_int = len(sys.argv)):

   # = = = = = = = = = = = = = = = = = = positional parameters = = = = = = = = = = = = = = = = = =
   LOCUS_FILE = ARGUMENTS_DICT['LOCUS']
   BAM_FILE = ARGUMENTS_DICT['BAM']

   # = = = = = = = = = = = = = = = = = = optional parameters = = = = = = = = = = = = = = = = = =
   OUTPUT = ARGUMENTS_DICT['OUTPUT']
   if OUTPUT == '':
      locus_basename_str = path.splitext(path.split(LOCUS_FILE)[-1])[0]  # LOCUS_FILE 如果是'~/locus/locus.vcf'，locus_basename_str就是locus
      bam_basename_str = path.splitext(path.split(BAM_FILE)[-1])[0]  # BAM_FILE的文件名
      output_str = bam_basename_str + '_' + locus_basename_str + '.tsv'
   else:
      output_str = path.realpath(path.expanduser(OUTPUT))

   REFERENCE_FILE = ARGUMENTS_DICT['REFERENCE'] # /share/data/reference/human/b37/Homo_sapiens_assembly19.fasta
   GOLDEN_FILE = ARGUMENTS_DICT['VCF_FILE']

   CONTEXT_FLANK = ARGUMENTS_DICT['CONTEXT_FLANK']
   if CONTEXT_FLANK < 0:
      CONTEXT_FLANK = 0

   FORAMT_STRING = ARGUMENTS_DICT['FORAMT_STRING']
   IS_NO_HEADER = ARGUMENTS_DICT['IS_NO_HEADER']
   LOCUS_AS_STANDARD = ARGUMENTS_DICT['LOCUS_AS_STANDARD']
   PROCESS = ARGUMENTS_DICT['PROCESS']

   # = = = = = = = = = = = = = = = = = = analysis = = = = = = = = = = = = = = = = = =
   manager = mp.Manager()
   q = manager.Queue()
   try:
      os.remove(output_str)
   except:
      pass


   print('读取位点...')
   if LOCUS_AS_STANDARD and (LOCUS_FILE.endswith('.vcf.gz') or LOCUS_FILE.endswith('.vcf')):
      GOLDEN_FILE = LOCUS_FILE

   if GOLDEN_FILE != '':
      format_list = ['chrom', 'pos', 'reference', 'context', 'coverage', 'A_count', 'T_count', 'C_count', 'G_count', 'N_count', 'miss_count', 'background_count', 'query_snp_counter', 'real_allele_snp', 'matched_snp_count', 'unmatched_snp_count', 'query_indel_counter', 'real_allele_indel', 'matched_indel_count', 'unmatched_indel_count']
      if GOLDEN_FILE.endswith('.realsite') or os.access(GOLDEN_FILE.split('.gz')[0] + '.realsite', os.R_OK):
         real_site_dict = vcf.get_real_variants_from_realsite(GOLDEN_FILE)
      elif GOLDEN_FILE.endswith('.vcf') or GOLDEN_FILE.endswith('.vcf.gz'):
         real_site_dict = vcf.get_real_variants_from_vcf(GOLDEN_FILE)
      else:
         message = f'标准位点必须是vcf文件或者是realsite文件，输入为{GOLDEN_FILE}'
         sys.exit(message)
   else:
      real_site_dict = None
      format_list = ['chrom', 'pos', 'reference', 'context', 'coverage', 'A_count', 'T_count', 'C_count', 'G_count', 'N_count', 'miss_count', 'background_count', 'query_snp_counter', 'query_indel_counter']

   if FORAMT_STRING != '':
      format_list.extend(FORAMT_STRING.split(','))

   print('读取位置...')
   locus_iter = utils.parse_locus(LOCUS_FILE, ARGUMENTS_DICT['LOCUS_FORMAT'])  #
   locus_lst = list(locus_iter)  # [[chrom, int, other], [chrom, int, other], ...]
   chunk_int = min([len(locus_lst), PROCESS])
   item_lst = utils.slice_list(locus_lst, chunk_int)

   file_pool = mp.Pool(1)
   file_pool.apply_async(write_file, (q, output_str, ))


   header_str = '\t'.join(format_list)
   if not IS_NO_HEADER:
      q.put(header_str + '\n')

   pool = mp.Pool(chunk_int)
   counter = manager.Value('i', 0)
   counter_lock = manager.Lock()
   jobs = []
   for loci_lst in item_lst:
      # real_site_dict
      job = pool.apply_async(multiple_process_helper, (BAM_FILE, loci_lst, format_list, q, REFERENCE_FILE, real_site_dict, CONTEXT_FLANK, counter, counter_lock, ))
      jobs.append(job)

   for job in jobs:
      job.get()

   q.put('#done#')  # all workers are done, we close the output file
   pool.close()
   pool.join()

   print(counter.value, 'loci Done', output_str)
   return


if __name__ == '__main__':

   r = get_arguments()
   main()

