# 其他函数
import os
import sys
import os.path as path
import gzip
import pysam
import math
from collections.abc import Iterator


BASES = ['A', 'T', 'C', 'G']

# 所有的__parse_*函数都接受一个文件，并且返回一个包含染色体，位置和本行其他信息的Iterator
def __parse_vcf(vcf_file: str, pass_only: bool = True, qual = 0) -> Iterator[str, int, str]:
   '''
   解析vcf位置文件，返回一个包含染色体和位置的Iterator
   '''

   vcf_file_str = path.realpath(path.expanduser(vcf_file))
   with open(vcf_file_str) if not vcf_file_str.endswith('.gz') else gzip.open(vcf_file_str) as in_f:
      for line in in_f:

         if isinstance(line, bytes):
            line_str = line.decode().strip()
         else:
            line_str = line.strip()

         if line_str == '' or line.startswith('#'):
            continue

         try:
            line_lst = line_str.split()
            chrom = line_lst[0]
            pos = int(line_lst[1])

            qual_int = int(line_lst[5])
            filter_str = line_lst[6]

            if pass_only and filter_str != 'PASS':
               continue

            if qual_int < qual:
               continue

         except Exception as ex:
            message = f'__parse_vcf: {ex} {line_lst} 格式错误。 跳过'
            print(message)
            continue

         other = '\t'.join(line_lst[2:])

         yield chrom, pos, other


   return None


def __parse_bed(bed_file: str) -> Iterator[str, int, str]:
   '''
   解析bed位置文件，返回一个包含染色体和位置的Iterator
   '''

   bed_file_str = path.realpath(path.expanduser(bed_file))
   with open(bed_file_str) as in_f:
      for line in in_f:
         if line.strip() == '' or line.startswith('#'):
            continue

         try:
            line_lst = line.split()
            chrom = line_lst[0]
            start = int(line_lst[1])
            end = int(line_lst[2])

            other = '\t'.join(line_lst[3:])
         except Exception as ex:
            message = f'__parse_pos: {ex} {line_lst} 格式错误。 跳过'
            print(message)
            continue

         for pos in range(start, end + 1, 1 if end + 1 > start else -1):
            yield chrom, pos, other

   return None

def __parse_pos(pos_file: str) -> Iterator[str, int, str]:
   '''
   解析POS位置文件，返回一个包含染色体和位置的Iterator
   '''
   pos_file_str = path.realpath(path.expanduser(pos_file))
   with open(pos_file_str) as in_f:
      for line in in_f:
         if line.strip() == '' or line.startswith('#'):
            continue

         try:
            line_lst = line.split()
            chrom = line_lst[0]
            pos = int(line_lst[1])
            other = '\t'.join(line_lst[2:])
         except Exception as ex:
            message = f'__parse_pos: {ex} {line_lst} 格式错误。 跳过'
            print(message)
            continue

         yield chrom, pos, other

   return None

# 解析位置文件（POS格式，BED格式，或VCF格式）
# locus_iter = parse_locus(locus_file, format = 'POS')
def parse_locus(locus_file: str, file_format: str) -> Iterator[str, int, str]:
   '''
   根据输入文件的格式，返回一个包含染色体和位置的Iterator

   Parameter:
      **locus_file**: str
         位置文件, 可以是POS格式，BED格式，或VCF格式

      **file_format**: str
         位置文件的格式, 'POS'，'BED'，'VCF'

   Return: locus_iter：Iterator
            返回Iterator[str, int]
            str 为染色体名称，
            int 为位置，
   '''

   if file_format == 'VCF':
      iterator = __parse_vcf(locus_file)
      return iterator

   if file_format == 'POS':
      iterator = __parse_pos(locus_file)
      return iterator

   if file_format == 'BED':
      iterator = __parse_bed(locus_file)
      return iterator

   message = f'parse_locus: 文件格式错误{file_format}，文件格式必须为POS，BED，VCF之一'
   sys.exit(message)

# 根据segment的flag （f1， f2, r1, r2）返回0,1,2,3. 或者返回None（任何错误或者无法判断）
def get_index(segment: pysam.AlignedSegment):

   if segment.is_read1 or segment.is_read2:  #  双端测序
      if segment.is_forward and segment.is_read1:
         return 0

      if segment.is_forward and segment.is_read2:
         return 1

      if segment.is_reverse and segment.is_read1:
         return 2

      if segment.is_reverse and segment.is_read2:
         return 3

   else: #  单端测序
      if segment.is_forward:
         return 0
      elif segment.is_reverse:
         return 2
      else:
         return None

   return None


# 提取基因组特定序列
# genome_reference_file_handle, index_dict = read_reference('/share/home/yaotianran/data/ecoli_genome/ecoli_total_2982.fasta')
# seq = get_base_fast(genome_reference_file_handle, index_dict, '1', 60000, 65000)
def read_reference(reference: str) -> tuple:
   '''
   This function is the helper of function get_base_fast. Input a genome fasta reference, return a handle and and index dictionary

       Parameter:
       **reference**:
           a genome reference fasta file

   Return:
       genome_reference_file_handle: File IO
       index_dict: Dictionary

   '''
   import os.path as path

   refer_str = path.realpath(path.expanduser(reference))
   indexfile = reference + '.fai'

   if not os.access(refer_str, os.R_OK):
      message = 'Reference file {} not found or not readable'.format(refer_str)
      raise FileNotFoundError(message)

   if not os.access(indexfile, os.R_OK):
      message = 'Index file {} not found or not readable'.format(indexfile)
      raise FileNotFoundError(message)

   hanel = open(refer_str, 'rt')

   index_dict = {}
   with open(file = indexfile, mode = 'rt') as index_f:
      for i, line in enumerate(index_f):
         line_lst = line.strip().split()
         try:
            index_dict[line_lst[0]] = list(map(int, line_lst[1:]))
         except Exception:
            message = 'File {0} is corrupted at line {1}. Position query may raise an IndexError.'.format(indexfile, i + 1)
            print(message)
            continue
   print('Read', len(index_dict), 'chromosomes')

   return hanel, index_dict


def get_base_fast(genome_reference_file_handle, index_dict, chrom, start, end=0) -> str:
   '''
   Get the specific bases from an index genome file. (samtools faidx <ref.fasta>)
   This function is fast enough to be used in millions of iterations
   For the speed we need to pass the genome_reference_file_handle and index dictionary

   Parameter:
       **genome_reference_file_handle**: file IO
           a genome reference fasta file open handle

       **index_dict**: dict
           index is dictionary in faidx format, can be read from an index genome file (samtools faidx <ref.fasta>)
           {'1': [249250621, 52, 60, 61],
           '2': [243199373, 253404903, 60, 61],
           '3': [198022430, 500657651, 60, 61], .... }

       **chrom**: str
           chromosome name

       **start**: int
           start position (1-based)

       **end**: optional, int
           end position (1-based)

   Return:
       Bases from chromosome
   '''

   if start < 1:
      start = 1

   start_offset_int = -1
   end_offset_int = -1

   # in-chrome offset + in-line offset + in-file offset
   if start % index_dict[chrom][2] != 0:
      start_offset_int = math.floor(start / index_dict[chrom][2]) * index_dict[chrom][3] + start % index_dict[chrom][2] + index_dict[chrom][1]
   else:
      start_offset_int = math.floor((start / index_dict[chrom][2] - 1)) * index_dict[chrom][3] + index_dict[chrom][2] + index_dict[chrom][1]

   genome_reference_file_handle.seek(start_offset_int - 1)

   if end == 0:
      end_offset_int = start_offset_int
   else:
      if end % index_dict[chrom][2] != 0:
         end_offset_int = math.floor(end / index_dict[chrom][2]) * index_dict[chrom][3] + end % index_dict[chrom][2] + index_dict[chrom][1]
      else:
         end_offset_int = math.floor((end / index_dict[chrom][2] - 1)) * index_dict[chrom][3] + index_dict[chrom][2] + index_dict[chrom][1]

   try:
      seq_str = genome_reference_file_handle.read(end_offset_int - start_offset_int + 1).replace('\n', '')
   except Exception as ex:
      print(ex)
      return('')

   return(seq_str.upper())


# 将数组分成N等份
def slice_list(in_list: list, chunk_num: int):
   return [in_list[x::chunk_num] for x in range(chunk_num)]

