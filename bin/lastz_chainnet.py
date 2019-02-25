#!/usr/bin/env python
import sys
import os
import getopt
from os import path
import subprocess
from shutil import rmtree

#some basic directories
#comparative_dir = os.path.abspath(os.path.dirname(sys.argv[0]) + '/../../../../../')
chainNet_bin_dir ='/home/yliao/prog/kentUtils/bin/linux.x86_64' 
lastz_prg = '/home/yliao/lastz-distrib/bin/lastz'
split_fasta_prg = '/home/yliao/yliao/2018_SV/bin/split_fasta.py'

#chainNet tools
axtChain_prg = chainNet_bin_dir + '/axtChain'
chainMergeSort_prg = chainNet_bin_dir + '/chainMergeSort'
chainPreNet_prg = chainNet_bin_dir + '/chainPreNet'
chainNet_prg = chainNet_bin_dir + '/chainNet'
netSyntenic_prg = chainNet_bin_dir + '/netSyntenic'
netToAxt_prg = chainNet_bin_dir + '/netToAxt'
axtSort_prg = chainNet_bin_dir + '/axtSort'
axtToMaf_prg = chainNet_bin_dir + '/axtToMaf'
faToTwoBit_prg = chainNet_bin_dir + '/faToTwoBit'
faSize_prg = chainNet_bin_dir + '/faSize'

#subdirectories for chainNet ouput
chainnet_subdir = ['0_scirpts', '1_axtChain_out', '2_chainPreNet_out', '3_chainNet_out', '4_netSyntenic_out',\
			'5_netToAxt_out', '6_axtSort_out', '7_axtToMaf_out']

#check and make a new directory
def check_mkdir(new_dir):
	if not os.path.exists(new_dir):
		try:
			os.makedirs(new_dir)
		except OSError, err:
			print >>sys.stderr, 'Error: %s' %(str(err))
			sys.exit(1)
	elif os.path.isdir(new_dir):
		print >>sys.stderr, 'Warning: %s exists' %(new_dir)
	elif os.path.isfile(new_dir):
		print >>sys.stderr, 'Error: a file named %s exists, can not mkdir %s' %(new_dir, new_dir)
		sys.exit(1)

#convert the target and query files in FASTA format to .2bit files
def convert_2bit(args):
	check_mkdir(args.two_bit_out)
	
	args.tseqs_2bit = args.two_bit_out + '/' + os.path.basename(args.tseqs) + '.2bit'
	popen_ret = subprocess.Popen([faToTwoBit_prg, args.tseqs, args.tseqs_2bit])
	popen_ret.wait()
	
	args.qseqs_2bit = args.two_bit_out + '/' + os.path.basename(args.qseqs) + '.2bit'
	popen_ret = subprocess.Popen([faToTwoBit_prg, args.qseqs, args.qseqs_2bit])
	popen_ret.wait()
	
	print >>sys.stderr, 'Convert the input files to .2bit format...'
	
#split the input when the files are large, note the output of splitting is in .2bit format
def split_process(args):
	print >>sys.stderr, 'Split the input files...'
		
	if args.tsplitm == '3':
		if os.path.exists(args.tsplitd):
			rmtree(args.tsplitd)
			print >>sys.stderr, 'Warning: %s exists, removed for new output of splitting!' %(args.tsplitd)
			
		split_cmd = ['python', split_fasta_prg, '--in', os.path.abspath(args.tseqs), '--mode', args.tsplitm,'--out', args.tsplitd, '--prefix', '','--2bit', '--remove_fasta']
		popen_ret = subprocess.Popen(split_cmd)
		popen_ret.wait()
	elif args.tsplitn > 0:
		if os.path.exists(args.tsplitd):
			rmtree(args.tsplitd)
			print >>sys.stderr, 'Warning: %s exists, removed for new output of splitting!' %(args.tsplitd)
		split_cmd = ['python', split_fasta_prg, '--in', os.path.abspath(args.tseqs), '--mode', args.tsplitm, '--nf', str(args.tsplitn), '--out', args.tsplitd, '--2bit', '--remove_fasta']
		#print '#split_cmd: ', ' '.join(split_cmd)
		popen_ret = subprocess.Popen(split_cmd)
		popen_ret.wait()
	
	if args.qsplitm == '3':
		if os.path.exists(args.qsplitd):
			rmtree(args.qsplitd)
			print >>sys.stderr, 'Warning: %s exists, removed for new output of splitting!' %(args.qsplitd)
		split_cmd = ['python', split_fasta_prg, '--in', os.path.abspath(args.qseqs), '--mode', args.qsplitm,'--out', args.qsplitd, '--prefix', '', '--2bit', '--remove_fasta']
		popen_ret = subprocess.Popen(split_cmd)
		popen_ret.wait()
	elif args.qsplitn > 0:
		if os.path.exists(args.qsplitd):
			rmtree(args.qsplitd)
			print >>sys.stderr, 'Warning: %s exists, removed for new output of splitting!' %(args.qsplitd)
		split_cmd = ['python', split_fasta_prg, '--in', os.path.abspath(args.qseqs), '--mode', args.qsplitm, '--nf', str(args.qsplitn), '--out', args.qsplitd, '--2bit', '--remove_fasta']
		#print '#split_cmd: ', ' '.join(split_cmd)
		popen_ret = subprocess.Popen(split_cmd)
		popen_ret.wait()

#gather the input information for following steps
def assign_input(args):
	if args.tsplitn > 0 or args.tsplitm == '3':
		args.tseq_dir = os.path.abspath(args.tsplitd)
		args.tseq_list = os.listdir(args.tsplitd)
	else:
		args.tseq_list.append(os.path.basename(args.tseqs) + '.2bit')
		args.tseq_dir = os.path.abspath(args.two_bit_out)
	if args.qsplitn > 0 or args.qsplitm == '3':
		args.qseq_dir = os.path.abspath(args.qsplitd)
		args.qseq_list = os.listdir(args.qsplitd)
	else:
		args.qseq_list.append(os.path.basename(args.qseqs) + '.2bit')
		args.qseq_dir = os.path.abspath(args.two_bit_out)
	
	args.qseqs_2bit = args.two_bit_out + '/' + os.path.basename(args.qseqs) + '.2bit'
	args.tseqs_2bit = args.two_bit_out + '/' + os.path.basename(args.tseqs) + '.2bit'
	
#lastz process
def lastz_process(args):
	
	if not os.path.exists(args.lastz_out):
		try:
			os.makedirs(args.lastz_out)
		except OSError:
			print >>sys.stderr, 'Error: Can not mkdir %s!' %(args.lastz_out)
			sys.exit(1)
	
	if len(args.tseq_list) > 1:
		for tfile in args.tseq_list:
			new_dir = args.lastz_out + '/' + tfile 
			if not os.path.exists(new_dir):
				os.makedirs(new_dir)
			elif os.path.isdir(new_dir):
				print >>sys.stderr, 'Warning: %s exists' %(new_dir)
			elif os.path.isfile(new_dir):
				print >>sys.stderr, 'Warning: a file named %s exists, can not mkdir %s' %(new_dir, new_dir)
				sys.exit(1)
		
	#get the file extension of lastz output
	lastz_out_ext = '.axt'
	format_pos = args.lastz_para.find('format=')
	if format_pos >= 0:
		lastz_out_ext = '.' + args.lastz_para[format_pos + len('format=') : ].split()[0]
	else:
		print >>sys.stderr, 'Warning: output format not given, use AXT by default'
		args.lastz_para += ' --format=axt'
		lastz_out_ext = '.axt'
	
	lastz_sh_file = 'step1_' + 'lastz_' + args.tname + '_' + args.qname + '.sh'
	if os.path.exists(lastz_sh_file):
		print >>sys.stderr, 'Warning: %s exists and will be renamed as %s.old' %(lastz_sh_file, lastz_sh_file)
		try:
			os.rename(lastz_sh_file, lastz_sh_file + '.old')
		except OSError:
			print >>sys.stderr, 'Error: "%s" and "%s" both exist, you have to remove or rename them.' %(lastz_sh_file, lastz_sh_file + '.old')
			sys.exit(1)

	flastz_sh = open(lastz_sh_file, 'w')
	for tfile in args.tseq_list:
		for qfile in args.qseq_list:
			if len(args.tseq_list) > 1:
				lastz_out_file = os.path.join(args.lastz_out, tfile, tfile + '_' + qfile + lastz_out_ext)
			else:
				lastz_out_file = os.path.join(args.lastz_out, tfile + '_' + qfile + lastz_out_ext)
			
			print >>flastz_sh, ' '.join([lastz_prg, args.tseq_dir + '/' + tfile + '[multi]', \
					args.qseq_dir + '/' + qfile + '[multi]',  args.lastz_para, '>', lastz_out_file])

#chain process
def chain_process(args):
	if args.tsizes == '':
		args.tsizes = os.getcwd() + '/' + os.path.basename(args.tseqs) + '.sizes'
		faSize_cmd = ' '.join([faSize_prg, '-detailed', args.tseqs, '>', args.tsizes])
		subprocess.call(faSize_cmd, shell=True)
	if args.qsizes == '':
		
		args.qsizes = os.getcwd() +  '/' + os.path.basename(args.qseqs) + '.sizes'
		faSize_cmd = ' '.join([faSize_prg, '-detailed', args.qseqs, '>', args.qsizes])
		subprocess.call(faSize_cmd, shell=True)
		#print >>sys.stderr, 'Error: you must provide the file containing information of sequence size(--tsizes & --qsizes)!'
		#sys.exit(1)
	
	#print '#args.tsizes, args.qsizes', args.tsizes, args.qsizes
	if not os.path.exists(args.chainnet_out):
		try:
			os.makedirs(args.chainnet_out)
		except OSError:
			print >>sys.stderr, 'Error: Can not mkdir %s!' %(args.chainnet_out)
			sys.exit(1)
	
	check_mkdir( '/'.join([args.chainnet_out, chainnet_subdir[0]] ))
	check_mkdir( '/'.join([args.chainnet_out, chainnet_subdir[1]] ))
	
	chainnet_sh_file = 'step2_' + 'chain_' + args.tname + '_' + args.qname + '.sh'
	if os.path.exists(chainnet_sh_file):
		print >>sys.stderr, 'Warning: %s exists and will be renamed as %s.old' %(chainnet_sh_file, chainnet_sh_file)
		try:
			os.rename(chainnet_sh_file, chainnet_sh_file + '.old')
		except OSError:
			print >>sys.stderr, 'Error: "%s" and "%s" both exist, you have to remove or rename them.' %(chainnet_sh_file, chainnet_sh_file + '.old')
			sys.exit(1)

	fchainnet = open(chainnet_sh_file, 'w')
	
	cat_taxt_list = []
	for tfile in args.tseq_list:
		if len(args.tseq_list) > 1:
			curr_lastz_out = args.lastz_out + '/' + tfile
		else:
			curr_lastz_out = args.lastz_out
		cat_tfile = args.lastz_out + '/' + tfile + '.cat.axt'
		
		script_file = args.chainnet_out + '/' + chainnet_subdir[0] + '/' + tfile + '.sh'
		fscript = open(script_file, 'w')
		
		if len(args.qseq_list) > 1:
			print >>fscript, 'cat ' + curr_lastz_out + '/* > ' + cat_tfile
			axt_file_left = tfile + '.cat.axt'
		else:
			axt_file_left = os.path.join(tfile + '_' + args.qseq_list[0] + '.axt')
			if len(args.tseq_list) > 1:
				axt_file_left = tfile + '/' + axt_file_left
			
		chain_file = args.chainnet_out + '/' + chainnet_subdir[1] + '/' + axt_file_left + '.chain'
				
		axt_file_full = args.lastz_out + '/' + axt_file_left
		tfile_file_full = args.tseq_dir + '/' + tfile
		
		print >>fscript, '%s -linearGap=%s %s %s %s %s' %(axtChain_prg, args.linearGap, axt_file_full, tfile_file_full, args.qseqs_2bit, chain_file)
				
		print >>fchainnet, 'sh %s' %(script_file) 
	fchainnet.close()

#net process
def net_process(args):
	for i in range(2, len(chainnet_subdir)):
		check_mkdir('/'.join([args.chainnet_out, chainnet_subdir[i]]))
	
	net_sh_file = 'step3_' + 'net_' + args.tname + '_' + args.qname + '.sh'
	fnet = open(net_sh_file, 'w')
	
	all_chain_file = args.chainnet_out + '/' + chainnet_subdir[2] + '/' + 'all.chain'
	chain_filter_file = all_chain_file + '.filter'
	
	qnet_file = args.chainnet_out + '/' + chainnet_subdir[3] + '/' + 'all.chain.filter.qnet'
	tnet_file = args.chainnet_out + '/' + chainnet_subdir[3] + '/' + 'all.chain.filter.tnet'
	synnet_file = args.chainnet_out + '/' + chainnet_subdir[4] + '/' + 'all.chain.filter.tnet.synnet'
	net2axt_file = args.chainnet_out + '/' + chainnet_subdir[5] + '/' + 'all.chain.filter.tnet.synnet.axt'
	axtsort_file = args.chainnet_out + '/' + chainnet_subdir[6] + '/' + 'all.chain.filter.tnet.synnet.axt'
	axt2maf_file = args.chainnet_out + '/' + chainnet_subdir[7] + '/' + 'all.chain.filter.tnet.synnet.axt.maf'

	print >>fnet, chainMergeSort_prg + ' ' + args.chainnet_out + '/' + chainnet_subdir[1] + '/*chain > '\
			+ all_chain_file
	print >>fnet, '%s %s %s %s %s' %(chainPreNet_prg, all_chain_file, args.tsizes, args.qsizes, chain_filter_file)
	print >>fnet, '%s %s %s %s %s %s' %(chainNet_prg, chain_filter_file, args.tsizes, args.qsizes, tnet_file, qnet_file)
	print >>fnet, '%s %s %s' %(netSyntenic_prg, tnet_file, synnet_file)
		#nib or 2bit file
	print >>fnet, '%s %s %s %s %s %s' %(netToAxt_prg, synnet_file, chain_filter_file, args.tseqs_2bit, args.qseqs_2bit, net2axt_file)
	print >>fnet, '%s %s %s' %(axtSort_prg, net2axt_file, axtsort_file)
	print >>fnet, '%s %s %s %s -tPrefix=%s. -qPrefix=%s. %s' %(axtToMaf_prg, net2axt_file, args.tsizes, args.qsizes, args.tname, args.qname, axt2maf_file)

	fnet.close()	
	
#the class is designe
class Args:
	def __init__(self):
		self.qseqs = ''
		self.tseqs = ''
		self.tname = 'target'
		self.qname = 'query'
		self.tseq_list = []
		self.tseq_dir = ''
		self.qseq_list = []
		self.qseq_dir = ''
		
		self.tsplitm = '1'
		self.qsplitm = '1'
		self.tsplitn = 0
		self.qsplitn = 0
		self.tsplitd = os.getcwd() + '/tsplit_out'
		self.qsplitd = os.getcwd() + '/qsplit_out'
		self.two_bit_out = os.getcwd() + '/2bit_out'
		self.qseqs_2bit = ''
		self.tseqs_2bit = ''
		
		self.lastz_para = '--format=axt'
		self.lastz_out = os.getcwd() + '/lastz_out'
		
		self.chainnet_out = os.getcwd() + '/chainnet_out'
		self.tsizes = ''
		self.qsizes = ''
		self.linearGap = chainNet_bin_dir + '/medium'
		
		self.step = '0123'
		#self.run_flag = '0'	
#usage information
def print_usage():
	#print >>stderr, "TODO"
	usage = '''
Description:
    This program is designed for lastz-chain-net whole genome alignment pipeline.   

Usage:
    python lastz_chainnet.py [options]

Options:
    --tseqs, --qseqs      <file> target file and query file <required>
    --tname, --qname      <str>  target name and query name [target, query]
    --tsplitm, --qsplitm  <int>  split mode:1, sequentail; 2, random; 3, 1 sequence in 1 file [1]
    --tsplitn, --qsplitn  <int>  number of subfiles splitting into
    --tsplitd, --qsplitd  <dir>  output directory for splitting results [`pwd`/tsplit_out, `pwd`/qsplit_out]
    --2bit_out            <dir>  output directory for .2bit conversion results [`pwd`/2bit_out]
    --lastz_para          <str>  lastz parameters, with double quotation marks(") outside ["--format=axt"]
    --lastz_out           <dir>  lastz output directory [`pwd`/lastz_out]
    --chainnet_out        <dir>  chainnet output directory [[`pwd`/chainnet_out]]
    --tsizes, --qsizes    <file> files containing the size information of target and query sequences
    --linearGap           <str>  gap scoring system for chaining("axtChain"):loose | medium | filename [medium]
    --step                <str>  which step or steps to run: 0, split; 1, lastz; 2, chain; 3, net [0123]
    -h, --help                   print this help
    
Note:
    1. This program is just used to GENERATE jobs, NOT to RUN the alignment jobs. You have to run the 
    jobs MANUALLY and ORDERLY.
    
    2. Usually, repeats of input sequences should be masked to lower case(soft masking), NOT masked as Ns. 
    If masked as Ns, the result of alignment will be very short. Soft masking can accelerate the aligning
    and repeats can be used after seeding to make the blocks longer.
	
    3. The input sequences are always converted to .2bit format for LASTZ and chainNet, because the .2bit 
    files cost less disk space and are easier to process.
    
    4. For LASTZ/chainNet option setting, you should refer to their own documentation.

Example:
    1. If the target and the query are both small files, so they are not split.
    
      python lastz_chainnet.py --tseqs target.fa --qseqs query.fa --lastz_para "C=2 format=axt" --step 012 > gen.out 2> gen.log
    
    When the command is run over without any error, two shell scripts should be generated:
      step1_lastz_target_query.sh
      step2_chain_target_query.sh
    You can ORDERLY run the jobs in the cluster or in the current node. 
  
    2. If the target and the query are large files and each contains <100 chromosome sequences.
    
      python lastz_chainnet.py --tseqs target.fa --tsplitm 3 --qseqs query.fa --qsplitm 3 --step 012 > gen.out 2> gen.log 

    3. If you want to get the NET results of Example 2,  you can type:
    
      python lastz_chainnet.py --tseqs target.fa --tsplitm 3  --qseqs query.fa --qsplitm 3 --step 3 > gen.out 2> gen.log
   
    4. If the target and the query are large files and both contain many sequeces(much more than 100).
    
      python lastz_chainnet.py --tseqs target.fa --tsplitn 50 --tsplitm 2 --qseqs query.fa --qsplitn 50 --tsplitm 2 \\
      --tsizes target.fa.sizes --qsizes query.fa.sizes --tname species1 --qname species2 > gen.out 2> gen.log
    
    When the command is run over without any error, three shell scripts should be generated:
      step1_lastz_species1_species2.sh
      step2_chain_species1_species2.sh
      step3_net_species1_species2.sh
    You can ORDERLY run the jobs in the cluster or in the current node. Note that 'step3_net_species1_species2.sh' 
    is just one single job and could not be split.
    
	'''
	
	print >>sys.stderr, usage

def main():
	
	long_opts = ['help', 'tseqs=', 'qseqs=', 'tname=', 'qname=', 'lastz_para=', 'lastz_out=',\
			'chainnet_out=', 'tsizes=', 'qsizes=',  'linearGap=', 'step=', \
			'tsplitn=', 'qsplitn=', 'tsplitd=', 'qsplitd=', 'tsplitm=', 'qsplitm=', '2bit_out=']
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'h', long_opts)
	except getopt.GetoptError, err:
		print >>sys.stderr, 'Error: ', str(err)
		print_usage()
		sys.exit(1)

	args = Args()
	for o, a in opts:
		#print o, a
		if o in ('--qseqs'):
			args.qseqs = a
		elif o in ('--tseqs'):
			args.tseqs = a
		elif o in ('--tname'):
			args.tname = a
		elif o in ('--qname'):
			args.qname = a
		elif o in ('--lastz_para'):
			args.lastz_para = a
		elif o in ('--lastz_out'):
			args.lastz_out = a
		elif o in ('--tsplitm'):
			args.tsplitm = a
		elif o in ('--qsplitm'):
			args.qsplitm = a
		elif o in ('--tsplitn'):
			args.tsplitn = int(a)
		elif o in ('--qsplitn'):
			args.qsplitn = int(a)
		elif o in ('--tsplitd'):
			args.tsplitd = a
		elif o in ('--qsplitd'):
			args.qsplitd = a
		elif o in '--2bit_out':
			args.two_bit_out = a
		elif o in ('--chainnet_out'):
			args.chainnet_out = a
		elif o in ('--tsizes'):
			args.tsizes = a
		elif o in ('--qsizes'):
			args.qsizes = a
		elif o in ('--linearGap'):
			args.linearGap = chainNet_bin_dir + '/' + a
		elif o in ('--step'):
			args.step = list(a)
		elif o in ('--help', '-h'):
			print_usage()
			sys.exit(1)
		else:
			print >>sys.stderr, 'Error: Unknow option "\s" ' %(o)
			print_usage()
			sys.exit(1)
	if args.qseqs == '' or args.tseqs == '' and '1' in args.step:
		print >>sys.stderr, 'Error: You must specify the target file and the query file'
		print_usage()
		sys.exit(1)
	
	#The input sequences are always converted to .2bit format for LASTZ and chainNet,
	#because the .2bit files cost less disk space and are easier to process.
	convert_2bit(args)

	if '0' in args.step:
		print >>sys.stderr, 'step 0: Preparing the required input...'
		#convert_2bit(args)
		split_process(args)
	
	assign_input(args)
	
	if '1' in args.step:	
		print >>sys.stderr, 'step 1: Generating LASTZ jobs...'
		lastz_process(args)
	if '2' in args.step:
		print >>sys.stderr, 'step 2: Generating CHAIN jobs...'
		chain_process(args)
	if '3' in args.step:
		print >>sys.stderr, 'step 3: Generating NET jobs...'
		net_process(args)
if __name__ == '__main__':
	main()


