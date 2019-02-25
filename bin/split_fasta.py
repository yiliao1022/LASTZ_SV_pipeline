import getopt
import sys
import os
import subprocess
import random

#comparative_dir = '/share/project002/licai/develop/comparative/trunk'
#comparative_dir = os.path.abspath(os.path.dirname(sys.argv[0]) + '/../../../../../../../')
chainNet_bin_dir = '/home/yliao/prog/kentUtils/bin/linux.x86_64' 
faToTwoBit_prg = '/home/yliao/prog/kentUtils/bin/linux.x86_64/faToTwoBit'

def read_fasta_to_dict(fasta_file):
	try:
		fseq = open(fasta_file, 'r')
	except IOError:
		print >>sys.stderr, 'Error: %s does not exist or is not readable!'
		sys.exit(1)
	
	seq_dict = {}
	curr_name = ''
	curr_seq = []
	for line in fseq:
		if line.find('>') >= 0:
			#print '#line: ', line
			if curr_name != '' and curr_seq != []:
				seq_dict[curr_name] = ''.join(curr_seq)
				curr_seq = []
			curr_name = line[(line.find('>') + 1) : -1]
			#print '#curr_name: ', curr_name
		else:
			if len(line) > 1:
				curr_seq.append(line[:-1])
	if curr_name != '' and curr_seq != []:
		seq_dict[curr_name] = ''.join(curr_seq)
	if len(seq_dict) == 0:
		print >>sys.stderr, 'Error: %s has no sequence or is not in FASTA format!' %(fasta_file)
		sys.exit(1)
	
	return seq_dict

def write_seqs_in_fasta(seq_list, out_file):
	#print 'in write_seqs_in_fasta...'

	try:
		fout = open(out_file, 'w')
	except IOError:
		print >>sys.stderr, 'Error: Can not open %s to write sequences!' %(out_file)
		sys.exit(1)
	for id, seq in seq_list:
		print >>fout, '>%s' %(id)
		for i in range(len(seq) / 50):
			print >>fout, seq[i*50: (i+1)*50]
		if len(seq) % 50 > 0:
			print >>fout, seq[-(len(seq) % 50):]
	fout.close()
	
#def get_next_seq(fin):
	

#TODO
#1 to 1 split
#not split but 2bit
#2bit directory
#only chainNet
#TODO

########## split_1to1	

def split_1to1(args):
	seqs_dict = read_fasta_to_dict(args.input)
	
	for id in seqs_dict:
		fa_file = args.output+ '/' + args.prefix + id.split()[0] + '.fa'
		write_seqs_in_fasta([[id, seqs_dict[id]]], fa_file)
		
		if args.two_bit_flag:
				two_bit_file = args.output+ '/' + args.prefix + id.split()[0] + '.2bit'
				fa_to_2bit(fa_file, two_bit_file, args.remove_flag)
		if args.nib_flag:
				nib_file = args.output + '/' + id.split()[0] + '.nib'
				fa_to_nib(fa_file, nib_file, args.remove_flag)

########## split_1to1

def fa_to_2bit(fa_file, two_bit_file, remove_flag):
	#faToTwoBit_prg = chainNet_bin_dir + '/' + 'faToTwoBit'
	faToTwoBit_cmd = ' '.join([faToTwoBit_prg, fa_file, two_bit_file])
	popen_ret = subprocess.Popen([faToTwoBit_prg, fa_file, two_bit_file])
	popen_ret.wait()
	if remove_flag:
		os.remove(fa_file)
def fa_to_nib(fa_file, nib_file, remove_flag):
	faToNib_path = chainNet_bin_dir + '/' + 'faToNib'
	faToNib_cmd = ' '.join([faToNib_path, fa_file, nib_file])
	popen_ret = subprocess.Popen([faToNib_path, fa_file, nib_file])
	popen_ret.wait()
	if remove_flag:
		os.remove(fa_file)

def split_sequential_ns(args):
	finput = open(args.input, 'r')
	
	seq_cnt = 0
	name = 1
	curr_fh = None
	fa_file = args.output+ '/' + args.prefix  + str(name) + '.fa'
	curr_fh = open(fa_file, 'w')
	for line in finput:
		if line.find('>') >=0:
			seq_cnt += 1
			if seq_cnt > args.ns:
				if curr_fh:
					curr_fh.close()
					if args.two_bit_flag:
						two_bit_file = args.output+ '/' + args.prefix + str(name) + '.2bit'
						fa_to_2bit(fa_file, two_bit_file, args.remove_flag)
					curr_fh = None
				
				seq_cnt = 1	
				name += 1
				fa_file = args.output+ '/' + args.prefix  + str(name) + '.fa'
				curr_fh = open(fa_file, 'w')
		
		print >>curr_fh, line[:-1]
			
	if curr_fh:
		curr_fh.close()
		if args.two_bit_flag:
			two_bit_file = args.output+ '/' + args.prefix + str(name) + '.2bit'
			fa_to_2bit(fa_file, two_bit_file, args.remove_flag)

def split_sequential_nf(args):
	finput = open(args.input, 'r')
	
	file_size = os.path.getsize(args.input)
	mean_size = file_size / float(args.nf)
	
	size_cnt = 0
	name = 1
	curr_fh = None
	fa_file = args.output+ '/' + args.prefix  + str(name) + '.fa'
	curr_fh = open(fa_file, 'w')
	for line in finput:
		if line.find('>') >=0:
			if size_cnt > mean_size:
				if curr_fh:
					curr_fh.close()
					if args.two_bit_flag:
						two_bit_file = args.output+ '/' + args.prefix + str(name) + '.2bit'
						fa_to_2bit(fa_file, two_bit_file, args.remove_flag)
					curr_fh = None
				name += 1
				fa_file = args.output+ '/' + args.prefix + str(name) + '.fa'
				curr_fh = open(fa_file, 'w')
				size_cnt = 0
		size_cnt += len(line)
		print >>curr_fh, line[:-1]
			
	if curr_fh:
		curr_fh.close()
		if args.two_bit_flag:
			two_bit_file = args.output+ '/' + args.prefix + str(name) + '.2bit'
			fa_to_2bit(fa_file, two_bit_file, args.remove_flag)
			
def split_random_ns(args):
	seqs_dict = read_fasta_to_dict(args.input)
	id_list = seqs_dict.keys()
	random.shuffle(id_list)
	
	seq_list = []
	name = 1
	count = 0
	total_len = 0
	
	for i in range(len(id_list)):
		seq_list.append([id_list[i], seqs_dict[id_list[i]]])
		count += 1
		if args.ns == count or i == (len(id_list) - 1):
			fa_file = args.output+ '/' + args.prefix  + str(name) + '.fa'
			write_seqs_in_fasta(seq_list, fa_file)
			if args.two_bit_flag:
				two_bit_file = args.output+ '/' + args.prefix + str(name) + '.2bit'
				fa_to_2bit(fa_file, two_bit_file, args.remove_flag)
			name += 1
			seq_list = []
			count = 0

def split_random_nf(args):
	seqs_dict = read_fasta_to_dict(args.input)
	id_list = seqs_dict.keys()
	random.shuffle(id_list)
	
	seq_list = []
	name = 1
	count = 0
	total_len = 0
	
	for id in id_list:
		total_len = total_len + len(seqs_dict[id])
	average_len = total_len / args.nf
	
	print "average_len: ", average_len
	for i in range(len(id_list)):
		seq = seqs_dict[id_list[i]]
		seq_list.append([id_list[i], seq])
		count += len(seq)
		#print count
		if count >= average_len or i == (len(id_list) - 1):
			#print 'part %d is finished!' %(name)
			fa_file = args.output+ '/' + args.prefix + str(name) + '.fa'
			write_seqs_in_fasta(seq_list, fa_file)
			if args.two_bit_flag:
				two_bit_file = args.output + '/' + args.prefix + str(name) + '.2bit'
				fa_to_2bit(fa_file, two_bit_file, args.remove_flag)
				
			name += 1
			count = 0
			seq_list=[]
			
	print "split completed!"

class Args:
	def __init__(self):
		self.input = ''
		self.mode = '1'
		self.nf = 0
		self.ns = 0
		self.prefix = 'split_'
		self.output = 'output'
		self.two_bit_flag = False
		self.nib_flag = False
		self.remove_flag = False



def print_usage():
	usage = '''
Usage: 
	python split_fasta.py [Options]

Opthion:
	-i, --in <path>    input file in FASTA format
	--mode  [1|2|3]    split mode:1, sequentail; 2, random; 3, 1 sequence in 1 file
	--nf <int>         split into <int> file
	--ns <int>         <int> sequences in each file
	--prefix <str>     args.prefix for output file
	--out <dir>     output directory for split files
	--2bit             convert the output files in FASTA format into .2bit format
	--nib              convert the output files in FASTA format into .nib format
	--remove_fasta     remove the output files in FASTA, work with '--2bit' or '--nib'
	-h, --help         print this help
	'''
	print >>sys.stderr, usage	
def main():
	
	long_opts = ['help', 'in=', 'mode=', 'nf=', 'ns=', 'prefix=', 'out=', '2bit', 'nib', 'remove_fasta']
	opts, args = getopt.getopt(sys.argv[1:], 'hi:n:o:', long_opts)
	
	args = Args()
	for o, a in opts:
		#print o, a
		if o in ('-i', '--in'):
			args.input = a
		elif o in ('--mode'):
			args.mode = a
		elif o in ('--nf'):
			args.nf = int(a)
		elif o in ('--ns'):
			args.ns = int(a)
		elif o in ('--prefix'):
			args.prefix = a
		elif o in ('-o', '--out'):
			args.output = a
			if not os.path.exists(args.output):
				os.mkdir(args.output)
			elif os.path.isfile(args.output):
				print >>sys.stderr, 'Error: file %s exists, can not mkdir %s!' %(args.output, args.output)
			
		elif o in ('--2bit'):
			args.two_bit_flag = True
		elif o in ('--nib'):
			args.nib_flag = True
		elif o in ('--remove_fasta'):
			args.remove_flag = True
		elif o in ('-h', '--help'):
			print_usage()
			sys.exit(1)
		else:
			print >>sys.stderr, 'Error: Unknow option "%s" ' %(o)
			print_usage()
			sys.exit(1)
			
	if args.input == '':
		print >>sys.stderr, 'Error: You must specify the input file!'
		print_usage()
		sys.exit(1)
	
	if args.mode == '3':
		split_1to1(args)
		if args.ns > 0 or args.ns > 0:
			print >>sys.stderr, 'Warning: "--ns" and "--nf" are ignored when using "--mode 3".'
	else:
		if args.nf <= 0 and args.ns <= 0:
			print >>sys.stderr, 'Error: You must set "--nf" or "--ns" with mode %s.' %(args.mode)
			sys.exit(1)
		if args.nf > 0 and args.ns > 0:
			print >>sys.stderr, 'Error: "--nf" and "--ns" are mutual exclusive, you can only set one of them.'
			sys.exit(1)
		if args.mode == '1':
			if args.ns > 0:
				split_sequential_ns(args)
			if args.nf > 0:
				split_sequential_nf(args)
		if args.mode == '2':
			if args.ns > 0:
				split_random_ns(args)
			elif args.nf > 0:
				split_random_nf(args)

		#split_sequential(args)
		
if __name__ == '__main__':
	main()


