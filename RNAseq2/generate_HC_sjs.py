"""Script for filtering splice junctions from STAR mapping round 1. STep 1: Filter splice junctions based on canonicity, read coverage and overhang etc. based on results from each sample. Then pool the filtered SJs from all samples into one file"""

import os
import argparse


class StarSjOutput:
	def __init__(self,line):
		self.line = line
		line1 = line.strip("\n").split("\t")
		self.chromosome = line1[0]
		self.start = line1[1]
		self.end = line1[2]
		self.strand = line1[3]
		self.intronmotif = int(line1[4]) #canonical or non-canonical
		self.is_annotated = int(line1[5])#0 means not annotated, 1 is annotated
		self.read_support_single = int(line1[6])
		self.read_support_multi = int(line1[7])
		self.max_overhang = int(line1[8])
	def filter_sjs(self, filter_non_canonical, readsupportcutoff, min_overhang, outfile):
		if self.is_annotated:
			return
		if filter_non_canonical and not self.intronmotif:
			return
		if self.max_overhang < min_overhang:
			return
		if self.read_support_single < readsupportcutoff and self.read_support_multi < readsupportcutoff:
			return
		outfile.write(self.line)
		return
	def return_id(self): #Returns splice junction unique identifier
		return self.chromosome + "_" + self.start + "_" + self.end + "_" + str(self.strand)


def main():
	parser = argparse.ArgumentParser(description='filter and concatonate STAR SJ output for two pass mapping')
	parser.add_argument('-i', dest = 'inpath', type = str, help = 'The directory containing the input SJ files')
	parser.add_argument('--s', dest = 'readsupportcutoff', type = int, help = 'Minimum number of reads required for a splice junction to be supported. Splice junctions with read support below this number will be removed', default = 5)
	parser.add_argument('--oh', dest = 'min_overhang', type = int, help = 'Minimum overhang (bp) for SJs. If below this number the SJ will be removed', default = 10)
	parser.add_argument('-o', dest = 'outfile',type = str, help = 'output file name')
	args = parser.parse_args()
	input_count = 0
	input_sjs = set()
	with open(args.outfile,"w") as outfile:
		for file_name in os.listdir(args.inpath):
			if file_name.endswith("_SJ.out.tab"):
				print(f"Filtering {file_name}")
				for line in open(args.inpath + "/" + file_name):
					input_count += 1
					starline = StarSjOutput(line)
					input_sjs.add(starline.return_id())
					starline.filter_sjs(True, args.readsupportcutoff, args.min_overhang, outfile)
	output_sjs = set()
	for line in open(args.outfile):
		output_sjs.add(StarSjOutput(line).return_id())
	output_count = len(open(args.outfile).readlines())
	print(f"NUmber of sjs before filtering (total): {input_count}")
	print(f"NUmber of sjs before filtering (unique): {len(input_sjs)}")
	print(f"Number of sjs after filtering (total): {output_count}")
	print(f"Number of sjs after filtering (unique): {len(output_sjs)}")
	print(f"Total number of sjs removed (total): {input_count - output_count}")
	print(f"Total number of sjs removed (unique): {len(input_sjs - output_sjs)}")

if __name__ == "__main__":
	main()

