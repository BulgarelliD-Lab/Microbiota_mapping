"""Mummerplot not working properly. this should fix it"""

import argparse


parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', dest = 'input_file', type = str, help = 'input file')
parser.add_argument('-o', dest = 'output', type = str, help = 'output file')
args = parser.parse_args()

print("Fixing gnu plot script...so probably ignore error message above")
with open(args.output,"w") as out:
	for n,line in enumerate(open(args.input_file)):
		if n == 10 or n == 11:
			out.write("#" + line)
		else:
			out.write(line)