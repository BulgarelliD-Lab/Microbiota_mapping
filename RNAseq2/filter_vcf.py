"""Filter vcf file based on quality and samples. All samples must carry same allele for variant to be considered"""

import sys

input_file = sys.argv[1]
out_file = sys.argv[2]

total_calls = 0
passed_multiallelic = 0
passed_monoallelic = 0
with open(out_file,"w") as out:
	for line in open(input_file):
		if line.startswith("#"):
			out.write(line)
			continue
		fields = line.rstrip().split("\t")
		total_calls += 1
		if fields[6] == "PASS":
			genotypes = fields[9:]
			alleles = []
			for genotype in genotypes:
				allele = genotype.split(":")[0]
				alleles.append(allele)
			if alleles[0] == "./." and alleles[1] == "./." and alleles[2] == "1/1":#Take only calls where Barke and 124_17 are reference and 124_52 is homozygous alternative
				out.write(line)
				passed_monoallelic += 1


print(f"Total calls in input: {total_calls}")
print(f"Calls that passed filtering and are monoallelic - these are kept: {passed_monoallelic}")

print("Script finished")




