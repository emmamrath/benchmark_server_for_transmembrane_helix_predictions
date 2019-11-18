#!/usr/bin/env python

# python extract_superfamily_from_pdbe_taxonomy.py taxonomy_temp4 taxonomy_temp4.superfamily

# python extract_superfamily_from_pdbe_taxonomy.py [input-file] [output-file]

import sys
import re

#1gwe_A:::::Genus: Micrococcus;;;;;Family: Micrococcaceae;;;;;Suborder: Micrococcineae;;;;;Order: Actinomycetales;;;;;Subclass: Actinobacteridae;;;;;Class: Actinobacteria (class);;;;;Phylum: Actinobacteria;;;;;Superkingdom: Bacteria;;;;;cellular organisms;;;;;:::::
#1i1w_A:::::Genus: Thermoascus;;;;;Family: Trichocomaceae;;;;;Order: Eurotiales;;;;;Subclass: Eurotiomycetidae;;;;;Class: Eurotiomycetes;;;;;leotiomyceta;;;;;Subphylum: Pezizomycotina;;;;;saccharomyceta;;;;;Phylum: Ascomycota;;;;;Subkingdom: Dikarya;;;;;Kingdom: Fungi;;;;;Opisthokonta;;;;;Superkingdom: Eukaryota;;;;;cellular organisms;;;;;:::::
#2o9s_A:::::Genus: Homo;;;;;Subfamily: Homininae;;;;;Family: Hominidae;;;;;Superfamily: Hominoidea;;;;;Parvorder: Catarrhini;;;;;Infraorder: Simiiformes;;;;;Suborder: Haplorrhini;;;;;Order: Primates;;;;;Superorder: Euarchontoglires;;;;;Eutheria;;;;;Theria;;;;;Class: Mammalia;;;;;Amniota;;;;;Tetrapoda;;;;;Sarcopterygii;;;;;Euteleostomi;;;;;Teleostomi;;;;;Superclass: Gnathostomata;;;;;Vertebrata;;;;;Subphylum: Craniata;;;;;Phylum: Chordata;;;;;Deuterostomia;;;;;Coelomata;;;;;Bilateria;;;;;Eumetazoa;;;;;Kingdom: Metazoa;;;;;Opisthokonta;;;;;Superkingdom: Eukaryota;;;;;cellular organisms;;;;;:::::



# global variables
save_superfamily = []


######################################################
def main():

	global outfile
	input_file = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	count_Eukaryota = 0
	count_Bacteria = 0
	count_Archaea = 0
	count_Viruses = 0
	count_other = 0
	count_blank = 0
	count_unclassified = 0
	count_Vertebrata = 0

	count_array_names = [ 'Superkingdom: Eukaryota','Superkingdom: Bacteria','Superkingdom: Archaea','Superkingdom: Viruses','other sequences','-','unclassified sequences','Vertebrata' ]
	count_array_counts = [ 0, 0, 0, 0, 0, 0, 0, 0 ]

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			bits1 = inline.split( ':::::' )
			this_id = bits1[0]
			this_pdbe = bits1[1]
			bits2 = this_pdbe.split( ';;;;;' )
			num_class = len(bits2)
			last_class = ''
			last_2nd_class = ''
			ix = num_class - 2
			if (num_class >= 1):
				last_class = bits2[ ix ]
			ix = num_class - 3
			if (num_class >= 2):
				last_2nd_class = bits2[ ix ]
			this_superfamily = last_class
			if (last_class == 'cellular organisms'):
				this_superfamily = last_2nd_class
			last_class = 'xxx' + last_class + 'yyy'
			if (this_superfamily == ''):
				this_superfamily = '-'
			output_line =  this_id + ':::::' + this_superfamily + ':::::' + "\r\n"
			outfile.write( output_line )

			found_it = 0
			for this_save_superfamily in save_superfamily:
				if (this_save_superfamily == this_superfamily):
					found_it = 1
			if (found_it == 0):
				save_superfamily.append( this_superfamily )
			for i in range( 0, len(count_array_names) ):
				if (count_array_names[i] == this_superfamily):
					count_array_counts[i] = count_array_counts[i] + 1
	outfile.close()

	for this_save_superfamily in save_superfamily:
		this_count = ''
		for i in range( 0, len(count_array_names) ):
			if (count_array_names[i] == this_save_superfamily):
				this_count = count_array_counts[i]
		print this_save_superfamily, this_count

	# Superkingdom: Eukaryota 4002
	# Superkingdom: Bacteria 4179
	# Superkingdom: Archaea 432
	# Superkingdom: Viruses 617
	# other sequences 21
	# - 11
	# unclassified sequences 21


if __name__ == "__main__":
	# Someone is launching this directly
	main()

