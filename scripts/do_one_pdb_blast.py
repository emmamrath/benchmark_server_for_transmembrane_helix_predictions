#!/usr/bin/env python

# python do_one_pdb_blast.py 1zoy_A KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL 2lyz.blastpdb

# python do_one_pdb_blast.py [input-pdbid_chain] [input-sequence] [output-append-file-for-sequences]

# This program makes a call to PDB to get the sequence for a PDBID+chain.

import sys
import os
import re
import SOAPpy


#Examples of output formats from blastall flags
#    pairwise (flag:-m 0)
#    master-slave showing identities (flag:-m 1)
#    master-slave no identities (flag:-m 2)
#    flat master-slave, show identities (flag:-m 3)
#    flat master-slave, no identities (flag:-m 4)
#    master-slave no identities and blunt ends (flag:-m 5)
#    flat master-slave, no identities and blunt ends (flag:-m 6)
#    XML Blast output (flag:-m 7)
#    tabular (flag:-m 8)
#    tabular with comment lines (flag:-m 9)


#<PRE>
#><a name = 5612></a>1E39:1:A|pdbid|entity|chain(s)|sequence
#          Length = 571
#
# Score =  114 bits (284), Expect = 3e-25,   Method: Composition-based stats.
# Identities = 108/411 (26%), Positives = 186/411 (45%), Gaps = 68/411 (16%)
#
#Query: 61  AQGGINAALGNMEE-----DNWRWHFYDTVKGSDWLGDQDAIHYMTEQAPASVVELENYG 115
#           A GG+NAA  + ++     D+    F DT+KG   + D   +  ++  +  SV  +   G
#Sbjct: 168 AAGGMNAAWTDQQKAKKITDSPELMFEDTMKGGQNINDPALVKVLSSHSKDSVDWMTAMG 227
#
#Query: 116 MPFSRTEDGKIYQRAFGGQSLKFGKGGQAHRCCCVADRTGHSLLHTLYGRSLRYDTSYFV 175
#                T+ G +     GG S+      +AHR    A    H ++  LY  +++ +    +
#Sbjct: 228 ADL--TDVGMM-----GGASVN-----RAHRPTGGAGVGAH-VVQVLYDNAVKRNIDLRM 274
#
#Query: 176 EYFALDLLMEN-GECRGVIALCIEDGSIHRIRARNTVVATGGYGRTY------------F 222
#               +++L ++ G  +G++   +  G  + ++A   ++ATGG+ +              F
#Sbjct: 275 NTRGIEVLKDDKGTVKGILVKGMYKG-YYWVKADAVILATGGFAKNNERVAKLDPSLKGF 333
#
#Query: 223 SCTSAHTSTGDGTAMVTRAGLPCQDLEFVQFHPTGIYGAGCLITEGCRGEGGILINSQGE 282
#             T+   + GDG  +   AG   +D++++Q  PT     G ++TE  RG G IL+N +G+
#Sbjct: 334 ISTNQPGAVGDGLDVAENAGGALKDMQYIQAAPTLSVKGGVMVTEAVRGNGAILVNREGK 393
#
#Query: 283 RFMERYAPVAKDLASRDVVS---RSMTLEIREGRGCGPEKDHVYLQLHHLPPEQLAVRL- 338
#           RF+       +D AS  +++   +S  L   +       K   Y+ L   P     V+L 
#Sbjct: 394 RFVNEIT--TRDKASAAILAQTGKSAYLIFDDSVRKSLSKIDKYIGLGVAPTADSLVKLG 451
#
#Query: 339 -------PGISET-----AMIFAGVDVTKE--------------PIPVLPTVHYNMGGIP 372
#                    ++ET     +++ +G D   E               I V P VH+ MGG+ 
#Sbjct: 452 KMEGIDGKALTETVARYNSLVSSGKDTDFERPNLPRALNEGNYYAIEVTPGVHHTMGGVM 511
#
#Query: 373 TNYKGQVLRHVNGQDQVVPGLYACGEAACASVHGANRLGANSLLDLVVFGR 423
#            + K +V+   N + QV+PGLY  GE     VHGANRLG N++ D++ FGR
#Sbjct: 512 IDTKAEVM---NAKKQVIPGLYGAGEVT-GGVHGANRLGGNAISDIITFGR 558
#</PRE>



######################################################
# global variables
save_date = []

######################################################
def main():

	global outfile
	input_id = sys.argv[1]
	input_seq = sys.argv[2]
	output_file_name = sys.argv[3]

	server = SOAPpy.SOAPProxy("http://www.pdb.org/pdb/services/pdbws")

	this_cutoff = float(0.005)
	this_matrix = 'BLOSUM62'
	this_output_format = '1'

	# String blastPDB (String structureId, String chainId, double eCutOff, String matrix, String outputFormat) 
	# soappy_output = server.blastPDB ( this_pdbid, this_chain, this_cutoff, this_matrix, this_output_format );

	# String blastPDB (String sequence, double eCutOff, String matrix, String outputFormat) 
	soappy_output = server.blastPDB ( input_seq, this_cutoff, this_matrix, this_output_format )

	soappy_output_lines = soappy_output.rsplit( "\n" )

	i = 0
	for inline in soappy_output_lines:

		#><a name = 5612></a>1E39:1:A|pdbid|entity|chain(s)|sequence

		if ( inline.find('><a name = ') > -1 ):
			bits = inline.split( '></a>' )
			bit1 = bits[1]
			bits2 = bit1.split( '|' )
			bit2 = bits2[0]
			this_pdbid = bit2[0:4]
			this_pdbid = this_pdbid.lower()
			this_chain = bit2[7:8]
			# Identities = 108/411 (26%), Positives = 186/411 (45%), Gaps = 68/411 (16%)
			j = i + 4
			inline2 = soappy_output_lines[j]
			if ( inline2.find('Identities = ') == -1 ):
				print 'Did not find the identities line where expected for sequence', input_seq[0:10], '... :', inline
			else:
				bits3 = inline2.split( ' ' )
				identity_string = bits3[4]
				identity = identity_string
				identity = identity.replace( ',', '' )
				identity = identity.replace( '(', '' )
				identity = identity.replace( ')', '' )
				identity = identity.replace( '%', '' )

				identity = float(identity)
				if (identity >= 30):

					soappy_output2 = server.getReleaseDates ( this_pdbid );

					bits1 = soappy_output2[0]
					bits2 = bits1.split( "\t" )
					this_date = bits2[1]

					if (this_date == 'N/A'):
						print 'PROBLEM :',this_pdbid, 'returned', this_date, 'when processing', input_id
					else:
						#print input_id, ':', this_pdbid, ':', this_date
						save_date.append( this_date )
		i = i + 1

	return_date = '9999-99-99'
	if ( len(save_date) > 0 ):
		for this_save_date in save_date:
			if (this_save_date < return_date):
				return_date = this_save_date

	outfile = open( output_file_name, "a" ) # open output file for writing in append mode
	output_line = input_id + ':::::' + return_date + ':::::' + "\r\n"
	outfile.write( output_line )
	outfile.close()

	return





if __name__ == "__main__":
	# Someone is launching this directly
	main()

