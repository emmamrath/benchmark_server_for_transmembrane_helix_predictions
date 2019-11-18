#!/usr/bin/env python

# python read_pdbtm_xml_all_types.py PDBTM__all_helical_transmembrane_proteins pdbtm_feb2012.all_pdbtm.TMH

# python read_pdbtm_xml_all_types.py [input-file] [output-file]

import sys
from xml.etree import ElementTree as ET
from xml.dom.minidom import parse, parseString
#import xml


# global variables
outfile = ''


######################################################
def getText(nodelist):
	rc = []
	for node in nodelist:
		if node.nodeType == node.TEXT_NODE:
			rc.append(node.data)
	return ''.join(rc)


######################################################
def handle_all_PDBTM(all_PDBTM):
	PDBTM = all_PDBTM.getElementsByTagName("PDBTM")[0]
	each_pdbtm = PDBTM.getElementsByTagName("pdbtm")
	for pdbtm in each_pdbtm:
		handle_pdbtm(pdbtm)


######################################################
def handle_pdbtm(pdbtm):
	#print pdbtm.getElementsByTagName("ID")
	pdbtm_ID = pdbtm.getAttribute("ID")
	pdbtm_chains = pdbtm.getElementsByTagName("CHAIN")
	#xml.etree.ElementTree.dump(pdbtm_chains)
	for chain in pdbtm_chains:
		handle_chain( pdbtm_ID, chain )


######################################################
def handle_chain( pdbtm_ID, chain ):

	# alpha:::::1,2,H,I,L,U
	# beta:::::1,2,B,I,L,U
	# coil:::::1,2,C,U
	# non_tm:::::1,2,H,I,L,U
	# tm_part:::::1:2,U
	# 1 = Side1
	# 2 = Side2
	# B = Beta-strand
	# H = alpha-helix
	# C = coil
	# I = membrane-inside
	# L = membrane-loop
	# U = unknown

	global output_file_name
	pdbtm_chain_CHAINID = chain.getAttribute("CHAINID")
	pdbtm_chain_NUM_TM = chain.getAttribute("NUM_TM")
	pdbtm_chain_TYPE = chain.getAttribute("TYPE")
	pdbtm_chain_TYPE_upper = pdbtm_chain_TYPE.upper()

	# pdbtm_chain_TYPE = alpha, beta, coil, non_tm, tm_part

	output_side_1 = ''
	output_side_2 = ''
	output_alpha_helix = ''
	output_beta_strand = ''
	output_coil = ''
	output_membrane_inside = ''
	output_membrane_loop = ''
	output_unknown = ''

	pdbtm_chain_regions = chain.getElementsByTagName("REGION")
	for region in pdbtm_chain_regions:

		segment = handle_region( pdbtm_ID, pdbtm_chain_CHAINID, region )
		if (segment != ''):
			segment_type = segment[0:1]
			segment_location = segment[1:]

			# segment_type = 1, 2, H, B, C, I, L, U

			if (segment_type == '1'):
				output_side_1 = output_side_1 + segment_location + ','
			if (segment_type == '2'):
				output_side_2 = output_side_2 + segment_location + ','
			if (segment_type == 'H'):
				output_alpha_helix = output_alpha_helix + segment_location + ','
			if (segment_type == 'B'):
				output_beta_strand = output_beta_strand + segment_location + ','
			if (segment_type == 'C'):
				output_coil = output_coil + segment_location + ','
			if (segment_type == 'I'):
				output_membrane_inside = output_membrane_inside + segment_location + ','
			if (segment_type == 'L'):
				output_membrane_loop = output_membrane_loop + segment_location + ','
			if (segment_type == 'U'):
				output_unknown = output_unknown + segment_location + ','

	output_segments = ''
	if (output_side_1 != ''):
		output_segments = output_segments + pdbtm_chain_TYPE_upper + '__side_1,,,,,' + output_side_1 + ';;;;;'
	if (output_side_2 != ''):
		output_segments = output_segments + pdbtm_chain_TYPE_upper + '__side_2,,,,,' + output_side_2 + ';;;;;'
	if (output_alpha_helix != ''):
		output_segments = output_segments + pdbtm_chain_TYPE_upper + '__alpha_helix,,,,,' + output_alpha_helix + ';;;;;'
	if (output_beta_strand != ''):
		output_segments = output_segments + pdbtm_chain_TYPE_upper + '__beta_strand,,,,,' + output_beta_strand + ';;;;;'
	if (output_coil != ''):
		output_segments = output_segments + pdbtm_chain_TYPE_upper + '__coil,,,,,' + output_coil + ';;;;;'
	if (output_membrane_inside != ''):
		output_segments = output_segments + pdbtm_chain_TYPE_upper + '__membrane_inside,,,,,' + output_membrane_inside + ';;;;;'
	if (output_membrane_loop != ''):
		output_segments = output_segments + pdbtm_chain_TYPE_upper + '__membrane_loop,,,,,' + output_membrane_loop + ';;;;;'
	if (output_unknown != ''):
		output_segments = output_segments + pdbtm_chain_TYPE_upper + '__unknown,,,,,' + output_unknown + ';;;;;'
	if (output_segments == ''):
		output_segments = '-'
	output_line = pdbtm_ID + '_' + pdbtm_chain_CHAINID + ':::::' + output_segments + ':::::' + "\r\n"
	outfile.write( output_line )

	return


######################################################
def handle_region( pdbtm_ID, pdbtm_chain_CHAINID, region ):

	pdbtm_chain_region_seq_beg = region.getAttribute("seq_beg")
	pdbtm_chain_region_pdb_beg = region.getAttribute("pdb_beg")
	pdbtm_chain_region_seq_end = region.getAttribute("seq_end")
	pdbtm_chain_region_pdb_end = region.getAttribute("pdb_end")
	pdbtm_chain_region_type = region.getAttribute("type")
	# print pdbtm_ID, pdbtm_chain_CHAINID, pdbtm_chain_region_pdb_beg, pdbtm_chain_region_pdb_end, pdbtm_chain_region_type #####

	segment = pdbtm_chain_region_type + pdbtm_chain_region_pdb_beg + '-' + pdbtm_chain_region_pdb_end

	return segment

#    <REGION seq_beg="1" pdb_beg="1"  seq_end="201" pdb_end="201"  type="1"/>		# 2agv_A
#    <REGION seq_beg="202" pdb_beg="202"  seq_end="223" pdb_end="223"  type="H"/>	# 2agv_A
#    <REGION seq_beg="224" pdb_beg="224"  seq_end="241" pdb_end="241"  type="2"/>	# 2agv_A
#    <REGION seq_beg="1" pdb_beg="1"  seq_end="2" pdb_end="2"  type="U"/>		# 3ag1_C
#    <REGION seq_beg="51" pdb_beg="51"  seq_end="67" pdb_end="67"  type="L"/>		# 2ahy_A	a loop, only half of it is a helix
#    <REGION seq_beg="1" pdb_beg="1"  seq_end="6" pdb_end="6"  type="I"/>		# 3abk_G	not a loop, not a helix

# Each CHAIN record contains one ore more REGION records, which locates the chain segment in the space relative to the membrane. 
# The type of REGION can be 1, 2, B, H, C, I, L and U for Side1, Side2, Beta-strand, alpha-helix, coil, membrane-inside, membrane-loop and unknown localizations, respectively. 
# Side1 and Side2 refers to the two sides of the membrane (we do not know which is inside or outside using the information only from the pdb file). 
# Membrane-inside is the inside part of a beta barrel. 
# Membrane- loop correspond to a region of the polypeptide chain which does not cross the membrane, just dips into the membrane (for example in aquaporins or potassium-channels). 

# The pdb_beg and pdb_end attributes contain the segment localization using the pdb numbering, 
# while the seq_beg and seq_end use the numbering in the sequence given in the SEQ record. 
# The sequence in SEQ record is generated by the alignment of sequences found in the SEQRES and ATOM records of the pdb file. 


######################################################
def main():

	global outfile
	input_xml_file = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )

	#pdbtm_xml_file = 'all_pdb_proteins.xml'
	#pdbtm_xml_file = 'all_transmembrane_proteins.xml'
	#pdbtm_xml_file = 'all_helical_transmembrane_proteins.xml'
	#pdbtm_xml_file = 'all_beta_barrel_transmembrane_proteins.xml'
	pdbtm_xml_file = input_xml_file

	datasource = open(pdbtm_xml_file)
	dom = parse(datasource)   # parse an open file

	handle_all_PDBTM(dom)

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

