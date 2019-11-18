#!/usr/bin/env python

# python parse_philius_xml.py PDBTM__all_helical_transmembrane_proteins pdbtm_feb2012.all_pdbtm.TMH
# python parse_philius_xml.py alpha_polytopic_biotopic_noblanks.fasta.philius.xml alpha_polytopic_biotopic_noblanks.philius

# python parse_philius_xml.py [input-file] [output-file]

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
def handle_all_philiusProteinSequenceResponse(all_philiusProteinSequenceResponse):

	each_philiusProteinSequenceResponse = all_philiusProteinSequenceResponse.getElementsByTagName("philiusProteinSequenceResponse")
	for philiusProteinSequenceResponse in each_philiusProteinSequenceResponse:
		fastaHeader_element = philiusProteinSequenceResponse.getElementsByTagName("fastaHeader")[0]
		fastaHeader = getText(fastaHeader_element.childNodes)
		this_id = fastaHeader.replace( '>', '' )
		psa_element = philiusProteinSequenceResponse.getElementsByTagName("psa")[0]
		sequence = psa_element.getAttribute("sequence")
		segmentList = psa_element.getElementsByTagName("segmentList")[0]
		each_philiusSegment = segmentList.getElementsByTagName("philiusSegment")
		# <philiusSegment typeConfidence="0.93" typeString="Cytoplasmic" type="2" end="8" start="1" sequenceID="1483985" id="2056516"/>
		# <philiusSegment typeConfidence="0.92" typeString="Transmembrane Helix" type="3" end="30" start="9" sequenceID="1483985" id="2056517"/>
		# <philiusSegment typeConfidence="0.87" typeString="Non-Cytoplasmic" type="1" end="45" start="31" sequenceID="1483985" id="2056518"/>
		prev_end = 0
		segments = ''
		for philiusSegment in each_philiusSegment:
			typeString = philiusSegment.getAttribute("typeString")
			this_start = philiusSegment.getAttribute("start")
			this_end = philiusSegment.getAttribute("end")
			this_start = int(this_start)
			this_end = int(this_end)
			if ((prev_end + 1) != this_start):
				print 'ERROR :', this_id, 'Not reading sequential sequences.'
			this_char = '_'
			if (typeString == 'Cytoplasmic'):
				this_char = 'i'
			elif (typeString == 'Non-Cytoplasmic'):
				this_char = 'o'
			elif (typeString == 'Transmembrane Helix'):
				this_char = 'H'
			elif (typeString == 'Signal Peptide'):
				this_char = 'S'
			else:
				print 'Unknown typeString ', typeString
			for i in range( this_start, (this_end + 1) ):
				segments = segments + this_char
			prev_end = this_end
		if (prev_end != len(sequence)):
			print 'ERROR :', this_id, 'Sequence length is', len(sequence), 'last segments ends at', prev_end
		if (segments == ''):
			segments = '-'
		output_line = this_id + ':::::' + segments + ':::::' + "\r\n"
		outfile.write( output_line )

	return


######################################################
def main():

	global outfile
	input_xml_file = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )

	pdbtm_xml_file = input_xml_file

	datasource = open(pdbtm_xml_file)
	dom = parse(datasource)   # parse an open file

	handle_all_philiusProteinSequenceResponse(dom)

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

