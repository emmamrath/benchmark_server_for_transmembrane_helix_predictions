#!/usr/bin/env python

# python ../extract_mpstruc_category_for_pdbids_from_html.py mpstruc_expanded_list.html mpstruc_category.pdbid

# python extract_mpstruc_category_for_pdbids_from_html.py [input-file] [output-file]

import sys
import re


#    <a id="id_MONOTOPIC_MEMBRANE_PROTEINS:Cyclooxygenases"></a>
#      <div class="sectionNameBody">
#        Cyclooxygenases
#    </div>
#
#    <div class="pdbPage">
#      <strong><a href="http://www.pdb.org/pdb/search/structidSearch.do?structureId=1PTH" title="Visit the PDB entry page.">1PTH</a></strong>
#    </div>
#
#    <a id="id_TRANSMEMBRANE_PROTEINS_ALPHAHELICAL:Electron_Transport_Chain_Complexes_Complex_II"></a>
#      <div class="sectionNameBody">
#        Electron Transport Chain Complexes: Complex II
#    </div>
#
#    <div class="pdbPage">
#      <strong><a href="http://www.pdb.org/pdb/search/structidSearch.do?structureId=1FUM" title="Visit the PDB entry page.">1FUM</a></strong>
#    </div>
#
#      <strong>Ram Prostaglandin H2 synthase-1 (COX-1): <em>Ovis aries</em>,&nbsp;2.61&nbsp;A</strong>
#      <br>
#      1EQG is complex with ibuprofen. <br> Complex with flurbiprofen, 2.70&nbsp;A: <a href="http://www.pdb.org/pdb/search/structidSearch.do?structureId=1EQH">1EQH</a> <br> Complex with flurbiprofen methyl ester, 2.75&nbsp;A: <a href="http://www.pdb.org/pdb/search/structidSearch.do?structureId=1HT5">1HT5</a> <br> Complex with alclofenac, 2.69&nbsp;A: <a href="http://www.pdb.org/pdb/search/structidSearch.do?structureId=1HT8">1HT8</a>
#    </div>
#  </td>
#  <td>
#    <div class="pdbPage">
#      <strong><a href="http://www.pdb.org/pdb/search/structidSearch.do?structureId=1EQG" title="Visit the PDB entry page.">1EQG</a></strong>
#    </div>
#
#        <a id="tableDescription"></a>
#        <h3>
#          Description of Table
#        </h3>


# global variables
save_pdbid = []

######################################################
def main():

	global outfile
	input_file = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	current_category = ''

	i = 0
	reached_end = 0
	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):

			if (reached_end == 0):

				posn = inline.find( 'Description of Table' )
				if (posn > -1):
					reached_end = 1
				else:

					posn = inline.find( '<div class="sectionNameBody">' )
					if (posn > -1):
						current_category = inlines[i+1]
						current_category = current_category.strip()
						current_category = current_category.replace( '<sup>', '' )
						current_category = current_category.replace( '</sup>', '' )
						current_category = current_category.replace( '<sub>', '' )
						current_category = current_category.replace( '</sub>', '' )
						current_category = current_category.replace( '<em>', '' )
						current_category = current_category.replace( '</em>', '' )

					posn = inline.find( 'http://www.pdb.org/pdb/search/structidSearch.do?structureId=' )
					if (posn > -1):
						posn2 = inline.find( 'title="Visit the PDB entry page.">' )
						if (posn2 > -1):
							bits = inline.rsplit( 'title="Visit the PDB entry page.">' )
							bit1 = bits[1]
							bits2 = bit1.rsplit( '</a>' )
							current_pdbid = bits2[0]
							current_pdbid = current_pdbid.lower()
							output_line = current_pdbid + ':::::' + current_category + ':::::' + "\r\n"
							outfile.write( output_line )
						else:
							bits = inline.rsplit( 'http://www.pdb.org/pdb/search/structidSearch.do?structureId=' )
							for bit1 in bits:
								current_pdbid = bit1[0:4]
								bit3 = bit1[4:6]
								if (bit3 == '">'):
									current_pdbid = current_pdbid.lower()
									output_line = current_pdbid + ':::::' + current_category + ':::::' + "\r\n"
									outfile.write( output_line )
		i = i + 1
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

