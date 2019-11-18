#!/usr/bin/env python

# python read_mempype_soap_for_list.py 3o7q_A.1line 3o7q_A.1line.mempype

# python read_mempype_soap_for_list.py [input-file] [output-file]

# http://mu2py.biocomp.unibo.it/mempype/rpc/call/soap?WSDL

import sys
import re
import hashlib
from pysimplesoap.client import SoapClient, SoapFault
#import pysimplesoap 


# global variables


######################################################
def main():

	global outfile
	input_file = sys.argv[1]
	output_file_name = sys.argv[2]

	# create a SOAP client
	client = SoapClient(wsdl="http://mu2py.biocomp.unibo.it/mempype/rpc/call/soap?WSDL")

	outfile = open( output_file_name, "w" )

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.split( ":::::" )
			this_id = bits[0]
			this_fastaid = bits[1]
			this_seq = bits[2]

			this_seq_hashlib = hashlib.sha1()
			this_seq_hashlib.update( this_seq )
			this_seq_hashlib_string = this_seq_hashlib.hexdigest()
			print this_seq_hashlib_string
			print

			# call SOAP method

			#response = client.get_prediction_results( seqhash='ce90a8c334a89b821ca6da78ce089eb88b4ec842' ) # <== this works, is the example given by andrea@biocomp.unibo.it
			#response = client.get_prediction_results( seqhash=this_seq_hashlib_string )
			#response = client.get_prediction_results( seqhash='4eaf2077845b24165649d6f5ad1bcf8c23738ac5' ) # <== this doesn't work, is the sha1 hash for my sequence, is not in their database?
			this_seq_hashlib_string = 'ce90a8c334a89b821ca6da78ce089eb88b4ec842'
			response = client.get_prediction_results( seqhash=this_seq_hashlib_string )

			print response

			# {'signal_peptide': u'1-32', 'prediction_status': u'Done', 'prediction_summary': u'Cell Membrane, Globular', 'cytoplasmic_region': u'NO', 'memloci_cell_membrane_score': u'99%', 'memloci_organellar_membrane_score': u'-88%', 'non_cytoplasmic_region': u'NO', 'memloci_internal_membrane_score': u'-21%', 'transmembrane_region': u'NO', 'GPI-anchor': u'NO'}

			#try:
			#	status = response['status']
			#	print dict(xml_request=client.xml_request, xml_response=client.xml_response, status = status )
			#except SoapFault, e:
			#	print dict(xml_request=client.xml_request, xml_response=client.xml_response, error=str(e))

			#output_line = output_line + "\r\n"
			#outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

