#!/usr/bin/env python

# python calc_avg_helix_end_difference.py /var/www/benchmark/data/alpha_polytopic_bitopic.opm_tm_subunits /var/www/benchmark/data/alpha_polytopic_bitopic.pdbtm_segments_topology

# python calc_avg_helix_end_difference.py [input-file-1] [input-file-2]

import sys
import re

# file 1 :
# 3tx3_A:::::OPM-TRANSMEMBRANE-SUBUNIT,,,,,30-52,80-102,143-156,208-218,;;;;;:::::
# 3ug9_A:::::OPM-TRANSMEMBRANE-SUBUNIT,,,,,89-109,121-138,157-175,188-207,212-232,252-274,281-304,;;;;;:::::
# 3ukm_A:::::OPM-TRANSMEMBRANE-SUBUNIT,,,,,20-44,102-116,129-157,181-206,211-224,243-268,;;;;;:::::

# file 2 :
#3tx3_A:::::_UUUiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiHHHHHHHHHHHHHHHHHHooUUUUUUUUUUUUUUUoooooHHHHHHHHHHHHHHHHHHHHHHiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii:::::
#3ug9_A:::::_______________________UUUUUUUUUUUUUUUUUUUUUUUUUooooooooooooooooooooooooooooooooooooooooooHHHHHHHHHHHHHHHHHHHUUUUUUUUiiiHHHHHHHHHHHHHHHHHHoooooooooooooooooHHHHHHHHHHHHHHHHHiiiiiiiiiiiiiiiiiHHHHHHHHHHHHHHHHHHooooHHHHHHHHHHHHHHHHHHHHiiiiiiiiiiiiiiiiiiiiiiHHHHHHHHHHHHHHHHHHHoooooooooooHHHHHHHHHHHHHHHHHHHiiiiiiiiiiiiiiiiiiiiiiiiUUUUUUiiiiiiiiiiUUUUUUUUUUUUUU:::::
#3ukm_A:::::__________UUUUUUUiiiiHHHHHHHHHHHHHHHHHHHHHHHoooooooooooooooooooooooooooooooooooooooooooooooooUUUUUUoooooLLLLLLLLLLLLLLLLLoooooooHHHHHHHHHHHHHHHHHHHHHHHHHHHHHiiiiiiiiiiiUUUUUUiiiiiiiiHHHHHHHHHHHHHHHHHHHHHHHHooooooLLLLLLLLLLLLLLLLLooooooooooooooHHHHHHHHHHHHHHHHHHHHHHHHHiiiiiiiiiiiiiUUUUUUUUU:::::
#3um7_A:::::UUUUUUUUUUUUUUUUUUUUUUUUiiiiiiHHHHHHHHHHHHHHHHHHHHHHHoooooooooooooooooooooooooooooooooooooooooooooooooooooUUUUUoooooLLLLLLLLLLLLLLLLLoooooooHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHiiiiiiiiiiiiiiUUiiiiiiiHHHHHHHHHHHHHHHHHHHHHHHHoooooLLLLLLLLLLLLLooooooooooooooooHHHHHHHHHHHHHHHHHHHHHHHHHHiiiiiiiiiiUUUUUUUUUUUUUUUUUUU:::::
#3zrs_A:::::UUUUUUUUUUUiiiiiiiiiiiiiiiUUUUUUiiiiiiiiiiiiiiHHHHHHHHHHHHHHHHHHHHHHooooooooooooooooooLLLLLLLLLLLLooooooooooHHHHHHHHHHHHHHHHHHHHHHHHiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiUUUiiiiiiiiiiiiiii:::::


# global variables
save_id = []
save_seq = []
found_seq = ''


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()

	# 3ukm_A:::::__________UUUUUUUiiiiHHHHHHHHHHHHHHHHHHHHHHHoooooooooooooooooooooooooooooooooooooooooooooooooUUUUUUoooooLLLLLLLLLLLLLLLLLoooooooHHHHHHHHHHHHHHHHHHHHHHHHHHHHHiiiiiiiiiiiUUUUUU
	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			save_id.append( this_id )
			save_seq.append( this_seq )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	num_ends = 0
	num_diff_residues = 0

	# 3tx3_A:::::OPM-TRANSMEMBRANE-SUBUNIT,,,,,30-52,80-102,143-156,208-218,;;;;;:::::
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_subunits = bits[1]
			this_subunits = this_subunits.replace( ';;;;;', '' )
			bits2 =  this_subunits.rsplit( ',,,,,' )
			this_pairs = bits2[1]
			bits3 = this_pairs.rsplit( ',' )
			for this_pair in bits3:
				if ((this_pair.find('--') == -1) and (this_pair != '')):
					if ( this_pair[0:1] == '-' ):
						this_pair = this_pair[1:]
					bits4 = this_pair.split( '-' )
					this_from = int(bits4[0])
					this_to = int(bits4[1])

					found_it = 0
					i = 0
					for this_save_id in save_id:
						if (this_save_id == this_id):
							this_save_seq = save_seq[i]
							found_it = 1
						i = i + 1

					if (found_it == 1):
						midpoint = round((this_to - this_from - 1) / 2) + this_from - 1
						midpoint = int(midpoint)
						if ((this_save_seq[midpoint] == 'H') or (this_save_seq[midpoint] == 'L')):

							# it's a matching membrane helix, now find how different the ends are
							j = midpoint - 1 # start from the midpoint, but pairs are counted from 1 whereas string numbering is accessed as if from zero
							end_search = 0
							while (end_search == 0):
								j = j - 1
								k = j + 1
								if ((this_save_seq[j:k] == 'H') or (this_save_seq[j:k] == 'L')):
									end_search = 0
								else:
									end_search = 1
									num_ends = num_ends + 1
									this_diff = abs(k - this_from)
									num_diff_residues = num_diff_residues + this_diff

							j = midpoint - 1
							end_search = 0
							while (end_search == 0):
								j = j + 1
								k = j + 1
								if ((this_save_seq[j:k] == 'H') or (this_save_seq[j:k] == 'L')):
									end_search = 0
								else:
									end_search = 1
									num_ends = num_ends + 1
									this_diff = abs(k - this_to)
									num_diff_residues = num_diff_residues + this_diff

	print '========================================'
	print 'Number of helix ends =', num_ends
	avg_diff = num_diff_residues / num_ends
	print 'Average difference in helix ends =', avg_diff
	print '========================================'

if __name__ == "__main__":
	# Someone is launching this directly
	main()

# ========================================
# Number of helix ends = 4438
# Average difference in helix ends = 2
# ========================================

