#!/usr/bin/env python

# python parse_dastmfilter_results.py dastmfilter_results2.txt dastmfilter_noblanks_input2.dastmfilter

# python parse_dastmfilter_results.py [input-file] [output-file]

import sys
import re


# input file :
# <html><head>
# <meta http-equiv="content-type" content="text/html; charset=UTF-8">
# <title> DAS-TMfilter prediction results </title>
#
# </head><body bgcolor="#FFFFFF">
#
# <xmp>
# Calculating prediction for the following proteins
# with reference library 08:
#
# &gt;2fyn_B
# &gt;2vt4_B
# &gt;3eff_K
# &gt;3e86_A
#
# === Result of the prediction ===
#
# &gt;1q16_C
# # TMH:  5 Q: trusted
# @   19   3.654 core:   11 ..   23 8.134e-04 !!! Warning! Potential signal peptide.
# @   62   4.625 core:   54 ..   70 2.640e-05
# @  102   4.160 core:   97 ..  109 1.360e-04
# @  133   5.391 core:  125 ..  144 1.765e-06
# @  194   5.774 core:  181 ..  205 4.568e-07
#
# &lt;-------- end of list --------&gt;
#
#
# &gt;1h6i_A
# # TMH:  7 Q: trusted
# @   25   4.679 core:   15 ..   32 2.603e-05
# @   60   2.862 core:   55 ..   62 1.590e-02
# @   84   3.530 core:   78 ..  115 1.503e-03 Twin peaks - two TMH with a short linker
# @  107   3.885 core:   78 ..  115 4.298e-04
# @  149   3.482 core:  143 ..  153 1.780e-03
# @  176   3.886 core:  170 ..  184 4.289e-04
# @  222   4.261 core:  214 ..  230 1.138e-04
#
# &lt;-------- end of list --------&gt;
#
#
# &gt;2hyd_A
# # TMH:  6 Q: trusted
# @   25   3.499 core:   18 ..   37 3.601e-03
# @   68   4.690 core:   60 ..   74 5.386e-05
# @  153   4.425 core:  142 ..  182 1.371e-04 Twin peaks - two TMH with a short linker
# @  173   4.879 core:  142 ..  182 2.762e-05
# @  265   3.368 core:  259 ..  291 5.718e-03
# !!! Warning! Weak twin peak. Ignored. @  283   3.145
# !!! Warning! Peak assignment overruled. Corrected no. of TMH:  5
#
# &lt;-------- end of list --------&gt;
#
#
# &gt;3arc_K
# # TMH:  1 Q:  0.93
# @   23   6.120 core:   12 ..   30 2.213e-08
#
# &lt;-------- end of list --------&gt;
#
#
# &gt;1okc_A
# # TMH:  0 Q: trusted !!! Warning! Non-TM protein!
#
# &lt;-------- end of list --------&gt;




# global variables

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	output_file_name = sys.argv[2]

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	outfile = open( output_file_name, "w" )

	seen_start_of_segments_info = 0
	curr_id = ''
	curr_seq = ''
	curr_signalpeptide = ''
	prev_was_twin_peaks = 0
	prev_from = 0
	prev_to = 0
	prev_peak = 0
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			if (seen_start_of_segments_info == 0):
				if (inline == '=== Result of the prediction ==='):
					seen_start_of_segments_info = 1
			else:
				# &gt;1q16_C
				is_a_new_id = 0
				if (len(inline) >= 10):
					if (inline[0:4] == '&gt;'):
						is_a_new_id = 1

				if (is_a_new_id == 1):

					# finish off previous id

					if (curr_id != ''):
						if (curr_seq == ''):
							curr_seq = '-' 
						if (curr_signalpeptide == ''):
							curr_signalpeptide = '-'
						output_line = curr_id + ':::::SIGNALPEPTIDE,,,,,' + curr_signalpeptide + ';;;;;TMH,,,,,' + curr_seq + ';;;;;:::::' + "\r\n"
						outfile.write( output_line )

					# initialise this id
					curr_id = inline[4:10]
					curr_seq = ''
					curr_signalpeptide = ''
					prev_was_twin_peaks = 0
					prev_from = 0
					prev_to = 0
					prev_peak = 0

				# @   19   3.654 core:   11 ..   23 8.134e-04 !!! Warning! Potential signal peptide.
				# @   68   4.690 core:   60 ..   74 5.386e-05
				# @  153   4.425 core:  142 ..  182 1.371e-04 Twin peaks - two TMH with a short linker
				# @  173   4.879 core:  142 ..  182 2.762e-05

				if (inline[0:1] == '@'):
					this_from = inline[20:25]
					this_from = this_from.strip()
					this_to = inline[28:33]
					this_to = this_to.strip()
					this_peak = inline[1:6]
					this_peak = this_peak.strip()

					if (inline.find( 'Warning! Potential signal peptide' ) > -1):
						curr_signalpeptide = curr_signalpeptide + this_from + '-' + this_to + ','
					else:

						if (prev_was_twin_peaks == 1):
							prev_to = str(int((int(this_peak) - int(prev_peak)) / 2) + int(prev_peak))
							this_from = str(int(prev_to) + 2)
							curr_seq = curr_seq + prev_from + '-' + prev_to + ',' + this_from + '-' + this_to + ','
							prev_was_twin_peaks = 0
						else:
							if (inline.find( 'Twin peaks - two TMH with a short linker' ) > -1):
								prev_was_twin_peaks = 1
							else:
								curr_seq = curr_seq + this_from + '-' + this_to + ','
								prev_was_twin_peaks = 0

					prev_from = this_from
					prev_to = this_to
					prev_peak = this_peak

	if (curr_id != ''):
		if (curr_seq == ''):
			curr_seq = '-' 
		if (curr_signalpeptide == ''):
			curr_signalpeptide = '-'
		output_line = curr_id + ':::::SIGNALPEPTIDE,,,,,' + curr_signalpeptide + ';;;;;TMH,,,,,' + curr_seq + ';;;;;:::::' + "\r\n"
		outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

