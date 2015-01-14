# McFee
Python for Richard McFee's 1959 paper 'Optimum Input Leads for Cryogenic Apparatus'
In McFee's 1959 paper on current lead optimization in cryogenic applications, heat load and length-to-area ratio of conductor
are found for an optimized current carrying lead that goes from a warm end to a cold end.

Summary of Files:

burnout_analysis.py finds the minimum conductor cross-sectional area needed for given parameters.

mcfee_analysis.py finds conductor heat load at the cold end along with the optimal length to area ratio for given parameters.

McFee_Functions.py holds the functions used by mcfee_analysis.py

heat_xfer.py holds functions and constants relating to physical properties of OFHC copper.

All software files are licensed by the MIT License.