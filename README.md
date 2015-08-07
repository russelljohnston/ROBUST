# ROBUST
Module to compute the Tc and Tv completeness statistics

Based on the following publications:

http://arxiv.org/abs/astro-ph/0010357  (R01)

http://arxiv.org/abs/astro-ph/0703040 (JTH07)

The module ROBUST.F90 contains the following routines:

tctv_R01: Computes the Rauzy (2001) Tc completeness statistic as well as the 
	  JTH07 variant, Tv. 
	  
tctv_JTH07: Computes the Tc and Tv statistics taking into about doubly truncated data
	   in the form of bright and faint apparemnt magniutde limits.
	  
SORT3.F90:  based on  Numerical Recipes, Created by  http://www.science-softcon.de/autochem/box-dir/box_html/998.html


