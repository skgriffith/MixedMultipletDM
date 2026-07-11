Mathematica files necessary to reproduce the results of "Minimal Dark Matter: Generalized Framework and Direct-Detection Sensitivity" (arxiv #2602.17764).

The files must be placed in the same folder to run correctly.

mixed_multiplet.m contains the basic functions necessary to compute the effective cross section (see paper for details).

overlap_integrals.m contains pre-computed overlap integrals (eq. D4 in the paper) for the relevant bound states.

evaluate_CrossSections.m computes the effective cross section for a given y value (small, mid, large) and given combination of multiplets.

analyze.nb then uses the output of evaluate_CrossSections.m to compute the necessary particle mass to match the observed dark matter abundance.

Please contact me at griffith.1037@osu.edu for more detailed instructions on how to use the code.
