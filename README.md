# lampinen_mrm_2019

Supplementary code to the manuscript: "Towards unconstrained compartment modeling in white matter using diffusion-relaxation correlation MRI with tensor-valued diffusion encoding" accepted in MRM 2020.

---Instructions for use---

1) Clone this repository
2) Clone the forked repository https://github.com/belampinen/md-dmri into a folder next to this repository.
3) Run the script setup.m for setting up paths.

To process data according to in Lampinen et al (2020) MRM, run the script smr_process_example after changing the content under 'User provided details'.
NB: the processing requires that Elastix (www.elastix.org) is installed.



---FOLDER CONTENT----
1) estimates: 		Parameter estimates for 15 subjects. 

2) figures: 		Code for creating Figures 2, 3 and 4, as well as Supporting
			Information Figures S1, S4 and S5

3) model: 		Code for the two-compartment diffusion-relaxation model.
			Code for model fitting using the MDM framework (dtd_smr*)

4) optimization: 	Protocol optimization based on the CRLB using SOMA.
		        Run from smr_optimize_master.

5) tools: 		Supporting functions used in various scripts.

6) waveforms: 		The diffusion encoding waveforms used for the 'in vivo protocol'.
7) dependencies:	Additional necessary scripts