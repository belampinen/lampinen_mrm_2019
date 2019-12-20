# lampinen_mrm_2019

Supplementary code to the manuscript: "Towards unconstrained compartment modeling in white matter using diffusion-relaxation correlation MRI with tensor-valued diffusion encoding" submitted to MRM in 2019.

Run the script setup.m for setting up paths.

-Contains 6 folders-
1) estimates: 		Parameter estimates for 15 subjects. 

2) figures: 		Code for creating Figures 2, 3 and 4, as well as Supporting
			Information Figures S1, S4 and S5

3) model: 		Code for the two-compartment diffusion-relaxation model.
			Code for model fitting using the MDM framework (dtd_smr*)

4) optimization: 	Protocol optimization based on the CRLB using SOMA.
		        Run from smr_optimize_master.

5) tools: 		Supporting functions used in various scripts.

6) waveforms: 		The diffusion encoding waveforms used for the 'in vivo protocol'.
