To run Cellrank using a multivelo results object, run these scripts in this order:
1. compute_vk.py
2a. set_macrostates.py #if using previous knowledge to set macrostates, in this case, latent time (see how this was done in run_multivelo_pt2.py)
2b. infer_fate.py #if want to use Cellrank to predict macrostates
3. get_fate_probabilities.py
4. uncover_driver_genes.py
