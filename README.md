# Spheroidal Harmonics Decomposition with Fast Ellipsoidal Conformal, Quasi-Conformal Map, and Density-Equalizing Parameterizations

<img src = "https://github.com/garyptchoi/ellipsoidal-map/blob/main/cover.png" height="360" />

This code combines the novel spheroidal harmonics (SOH) [3] for analyzing 3D parameteric surfaces with the recent ellipsoidal parameterization techniques.

* **Fast Ellipsoidal Conformal Map (FECM)**: Compute an ellipsoidal conformal parameterization of a genus-0 closed surface using the method in [1].

* **Fast Ellipsoidal Quasi-Conformal Map (FEQCM)**: Compute an ellipsoidal quasi-conformal parameterization of a genus-0 closed surface with prescribed landmark constraints using the method in [1].

* **Ellipsoidal density-equalizing map (EDEM)**: Compute an ellipsoidal density-equalizing map of a genus-0 closed surface onto a prescribed ellipsoid [2].

* **Ellipsoidal density-equalizing quasi-conformal map (EDEQ)**: Compute an ellipsoidal density-equalizing quasi-conformal map of a genus-0 closed surface, with the ellipsoidal shape automatically determined [2].

If you use this code in your work, please cite the following papers:

[1] G. P. T. Choi, 
    "[Fast ellipsoidal conformal and quasi-conformal parameterization of genus-0 closed surfaces.](https://doi.org/10.1016/j.cam.2024.115888)"
    Journal of Computational and Applied Mathematics, 447, 115888, 2024.

[2] Z. Lyu, L. M. Lui, and G. P. T. Choi, "[Ellipsoidal Density-Equalizing Map for Genus-0 Closed Surfaces.](https://arxiv.org/pdf/2410.12331)" Preprint, arXiv:2410.12331, 2024.

[3] M. Shaqfa, and W. M. van Rees "[Spheroidal harmonics for generalizing the morphological decomposition of closed parametric surfaces.](https://arxiv.org/abs/2407.03350)"
    Journal of Construction and Building Materials, 2024.



Copyright (c) 2023-2024, Mahmoud Shaqfa, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi, Wim M. van Rees

===============================================================

Usage:
* To run the demos, follow the instructions provided per the `AP_SOH_main.m` file.

This work is complementary to the SOH [3] paper, and their Python3 codes can be found on [this link](https://github.com/msshaqfa/spheroidal-harmonics). For open surfaces, we also proposed the hemispheroidal harmonics (HSOH), where you can find our codes online on [this link.](https://github.com/msshaqfa/hemispheroidal-harmonics)
