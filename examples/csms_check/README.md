# CSMS Verification Example

The configurations in this directory are intended to reproduce (a subset of) the results published
in tables I and II of 
[The high energy neutrino cross-section in the Standard Model and its uncertainty](https://doi.org/10.1007/JHEP08%282011%29042)
(Also [on the arXiv](https://arxiv.org/abs/1106.3723v2)).  

It should be noted that all of the entries in the published tables are rounded to two significant
figures (reasonably, given the argument uncertainties presented there), so the un-rounded data 
produced directly by running these examples will not be identical, but should agree when 
rounding is applied.

The `generate_all.sh` script in this directory will generate all four of the particle 
type/interaction type combinations. It can accept an option, `--parallel`, which will cause it to
(very simplistically) run them all concurrently. This has no regard for available CPU and memory
capacity, and may not be desirable on less powerful computers. 
