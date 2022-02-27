# DISintigrate

This is a refactoring of [DISPred](https://dispred.hepforge.org) by J. Ferrando with simpler 
required dependencies, a different configuration mechanism, and more options for storing 
calculated results.

Not all features of the original code have been successfully preserved, as this version has been
used primarily for producing neutrino interaction cross sections for use in experiment simulations,
so features not directly related to such may have been been omitted or not be functional.

Furthermore, this code is currently in a somewhat untidy state internally, as it was composed
hastily as needed, and has not yet undergone complete clean up. 

## Dependencies

This rewrite removes the need for ROOT. Like the original, it does require LHAPDF 5. For 
perfoming integrals it can use either the GNU Scientific Library, or Cubpack++, although the
required version of the latter is not yet publicly available. For writing results, it can use the
Photospline library to fit spline surfaces. A pair of small, header-only libraries are used for
one-dimensional integration and command line option parsing.

The configuration script for this project has a best-effort ability to download (and when
necessary, build) the required dependencies itself. 

Summary of dependencies:

- REQUIRED: [LHAPDF 5](https://lhapdf.hepforge.org/lhapdf5/)
  ([source code](https://lhapdf.hepforge.org/downloads?f=old/lhapdf-5.9.1.tar.gz))
  This is the underlying source used for physics data, specifically the parton distribution functions. 
  This dependency can be fetched and compiled automatically, but all of its dependencies
  (primarily, a Fortran compiler) must already be satisfied for this to work.  
- REQUIRED: [cl_options](https://github.com/cnweaver/cl_options)
  This library is used for parsing command line arguments. 
  Because it is uncommon, but simple, the configuration script can fetch it automatically. 
- REQUIRED: [AdaptiveQuad](https://github.com/cnweaver/AdaptiveQuad)
  This library is used for one dimensional integrals, as it tends to be much more effiicent than the
  GSL one-dimensional algorithms, for this program's needs.  
  Because it is uncommon, but simple, the configuration script can fetch it automatically. 
- OPTIONAL: The [GNU Scientific Library](https://www.gnu.org/software/gsl/)
  This is the default choice for performing two-dimensional integrals. 
- OPTIONAL: Cubpack++ 
  This is an alternate choice for performing two-dimensional integrals, which has far better performance. 
  A version of this library sufficiently modernized to be used by this program is not currently publicly available.
- OPTIONAL: [Photospline 2](https://github.com/icecube/photospline)
  This provides an alternate output format, of multi-dimensional tensor-product b-spline surfaces. 
  These provide interpolation, gradients, ease of sampling, and may be smaller than the tabulated data from which they are produced.  

## Compilation

First, run the provided `configure` script. This should detect or automatically handle dependencies. 
It will fail if a depenency cannot be located, in which case manual specification may be required. 
See the output of `configure --help` for a full listing of its options, particularly for specifying 
depenencies. 

## Usage

By default the program does nothing except initiallize the PDFs and exit. To make it do anything 
productive, one must at a minimum specify whether total, singly-differential, or doubly-differential
cross sections should be produced. All other required parameters have default options, but the 
user is cautioned to check them carefully in order to compute results which are specificly desired.

The `--total` (`-t`), `--singly-differential` (`-s`), and `--doubly-differential` (`-d`)
options are used to select the form of cross sectio to compute. Any combination of these may be
used together. The corresponding `--totalXSOutput`, `--singleXSOutput`, and
`--doubleXSOutput` options should be used to set the names of the ouput files which will be
written. The `--format` (`-f`) option can be used to set the output format used, and the
`--unit` (`-u`) sets the units of the output. 

The main physics options which should be set are as follows:

- `--interaction` (`-i`): The interaction process to computed, either `CC` for charged-current or `NC` for neutral-current. 
- `--particle` (`-p`): The particle or anti-particle type to use, `nu` or `nubar`
- `--pdf`: The parton distibution function data set from LHAPDF to use for the calculation
- `--protonFraction`: The composition of the target. A value of 1 corresponds to pure protons, 0 to pure neutrons, and 0.5 to an 'isoscalar' mixture. 
- `--leptonFlavor`: The flavor of the incoming particle. This is only important if `--enforceKinematics` is also enabled an the interaction is charged-current, such that the final-state lepton mass affects the result.
- `--Q2min`: The minimum value of `Q^2` which contributes to the calculation. 

Also important are the options for the tabulation:

- `--Emin`: The minimum incoming particle energy for which to tabulate.
- `--Emax`: The maximum incoming particle energy for which to tabulate.
- `--Esteps`: The number of incoming energy values to tabulate, including the minimum and maximum. 
- `--xysteps`: The number of values to tabulate in both `x` and `y` for doubly-differential cross sections.
- `--squareSingleTable`: This option is relevant only to tabulation of singly-differential cross sections, 
  and switches tabulation to be in the 'square' form described in 
  [https://arxiv.org/abs/2112.13804 appendix D](https://arxiv.org/abs/2112.13804), which is more
  useful for interpolation.
- `--singleDiffPeakHack`: This option is relevant only to tabulation of singly-differential cross sections,
  and performs a replacement of the E_out == E_in entries tabulated to better capture the sharp 
  peak of the cross section at hig energies, as described in 
  [https://arxiv.org/abs/2112.13804 appendix D](https://arxiv.org/abs/2112.13804)

Finally, some of the tabulation options are only relevant if the program is compiled with the optional photospline support:

- `--Eknots`: Sets the number of knots used in the incoming particle energy dimension of the spline surface. 
- `--xyknots`: Sets the number of knots used in the `x` and `y` dimensions of the doubly-differential spline surface. 
- `--zknots`: Sets the number of knots used in the `z` dimensions of the singly-differential spline surface when `--squareSingleTable` is enabled. 
- `--single-regularization` (`--rs`): Sets the strength of regularization for fitting the singly-differential spline surface. 

See the `examples/csms_check` directory for a set of configurations which can be used to
reproduce previously published results from the original DISPred code. These serve both to
provide some minimal verification that this code is functioning correctly, and to provide a starting
point for users to produce additional tables. 
