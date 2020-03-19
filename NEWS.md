# nmathresh 0.1.5

* Fix: `recon_vcov` now allows `prior.prec = 0`, reconstructing the covariance matrix from a frequentist analysis of the contrast data (or a Bayesian analysis with improper flat priors).
* Fix: Bug when excluding treatment 1 from the decision set where incorrect base-case optimal treatment could be chosen.

# nmathresh 0.1.4

* Feature: Diagnostics from `recon_vcov` when using NNLS are now returned as additional attributes to the matrix (as well as being printed).
* Updated references to Phillippo et al. JRSS:A paper, which is now in print.

# nmathresh 0.1.3

* Fix: Previous release introduced a bug where a single column passed to `add.columns` resulted in an error.
* Feature: Warn when arguments passed to `nma_thresh` are mismatched (`nmatype` options with `X`, `mu.design`, and `delta.design`), indicating possible user error.

# nmathresh 0.1.2

* Fix: Additional columns specified with `add.columns` in `thresh_forest` are now sorted along with the rest of the table when `orderby` is given.
* Fix: Bugs in `thresh_forest` with infinite or undefined data estimates and variances.
* Fix: Covariance matrices created using `Matrix` package are acceptable always (not just in the block diagonal case).

# nmathresh 0.1.1

* Fix DOI brackets in DESCRIPTION for CRAN.

# nmathresh 0.1.0

* Feature: Provide additional columns in the table produced by `thresh_forest` using argument `add.columns`.
* Feature: Thresholds for decisions based on minimal clinically important difference. Some internal API changes required, but these should not affect the end user.

# nmathresh 0.0.5

* Update citation with DOI
* Minor display fix in vignette

# nmathresh 0.0.4

Version accompanying Phillippo et al. (2017), available from JRSS:A code archive.
