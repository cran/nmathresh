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
