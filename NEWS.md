### Package changes from previous ecolRxC version 0.1.1-10

Compared to version 0.1.1-10 of ecolRxC, the new version:

* fixes a bug in the estimation of confidence transition probabilities when `method = "Thomsen"`. As a result, more accurate estimates are obtained on average.
* fixes a bug in the estimation of confidence intervals for the transfer probabilities when `method = "IPF"`.
* includes a new argument (`regions`) to split the set of units into subsets, applying the selected ecolRxC specification separately to each subset.  
* includes a new argument (`ref.combination`) to specify how unit table estimates should be combined to generate the default output  when `method = "Thomsen"` and `reference = NULL`.
* includes, as new outputs, the bounds implied by the observed data for the transfer probabilities.
* offers a new method to estimate confidence intervals for the transfer probabilities adapted to this setting, based on Pavía-Miralles (2005) <https://www.jstor.org/stable/27590658>.
* computes bivariate Gaussian probabilities using `mvtnorm::pmvnorm()` to acelerate the process.
* translates some internal functions to `C++` to accelerate computations. 