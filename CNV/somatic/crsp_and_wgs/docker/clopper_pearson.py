import numpy as np
import scipy
import scipy.stats


def clopper_pearson(k, n, alpha=0.05):
    """
    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    alpha confidence intervals for a binomial distribution of k expected successes on n trials
    Clopper Pearson intervals are a conservative estimate.
    """
    lo = scipy.stats.beta.ppf(alpha / 2, k, n - k + 1)
    hi = scipy.stats.beta.ppf(1 - alpha / 2, k + 1, n - k)
    if (n - k) == 0:
        return lo, 1.0
    if np.isnan(lo):
        lo = 0.0
    return lo, hi
