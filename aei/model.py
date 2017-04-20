#####
# functions to run and fit statistical functions and models
#
# c. 2016-2017 Christopher Anderson
#####

# method for fitting a pre-defined model and returning the fit data and
#  r-squared/rmse
def fit(function, x, y, *args, **kwargs):
    import scipy
    import numpy as np
    
    opt, cov = scipy.optimize.curve_fit(function, x, y, *args, **kwargs)
    y_fit = function(x, *opt)
    rsq = metrics.r2_score(y, y_fit)
    rms = np.sqrt(metrics.mean_squared_error(y, y_fit))
    return [y_fit, rsq, rms]
    
# do the same, but for an n-dimensional polynomial
def fit_poly(x, y, degree):
    import numpy as np
    coeffs = np.polyfit(x, y, degree)
    y_fit = np.polyval(coeffs, x)
    rsq = metrics.r2_score(y, y_fit)
    rms = np.sqrt(metrics.mean_squared_error(y, y_fit))
    return [y_fit, rsq, rms]

# linear model
def linear(x, a, b):
    return a * x + b
    
# gaussian model
def gaus(x, a, x0, sigma):
    return a * exp(-(x - x0) ** 2 / (2 * sigma ** 2))
    
# logistic sigmoid model
def sigmoid(x, x0, k, a, c):
    return a / (1 + np.exp(-k * (x - x0))) + c
    
