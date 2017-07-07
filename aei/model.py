#####
# functions to run and fit statistical functions and models
#
# c. 2016-2017 Christopher Anderson
#####
import numpy as np

# method for fitting a pre-defined model and returning the fit data and
#  r-squared/rmse
def fit(function, x, y, *args, **kwargs):
    import scipy as sp
    import numpy as np
    from sklearn import metrics
    
    opt, cov = sp.optimize.curve_fit(function, x, y, *args, **kwargs)
    y_fit = function(x, *opt)
    rsq = metrics.r2_score(y, y_fit)
    rms = np.sqrt(metrics.mean_squared_error(y, y_fit))
    return [y_fit, rsq, rms]
    
# do the same, but for an n-dimensional polynomial
def fit_poly(x, y, degree):
    import numpy as np
    from sklearn import metrics
    
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
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
    
# logistic sigmoid model
def sigmoid(x, x0, k, a, c):
    return a / (1 + np.exp(-k * (x - x0))) + c
    
# decay (?) model
def decay(x, a, b, c):
    return c + 1 / (a * x + b)
    
# exponential decay (positive)
def exp_decay_pos(x, a, b, c):
    return a * (1 - np.exp(-b * t)) + c
    
# exponential decay (negative)
def exp_decay_neg(x, a, b, c):
    return a * np.exp(-b * x) + c
    
# function to normalize data
def normalize(x, axis = 0):
    import numpy as np
    
    # get the mean/stdev to normalize with
    mn = np.mean(x, axis = axis)
    sd = np.std(x, axis = axis)
    
    # return the normalized  
    return (x - mn) / sd
    
# function to linearly scale data to new range
def scale(x, xlim = [0,1]):
    import numpy as np
    

# function to test normality of a data set
def test_normal(x, method = 'shapiro-wilk', *args, **kwargs):
    """
    """
    import scipy as sp
    
    methods = ['shapiro-wilk', 'anderson', 'kstest']
    
    if method.lower() not in method:
        print("[ ERROR ]: method must be one of: {}".format(aei.fn.strJoin(methods, ', ')))
        return -1
        
    if method.lower() == 'shapiro-wilk':
        return sp.stats.shapiro(x, *args, **kwargs)
        
    elif method.lower() == 'anderson':
        return sp.stats.anderson(x, *args, **kwargs)
        
    elif method.lower() == 'kstest':
        return sp.stats.kstest(x, *args, **kwargs)
        
# function to calculate sum of squares for each kmeans cluster
def kmeans_ss(data, km_model):
    """
    """
    
    # create output array based on number of clusters in model
    n_clusters = km_model.cluster_centers_.shape[0]
    out_ss = np.zeros(n_clusters)
    
    # predict the clusters for the data presented
    y_pred = km_model.predict(data)
    
    # loop through each cluster and calculate the sum of
    #  squared residuals 
    for i in range(n_clusters):
        y_hat = km_model.cluster_centers_[i,:]
        y_cls = data[y_pred == i]
        out_ss[i] = ((y_cls - y_hat) ** 2).sum()
        
    # finally, sum and return the within-cluster residuals
    return out_ss.sum()
    
# function to calculate multiple kmeans cluster centers and return a vector of
#  the sum of squared residuals for each cluster setup
def kmeans_fit_multi(data, n_clusters, n_jobs = -1, **kwargs):
    """
    """
    from sklearn import cluster
    
    # create an output array based on the number of output clusters
    out_ss = np.zeros(n_clusters)
    
    # loop through each clustering iteration, fit the model, then
    #  record the sum of squares fits for each
    for i in range(n_clusters):
        km_model = cluster.KMeans(n_clusters = i+1, n_jobs = n_jobs,
            **kwargs)
        km_model.fit(data)
        out_ss[i] = kmeans_ss(data, km_model)
        
    # return the output sum of squares array
    return out_ss
    
# function to find the minimum rate of change for multi kmeans clusters
def kmeans_find_min(data, n_clusters, start_cluster = 1, 
        function = exp_decay_neg, **kwargs):
    """
    """
    import scipy as sp
    import numpy as np
    from sklearn import metrics
    from matplotlib import pyplot as plt
        
    # run the multi-fit function
    multi_ss = kmeans_fit_multi(data, n_clusters, **kwargs)
    
    # exclude points prior to the start cluster
    ydata = multi_ss[start_cluster:]
    xdata = np.arange(n_clusters-start_cluster)
    
    # linearly scale the data to simplify curve fitting
    ymax = ydata.max()
    ymin = ydata.min()
    ydata = (ydata - ymin) / (ymax - ymin)
    
    # fit an exponential decay function, and do not include
    #  the values from start_cluster on
    opt, cov = sp.optimize.curve_fit(function, xdata, ydata)
        
    # calculate the fit from the calibration data
    yfit = function(xdata, *opt)
    
    # plot the original and the fit data
    plt.plot(xdata, ydata, color = 'purple', linewidth = 2, label = 'orig. data')
    plt.plot(xdata, yfit, color = 'orange', linewidth = 2, label = 'fit data')
    plt.xlabel("Number of clusters starting from cluster {:02d}".format(
        start_cluster+1))
    plt.ylabel("Scaled sum of squares within each cluster group")
    plt.legend()
    
    # calculate the second derivative from the fit function
    n_derivs = n_clusters-start_cluster-2
    deriv = np.zeros(n_derivs)
    for i in range(start_cluster, n_derivs + start_cluster):
        deriv[i - start_cluster] = sp.misc.derivative(function, 
            xdata[i], args = opt, n = 2)
            
    # return the plot for the user to manipulate
    return plt
        