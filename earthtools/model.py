#####
# functions to run and fit statistical functions and models
#
# c. 2016-2017 Christopher Anderson
#####
import numpy as np


# method for fitting a pre-defined model and returning the fit data and
#  r-squared/rmse
def fit(func, x, y, *args, **kwargs):
    import scipy as sp
    import numpy as np
    from sklearn import metrics

    opt, cov = sp.optimize.curve_fit(func, x, y, *args, **kwargs)
    y_fit = func(x, *opt)
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
def exp_decay_pos(x, a, b, c, t):
    return a * (1 - np.exp(-b * t)) + c


# exponential decay (negative)
def exp_decay_neg(x, a, b, c):
    return a * np.exp(-b * x) + c


# function to normalize data
def normalize(x, axis=0):
    import numpy as np

    # get the mean/stdev to normalize with
    mn = np.mean(x, axis=axis)
    sd = np.std(x, axis=axis)

    # return the normalized
    return (x - mn) / sd


# function to linearly scale data to new range
def scale(x, xlim=[0, 1], axis=0):
    import numpy as np


# function to test normality of a data set
def test_normal(x, method='shapiro-wilk', *args, **kwargs):
    """
    """
    import scipy as sp
    import earthtools as et

    methods = ['shapiro-wilk', 'anderson', 'kstest']

    if method.lower() not in method:
        print("[ ERROR ]: method must be one of: {}".format(et.fn.strJoin(methods, ', ')))
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
        y_hat = km_model.cluster_centers_[i, :]
        y_cls = data[y_pred == i]
        out_ss[i] = ((y_cls - y_hat) ** 2).sum()

    # finally, sum and return the within-cluster residuals
    return out_ss.sum()


# function to calculate multiple kmeans cluster centers and return a vector of
#  the sum of squared residuals for each cluster setup
def kmeans_fit_multi(data, n_clusters, n_jobs=-1, **kwargs):
    """
    """
    from sklearn import cluster

    # create an output array based on the number of output clusters
    out_ss = np.zeros(n_clusters)

    # loop through each clustering iteration, fit the model, then
    #  record the sum of squares fits for each
    for i in range(n_clusters):
        km_model = cluster.KMeans(n_clusters=i + 1, n_jobs=n_jobs,
                                  **kwargs)
        km_model.fit(data)
        out_ss[i] = kmeans_ss(data, km_model)

    # return the output sum of squares array
    return out_ss


# function to find the minimum rate of change for multi kmeans clusters
def kmeans_find_min(data, n_clusters, start_cluster=1,
                    function=exp_decay_neg, **kwargs):
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
    xdata = np.arange(n_clusters - start_cluster + 1)

    # linearly scale the data to simplify curve fitting
    ymax = ydata.max()
    ymin = ydata.min()
    ydata = (ydata - ymin) / (ymax - ymin)

    # fit an exponential decay function, and do not include
    #  the values from start_cluster on
    opt, cov = sp.optimize.curve_fit(function, xdata, ydata)

    # calculate the fit from the calibration data
    yfit = function(xdata, *opt)

    # re-scale the xdata to include the starting cluster
    xplot = xdata + start_cluster + 1

    # plot the original and the fit data
    plt.figure(np.random.randint(256))
    plt.plot(xplot, ydata, color='purple', linewidth=2, label='orig. data')
    plt.plot(xplot, yfit, color='orange', linewidth=2, label='fit data')
    plt.xlabel("Number of clusters")
    plt.ylabel("Scaled sum of squares within each cluster group")
    plt.tight_layout()
    plt.legend()

    # calculate the second derivative from the fit function
    n_derivs = n_clusters - start_cluster - 2
    deriv = np.zeros(n_derivs)
    for i in range(start_cluster, n_derivs + start_cluster):
        deriv[i - start_cluster] = sp.misc.derivative(function,
                                                      xdata[i], args=opt, n=2)

    # return the plot for the user to manipulate
    return plt


# functions to list the available options for various sklearn parameters
def list_scoring():
    options = {
        'Classification': [
            'accuracy',
            'average_precision',
            'f1',
            'f1_micro',
            'f1_macro',
            'f1_weighted',
            'f1_samples',
            'neg_log_loss',
            'precision',
            'recall',
            'roc_auc'
        ],
        'Clustering': ['adjusted_rand_score'],
        'Regression': [
            'neg_mean_absolute_error',
            'neg_mean_squared_error',
            'neg_median_absolute_error',
            'r2'
        ]
    }


# a class for model tuning that has the default parameters pre-set
class tune:
    def __init__(self, x, y, optimizer=None, param_grid=None,
                 scoring=None, fit_params=None, n_jobs=None, refit=True,
                 verbose=1, error_score=0, return_train_score=True,
                 cv=None, n_splits=5):
        from earthtools import params
        from sklearn import model_selection

        # set defaults for each of the potential parameters
        if optimizer is None:
            optimizer = model_selection.GridSearchCV
        if n_jobs is None:
            n_jobs = params.cores - 1

        # set the variables to the tune class
        self.x = x
        self.y = y
        self.optimizer = optimizer
        self.param_grid = param_grid
        self.scoring = scoring
        self.fit_params = fit_params
        self.n_jobs = n_jobs
        self.refit = refit
        self.verbose = verbose
        self.error_score = error_score
        self.return_train_score = return_train_score
        self.cv = cv
        self.n_splits = n_splits

    # function to actually run the tuning process and report outputs
    def run_gs(self, estimator):

        # create the grid search
        gs = self.optimizer(estimator, param_grid=self.param_grid,
                            scoring=self.scoring, fit_params=self.fit_params,
                            n_jobs=self.n_jobs, cv=self.cv, refit=self.refit,
                            verbose=self.verbose, error_score=self.error_score,
                            return_train_score=self.return_train_score)

        # begin fitting the grid search
        gs.fit(self.x, self.y)

        # update the tuning object with the outputs of the tuning
        self.cv_results = gs.cv_results_
        self.best_estimator = gs.best_estimator_
        self.best_score = gs.best_score_
        self.best_params = gs.best_params_
        self.best_index = gs.best_index_
        self.scorer = gs.scorer_
        self.n_splits = gs.n_splits_
        self.gs = gs

    #####
    # create functions to tune each different model type based on
    #  a set of default parameter grids

    #####
    # LINEAR MODELS

    # OLS linear regression
    def LinearRegression(self, optimizer=None, param_grid=None,
                         scoring=None, fit_params=None, cv=None):
        from sklearn import linear_model, model_selection, metrics

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'normalize': (True, False),
                    'fit_intercept': (True, False)
                }
                self.param_grid = param_grid
        else:
            self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = metrics.explained_variance_score
                self.scoring = scoring
        else:
            self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
                self.cv = cv
        else:
            self.cv = cv

        # create the estimator and run the grid search
        estimator = linear_model.LinearRegression()
        self.run_gs(estimator)

    # the logistic (i.e. MaxEnt) model
    def LogisticRegression(self, optimizer=None, param_grid=None,
                           scoring=None, fit_params=None, cv=None):
        from sklearn import linear_model, model_selection

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'C': (1e-2, 1e-1, 1e0, 1e1),
                    'tol': (1e-3, 1e-4, 1e-5),
                    'fit_intercept': (True, False)
                }
                self.param_grid = param_grid
        else:
            self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = 'roc_auc'
                self.scoring = scoring
        else:
            self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
                self.cv = cv
        else:
            self.cv = cv

        # create the estimator and run the grid search
        estimator = linear_model.LogisticRegression()
        self.run_gs(estimator)

    #####
    # DECISION TREE FUNCTIONS

    # function for decision tree classifier
    def DecisionTreeClassifier(self, optimizer=None, param_grid=None,
                               scoring=None, fit_params=None, cv=None):
        from sklearn import tree, model_selection

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'criterion': ('gini', 'entropy'),
                    'splitter': ('best', 'random'),
                    'max_features': ('sqrt', 'log2', None),
                    'max_depth': (2, 5, 10, None),
                    'min_samples_split': (2, 0.01, 0.1),
                    'min_impurity_split': (1e-7, 1e-6)
                }
                self.param_grid = param_grid
        else:
            self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = 'roc_auc'
                self.scoring = scoring
        else:
            self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
                self.cv = cv
        else:
            self.cv = cv

        # create the estimator and run the grid search
        estimator = tree.DecisionTreeClassifier()
        self.run_gs(estimator)

    #####
    # SVM FUNCTIONS

    # function for Support Vector Classifier (SVC)
    def SVC(self, optimizer=None, param_grid=None,
            scoring=None, fit_params=None, cv=None,
            class_weight=None):
        from sklearn import svm, model_selection

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'C': (1e-3, 1e-2, 1e-1, 1e0, 1e1),
                    'kernel': ('rbf', 'linear'),
                    'gamma': (1e-3, 1e-4, 1e-5, 1e-6, 1e-7)
                }
        self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = 'roc_auc'
        self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
        self.cv = cv

        # set the class weight
        if class_weght is not None:
            if self.class_weight is None:
                class_weight = 'balanced'
        self.param_grid['class_weight'] = class_weight

        # create the estimator and run the grid search
        estimator = svm.SVC()
        self.run_gs(estimator)

    # function for Support Vector Regression (SVR)
    def SVR(self, optimizer=None, param_grid=None,
            scoring=None, fit_params=None, cv=None):
        from sklearn import svm, model_selection, metrics

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'C': (1e-2, 1e-1, 1e0, 1e1),
                    'epsilon': (0.01, 0.1, 1),
                    'kernel': ('rbf', 'linear', 'poly', 'sigmoid'),
                    'gamma': (1e-2, 1e-3, 1e-4)
                }
                self.param_grid = param_grid
        else:
            self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = metrics.explained_variance_score
                self.scoring = scoring
        else:
            self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
                self.cv = cv
        else:
            self.cv = cv

        # create the estimator and run the grid search
        estimator = svm.SVR()
        self.run_gs(estimator)

    # function for Linear SVC
    def LinearSVC(self, optimizer=None, param_grid=None,
                  scoring=None, fit_params=None, cv=None):
        from sklearn import svm, model_selection

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'C': (1e-2, 1e-1, 1e0, 1e1),
                    'loss': ('hinge', 'squared_hinge'),
                    'tol': (1e-3, 1e-4, 1e-5),
                    'fit_intercept': (True, False)
                }
                self.param_grid = param_grid
        else:
            self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = 'roc_auc'
                self.scoring = scoring
        else:
            self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
                self.cv = cv
        else:
            self.cv = cv

        # create the estimator and run the grid search
        estimator = svm.LinearSVC()
        self.run_gs(estimator)

    # function for Linear SVR
    def LinearSVR(self, optimizer=None, param_grid=None,
                  scoring=None, fit_params=None, cv=None):
        from sklearn import svm, model_selection, metrics

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'C': (1e-2, 1e-1, 1e0, 1e1),
                    'loss': ('epsilon_insensitive', 'squared_epsilon_insensitive'),
                    'epsilon': (0, 0.01, 0.1),
                    'dual': (False),
                    'tol': (1e-3, 1e-4, 1e-5),
                    'fit_intercept': (True, False)
                }
                self.param_grid = param_grid
        else:
            self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = metrics.explained_variance_score
                self.scoring = scoring
        else:
            self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
                self.cv = cv
        else:
            self.cv = cv

        # create the estimator and run the grid search
        estimator = svm.LinearSVR()
        self.run_gs(estimator)

    #####
    # ENSEMBLE METHODS

    # Ada boosting classifier
    def AdaBoostClassifier(self, optimizer=None, param_grid=None,
                           scoring=None, fit_params=None, cv=None):
        from sklearn import ensemble, model_selection

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'n_estimators': (25, 50, 75, 100),
                    'learning_rate': (0.1, 0.5, 1.)
                }
                self.param_grid = param_grid
        else:
            self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = 'roc_auc'
                self.scoring = scoring
        else:
            self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
                self.cv = cv
        else:
            self.cv = cv

        # create the estimator and run the grid search
        estimator = ensemble.AdaBoostClassifier()
        self.run_gs(estimator)

    # Ada boosting regressor
    def AdaBoostRegressor(self, optimizer=None, param_grid=None,
                          scoring=None, fit_params=None, cv=None):
        from sklearn import ensemble, model_selection, metrics

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'n_estimators': (25, 50, 75, 100),
                    'learning_rate': (0.1, 0.5, 1.),
                    'loss': ('linear', 'exponential', 'square')
                }
                self.param_grid = param_grid
        else:
            self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = metrics.explained_variance_score
                self.scoring = scoring
        else:
            self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
                self.cv = cv
        else:
            self.cv = cv

        # create the estimator and run the grid search
        estimator = ensemble.AdaBoostRegressor()
        self.run_gs(estimator)

    # Gradient boosting classifier
    def GradientBoostClassifier(self, optimizer=None, param_grid=None,
                                scoring=None, fit_params=None, cv=None):
        from sklearn import ensemble, model_selection

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'n_estimators': (10, 100, 500),
                    'learning_rate': (0.01, 0.1, 0.5),
                    'max_features': ('sqrt', 'log2', None),
                    'max_depth': (1, 10, None),
                    'min_samples_split': (2, 0.01, 0.1),
                    'min_impurity_split': (1e-7, 1e-6)
                }
                self.param_grid = param_grid
        else:
            self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = 'roc_auc'
                self.scoring = scoring
        else:
            self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
                self.cv = cv
        else:
            self.cv = cv

        # create the estimator and run the grid search
        estimator = ensemble.GradientBoostingClassifier()
        self.run_gs(estimator)

    # Gradient boosting regressor
    def GradientBoostRegressor(self, optimizer=None, param_grid=None,
                               scoring=None, fit_params=None, cv=None):
        from sklearn import ensemble, model_selection

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'n_estimators': (10, 100, 500),
                    'learning_rate': (0.01, 0.1, 0.5),
                    'max_features': ('sqrt', 'log2', None),
                    'max_depth': (1, 10, None),
                    'min_samples_split': (2, 0.01, 0.1),
                    'min_impurity_split': (1e-7, 1e-6)
                }
                self.param_grid = param_grid
        else:
            self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = 'neg_mean_absolute_error'
                self.scoring = scoring
        else:
            self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
                self.cv = cv
        else:
            self.cv = cv

        # create the estimator and run the grid search
        estimator = ensemble.GradientBoostingRegressor()
        self.run_gs(estimator)

    # Random Forest classifier
    def RandomForestClassifier(self, optimizer=None, param_grid=None,
                               scoring=None, fit_params=None, cv=None,
                               class_weight=None):
        from sklearn import ensemble, model_selection

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'criterion': ('gini', 'entropy'),
                    'n_estimators': (10, 100, 500),
                    'max_features': ('sqrt', 'log2', None),
                    'max_depth': (1, 10, None),
                    'min_samples_split': (2, 0.01, 0.1),
                    'min_impurity_split': (1e-7, 1e-6)
                }
                self.param_grid = param_grid
        else:
            self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = 'roc_auc'
                self.scoring = scoring
        else:
            self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
                self.cv = cv
        else:
            self.cv = cv

        # create the estimator and run the grid search
        estimator = ensemble.RandomForestClassifier()
        self.run_gs(estimator)

    # Random Forest regressor
    def RandomForestRegressor(self, optimizer=None, param_grid=None,
                              scoring=None, fit_params=None, cv=None):
        from sklearn import ensemble, model_selection

        # check if the optimizer has changed, otherwise use default
        if optimizer is not None:
            self.optimizer = optimizer

        # check if the parameter grid has been set, otherwise set defaults
        if param_grid is None:
            if self.param_grid is None:
                param_grid = {
                    'n_estimators': (10, 100, 500),
                    'max_features': ('sqrt', 'log2', None),
                    'max_depth': (1, 10, None),
                    'min_samples_split': (2, 0.01, 0.1),
                    'min_impurity_split': (1e-7, 1e-6)
                }
                self.param_grid = param_grid
        else:
            self.param_grid = param_grid

        # set the scoring function
        if scoring is None:
            if self.scoring is None:
                scoring = 'neg_mean_absolute_error'
                self.scoring = scoring
        else:
            self.scoring = scoring

        # set the default fit parameters
        if fit_params is not None:
            self.fit_params = fit_params

        # set the cross validation strategy
        if cv is None:
            if self.cv is None:
                cv = model_selection.StratifiedKFold(n_splits=self.n_splits)
                self.cv = cv
        else:
            self.cv = cv

        # create the estimator and run the grid search
        estimator = ensemble.RandomForestRegressor()
        self.run_gs(estimator)
