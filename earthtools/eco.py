#####
# functions to run ecological models
#
# c. 2016-2017 Christopher Anderson
#####

# a series of functions to calculate diversity metrics
class Diversity:

    # initialize the class
    def __init__(self, data, weights=None):
        """Calculates common diversity metrics

        Args:
            data: the data set to initialize. should be a 2-d np array
                  where nrows = nsites to asses, and ncol = nspecies in community.
                  the values in the array should be the observations
            weights: weights for each site. default is sum of

        Returns:
            an object with properties obj.shannon() and obj.simpson() to calculate
            those diversity metrics.
        """
        import numpy as np

        # set the data in this object to carry forward
        self.data = data

        # calculate relative abundance
        self.rel_a = data / np.expand_dims(data.sum(axis=1), 1)

        # calculate relative weights for all sites
        if weights is None:
            weights = data.sum(axis=1) / data.sum()
        self.weights = weights / weights.sum()
        self.total_list = data.sum(axis=0)
        self.total_rel = self.total_list / self.total_list.sum()

    # set function to calculate simpsons diversity
    def calc_simpson(self):
        """Calculates Simpson's D

        Eqn:
            D = sum(rel_a^2)
            where rel_a = species relative abundance in a site.

        Args:
            None.

        Returns:
            object updated with obj.alpha_simpson, obj.beta_simpson, obj.gamma_simpson
        """
        import numpy as np

        # calculate raw alpha
        self.alpha_simpson = np.sum(self.rel_a ** 2, axis=1)

    # set function to calculate gini-simpson
    def calc_gini_simpson(self):
        """Calculates gini-simpson diversity coeff

        Eqn:
            D = 1 - sum(rel_a^2)

        Args:
            None.

        Returns:
            object updated with obj.alpha_gini, obj.beta_gini, obj.gamma_gini
        """

        # calculate raw alpha
        self.alpha_gini = 1 - np.sum(rel_a ** 2, axis=1)

    # set function to calculate shannon diversity
    def calc_shannon(self):
        """Calculates Shannon index

        Eqn:
            D = -sum(rel_a * log(rel_a))

        Args:
            None.

        Returns:
            object updated with obj.alpha_shannon, obj.beta_shannon, obj.gamma_shannon
        """

        # calculate raw alpha
        self.alpha_shannon = -sum(self.rel_a * log(rel_a), axis=1)

    # set function to calculate dominance index
    def calc_dominance(self):
        """Calculates dominance index

        Eqn:
            D = max(rel_a)

        Args:
            None.

        Returns:
            object updated with obj.dominance
        """

        # calculate dominance
        self.dominance = self.rel_a.max(axis=1)

    def calc_richness(self):
        """Calculates species richness

        Eqn:
            D = sum(rel_a)

        Args:
            None

        Returns:
            object updated with obj.richness
        """

        # calculate richness
        self.richness = np.sum(self.rel_a, axis=1)

    def calc_evenness(self):
        """Calculates species evenness

        Eqn:
            shannons_h / log(sp_richness)

        Args:
            None

        Returns:
            object updated with obj.evenness
        """

        if not hasattr(self, "alpha_shannon"):
            self.calc_shannon()

        if not hasattr(self, "richness"):
            self.calc_richness()

        self.evenness = self.alpha_shannon / np.log(self.richness)
