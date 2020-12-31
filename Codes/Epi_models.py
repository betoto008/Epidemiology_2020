from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import networkx as nx
import numpy as np
import scipy as scipy
import scipy.integrate
import matplotlib.pyplot as plt
import math
import sys
import random
import seaborn

from scipy.interpolate import interp1d
from models import *

class Ensemble_EpiModel():
    """docstring for Ensemble_EpiModel"""

    def __init__(self, G, beta, sigma, gamma, xi=0, mu_I=0, mu_0=0, nu=0, beta_local=None, p=0,
                    Q=None, beta_D=None, sigma_D=None, gamma_D=None, mu_D=None, beta_D_local=None,
                    theta_E=0, theta_I=0, phi_E=0, phi_I=0, psi_E=1, psi_I=1, q=0,
                    initE=0, initI=10, initD_E=0, initD_I=0, initR=0, initF=0,
                    node_groups=None, store_Xseries=False):

    	self.T_regular = np.linspace(0, Total_t, intervals)
    	self.gamma = gamma



    def My_print(self, beta):
    	print('This is self.beta:',self.beta,' This is self.gamma:', self.gamma, ' and this is beta', beta)
