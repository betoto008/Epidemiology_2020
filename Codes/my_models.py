from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import networkx as networkx
import numpy as np
import scipy as scipy
import scipy.integrate
import matplotlib.pyplot as plt

class EpiModel():
    """docstring for EpiModel"""

    def __init__(self, N, I0, beta, gamma):

        #super(EpiModel, self).__init__()
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Model Parameters:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.N   = N
        self.I0   = I0
        self.beta = beta
        self.gamma  = gamma
        self.R0 = beta/gamma

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Initialize Timekeeping:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.t       = 0
        self.tmax    = 0 # will be set when run() is called
        self.tseries = np.array([0])

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Initialize Counts of inidividuals with each state:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.I = np.array([int(I0)])
        self.S = np.array([int(N)-self.I[-1]])
        self.R = np.array([0])

        assert(self.I[0] > 0), "The specified initial infected population size of n0 must be greater than zero."

    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    @staticmethod
    def ODE_system(t, variables, N, I0, beta, gamma ):

        S, I, R = variables

        dS = -beta*I*S/N

        dI = beta*I*S/N - gamma*I

        dR = gamma*I

        

        return [dS, dI, dR]

    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    def run(self, runtime, dt = 0.1):


        if(runtime>0):
            self.tmax = runtime
        else:
            return False
        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        t_eval    = np.arange(start=self.t, stop=self.t+runtime, step=dt)

        # Define the range of time values for the integration:
        t_span          = (self.t, self.t+runtime)

        # Define the initial conditions as the system's current state:
        # (which will be the t=0 condition if this is the first run of this model, 
        # else where the last sim left off)

        init_cond       = [self.S[-1], self.I[-1], self.R[-1]]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Solve the system of differential eqns:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        solution        = scipy.integrate.solve_ivp(lambda t, X: EpiModel.ODE_system(t, X, self.N, self.I0, self.beta, self.gamma), t_span=[self.t, self.tmax], y0=init_cond, t_eval=t_eval)

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Store the solution output as the model's time series and data series:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.tseries    = np.append(self.tseries, solution['t'])
        self.S       = np.append(self.S, solution['y'][0])
        self.I       = np.append(self.I, solution['y'][1])
        self.R       = np.append(self.R, solution['y'][2])

        self.t = self.tseries[-1]

        return self.tseries, self.S, self.I, self.R


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%