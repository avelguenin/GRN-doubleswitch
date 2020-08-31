import numpy as np


class GeneNetwork:

    def __init__(self,Theta,S,M,sigma):
        
        self.size = 3
        self.s = [1000,10]  ## Production rate vector
        self.d = [0.5,0.1] ## Degradation rate vector
        self.k = [[0.34,2.15],10] ## Basal switch rate vector : [k_on = [k_0,k_1] , k_off]
        
        self.Theta = Theta ## Interaction matrix
        self.S = S ## Effect threshold matrix
        self.M = M ## Interaction type matrix
        self.sigma = sigma ## Basal action potential

        self.t_simu = 0        



    def simulate(self,T_simu,Delta = [0,0,0]):
        
        if not self.t_simu:
            self.g = np.zeros(self.size)
            self.m = np.zeros(self.size)
            self.p = np.zeros(self.size)
            
            self.plots = dict()         
            self.plots["time"] = []
            for i in range(self.size) :
                self.plots["gene " + str(i+1)] = []
                self.plots["mARN " + str(i+1)] = []
                self.plots["protein " + str(i+1)] = []
                self.plots["phi " + str(i+1)] = []
        
        T_simu += self.t_simu
        while self.t_simu < T_simu:

            ## Printout
            print("___\nt = " + str(self.t_simu))
            print("e = ", self.g)
            print("p = ", self.p)

            ## Data writing
            self.plots["time"].append(self.t_simu)
            for i in range(self.size) :
                self.plots["gene " + str(i+1)].append(self.g[i])
                self.plots["mARN " + str(i+1)].append(self.m[i])
                self.plots["protein " + str(i+1)].append(self.p[i])
                

            ## Next step computation
            self.gillespieStep(Delta)
            
        
        
        
    def gillespieStep(self,Delta):
        
        ## Stochastic step forward
        w = self.computeSwitchRates(Delta)
        W = sum(w)
        dt = np.random.exponential(1/W)
        
        self.t_simu += dt
           
        ## Deterministic integration
        for i in range(self.size) :
            self.m[i] /= np.exp(self.d[0]*dt)
            self.m[i] += self.s[0]*self.g[i]*(1-np.exp(-self.d[0]*dt))/self.d[0]
            self.p[i] /= np.exp(self.d[1]*dt)
            self.p[i] += self.s[1]*self.m[i]*(1-np.exp(-self.d[1]*dt))/self.d[1]
        
        ## Stochastic state change 
        w = [i/W for i in w]
        e = np.random.choice(self.size,p=w)
        self.g[e] = 1 - self.g[e]
                

        
    def computeSwitchRates(self,Delta):
        
        phi = [1 for i in range(self.size)]
        for i in range(self.size):
            for j in range(self.size):
                phi[i] *= (1 + np.exp(self.Theta[i,j]) * ((self.p[j]/self.S[i,j]) ** self.M[i,j]))
                phi[i] /= (1 + ((self.p[j]/self.S[i,j]) ** self.M[i,j]))
        phi = np.array(phi) * self.sigma + np.array(Delta)
        print("phi = ", phi)
        for i in range(self.size):
            self.plots["phi " + str(i+1)].append(phi[i])

        
        k_on = [(self.k[0][0] + self.k[0][1] * phi[i])  for i in range(self.size)]         
        k_on = [(k_on[i] / (1 + phi[i] )) for i in range(self.size)]
        print("k_on = ",k_on)
        
        w = []
        for i in range(self.size) :
             if self.g[i] == 1:
                 w.append(self.k[1])
             elif self.g[i] == 0:
                w.append(k_on[i])
                
        return(w)
         



