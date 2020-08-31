import numpy as np


class GeneNetwork:

    def __init__(self,k,s,d,Theta):
        self.K = k ## switch rate matrix
        self.S = s ## production rate matrix
        self.D = d ## degradation rate matrix
        self.Theta = Theta ## interaction matrix
        self.size = len(Theta)
        self.t_simu = 0

    def simulate(self,T_simu,Delta):
        
        if not self.t_simu:
            self.g = np.zeros(self.size)
            self.m = np.zeros(self.size)
            self.p = np.zeros(self.size)
            
            self.plots = dict()         
            self.plots["time"] = [self.t_simu]
            for i in range(self.size) :
                self.plots["gene " + str(i+1)] = [self.g[i]]
                self.plots["mARN " + str(i+1)] = [self.m[i]]
                self.plots["protein " + str(i+1)] = [self.p[i]]
        
        T_simu += self.t_simu
        while self.t_simu < T_simu:
            print("___\nt = " + str(self.t_simu))
            self.gillespieStep(Delta)

            ## Data writing
            self.plots["time"].append(self.t_simu)
            for i in range(self.size) :
                self.plots["gene " + str(i+1)].append(self.g[i])
                self.plots["mARN " + str(i+1)].append(self.m[i])
                self.plots["protein " + str(i+1)].append(self.p[i])
        
        return()
        
        
    def gillespieStep(self,Delta):
        
        ## Stochastic step forward
        w = self.computeSwitchRates(Delta)
        W = sum(w)
        dt = np.random.exponential(1/W)
        
        if dt > 1/self.K[0]: ## a max time step must be introduced to ensure stability
            dt = 1/self.K[0]
            ## Deterministic integration
            for i in range(self.size) :            
                self.m[i] += self.S[0]*self.g[i]*dt
                self.m[i] /= np.exp(self.D[0]*dt)
                self.p[i] += self.S[1]*self.m[i]*dt
                self.p[i] /= np.exp(self.D[1]*dt)
            
        else:         
            ## Deterministic integration
            for i in range(self.size) :            
                self.m[i] += self.S[0]*self.g[i]*dt
                self.m[i] /= np.exp(self.D[0]*dt)
                self.p[i] += self.S[1]*self.m[i]*dt
                self.p[i] /= np.exp(self.D[1]*dt)
       
            ## Stochastic state change 
            w = [i/W for i in w]
            e = np.random.choice(self.size,p=w)
            self.g[e] = 1 - self.g[e]
        
        self.t_simu += dt
        
        

        
    def computeSwitchRates(self,Delta):
        w = []
        p_act = self.Theta.dot(np.array([x/(15000+x) for x in self.p]))
#        p_act = self.Theta.dot(self.p)/3
#        p_act = self.Theta.dot(self.p)/100000
        print("p_act = " + str(p_act))
        for i in range(self.size) :
            if self.g[i] == 1:
                w.append(self.K[1])
            elif self.g[i] == 0:
#                w.append(self.K[0] / (1 + 3*np.exp(- p_act[i] - Delta[i] - 5)))
                w.append(self.K[0] + p_act[i] + Delta[i])
        
        print("w = " + str(w))
        w = [max(x,self.K[0]/3) for x in w]
        return(w)
    



