import numpy as np

class LPSAlgo:
    def __init__(self, regions, gammas, H_0, H_4):
        self.R = 6370e3 #Radius of the earth
        self.omega = 2*np.pi/(24*60*60) #Rotational velocity of Earth
        self.H_0 = H_0 #Height of h_3 on eastern boundry
        self.H_4 = H_4
        
        self.fval = np.zeros([len(regions)]).astype(float)
        self.gamma = np.zeros([len(gammas)]).astype(float)
        self.gamma = gammas
        self.fval = regions
        self.fval = 2*self.omega*np.sin(np.radians(self.fval))
        
        
    def Depth(self, lamda_deg, phi_deg, TRegion, region, TLayer, EkmanIntegral):
        
        self.shadow = 0
              
        if region == 3: x=self.Region3(lamda_deg, phi_deg, TRegion, region, TLayer, EkmanIntegral)
        if region == 2: x=self.Region2(lamda_deg, phi_deg, TRegion, region, TLayer, EkmanIntegral)
        if region == 1: x=self.Region1(lamda_deg, phi_deg, TRegion, region, TLayer, EkmanIntegral)
            
        return x
        
    def Region3(self, lamda_deg, phi_deg,  TRegion, region, TLayer, EkmanIntegral):
        phi = np.radians(phi_deg)
        
        f = 2 * self.omega * np.sin(phi)
        beta = 2 * self.omega / self.R * np.cos(phi)

        self.D_0_squared = - 2 * f**2/(beta * self.gamma[TRegion-1]) * EkmanIntegral
        
        h3 = 1
        F = 1
        H = np.sqrt(self.D_0_squared + self.H_0**2)/ np.sqrt(F)
        
        Depths = [h3*H, self.H_4]

        self.Depths = -np.array(Depths)
        
        return self.Depths
        
    def Region2(self, lamda_deg, phi_deg,  TRegion, region, TLayer, EkmanIntegral):
        phi = np.radians(phi_deg)
        
        f = 2 * self.omega * np.sin(phi)
        beta = 2 * self.omega / self.R * np.cos(phi)

        self.D_0_squared = - 2 * f**2/(beta * self.gamma[TRegion-1]) * EkmanIntegral
        
        h3 = f/self.fval[TRegion-2]
        h2 = 1 - h3
        F = 1 + self.gamma[TRegion-2]/self.gamma[TRegion-1]*(1 - f/self.fval[TRegion-2])**2
        H = np.sqrt(self.D_0_squared + self.H_0**2)/ np.sqrt(F)
        Depths = [h2*H, (h2+h3)*H, self.H_4]
        
        Area_R = self.D_0_squared <= (F - 1) * self.H_0**2
        H2 = (self.gamma[TRegion-1]/self.gamma[TRegion-2]*self.D_0_squared)**0.5
        Depths[0] = np.where(Area_R, H2, Depths[0])
        Depths[1] = np.where(Area_R, self.H_0, Depths[1])
        self.shadow = np.where(Area_R, 1, self.shadow)
        
        self.Depths = -np.array(Depths)
        
        return self.Depths
        
    def Region1(self, lamda_deg, phi_deg, TRegion, region, TLayer, EkmanIntegral):
        phi = np.radians(phi_deg)
        
        f = 2 * self.omega * np.sin(phi)
        beta = 2 * self.omega / self.R * np.cos(phi)

        self.D_0_squared = - 2 * f**2/(beta * self.gamma[TRegion-1]) * EkmanIntegral
        
        h3 = (f/self.fval[TRegion-2])
        h2 = f/self.fval[TRegion-3]*(1-self.fval[TRegion-3]/self.fval[TRegion-2])*((1+self.gamma[TRegion-2]/self.gamma[TRegion-1]*(1-f/self.fval[TRegion-2]))/(1+self.gamma[TRegion-2]/self.gamma[TRegion-1]*(1-self.fval[TRegion-3]/self.fval[TRegion-2])))
        h1 = 1 - (h2+h3)
        F = (1 + self.gamma[TRegion-2]/self.gamma[TRegion-1]*(h1+h2)**2 + self.gamma[TRegion-3]/self.gamma[TRegion-1]*h1**2)
        H = np.sqrt(self.D_0_squared + self.H_0**2)/ np.sqrt(F)
        
        Depths = [h1*H, (h1+h2)*H, (h1+h2+h3)*H, self.H_4]
        
        Area_M = self.D_0_squared <= (F - 1) * self.H_0**2
        a = (1-f/self.fval[TRegion-3]*self.gamma[TRegion-2]/self.gamma[TRegion-1]*((1-self.fval[TRegion-3]/self.fval[TRegion-2])/(1+self.gamma[TRegion-2]/self.gamma[TRegion-1]*(1-self.fval[TRegion-3]/self.fval[TRegion-2]))))
        b = self.H_0*(f/self.fval[TRegion-3])*(1-self.fval[TRegion-3]/self.fval[TRegion-2])/(1 + self.gamma[TRegion-2]/self.gamma[TRegion-1]*(1-self.fval[TRegion-3]/self.fval[TRegion-2]))
        H2 = (self.gamma[TRegion-3]/self.gamma[TRegion-2]*a*b + (self.gamma[TRegion-1]/self.gamma[TRegion-2]*self.D_0_squared*(1+self.gamma[TRegion-3]/self.gamma[TRegion-2]*a**2)-self.gamma[TRegion-3]/self.gamma[TRegion-2]*b**2)**0.5)*(1+self.gamma[TRegion-3]/self.gamma[TRegion-2]*a**2)**-1
        h2 = f/self.fval[TRegion-3]*(self.H_0+self.gamma[TRegion-2]/self.gamma[TRegion-1]*H2)*(1-self.fval[TRegion-3]/self.fval[TRegion-2])/(1+self.gamma[TRegion-2]/self.gamma[TRegion-1]*(1-self.fval[TRegion-3]/self.fval[TRegion-2]))
        h1 = H2-h2
        Depths[0] = np.where(Area_M, h1, Depths[0])
        Depths[1] = np.where(Area_M, H2, Depths[1])
        Depths[2] = np.where(Area_M, self.H_0, Depths[2])
        self.shadow = np.where(Area_M, 2, self.shadow)
        
        Area_R = self.D_0_squared <= self.H_0**2 * self.gamma[TRegion-2]/self.gamma[TRegion-1] * (1-self.fval[TRegion-3]/self.fval[TRegion-2])**2 * (1 + self.gamma[TRegion-3]/ self.gamma[TRegion-2] * (1 - f/self.fval[TRegion-3])**2)
        H2 = (self.gamma[TRegion-1]/self.gamma[TRegion-2]*1/(1+self.gamma[TRegion-3]/self.gamma[TRegion-2]*(1-f/self.fval[TRegion-3])**2)*self.D_0_squared)**0.5
        h2 = f/self.fval[TRegion-3]*H2
        h1 = H2-h2
        Depths[0] = np.where(Area_R, h1, Depths[0])
        Depths[1] = np.where(Area_R, H2, Depths[1])
        Depths[2] = np.where(Area_R, self.H_0, Depths[2])
        self.shadow = np.where(Area_R, 3, self.shadow)
        
        self.Depths = -np.array(Depths)
        
        return self.Depths