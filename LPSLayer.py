import numpy as np
import Ekman as ek
import LPSAlgo as lps

class LPSLayer:
    def __init__(self, TRegion, regions, gammas, lamda_e, lamda_minmax, H_0, H_4, phi, we_max, ekmanpumping): 
        self.name = "LPS"
        
        self.ekmanpumping = ekmanpumping
        self.ekman = ek.EkmanActual(0)
                      
        self.TLayer = 3
        self.TRegion = TRegion
        
        if self.TRegion > self.TLayer: raise Exception("Layer needs to be greater or equal to regions")
    
        self.Points = 500
        
        self.regions = regions
        self.gammas = gammas
        self.lamda_minmax = lamda_minmax
        
        self.lamda_e = np.linspace(lamda_e, lamda_e, self.Points)
        self.lamda_e = np.rollaxis(self.lamda_e[np.newaxis],1)
        self.H_0 = H_0
        self.H_4 = H_4
        self.phi = phi
        self.we_max = we_max
        
        self.LPS = lps.LPSAlgo(self.regions, self.gammas, self.H_0, self.H_4)
        
    def Regions(self, Region_all):
        i=-1
        for i in range(self.TRegion-1):
            Region_all = np.where(np.logical_and(self.phi_all>=self.regions[-i-2], self.phi_all<=self.regions[-i-1]), self.TLayer-i, Region_all)
        i += 1
        Region_all = np.where(np.logical_and(self.phi_all>=self.regions[0], self.phi_all<=self.regions[-i-1]), self.TLayer-i, Region_all)
        return Region_all
        
    def EkmanIntegral(self, lamda_num, phi_num):
        if self.ekmanpumping == 'Quadratic':
            ekmanintegral = ek.EkmanQuadraticIntegral(self.lamda_all[phi_num][lamda_num], self.phi_all[phi_num][lamda_num], self.lamda_e[lamda_num][0], self.phi, self.we_max)
        elif self.ekmanpumping == 'Actual':
            ekmanintegral = ek.EkmanActualIntegral(self.lamda_all[phi_num][lamda_num], self.phi_all[phi_num][lamda_num], self.lamda_e[lamda_num][0], self.ekman)
        return ekmanintegral
        
    def Layers(self):
        self.lamda_range = np.linspace(self.lamda_minmax[0], self.lamda_minmax[1], self.Points)
        self.phi_range = np.linspace(self.regions[0], self.regions[-1], self.Points)
        self.lamda_all, self.phi_all = np.meshgrid(self.lamda_range, self.phi_range)
        self.phi_all = np.where(self.lamda_all<=self.lamda_e, self.phi_all, np.nan)
        self.lamda_all = np.where(self.lamda_all<=self.lamda_e, self.lamda_all, np.nan)
        self.Z_all = np.zeros((self.lamda_all.shape[0], self.lamda_all.shape[1], self.TRegion+1))
        self.Z_all[:] = np.nan
        self.shadow_all = np.zeros((self.lamda_all.shape[0], self.lamda_all.shape[1], self.TRegion+1))
        self.shadow_all[:] = np.nan
        self.Region_all = np.zeros((self.Z_all.shape[0], self.Z_all.shape[1]))
        self.Region_all[:] = np.nan
        self.Region_all = self.Regions(self.Region_all)
        
        for lamda_num in range(len(self.lamda_range)):
            for phi_num in range(len(self.phi_range)):
                ek = self.EkmanIntegral(lamda_num, phi_num)
                if not np.isnan(self.Region_all[phi_num][lamda_num]) and not np.isnan(ek):
                    holder = self.LPS.Depth(self.lamda_all[phi_num][lamda_num], self.phi_all[phi_num][lamda_num], TRegion=len(self.regions), region=self.Region_all[phi_num][lamda_num], TLayer=self.TLayer, EkmanIntegral=ek)
                    self.Z_all[phi_num][lamda_num][-holder.shape[0]:] = holder
                    shadow = self.LPS.shadow
                    self.shadow_all[phi_num][lamda_num][-holder.shape[0]:] = shadow
                    
        self.Z_all = np.rollaxis(self.Z_all, 2)
        self.shadow_all = np.rollaxis(self.shadow_all, 2)
        
    def Heights(self):
        self.plot = "Heights"
        fill = np.zeros((1, self.Z_all.shape[1], self.Z_all.shape[2]))
        H = np.concatenate((fill, self.Z_all), axis=0)
        H = np.where(np.isnan(H), 0, H)
        self.heights = H[:-1]-H[1:]
        return self.lamda_all, self.phi_all, self.heights, self.plot, self.name
        
    def Depths(self):
        self.plot = "Depths"
        self.depths = self.Z_all
        return self.lamda_all, self.phi_all, self.depths, self.plot, self.name
        
    def PV(self):
        self.plot = "Potential Vorticity"
        omega = 2*np.pi/(24*60*60)
        f = 2 * omega * np.sin(np.radians(self.phi_all))
        fill = np.zeros((1, self.Z_all.shape[1], self.Z_all.shape[2]))
        H = np.concatenate((fill, self.Z_all), axis=0)
        H = np.where(np.isnan(H), 0, H)
        h = H[:-1]-H[1:]
        self.pv = np.abs(np.where(f==0, np.nan, f)/np.where(h==0, np.nan, h))
        return self.lamda_all, self.phi_all, self.pv, self.plot, self.name
        
    def ShadowZone(self):
        self.plot = "Shadow Zone"
        self.shadowzones = self.shadow_all
        self.shadowzones2 = np.where((self.shadowzones[:,:,:-1]==1)!=(self.shadowzones[:,:,1:]==1), 1, 0)+ np.where((self.shadowzones[:,:,:-1]==2)!=(self.shadowzones[:,:,1:]==2), 2, 0)+np.where((self.shadowzones[:,:,:-1]==3)!=(self.shadowzones[:,:,1:]==3), 3, 0) 
        self.shadowzones2 = np.where(self.shadowzones2==0, np.nan, self.shadowzones2)
        return self.lamda_all, self.phi_all, self.shadowzones, self.plot, self.name