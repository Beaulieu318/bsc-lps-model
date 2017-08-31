import numpy as np
import ArgoData as ad

class ArgoLayer:
    def __init__(self):
        self.name = "Argo"
        
        data = ad.ArgoData(Type='Density', Month = 'Average')
        self.lamda_range, self.phi_range, self.press, self.press_units, self.t, self.t_units, self.t_title = data.AllData()
        
        self.t_lats = np.rollaxis(self.t, 1)
        self.t_lons = np.rollaxis(self.t, 2)
        self.press_grid = np.zeros(self.t.shape)    
        for i in range(len(self.press)):
            self.press_grid[i].fill(self.press[i])
            
        self.lamda_all, self.phi_all = np.meshgrid(self.lamda_range, self.phi_range)
        
    def Layers(self, density=[25.2, 26.2, 27.0, 27.4]):
        self.Z_all = np.zeros((len(density)-1, self.t.shape[1], self.t.shape[2]))
        for i in range(len(density)-1):
            self.Z_all[i] = -np.nanmax(np.where(np.logical_and(self.t >= density[i], self.t <= density[i+1]), self.press_grid, np.nan), axis=0)
            
        self.Z_all = np.concatenate((self.Z_all[:,:,340:], self.Z_all[:,:,:340]), axis=2)
        self.lamda_all = np.concatenate((self.lamda_all[:,340:]-360, self.lamda_all[:,:340]), axis=1)
        self.phi_all = np.concatenate((self.phi_all[:,340:], self.phi_all[:,:340]), axis=1)
            
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
        self.pv = np.abs(np.where(f==0, np.nan, f))/np.where(h==0, np.nan, h)
        return self.lamda_all, self.phi_all, self.pv, self.plot, self.name