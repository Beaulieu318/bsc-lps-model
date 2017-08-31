import numpy as np
import matplotlib.pyplot as plt
import Ekman as ek
from mpl_toolkits.mplot3d import Axes3D

class LPSAnalysis:
    def __init__(self, data):
        self.lamda_all, self.phi_all, self.Z_all, self.plot, self.name = data
        
        self.TLayer = len(self.Z_all)
        
        omega = 2*np.pi/(24*60*60)
        if self.plot == "Potential Vorticity": self.Z_all = self.Z_all*500/(2*omega)
        
    def RestrictArea(self, area):
        [lamda_min, lamda_max], [phi_min, phi_max] = area
        lamda_min_idx = (np.abs(self.lamda_all[0]-lamda_min)).argmin()
        lamda_max_idx = (np.abs(self.lamda_all[0]-lamda_max)).argmin()
        phi_min_idx = (np.abs(self.phi_all[:,0]-phi_min)).argmin()
        phi_max_idx = (np.abs(self.phi_all[:,0]-phi_max)).argmin()
        self.lamda_all = self.lamda_all[phi_min_idx:phi_max_idx, lamda_min_idx:lamda_max_idx]
        self.phi_all = self.phi_all[phi_min_idx:phi_max_idx, lamda_min_idx:lamda_max_idx]
        self.Z_all = self.Z_all[:,phi_min_idx:phi_max_idx, lamda_min_idx:lamda_max_idx]
    
    def Outcrops(self, Layers=[1,2], regions=[]):
        for Nlayer in Layers:
            layer = self.Z_all[Nlayer-(4-self.TLayer)-1]
            lat = np.where(np.isnan(layer), 0, self.phi_all)
            y_outcrop_args = np.nanargmax(lat, axis=0)
            x_args = np.arange(len(y_outcrop_args))
            plt.plot(self.lamda_all[y_outcrop_args,x_args][:-1], self.phi_all[y_outcrop_args,x_args][:-1], label=r'$y_%s$' % Nlayer)
            #x_args = np.arange(min(self.lamda_all[0]), max(self.lamda_all[0]))
            #y_args = [regions[Nlayer]] * len(x_args)
            #plt.plot(x_args, y_args, label=r'$y_%s$' % Nlayer)
            plt.legend(loc='best')
            plt.show()  
            
    def ShadowZone(self, data):
        lamda_all, phi_all, Z_all, plot, name = data
        shadow = Z_all[2]
        X = np.ma.masked_where(np.isnan(lamda_all),lamda_all)
        Y = np.ma.masked_where(np.isnan(phi_all),phi_all)
        Z = np.ma.masked_where(np.isnan(shadow),shadow)
        CS = plt.contour(X, Y, Z, [0, 2], colors='m')
        plt.clabel(CS, fontsize=20, inline=1, fmt={0:r'$\tilde{x}_3(y)$',2:r'$\tilde{x}_2(y)$'})
        plt.show()
    
    def Plot(self, Layers=[3], plot_type='contour', CS_longs=0, contoursave=False, color=True):
        CS_longs_idx = (np.abs(self.lamda_all[0]-CS_longs)).argmin()
        if plot_type=='surface': self.SurfacePlot_Init()
        for layer in Layers:
            if plot_type=='crosssection': self.CrossSectionPlot(self.phi_all[:,CS_longs_idx], self.Z_all[layer-(4-self.TLayer)-1][:,CS_longs_idx], layer, self.lamda_all[0][CS_longs_idx])
            if plot_type=='surface': self.SurfacePlot(self.lamda_all, self.phi_all, self.Z_all[layer-(4-self.TLayer)-1])
            if plot_type=='contour': self.ContourPlot(self.lamda_all, self.phi_all, self.Z_all[layer-(4-self.TLayer)-1], layer, addmap=False, save=contoursave, color=color)
            if plot_type=='contourmap': self.ContourPlot(self.lamda_all, self.phi_all, self.Z_all[layer-(4-self.TLayer)-1], layer, addmap=True, save=contoursave, color=color)

    def CrossSectionPlot(self, X, Y, layer, longs=""):
        fig = plt.figure('CrossSection ' + self.plot + ' ' + self.name + ' ' + str(self.TLayer))
        plt.plot(X, Y, label=layer)
        plt.title('Cross section of ' + self.plot + ' for ' + self.name + '\n across longitude = ' + str(int(longs)))
        plt.xlabel('Latitude (deg)')
        plt.ylabel(self.plot + ' (m)')
        plt.legend(loc='best')
        plt.xlim(max(X), min(X))
        plt.show()
        
    def ContourPlot(self, X, Y, Z, layer, addmap=False, save=False, color=True):
        figname = 'Contour ' + self.plot  + ' ' + self.name + ' ' + str(self.TLayer) + ' ' + str(layer)
        fig = plt.figure(figname)
        if addmap == True: ek.WorldMapMaskPlot(0)
        X = np.ma.masked_where(np.isnan(X),X)
        Y = np.ma.masked_where(np.isnan(Y),Y)
        Z = np.ma.masked_where(np.isnan(Z),Z)
        if self.plot == "Potential Vorticity":
            plt.title('Contour Plot of ' + self.plot + ' for ' + self.name + ' data \n ' + str(layer) + ' out of ' + str(self.TLayer) + ' layers' + r' ($f/h$ in units of $2\Omega / 500m$)')
            if layer == 3 or layer == 4:
                if color==True:plt.pcolor(X, Y, Z, vmin=0.0, vmax=1.0)
                if color==True:plt.colorbar()
                #CS = plt.contour(X, Y, Z, [0.1, 0.2, 0.3, 0.33, 0.4, 0.5], colors='k')
                CS = plt.contour(X, Y, Z, [0.2, 0.3, 0.33, 0.4, 0.5], colors='k')
            else:
                if color==True:plt.pcolor(X, Y, Z, vmin=0.0, vmax=3.5)
                if color==True:plt.colorbar()
                CS = plt.contour(X, Y, Z, [1., 1.5, 2.0, 2.5, 3.0], colors='k') 
        else:
            plt.title('Contour Plot of ' + self.plot + ' for ' + self.name + ' data \n ' + str(layer) + ' out of ' + str(self.TLayer) + ' layers' + r' (units of $m$)')
            if color==True:plt.pcolor(X, Y, Z)
            if color==True:plt.colorbar()
            #CS = plt.contour(X, Y, Z, colors='k')
            CS = plt.contour(X, Y, Z, [-800, -700, -650, -600, -500], colors='k')
        plt.xlabel('Longitude (deg)')
        plt.ylabel('Latitude (deg)')
        plt.clabel(CS, fontsize=9, inline=1)
        plt.show()
        location = "C:\\Users\\tvb\\OneDrive\\Documents\\Imperial\\Year 3\\Project\\Project - Shared\\Pictures\\" + self.plot + "/" + figname
        if save == True: fig.savefig(location)

    def SurfacePlot_Init(self):
        fig = plt.figure('Surface ' + self.plot  + ' ' + self.name + ' ' + str(self.TLayer))
        self.ax = fig.add_subplot(111, projection='3d')
        self.ax.set_xlabel('lamda - longitude')
        self.ax.set_ylabel('phi - latitude')
        self.ax.set_zlabel('Depth')
        plt.show()
        
    def SurfacePlot(self, X, Y, Z):
        self.ax.plot_surface(X, Y, Z)
        self.ax.set_xlabel('lamda - longitude')
        self.ax.set_ylabel('phi - latitude')
        self.ax.set_zlabel('Depth')
        plt.show()