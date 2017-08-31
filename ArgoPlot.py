import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap

class ArgoPlot:
    def __init__(self, ArgoData, plot='Temperature'):
        self.lons, self.lats, self.press, self.press_units, self.t, self.t_units, self.t_title= ArgoData
        self.t_lats = np.rollaxis(self.t, 1)
        self.t_lons = np.rollaxis(self.t, 2)
        self.plot = plot
        
    def Figure(self, Figure = None):
        if Figure != None:
            plt.figure(Figure)
            
    def Surface(self, pressure=0):
        t = self.t[pressure]
        m = Basemap(projection='mill',llcrnrlat=-70,urcrnrlat=80,\
                    llcrnrlon=10,urcrnrlon=360,resolution='c')
                    
        self.m = m
                    
        lon, lat = np.meshgrid(self.lons, self.lats)
        xi, yi = m(lon, lat)
        cs = m.pcolormesh(xi,yi,np.squeeze(t))
        m.drawparallels(np.arange(-60., 60.1, 30.), labels=[1,0,0,0], fontsize=10)
        m.drawmeridians(np.arange(0., 360., 40.), labels=[0,0,0,1], fontsize=10)
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        cbar.set_label(self.t_units)
        plt.title('ARGO ' + self.plot + ': ' + self.t_title + ' & ' + str(self.press[pressure]) + self.press_units)
        plt.show()
        
    def Coordinates(self, xpt, ypt):
        lonpt, latpt = self.m(xpt,ypt,inverse=True)
        return lonpt-self.lons[0], latpt-self.lats[0]
            
    def CrossSection(self, lats=None, lons=None):
        if lats != None:
            t = self.t_lats[lats]
            latlons = self.lons
            latlons_used = self.lats[lats]
            latlons_label = "Latitude"
            latlons_xaxis = "Longitude"
        elif lons != None:
            t = self.t_lons[lons]
            latlons = self.lats
            latlons_used = self.lons[lons]
            latlons_label = "Longitude"
            latlons_xaxis = "Latitude"
            
        latlons, press = np.meshgrid(latlons, -self.press)
        plt.pcolor(latlons, press, t)
        plt.colorbar()
        CS = plt.contour(latlons, press, t, [24.7, 25.5, 26.2, 27.0, 27.5], colors='k')
        plt.clabel(CS, fontsize=10, inline=1)
        plt.ylabel('Depth (m)')
        plt.xlabel(latlons_xaxis + ' (deg)')
        #plt.gca().invert_xaxis()
        plt.title('ARGO ' + self.plot + ': ' + self.t_title + ' & ' + str(latlons_used) + ' ' + latlons_label)
        plt.show()