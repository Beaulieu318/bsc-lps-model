from netCDF4 import Dataset
import numpy as np
import gsw as gsw
import os

class ArgoData:
    def __init__(self, Type = 'Temperature', Month = 'Average'):
        self.Type = Type
        self.Month = Month
        
    def AllData(self):
        if self.Type == 'Temperature' or self.Type == 'Salinity':
            self.Location = os.getcwd() + '/RG_ArgoClim_'+self.Type+'_2016.nc'
            self.__GetMonth()
            self.__GetData()
        elif self.Type == 'Density':
            self.__Density()
        else:
            raise Exception('Type needs to be Temperature, Salinity or Density')
        return self.lons, self.lats, self.press, self.press_units, self.t_all, self.t_units, self.t_title
        
    def __GetData(self):
        self.fh = Dataset(self.Location, mode='r')
        self.lons = self.fh.variables['LONGITUDE'][:]
        self.lats = self.fh.variables['LATITUDE'][:]
        self.press = self.fh.variables['PRESSURE'][:]
        self.time = self.fh.variables['TIME'][:]
        self.t_units = self.fh.variables['ARGO_'+self.Type.upper() +'_MEAN'].units
        self.press_units = self.fh.variables['PRESSURE'].units
        self.fh.close()
        
    def __GetMonth(self):
        self.fh = Dataset(self.Location, mode='r')
        if self.Month == 'Average':
            self.t_title = 'Average'
            self.t_all = self.fh.variables['ARGO_'+self.Type.upper() +'_MEAN'][:]
        elif self.Month >= 0 and self.Month < 144:
            self.t_title = str(self.Month)
            self.t_all = self.fh.variables['ARGO_'+self.Type.upper() +'_ANOMALY'][self.Month]
            return self.t_all
        else:
            raise Exception('Month between 0 and 143')
        self.fh.close()
        
    def __Density(self):
        temp = ArgoData(Type = 'Temperature', Month = self.Month)
        temp.AllData()
        sal = ArgoData(Type = 'Salinity', Month = self.Month)
        sal.AllData()
        
        self.lons = temp.lons
        self.lats = temp.lats
        self.press = temp.press
        self.press_units=temp.press_units
        self.t_units = 'kg/m^3'
        self.t_title = temp.t_title
        press = np.zeros(temp.t_all.shape)    
        for i in range(len(self.press)):
            press[i].fill(self.press[i])
        self.t_all = gsw.sigma0(sal.t_all, temp.t_all) #potential density (rho for density)