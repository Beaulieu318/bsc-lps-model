import LPSArgoAnalysis as laa
import numpy as np
import matplotlib.pyplot as plt

class LPSArgoScatter:
    def __init__(self, dataM, dataR):
        self.dataM = dataM
        self.dataR = dataR
        
    def SortData(self):
        analysis = laa.LPSAnalysis(self.dataR)
        analysis.RestrictArea([[min(self.dataM[0][0]), max(self.dataM[0][0])], [min(self.dataM[1][:,1]), max(self.dataM[1][:,1])]])
        self.dataR = np.array(self.dataR)
        self.dataR[:3] = [analysis.lamda_all, analysis.phi_all, analysis.Z_all]
        
        self.dataM_arg_x = np.zeros(self.dataR[0].shape, dtype=int)
        self.dataM_arg_y = np.zeros(self.dataR[0].shape, dtype=int)
        for y in range(self.dataR[0].shape[0]):
            for x in range(self.dataR[0].shape[1]):
                self.dataM_arg_y[y][x] = (np.abs(self.dataM[1][:,0] - self.dataR[1][:,0][y])).argmin()
                self.dataM_arg_x[y][x] = (np.abs(self.dataM[0][0] - self.dataR[0][0][x])).argmin()
                
    def ScatterDepthsPlot(self, layer=3, color='r'):
        plt.figure('Depths Scatter Plot')
        plt.scatter(-self.dataM[2][layer-1][self.dataM_arg_y,self.dataM_arg_x], -self.dataR[2][layer-1], color=color, label=layer)
        plt.xlabel('LPS (m)')
        plt.ylabel('ARGO (m)')
        plt.legend()
        plt.title(r'Depths of ARGO against LPS for different layers')
        x_range = np.linspace(0,1600,500)
        y_range = [x for x in x_range]
        plt.plot(x_range, y_range)
        plt.show()
        
    def ScatterHeightsPlot(self, layer=2, color='r'):
        plt.figure('Heights Scatter Plot')
        plt.scatter(self.dataM[2][layer-1][self.dataM_arg_y,self.dataM_arg_x], self.dataR[2][layer-1], color=color, label=layer)
        plt.xlabel('LPS (m)')
        plt.ylabel('ARGO (m)')
        plt.title('Heights of ARGO against LPS for different layers')
        plt.legend()
        x_range = np.linspace(0,800,500)
        y_range = [x for x in x_range]
        plt.plot(x_range, y_range)
        plt.show()