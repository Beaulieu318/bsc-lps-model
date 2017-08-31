import LPSLayer as ll
import LPSArgoAnalysis as laa
import LPSArgoScatter as las
import numpy as np
import ArgoLayer as al
import Ekman as ek

version = 'New'

if version == 'Paper':
    TRegion = 3
    
    density = np.array([25.2, 26.2, 27.0, 27.4, 28.0])
    regions = [10, 25, 37.5, 55]
    gammas = [0]+list((density[1:]-density[:-1])/np.mean(density+1000)*9.8)[:3]
    lamda_e = 343
    lamda_minmax = [280, lamda_e]
    H_0 = 800
    H_4 = 1400
    phi = [40, 25, 10]
    we_max = 30 / float(365*24*60*60)
    ekmanpumping = 'Quadratic'
    
    LdataA = ll.LPSLayer(TRegion, regions, gammas, lamda_e, lamda_minmax, H_0, H_4, phi, we_max, ekmanpumping)
    LdataA.Layers()
    
    #regions = [10, 35, 40, 45]
    #gammas = [0, 9.5e-3, 6.5e-3, 5.5e-3]
    #lamda_e = 230
    #lamda_minmax = [120, 230]
    #H_0 = 700
    #phi = [40, 25, 10]
    #we_max = 30 / float(365*24*60*60)
    #ekmanpumping = 'Actual'
    #
    #LdataP = ll.LPSLayer(TRegion, regions, gammas, lamda_e, lamda_minmax, H_0, phi, we_max, ekmanpumping)
    #LdataP.Layers()
    #
    #LdataAPDepths = (np.concatenate((LdataP.Depths()[0], LdataA.Depths()[0]), axis=1), np.concatenate((LdataP.Depths()[1], LdataA.Depths()[1]), axis=1),np.concatenate((LdataP.Depths()[2], LdataA.Depths()[2]), axis=2),LdataP.Depths()[3], LdataP.Depths()[4]) 
    
    Adata = al.ArgoLayer()
    Adata.Layers(density)

elif version == 'New':
    TRegion = 3
    
    #density = np.array([24.7, 25.5, 26.2, 27.0, 27.5])
    #regions = [10, 30, 45, 50]
    #gammas = [0]+list((density[1:]-density[:-1])/np.mean(density+1000)*9.8)[:3]
    #lamda_e = 343
    #lamda_minmax = [280, 343]
    #H_0 = 500
    #H_4 = 1100
    #phi = [45, 25, 10]
    #we_max = 30 / float(365*24*60*60)
    #ekmanpumping = 'Actual'
    #LdataA = ll.LPSLayer(TRegion, regions, gammas, lamda_e, lamda_minmax, H_0, H_4, phi, we_max, ekmanpumping)
    #LdataA.Layers()
    
    density = np.array([24.5, 25.2, 25.7, 27.0, 27.4])
    regions = [10, 32, 37, 45]
    gammas = [0]+list((density[1:]-density[:-1])/np.mean(density+1000)*9.8)[:3]
    lamda_e = 235
    lamda_minmax = [120, lamda_e]
    H_0 = 600
    H_4 = 1500
    phi = [40, 25, 10]
    we_max = 30 / float(365*24*60*60)
    ekmanpumping = 'Actual'
    LdataP = ll.LPSLayer(TRegion, regions, gammas, lamda_e, lamda_minmax, H_0, H_4, phi, we_max, ekmanpumping)
    LdataP.Layers()
    
    #LdataAPDepths = (np.concatenate((LdataP.Depths()[0], LdataA.Depths()[0]), axis=1), np.concatenate((LdataP.Depths()[1], LdataA.Depths()[1]), axis=1),np.concatenate((LdataP.Depths()[2], LdataA.Depths()[2]), axis=2),LdataP.Depths()[3], LdataP.Depths()[4]) 
    
    Adata = al.ArgoLayer()
    Adata.Layers(density)

#ek.OutcropPlot(density)
#ek.D0([lamda_minmax, [0,60]], gammas[-1])
#ek.EkmanPlot()
#
#analysis = laa.LPSAnalysis(LdataA.PV())
#analysis.Plot(Layers=[3], plot_type='contourmap', CS_longs=50)
#analysis.Outcrops(Layers=[4,3,2,1])
#analysis = laa.LPSAnalysis(Adata.PV())
##analysis.RestrictArea(area=[lamda_minmax, [regions[0],regions[-1]]])
#analysis.Plot(Layers=[3], plot_type='contourmap', CS_longs=50)
#
#scatter = las.LPSArgoScatter(LdataA.Depths(), Adata.Depths())
#scatter.SortData()
#scatter.ScatterDepthsPlot(layer=4, color='b')
#scatter.ScatterDepthsPlot(layer=3, color='r')
#scatter.ScatterDepthsPlot(layer=2, color='g')
#scatter.ScatterDepthsPlot(layer=1, color='k')
#
#scatter = las.LPSArgoScatter(LdataA.Heights(), Adata.Heights())
#scatter.SortData()
#scatter.ScatterHeightsPlot(layer=4, color='b')
#scatter.ScatterHeightsPlot(layer=3, color='r')
#scatter.ScatterHeightsPlot(layer=2, color='g')
#scatter.ScatterHeightsPlot(layer=1, color='k')