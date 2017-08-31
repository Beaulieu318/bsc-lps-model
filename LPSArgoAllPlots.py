import LPSLayer as ll
import LPSArgoAnalysis as laa
import numpy as np
import ArgoLayer as al

def LPSPlotAll():
    for r in range(1, 4):
        TRegion = r
        
        regions = [10, 25, 40, 50]
        gammas = [0, 9e-3, 6e-3, 5.5e-3]
        lamda_e = 343
        lamda_minmax = [280, 343]
        H_0 = 800
        phi = [40, 25, 10]
        we_max = 30 / float(365*24*60*60)
        ekmanpumping = 'Actual'
        LdataA = ll.LPSLayer(TRegion, regions, gammas, lamda_e, lamda_minmax, H_0, phi, we_max, ekmanpumping)
        LdataA.Layers()
        
        #regions = [10, 35, 40, 45]
        #gammas = [0, 9.5e-3, 6.5e-3, 5.5e-3]
        #lamda_e = 230
        #lamda_minmax = [120, 230]
        #H_0 = 700
        #phi = [40, 25, 10]
        #we_max = 30 / float(365*24*60*60)
        #ekmanpumping = 'Actual'
        #LdataP = ll.LPSLayer(TRegion, regions, gammas, lamda_e, lamda_minmax, H_0, phi, we_max, ekmanpumping)
        #LdataP.Layers()
        
        #LdataAPDepths = (np.concatenate((LdataP.Depths()[0], LdataA.Depths()[0]), axis=1), np.concatenate((LdataP.Depths()[1], LdataA.Depths()[1]), axis=1),np.concatenate((LdataP.Depths()[2], LdataA.Depths()[2]), axis=2),LdataP.Depths()[3], LdataP.Depths()[4]) 
        
        #Adata = al.ArgoLayer()
        #Adata.Layers()
        
        for i in range(TRegion):
            analysis = laa.LPSAnalysis(LdataA.PV())
            analysis.Plot(Layers=[3-i], plot_type='contourmap', CS_longs=50, contoursave=True)
            
def ArgoPlotAll():
    TRegion = 3
    
    Adata = al.ArgoLayer()
    Adata.Layers()
    
    for i in range(TRegion):
        analysis = laa.LPSAnalysis(Adata.PV())
        analysis.Plot(Layers=[3-i], plot_type='contourmap', CS_longs=50, contoursave=True)