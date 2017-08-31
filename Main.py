import ArgoData as ad
import ArgoPlot as ap

def LongitudePlot():
    Data = ad.ArgoData()
    Data.Month()
    Plot = ap.ArgoPlot(Data.AllData())
    Plot.Figure('Longitude')
    Plot.CrossSection(lats=60)
    
def SurfacePlot():
    Data = ad.ArgoData()
    Data.Month()
    Plot = ap.ArgoPlot(Data.AllData())
    Plot.Figure('Surface')
    Plot.Surface(pressure=0)
    
LongitudePlot()
SurfacePlot()

