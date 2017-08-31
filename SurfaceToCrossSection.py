import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import ArgoData as ad
import ArgoPlot as ap

def on_press(event):
    lon, lat = Plot.Coordinates(event.xdata, event.ydata)
    plt.clf()
    plt.subplot(221)
    Plot.Surface(pressure=0)
    fig = plt.subplot(222)
    Plot.CrossSection(lons=int(lon))
    fig = plt.subplot(212)
    Plot.CrossSection(lats=int(lat))
    event.canvas.figure.gca().grid()
    event.canvas.draw()

plot = 'Density'

figure = plt.figure()
Data = ad.ArgoData(Type=plot, Month = 'Average')
Plot = ap.ArgoPlot(Data.AllData(), plot)
#plt.subplot(221)
#Plot.Surface(pressure=0)
#fig = plt.subplot(222)
#Plot.CrossSection(lons=342.5-Plot.lons[0])
#fig = plt.subplot(212)
Plot.CrossSection(lats=50-Plot.lats[0])
plt.show()
figure.canvas.mpl_connect('button_press_event', on_press)