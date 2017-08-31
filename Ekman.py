import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import os

def EkmanActual(region=1):
    LPS_Wek_filename = os.getcwd() + "/LPS_Wek.mat"
    LPS_ARGOml_filename = os.getcwd() + "/LPS_ARGOml.mat"
    
    LPS_Wek = sio.loadmat(LPS_Wek_filename)
    LPS_ARGOml = sio.loadmat(LPS_ARGOml_filename)
    
    wek = LPS_Wek['wek']
    xw = LPS_Wek['xs']
    yw = LPS_Wek['ys']
    mask = LPS_Wek['mask']
    
    xm = LPS_ARGOml['xs']
    ym = LPS_ARGOml['ys']
    
    wek_new = np.nan * wek
    mask_new = np.nan * mask
    xw_new = xw
    if region==1:
        wek_new[0:240,:] = wek[240:480,:]
        wek_new[240:480,:] = wek[0:240,:]
        mask_new[0:240,:] = mask[240:480,:]
        mask_new[240:480,:] = mask[0:240,:]
        xw_new = xw-180
    else:
        wek_new = wek
        wek_new = wek
        mask_new = mask
        mask_new = mask
        xw_new = xw
    
    xw_new = np.squeeze(xw_new)
    yw = np.squeeze(yw)
        
    lon, lat = np.meshgrid(xw_new, yw)
    wek_new = np.rollaxis(wek_new,1)
    mask_new = np.rollaxis(mask_new,1)
    
    return lon, lat, wek_new, mask_new
    
def OutcropsActual(region=1):
    LPS_Wek_filename = os.getcwd() + "/LPS_Wek.mat"
    LPS_ARGOml_filename = os.getcwd() + "/LPS_ARGOml.mat"
    
    LPS_Wek = sio.loadmat(LPS_Wek_filename)
    LPS_ARGOml = sio.loadmat(LPS_ARGOml_filename)
    
    wek = LPS_Wek['wek']
    xw = LPS_Wek['xs']
    yw = LPS_Wek['ys']
    mask = LPS_Wek['mask']
    
    xm = LPS_ARGOml['xs']
    ym = LPS_ARGOml['ys']
    SSDsub = LPS_ARGOml['SSDsub']
    
    SSDsub = SSDsub
    mask_new = np.nan * mask
    xw_new = xm
    
    xw_new = np.squeeze(xw_new)
    ym = np.squeeze(ym)
        
    lon, lat = np.meshgrid(xw_new, ym)
    SSDsub = np.rollaxis(SSDsub,1)
    mask_new = np.rollaxis(mask_new,1)
    
    return lon, lat, SSDsub, mask_new
    
    
    
def EkmanActualIntegral(lamda_deg, phi_deg, lamda_e_deg, ekman):
    R = 6370e3 #Radius of the earth
    lon, lat, wek, mask = ekman
    wek = wek*mask
    lon_idx = (np.abs(lon[0]-lamda_deg)).argmin()
    lat_idx = (np.abs(lat[:,0]-phi_deg)).argmin()
    lon_e_idx = (np.abs(lon[0]-lamda_e_deg)).argmin()
    lon_integrate = lon[lat_idx, lon_idx:lon_e_idx]
    wek_integrate = wek[lat_idx, lon_idx:lon_e_idx]
    diff_lon_integrate = lon_integrate[1:] - lon_integrate[:-1]
    diff_x_integrate = R * np.cos(np.radians(lat[lat_idx][0])) * np.radians(diff_lon_integrate)
    ekmanintegral = np.trapz(y=wek_integrate, dx=diff_x_integrate)
    return ekmanintegral
    
def D0(area, gamma_3):
    omega = 2*np.pi/(24*60*60)
    R = 6370e3 #Radius of the earth
    [lamda_min, lamda_e_deg], [phi_min, phi_max] = area
    ekman = EkmanActual(0)
    lamda_range = np.linspace(lamda_min, lamda_e_deg, 100)
    phi_range = np.linspace(phi_min, phi_max, 100)
    area = np.zeros((len(phi_range), len(lamda_range)))
    for lamda in range(len(lamda_range)):
        for phi in range(len(phi_range)):
            f = 2 * omega * np.sin(np.radians(phi_range[phi]))
            beta = 2 * omega / R * np.cos(np.radians(phi_range[phi]))
            area[phi][lamda] = - 2 * f**2 / (beta * gamma_3) * EkmanActualIntegral(lamda_range[lamda], phi_range[phi], lamda_e_deg, ekman)
    area /= float(1000**2)
    lamda_all, phi_all = np.meshgrid(lamda_range, phi_range)
    plt.figure('D0Squared')
    plt.title(r'$D_0^2$ in units of km$^2$')
    plt.pcolor(lamda_all, phi_all, area, vmin=-0.5, vmax=0.5)
    plt.xlabel('Longitude (deg)')
    plt.ylabel('Latitude (deg)')
    plt.colorbar()
    plt.show()
    

def EkmanQuadratic(lamda_deg, phi_deg, phi = [50, 30, 0], we_max=9.51e-7):
    phi_n, phi_min, phi_s = phi
    we = -we_max #wind speed at eastern boundry
    a = we /  (phi_min**2 - phi_min*(phi_s + phi_n) - phi_n**2 + phi_n*(phi_s+phi_n))
    b = - a * (phi_s+phi_n)
    c = -phi_n**2*a- phi_n*b
    return a*phi_deg**2 + b*phi_deg + c
    
def EkmanQuadraticIntegral(lamda_deg, phi_deg, lamda_e_deg, phi = [50, 30, 0], we_max=9.51e-7):
    R = 6370e3 #Radius of the earth
    dist = R * np.cos(np.radians(phi_deg)) * (np.radians(lamda_e_deg) - np.radians(lamda_deg))
    ekmanintegral = EkmanQuadratic(lamda_deg, phi_deg, phi, we_max) * dist
    return ekmanintegral
    
    
    
def EkmanPlot(region=0):
    lon, lat, wek, mask = EkmanActual(region)
    day = 24*3600
    mperyr = 1./(365.*float(day))
    wek_peryear = wek/mperyr
    plt.figure('RealEkmanPlot')
    plt.pcolor(lon, lat, wek_peryear*mask, vmin=-50, vmax=50)
    plt.colorbar()
    plt.title('Ekman pumping in units of m/year')
    plt.xlabel('Longitude (deg)')
    plt.ylabel('Latitude (deg)')
    plt.show()
    
def EkmanQuadraticPlot():
    phi = [40, 25, 10]
    we_max = 30 / float(365*24*60*60)
    day = 24*3600
    mperyr = 1./(365.*float(day))
    lon, lat, wek, mask = EkmanActual(0)
    wek_quad = EkmanQuadratic(lon, lat, phi, we_max)
    plt.figure('QuadraticEkmanPlot')
    wek_peryear = wek_quad/mperyr
    plt.pcolor(lon, lat, wek_peryear*mask, vmin=-50, vmax=50)
    plt.colorbar()
    plt.title('Ekman pumping in units of m/year')
    plt.xlabel('Longitude (deg)')
    plt.ylabel('Latitude (deg)')
    plt.show()    
    
def WorldMapMaskPlot(region=1):
    lon, lat, wek, mask = EkmanActual(region)
    plt.pcolor(lon, lat, mask, vmin=-50, vmax=50)
    plt.title('Ekman')
    plt.show()
    return mask
    
def OutcropPlot(density):
    plt.figure('Outcrop')
    WorldMapMaskPlot(region=0)
    lon, lat, SSDsub, mask = OutcropsActual(region=0)
    lon = 180+lon
    SSDsub_new = SSDsub * np.nan
    SSDsub_new[:,0:180] = SSDsub[:,180:360]#[:,::-1]
    SSDsub_new[:,180:360] = SSDsub[:,0:180]
    plt.contour(lon, lat, SSDsub_new, [density[3]], colors='b')
    plt.contour(lon, lat, SSDsub_new, [density[2]], colors='m')
    plt.contour(lon, lat, SSDsub_new, [density[1]], colors='w')
    plt.title('Outcrop lines from Stommel')
    plt.xlabel('Longitude (deg)')
    plt.ylabel('Latitude (deg)')
    plt.show()
    