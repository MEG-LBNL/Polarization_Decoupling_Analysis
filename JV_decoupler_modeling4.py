#Physical Constants
Na = 6.022e23 # mol-1
F = 96485.3415 #C mol-1
kb = 1.38e-23 #J K-1
kbev = 8.614e-5 #eV K-1
n = 1
ne = 1
N = 1e14 #number catalytic sites/cm^2
Aaq = 1e-5 #solute concentration, M; [A]aq ~ [A]t - [B]aq
E0 = 1.23 #thermodynamic potential for OER-coupled HER
Vph = 0.0 #photovoltage

#Rate constants for anode-limited catalysis (light-coupled anode PEC case)
io = 1e-5 #C/s
ko = io*Na/(F*N) #s-1 per active site
k1 = 6e-5 #cm2/s (CO2)
k1b = k1*Aaq # s-1 
#k2 = 1e3  #s-1

Rs = 57 #ohms (57)
Rshunt = 1.5e4 #ohms (1e4)

#Libraries
import matplotlib.pyplot as mt; import numpy as np
from extractor_v2 import *
mt.rcParams.update({'font.size': 20, 'lines.markersize': 20, 'font.weight': 'bold'})
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import griddata

#BV description of interfacial CT
def ket(alpha, V, T):
    kcathode = ko*np.exp(-(n*F*V*alpha)/(Na*kb*T))
    kanode = ko*np.exp(n*F*V*(1-alpha)/(Na*kb*T))
    return [kcathode, kanode, kanode-kcathode]

def darkec(kcat, alpha, V, T):# dark electrocatalysis - values here assume reactions are cathode-limited
    N = 1e17    
    S = 5e11 #substrate concentration; molecules CO2/cm2
    k2 = kcat
    k1 = 0.158 #cm2 s-1
    k1b = k1*S #s-1
    #k2 = 1e3  #s-1
    ketf, ketb = ket(alpha, V, T)[0:2] #kf = cathode; kb = anode; using generic kamitaka equation for catalysis
    dark = ne*ketf*N*(1 - ketb/(k2 + ketb)) / ( ketf/(ketb + k2) + (1 + (1/(k1*S))*(k1b + k2*ketf/(ketb + k2)) ) ) #derived rate eq.
    return kcat / (1 + kcat/ketf + ketb/ketf), n*N*F*dark/Na

#PV current as a rate constant
def PV(Irs, Vd, T, Iphoto): #Irs = reverse saturation current, Vd = diode voltage, T = temp/K
    VT = kb*T/(1.602e-19) #Thermal voltage
    return (Na/(n*F))*((Vd)/Rshunt + Irs*(np.exp((Vd)/(n*VT) - 1)) - Iphoto)

#PV-E/PEC polarization dependence. Irs = reverse saturation current (Amps); Iphoto = short-circuit photocurrent (Amps).
def PEC(kcat, Voltage, alpha, T, Irs, Iphoto):
    ketf, ketb, ketsum = ket(alpha, Voltage, T)#forward and reverse interfacial ET rates
    k2 = kcat #0.1
    kpv = PV(Irs, Voltage, T, Iphoto)
    f = abs(N*ketsum / (N*ketsum + kpv))
    #ret = k1*ketf*(Aaq + kpv)*(1 - ketf*ketb/(k2 + ketf)) / (1 + k1b + ketf - (ketf*ketb)/(k2+ketf) + ketf )
    #ret = ketf*(Aaq + f*kpv/ketb)*(1 - ketb/(k2 + ketb)) / (1 + (1/k1)*(k1b + ketf - (ketf*ketb)/(k2+ketf)) + ketf/ketb )
    #Model 4
    ret = ne*ketf*(N + f*kpv/ketb)*(1 - ketb/(k2 + ketb)) /  ( (1/(k1*Aaq))*(k1b + k2*ketf/(ketb+k2)) + 1 + ketf/ketb )
    #return [ PEC device, PV current, sum electrode current, cathodic (f) current, anodic (b) current ]
    return [ n*F*ret/Na, n*F*kpv/Na, N*n*F*ketsum/Na, N*n*F*ketf/Na, N*n*F*ketb/Na  ]
   
#e.g. In [2]: plotPEC2(200, [-3, 1.5], 0.2, 293, 1e-10, 5e-3)
def plotPEC2(kcat, Voltages, alpha, T, Irs, Iphoto):#voltage range input as list, e.g. for 0 to 3 volt sweep --> [0, 3]
    V = np.arange(Voltages[0], Voltages[1], 0.01)
    Vdark = np.arange(-2.5, 1.5, 0.01)
    #mt.style.use('bmh')
    ax = mt.subplot(1, 3, 3)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    mt.fill_between(V-E0+Vph, PEC(kcat, V, alpha, T, Irs, Iphoto)[0], color='c') #PEC polarization.
    mt.plot(V-E0+Vph, PEC(kcat, V, alpha, T, Irs, Iphoto)[0], 'c-', linewidth = 1.5, label = 'd') #PEC polarization.
    mt.xlabel('Applied Potential [V]'); mt.ylabel('Current [A]'); mt.title('Decoupled Region')
    
    #Plot PEC device performance, PV current, and Electrode currents
    ax = mt.subplot(1, 3, 2)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    for suns in np.arange(0.0, 1.1, 0.1):
        mt.plot(V-E0+Vph, PEC(kcat, V, alpha, T, Irs, suns*Iphoto)[0], 'c-', linewidth = 1.5, label = 'c') #PEC polarization.
    mt.xlabel('Applied Potential [V]'); mt.ylabel('Current [A]'); mt.title('PEC Polarization')

    """mt.subplot(2, 2, 2)
    for suns in np.arange(0.2, 1.2, 0.2):
        mt.plot(V, PEC(kcat, V, alpha, T, Irs, suns*Iphoto)[1], 'b-', linewidth = 1.5, label = 'b') #PV load curve.
    mt.ylabel('Current [A]'); mt.title('Photovoltaic (PV) only')"""

    mt.subplot(1, 3, 1)
    mt.plot(Vdark-E0, -darkec(kcat, Vdark+E0, alpha, T)[0], 'k-', linewidth = 2.0) #unbounded forward (cathodic) current (ketf).
    mt.ylabel('Current [A]'); mt.title("Dark Polarization"); mt.xlabel('Applied Potential [V]')
    mt.ylim(-0.04, 0.002)
    figure = mt.gcf(); figure.set_size_inches (15,4.5)
    mt.tight_layout(); mt.show()

#e.g. In [2]: plotPEC2(1e-3, [-0.5,0.5], 0.5, 293, 1e-10, 1e-2)
def plotPEC3(kcat, Voltages, alpha, T, Irs, Iphoto):#voltage range input as list, e.g. for 0 to 3 volt sweep --> [0, 3]
    V = np.arange(Voltages[0], Voltages[1], 0.01)
    #mt.style.use('bmh')
    """mt.subplot(1, 2, 2)
    mt.fill_between(V-1.23, PEC(kcat, V, alpha, T, Irs, Iphoto)[0], color='c') #PEC polarization.
    mt.plot(V-1.23, PEC(kcat, V, alpha, T, Irs, Iphoto)[0], 'c-', linewidth = 1.5, label = 'd') #PEC polarization.
    mt.xlabel('Applied Potential [V]'); mt.ylabel('Current [A]')#; mt.title('Decoupled Region')"""
    
    #Plot PEC device performance, PV current, and Electrode currents
    mt.subplot(1, 2, 1)
    #I_list = [0]
    light = np.array([0.06, 0.13, 0.25, 0.55, 0.95, 0.99, 1])
    I_list = [PEC(kcat, Voltages[-1], alpha, T, Irs, 0.01*Iphoto)[0]]
    Icounter = 0
    for Vd in V:
        for suns in light: #np.arange(0.0, 1.1, 0.1):
            I = I_list[Icounter]
            mt.plot(Vd-E0+Vph, PEC(kcat, Vd-I*Rs, alpha, T, Irs, suns*Iphoto)[0], 'c.', markersize = 1.5, label = 'c') #PEC polarization.
            I = PEC(kcat, Vd-I*Rs, alpha, T, Irs, suns*Iphoto)[0]
            I_list.append(I)
            Icounter += 1
        I_list = [0]
        Icounter = 0
    mt.xlabel('Applied Potential [V]'); mt.ylabel('Current [A]'); mt.title('')
    #mt.annotate('Rsh = ' + str(Rshunt) + ' ' + r'$\Omega', (0, 0))
    #mt.annotate('Rs = ' + str(Rs) + ' ' + r'$\Omega', (0, 1))
    mt.xlim(-4, 0)
    #mt.ylim(-0.010, -0.005)
    figure = mt.gcf(); figure.set_size_inches (10,4.5)
    mt.tight_layout()
    mt.show()

#example call: Marcus_contour(np.arange(0.0, 2.1, 0.1), np.arange(2, 30.1, 0.1), 1.0, 1.4, 293)
def Marcus_contour(dVrange, Rrange, lambd, beta, T):
    ratematrix = np.zeros((len(Rrange), len(dVrange)))
    row,col = 0,0
    for dV in dVrange:
        for rda in Rrange:
            k = 1e13 * np.exp(-0.5*beta*rda) * np.exp(-(-dV + lambd)**2 / (4*lambd*kbev*T))
            ratematrix[row,col] = F*k/Na
            row += 1
        col +=1; row = 0
    print (ratematrix)
    mt.style.use('bmh'); levels = np.linspace(np.log10(ratematrix).min(), np.log10(ratematrix).max(), 300)
    mt.contourf(dVrange, Rrange, np.log10(ratematrix), levels, cmap = 'jet')
    mt.colorbar(label = 'log(k$_{ET}$) [s$^-$$^1$]')
    mt.colorbar(label = 'log(i$_{ET}$) [A]', ticks = [-18,-16,-14,-12,-10,-8,-6])
    mt.xlabel('Driving Force [V]'); mt.ylabel('D-A separation [$\AA$]')
    mt.show()

#example call: JVD_contour_sim(1e-3, np.arange(-3, 1.5, 0.01), 0.2, 293, 1e-10, 5e-3)
def JVD_contour_sim(kcat, Voltage, alpha, T, Irs, Iphoto):
    Suns = np.arange(0, 1.01, 0.01)
    currentmatrix = np.zeros((len(Suns), len(Voltage)))
    row,col = 0,0
    for V in Voltage:
        for phi in Suns:
            i = PEC(kcat, V, alpha, T, Irs, phi*Iphoto)[0]
            currentmatrix[row,col] = i
            row += 1
        col +=1; row = 0
    currentmatrix = currentmatrix*1000 #rescale matrix values to mA.
    figure = mt.gcf(); figure.set_size_inches (6,4.8)
    mt.style.use('bmh')
    levels = np.linspace(currentmatrix.min(), currentmatrix.max(), 300)
    mt.contourf(Voltage-E0+Vph, Suns*100, currentmatrix, levels, cmap = 'jet')
    cb = mt.colorbar(label='Current [mA]', ticks=[0,-4,-8,-12,-16,-20]); cb.set_label(label='Current [mA]', weight='normal')
    #ax.set_xticks([2.6,2.8,3.0,3.2])
    #mt.xlabel('Applied Potential [V]', fontweight = 'bold'); mt.ylabel('Illumination [%]', fontweight = 'bold')
    mt.xlabel('Applied Potential [V]'); mt.ylabel('Illumination [%]')
    mt.savefig('Fig_1f_inset(JVD-contour).png', dpi=720, Transparency=True); mt.show()

def JVD_voltagecontour_sim(kcat, Voltage, alpha, T, Irs, Iphoto):
    Suns = np.arange(0, 1.01, 0.01)
    currentmatrix = np.zeros((len(Suns), len(Voltage)))
    row,col = 0,0
    for V in Voltage:
        for phi in Suns:
            i = PEC(kcat, V, alpha, T, Irs, phi*Iphoto)[0]
            currentmatrix[row,col] = i
            row += 1
        col +=1; row = 0
    currentmatrix = currentmatrix*1000 #rescale matrix values to mA.
    figure = mt.gcf(); figure.set_size_inches (6,4.8)
    mt.style.use('bmh')
    levels = np.linspace(currentmatrix.min(), currentmatrix.max(), 300)
    mt.contourf(Voltage-E0+Vph, Suns*100, currentmatrix, levels, cmap = 'jet')
    cb = mt.colorbar(label='Current [mA]', ticks=[0,-4,-8,-12,-16,-20]); cb.set_label(label='Current [mA]', weight='normal')
    #ax.set_xticks([2.6,2.8,3.0,3.2])
    #mt.xlabel('Applied Potential [V]', fontweight = 'bold'); mt.ylabel('Illumination [%]', fontweight = 'bold')
    mt.xlabel('Applied Potential [V]'); mt.ylabel('Illumination [%]')
    mt.savefig('Fig_1f_inset(JVD-contour).png', dpi=720, Transparency=True); mt.show()
    
def JVD_contours_real(header, delim):
    #read in all .DTA files (Gamry pstat format) in a single directory, and extract relevant data:
    directory = np.sort(glob.glob(os.getcwd() + '/PLOT_INPUT/*.DTA'))
    V_app = extract_auto(directory[0], header, delim, 2, 4)[0][:,0] #vector of applied potential values.
        
    #Generate x,y,z matrices for contour mapping
    counter = 0
    for every in directory:
        counter += len(V_app)
    x = np.zeros((len(V_app), len(directory))) #super array contains all values from all powers measured for map generation
    y = np.zeros((len(V_app), len(directory)))
    z = np.zeros((len(V_app), len(directory)))
    
    #Add all values from various power tests to superarrays x,y and z:
    counter = 0
    for every in directory:
        x[:,counter] = V_app #matrix of applied potentials.
        y[:,counter] = extract_auto(every, header, delim, 4, 4)[0][:,0] #matrix of current values.
        z[:,counter] = extract_auto(every, header, delim, 7, 7)[0][:,0] #matrix of floating potentials.
        counter += 1
    
    #A1_avg[:,1] = np.round(A1_avg[:,1], 5) #format array values to only 2 decimal places
    #z = np.around(z, 1)
    
    mt.subplot(1,2,1) #generate contour of applied voltage (x), current (y) and floating back contact voltage (z).     
    figure = mt.gcf(); figure.set_size_inches (10, 4)
    mt.style.use('bmh')
    levels = np.linspace(z.min(), z.max(), 300)
    mt.contourf(x,y,z, levels, cmap = 'jet')
    #mt.contourf(x,y,z, levels, colors = 'r') #for front contact contour
    mt.colorbar(label='Vfloat [V]', ticks = [2.5, 2.6, 2.7, 2.8, 2.9])
    #mt.colorbar(label='Vfloat [V]', ticks = [3.004, 3.005, 3.006, 3.007]) #front contact contour
    for every in directory:
        current = extract_auto(every, header, delim, 4, 4)[0][:,0]
        mt.plot(V_app, current, 'k-', linewidth = 1.5) #plot each separate polarization curve.
    #mt.xlim(-3, 0)
    mt.ylim(-0.0025, 0)
    mt.xlabel('Time [s]'); mt.ylabel('Current [A]')
    """
    mt.subplot(1,2,1) #plot polarization only.     
    figure = mt.gcf(); figure.set_size_inches (10, 4)
    #mt.style.use('bmh')
    for every in directory:
        current = extract_auto(every, header, delim, 4, 4)[0][:,0]
        mt.plot(V_app, current, 'c-', linewidth = 1.5) #plot each separate polarization curve.
    mt.xlim(-4, 0)
    mt.ylim(-0.0025, 0)"""

    mt.savefig('Fig_back-contact_voltage_CA_02152023.png', dpi=720, Transparency=True)
    mt.show() 