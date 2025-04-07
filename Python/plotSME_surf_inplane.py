#!/usr/bin/python
# -*- coding: utf-8 -*-

# Reads ensemble-averaged scattering grids from ADDA and plots a mosaic of the
# intensity and polarization elements to compare the effects of different physical 
# properties such as the number of faces or material (switch on lines 180-184). 

import numpy as np
from matplotlib import pyplot as plt

# maximum zenith angle
Na = 89         # 39 recommended for refraction         
iang = [2, 4, 6]
refr = False    # refraction
            
def draw(inc0, ax, mn, mx):
    # Draw reference lines 
    if refr == True: 
        rfr = 180/np.pi*np.arcsin(np.sin(inc0*np.pi/180.0)/1.55)  # refraction
        ax.plot([rfr,rfr], [mn,mx], '--', color=[0.4,0.4,0.4], lw=2)
    else:
        ax.plot([inc0,inc0], [mn,mx], '--', color=[0.4,0.4,0.4], lw=2)  
        ax.plot([-inc0,-inc0], [mn,mx], '--', color=[0.7,0.7,0.7], lw=2)

def simpson_nonuniform(x, f):
    """
    Simpson rule for irregularly spaced data.
    x: Sampling points for the function values
    f: Function values at the sampling points
    """
    N = len(x) - 1
    h = [x[i + 1] - x[i] for i in range(0, N)]
    assert N > 0

    result = 0.0
    for i in range(1, N, 2):
        h0, h1 = h[i - 1], h[i]
        hph, hdh, hmh = h1 + h0, h1 / h0, h1 * h0
        result += (hph / 6) * ((2 - hdh) * f[i - 1] + (hph**2 / hmh) * f[i] + (2 - 1 / hdh) * f[i + 1])

    if N % 2 == 1:
        h0, h1 = h[N - 2], h[N - 1]
        result += f[N]     * (2 * h1 ** 2 + 3 * h0 * h1) / (6 * (h0 + h1))
        result += f[N - 1] * (h1 ** 2 + 3 * h1 * h0)     / (6 * h0)
        result -= f[N - 2] * h1 ** 3                     / (6 * h0 * (h0 + h1))
    return result

def main(tp, col, ax2, zext):


    Nf = Na * 25
    Nb = 181*25 - Nf
#     angles = np.arange(-Na+1,Na)
    
    Nang = 2*Na - 1
    angles = np.linspace(-Na+1,Na-1,Nang) 
    

    # orientation indices and the orientations in degrees
    x = np.array([0.5,1,1.5,2,2.5,3,4,5,6,7,8])
    wx = x ** (-3.0)
    weights = simpson_nonuniform(x, wx)
    
    for ori in range(3):        
        ScM0 = np.zeros((len(x),Nf,6)) #13213
        P1144 = np.zeros((Nf,6))
        for ix, sp in enumerate(x):
            if tp[0] == 'F':
            # vary shape
                if sp == 0.5:
                    data = np.loadtxt('Ensemble10_%s_Sp02_m217i0004_x05/inc%d_mueller_scatgrid' % (tp,iang[ori]*10), skiprows=1)
                else:
                    try:
                        data = np.loadtxt('Ensemble16_%s_Sp02_m217i0004_x%d/inc%d_mueller_scatgrid' % (tp,sp*10,iang[ori]*10), skiprows=1)
                    except:
                        data = np.loadtxt('Ensemble12_%s_Sp02_m217i0004_x%d/inc%d_mueller_scatgrid' % (tp,sp*10,iang[ori]*10), skiprows=1)

            # vary material
            else:
                if sp == 0.5:
                    data = np.loadtxt('Ensemble10_F12_Sp02_m%s_x05/inc%d_mueller_scatgrid' % (tp,iang[ori]*10), skiprows=1)
                else:
                    try:
                        data = np.loadtxt('Ensemble16_F12_Sp02_m%s_x%d/inc%d_mueller_scatgrid' % (tp,sp*10,iang[ori]*10), skiprows=1)
                    except:
                        data = np.loadtxt('Ensemble12_F12_Sp02_m%s_x%d/inc%d_mueller_scatgrid' % (tp,sp*10,iang[ori]*10), skiprows=1)

            if refr == True:
                ScM0[ix,:,0] = data[Nb:,2]
                ScM0[ix,:,1] = data[Nb:,3]    # F12
                ScM0[ix,:,2] = data[Nb:,7]    # F22
                ScM0[ix,:,3] = data[Nb:,12]    # F33
                ScM0[ix,:,4] = data[Nb:,13]    # F34
                ScM0[ix,:,5] = data[Nb:,17]    # F44
            else:
                ScM0[ix,:,0] = data[:Nf,2]
                ScM0[ix,:,1] = data[:Nf,3]    # F12
                ScM0[ix,:,2] = data[:Nf,7]    # F22
                ScM0[ix,:,3] = data[:Nf,12]    # F33
                ScM0[ix,:,4] = data[:Nf,13]    # F34
                ScM0[ix,:,5] = data[:Nf,17]    # F44
        
        for i in range(Nf):
            for el in range(6):
                wP1144x = wx*ScM0[:,i,el]
                P1144[i,el] = simpson_nonuniform(x, wP1144x) / weights
    
        valuesz = np.reshape(P1144,(Na,25,6))
        
        if refr == True:
            valuesz = np.flip(valuesz, axis=0)

#     #-- Plot ------------------------------------------------
        # Polarizations
        for i in range(1,6):
            polz = valuesz[:,:,i] / valuesz[:,:,0]            
            zen = np.hstack([np.flip(polz[:,0]),polz[1:,12]])
            
            if i == 1:                
                ax2[i,ori].plot(angles, -zen, '-', color=col, lw=2)
            else:
                if tp[0] == 'F':
                    ax2[i,ori].plot(angles, zen, '-', color=col, lw=2, label='%sf' % (tp[1:3]))
                else:
                    ax2[i,ori].plot(angles, zen, '-', color=col, lw=2, label='%s.%s' % (tp[0],tp[1:3]))
            ax2[i,ori].set_ylim([-1.01,1.02])
        
        # Intensity
        zen = np.hstack([np.flip(valuesz[:,0,0]),valuesz[1:,12,0]])
        
        mn, mx = min(zen), max(zen)    
        zext.append(0.9*mn)
        zext.append(1.1*mx)
        
        ax2[0,ori].semilogy(angles, zen, '-', color=col, lw=2)
    
    return zext
	
if __name__ == '__main__':
    
    fig2, ax2 = plt.subplots(nrows=6,ncols=3,sharex=True)
#   For shape comparison, also normal incidence can be included     
#     fig2, ax2 = plt.subplots(nrows=6,ncols=4,sharex=True)
    fig2.subplots_adjust(top=0.96, right=0.96, left=0.08, bottom=0.1, wspace=0.26)

    #ax2[0,0].set_title(r'$\theta_i$: 0°', fontsize=16)
    ax2[0,0].set_title(r'$\theta_i$: 20°', fontsize=16) 
    ax2[0,1].set_title(r'$\theta_i$: 40°', fontsize=16)
    ax2[0,2].set_title(r'$\theta_i$: 60°', fontsize=16)

    ax2[5,0].set_xlabel(r'$\theta_z$ [°]', fontsize=16)
    ax2[5,1].set_xlabel(r'$\theta_z$ [°]', fontsize=16) 
    ax2[5,2].set_xlabel(r'$\theta_z$ [°]', fontsize=16)
    #ax2[5,3].set_xlabel(r'$\theta_z$ [°]', fontsize=14)

    ax2[0,0].set_ylabel(r'$F_{11}$', fontsize=16)
    ax2[1,0].set_ylabel(r'$-F_{12}/F_{11}$', fontsize=16) 
    ax2[2,0].set_ylabel(r'$F_{22}/F_{11}$', fontsize=16)
    ax2[3,0].set_ylabel(r'$F_{33}/F_{11}$', fontsize=16)
    ax2[4,0].set_ylabel(r'$F_{34}/F_{11}$', fontsize=16)
    ax2[5,0].set_ylabel(r'$F_{44}/F_{11}$', fontsize=16)

    for a in ax2.flat:
        a.tick_params(axis='both',labelsize=16)
        a.set_xlim([-Na,Na])
    
    zext = []
    zext = main('217i0004','k',ax2,zext)
    zext = main('279i00155','r',ax2,zext)
    # Optionally, uncomment below and comment the above to compare the shapes
    # zext = main('F12','k',ax2,zext)
    # zext = main('F20','r',ax2,zext)
    
    mn, mx = min(zext), max(zext)
    for k in range(3):
#     for k in range(4):     # if normal incidence is included           
        draw(iang[k]*10,ax2[0,k], mn, mx)
        for i in range(1,6):
            draw(iang[k]*10,ax2[i,k], -1.01, 1.1)
            ax2[0,k].set_xlim([-Na,Na])
        ax2[0,k].set_ylim([mn, mx])
#         ax2[2,k].set_ylim([-0.05, 1.05])
    ax2[2,0].legend(loc=0, prop={'size':11})
    
    plt.show() 