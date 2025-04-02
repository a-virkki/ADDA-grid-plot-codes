#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Reads ensemble-averaged scattering grids computed using ADDA, computes
# the size-frequency-distribution weighting for size parameters 0.5, 1, 1.5, 
# 2, 2.5, 3, 4, 5, 6, 7, and 8 using nonuniform Simpson's rule and plots the
# intensity in polar plots for the zenith and nadir views and one for both
# in the scattering plane along the horizontal axis for incidence angles 
# 0°, 20°, 40°, and 60°. Other elements can be plotted but the interpretation
# of some elements could be ambiguous out of the scattering plane.

import numpy as np
from matplotlib import pyplot as plt

Na = 81
particletype = 'F12_Sp02_m217i0004'
            
def draw(inc0, ax, mn, mx):
    # Draw reference lines 
    
    inc = np.pi*0.5 + inc0  # incident
    rfl = np.pi*0.5 - inc0  # reflection
    rfr = 1.5*np.pi + np.arcsin(np.sin(inc0)/1.55)  # refraction
        
    ax.plot([inc,inc], [mn,mx], 'k--', lw=1)  
    ax.plot([rfr,rfr], [mn,mx], 'k--', lw=1)
    ax.plot([rfl,rfl], [mn,mx], 'k--', lw=1)

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
        

def main():

    Nf = Na * 25
    Nb = 4525-Nf
#     angles = np.arange(-Na+1,Na)
    angles = np.linspace(0,6.283,360)
    
#     #-- Contour plot ------------------------------------------------
    zeniths = np.arange(Na)
    azimuths = np.radians(np.arange(0, 361, 15))
    rz, thetaz = np.meshgrid(zeniths, azimuths)
    rn, thetan = np.meshgrid(zeniths[:41], azimuths)
    
    fig2, ax2 = plt.subplots(nrows=3,ncols=4, subplot_kw=dict(projection='polar'),figsize=(7,6))
    fig2.subplots_adjust(top=0.92, right=0.96, left=0.06, bottom=0.06, hspace=0.36)
    ax2[0,0].set_title(r'$\theta_i$: 0°', fontsize=16)
    ax2[0,1].set_title(r'$\theta_i$: 20°', fontsize=16) 
    ax2[0,2].set_title(r'$\theta_i$: 40°', fontsize=16)
    ax2[0,3].set_title(r'$\theta_i$: 60°', fontsize=16)

    for a in ax2.flat:
        a.set_xticklabels([])
        a.set_rlabel_position(-30)

    # orientation indices and the orientations in radians
    iang = [0, np.pi/9.0, np.pi/4.5, np.pi/3.0]
    x = np.array([0.5,1,1.5,2,2.5,3,4,5,6,7,8])
    wx = x ** (-3.0)
    weights = simpson_nonuniform(x, wx)
    
    for ori in range(4):
        ScM0 = np.zeros((len(x),Nf,2)) #13213
        ScM180 = np.zeros((len(x),Nf,2))
        P11n = np.zeros((Nf,2)) #13213
        P11z = np.zeros((Nf,2))
        for ix, sp in enumerate(x):
            if sp == 0.5:
                data = np.loadtxt('Ensemble10_%s_x05/inc%d_mueller_scatgrid' % (particletype,20*ori), skiprows=1)
            else:
                try:
                    data = np.loadtxt('Ensemble16_%s_x%d/inc%d_mueller_scatgrid' % (particletype,10*sp,20*ori), skiprows=1)
                except:
                    data = np.loadtxt('Ensemble12_%s_x%d/inc%d_mueller_scatgrid' % (particletype,10*sp,20*ori), skiprows=1)
        
            ScM0[ix,:,0] = data[:Nf,2]
            ScM180[ix,:,0] = data[Nb:,2] 
        # Uncomment to display a different scattering matrix element:
#         ScM0[ix,:,1] = data[:Nf,17]   # e.g. F21: 6, F22: 7, F44: 17
#         ScM180[ix,:,1] = data[Nb:,17]
        for i in range(Nf):
            for el in range(2):
                wP11zx = wx * ScM0[:,i,el]
                wP11nx = wx * ScM180[:,i,el]
                P11z[i,el] = simpson_nonuniform(x, wP11zx) / weights
                P11n[i,el] = simpson_nonuniform(x, wP11nx) / weights
    
        # az, dec maps
        valuesz = np.reshape(P11z,(Na,25,2))
        valuesn = np.reshape(P11n,(Na,25,2))
        valuesn = np.flip(valuesn,axis=0)

#     #-- Plot ------------------------------------------------
        if np.sum(ScM0[:,:,1]) == 0:
            ax2[0,ori].contourf(thetaz+np.pi, rz, np.log10(valuesz[:,:,0].T), cmap=plt.cm.Greens_r)
            ax2[1,ori].contourf(thetan+np.pi, rn, np.log10(valuesn[:41,:,0].T), cmap=plt.cm.Blues_r)
            
            # In-plane scattering
            zen = np.hstack([np.flip(valuesz[:,0,0]),valuesz[1:,12,0]])   # 0 -> 18 and 12 -> 54 for Prop81
            nad = np.hstack([np.flip(valuesn[:41,0,0]),valuesn[1:41,12,0]]) 
            zen = np.log10(zen)
            nad = np.log10(nad)           

        else:
            polz = valuesz[:,:,1] / valuesz[:,:,0]
            poln = valuesn[:,:,1] / valuesn[:,:,0]
            ax2[0,ori].contourf(thetaz+np.pi, rz, polz.T, cmap=plt.cm.Greens_r)
            ax2[1,ori].contourf(thetan+np.pi, rn, poln[:41,:].T, cmap=plt.cm.Blues_r)
                        
            zen = np.hstack([np.flip(polz[:,0]),polz[1:,12]])
            nad = np.hstack([np.flip(poln[:41,0]),poln[1:41,12]])
            
        mn, mx = min([min(zen),min(nad)]), max([max(zen),max(nad)])
            
        draw(iang[ori],ax2[2,ori], mn, mx)
                        
        ax2[2,ori].plot(np.flip(angles[10:171]), zen, 'g-', lw=2)
        ax2[2,ori].plot(angles[230:311], nad, 'b-', lw=2)  
          
    
    plt.show()    
	
if __name__ == '__main__':
    main()