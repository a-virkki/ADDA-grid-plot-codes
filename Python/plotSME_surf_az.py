#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Reads ensemble-averaged scattering grids computed using ADDA and plots the
# intensity in polar plots for the zenith and nadir views and one for both
# in the scattering plane along the horizontal axis for incidence angles 
# 0°, 20°, 40°, and 60°. Other elements can be plotted but the interpretation
# of some elements could be ambiguous out of the scattering plane.

import numpy as np
from matplotlib import pyplot as plt

Na = 81
            
def draw(inc0, ax, mn, mx):
    # Draw reference lines 
    
    inc = np.pi*0.5 + inc0  # incident
    rfl = np.pi*0.5 - inc0  # reflection
    rfr = 1.5*np.pi + np.arcsin(np.sin(inc0)/1.55)  # refraction
        
    ax.plot([inc,inc], [mn,mx], 'k--', lw=1)  
    ax.plot([rfr,rfr], [mn,mx], 'k--', lw=1)
    ax.plot([rfl,rfl], [mn,mx], 'k--', lw=1)        

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
    
    for ori in range(4):
        ScM0 = np.zeros((Nf,2)) #13213
        ScM180 = np.zeros((Nf,2))
        data = np.loadtxt('Ensemble16_F12_Sp02_m217i0004_x60/inc%d_mueller_scatgrid' % (20*ori), skiprows=1)
        
        ScM0[:,0] = data[:Nf,2] 
        ScM180[:,0] = data[Nb:,2] 
        # Uncomment to display a different scattering matrix element:
#         ScM0[:,1] = data[:Nf,17]    # e.g. F21: 6, F22: 7, F44: 17
#         ScM180[:,1] = data[Nb:,17]
    
        # az, dec maps
        valuesz = np.reshape(ScM0,(Na,25,2))
        valuesn = np.reshape(ScM180,(Na,25,2))
        valuesn = np.flip(valuesn,axis=0)

#     #-- Plot ------------------------------------------------
        if sum(ScM0[:,1]) == 0:
            # 1-1 element
            ax2[0,ori].contourf(thetaz+np.pi, rz, np.log10(valuesz[:,:,0].T), cmap=plt.cm.Greens_r)
            ax2[1,ori].contourf(thetan+np.pi, rn, np.log10(valuesn[:41,:,0].T), cmap=plt.cm.Blues_r)
            
            # In-plane scattering
            zen = np.hstack([np.flip(valuesz[:,0,0]),valuesz[1:,12,0]])   # 0 -> 18 and 12 -> 54 for Prop81
            nad = np.hstack([np.flip(valuesn[:41,0,0]),valuesn[1:41,12,0]]) 
            zen = np.log10(zen)
            nad = np.log10(nad)           

        else:
            # Other elements
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