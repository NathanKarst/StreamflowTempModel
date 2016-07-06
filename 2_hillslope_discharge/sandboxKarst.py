from vadoseZone import LaioVadoseZone
from groundwaterZone import GroundwaterZone
from REW import REW
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
        
def main():
    memory_address = lambda input: hex(id(input))
    
    rew = REW([0],[0],[0],LaioVadoseZone,GroundwaterZone,aspect=90)
        
    dt = 1/24.
    s0 = .9
    rew.vz.s = s0
    rew.vz.asd = 1
    T = 0
    s = [s0]
    t = [0]
    storage = [0]
    while T < 5:
        T += dt
        vzFlux = rew.vz.computeFluxes(0)
        rew.vz.s += -sum(vzFlux)*dt
        rew.gz.s += vzFlux[2]*dt
        s.append(rew.vz.s)
        storage.append(rew.gz.s)
        t.append(T)  
        
    b = rew.vz.b
    m = rew.vz.m
    eta = rew.vz.emax
    etaw = rew.vz.ew
    sfc = rew.vz.sfc
    sstar = rew.vz.sstar
    sw = rew.vz.sw
    sh = rew.vz.sh

    sClosedForm = np.array([LaioClosedForm(tt,s0,b,m,eta,etaw,sfc,sstar,sw,sh) for tt in t])
    
    ## 
#     tsfc = 1/(b*(m-eta))*(b*(sfc-s0) + np.log((eta - m + m*np.exp(b*(s0-sfc)))/eta))
#     tsstar = (sfc - sstar)/eta + tsfc
#     tsw     = (sstar-sw)/(eta - etaw)*np.log(eta/etaw) + tsstar
#     
#     t1 = np.linspace(0,tsfc,1000)
#     s1 = s0 - (1/b)*np.log(((eta - m + m*np.exp(b*(s0 - sfc)))*np.exp(b*(eta-m)*t1) - m*np.exp(b*(s0-sfc)))/(eta-m))
#     
#     t2 = np.linspace(tsfc,tsstar,1000)
#     s2 = sfc - eta*(t2 - tsfc)
#     
#     t3 = np.linspace(tsstar,tsw,1000)
#     s3 = sw + (sstar - sw)*(eta/(eta - etaw)*np.exp(-(eta-etaw)/(sstar - sw)*(t3 - tsstar)) - etaw/(eta-etaw))
#     
#     t4 = np.linspace(tsw,np.max(t),1000)
#     s4 = sh + (sw - sh)*np.exp(-etaw/(sw-sh)*(t4-tsw))
#     
#     print((tsfc,tsstar,tsw))
    
    plt.plot(t,s)
    plt.plot(t,sClosedForm,'--')
#     plt.plot(np.hstack((t1,t2,t3,t4)),np.hstack((s1,s2,s3,s4)),'r--')
    plt.plot(t,storage)
#     plt.plot([tsfc, tsfc],[0, 1.5],'k')
    #plt.plot([tsstar, tsstar],[0, 1.5],'k')
    #plt.plot([tsw, tsw],[0, 1.5],'k')
    plt.legend(('Soil moisture: Newton','Soil moisture: closed form','Groundwater zone storage','Field capacity time'))
    plt.xlabel('Time [d]')
    plt.ylabel('Soil moisture []')
    plt.title('Laio model validation')
    
    plt.figure()
    plt.plot(t,np.abs(s-sClosedForm)/s*100)
    plt.title('Validation: Laio simulation vs. Laio closed form')
    plt.xlabel('Time [d]')
    plt.ylabel('Absolute Percent Error')
    
    plt.show()
    

def LaioClosedForm(t,s0,b,m,eta,etaw,sfc,sstar,sw,sh): 
    tsfc = 1/(b*(m-eta))*(b*(sfc-s0) + np.log((eta - m + m*np.exp(b*(s0-sfc)))/eta))
    tsstar = (sfc - sstar)/eta + tsfc
    tsw     = (sstar-sw)/(eta - etaw)*np.log(eta/etaw) + tsstar

    if t < tsfc: return s0 - (1/b)*np.log(((eta - m + m*np.exp(b*(s0 - sfc)))*np.exp(b*(eta-m)*t) - m*np.exp(b*(s0-sfc)))/(eta-m))
    if t < tsstar and t > tsfc: return sfc - eta*(t - tsfc)
    if t < tsw and t > tsstar: return sw + (sstar - sw)*(eta/(eta - etaw)*np.exp(-(eta-etaw)/(sstar - sw)*(t - tsstar)) - etaw/(eta-etaw))
    if t > tsw: return sh + (sw - sh)*np.exp(-etaw/(sw-sh)*(t-tsw))

if __name__ == '__main__': main()