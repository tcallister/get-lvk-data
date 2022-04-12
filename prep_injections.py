import numpy as np
import h5py
import json
import astropy.cosmology as cosmo
import astropy.units as u
from astropy.cosmology import Planck15
import sys

sys.path.append('/Users/tcallister/Documents/Repositories/effective-spin-priors/')
from priors import chi_effective_prior_from_isotropic_spins
from priors import joint_prior_from_isotropic_spins

def loadInjections(ifar_threshold,snr_threshold):

    # Read injection file
    mockDetections = h5py.File('./injections/o1+o2+o3_bbhpop_real+semianalytic-LIGO-T2100377-v2.hdf5','r')

    # Total number of trial injections (detected or not)
    nTrials = mockDetections.attrs['total_generated']

    # Read out IFARs and SNRs from search pipelines
    ifar_1 = mockDetections['injections']['ifar_gstlal'][()]
    ifar_2 = mockDetections['injections']['ifar_pycbc_bbh'][()]
    ifar_3 = mockDetections['injections']['ifar_pycbc_hyperbank'][()]
    snr = mockDetections['injections']['optimal_snr_net'][()]

    # Determine which events pass IFAR threshold (O3) or SNR threshold (O1/O2)
    detected_O3 = np.where((ifar_1>ifar_threshold) | (ifar_2>ifar_threshold) | (ifar_3>ifar_threshold))[0]
    detected_O1O2 = np.where((mockDetections['injections']['name'][()]!=b'o3') & (snr>snr_threshold))[0]
    detected_full = np.concatenate([detected_O1O2,detected_O3])
    print(detected_O3.size,detected_O1O2.size,detected_full.size)

    # Get properties of detected sources
    m1_det = np.array(mockDetections['injections']['mass1_source'][()])[detected_full]
    m2_det = np.array(mockDetections['injections']['mass2_source'][()])[detected_full]
    s1x_det = np.array(mockDetections['injections']['spin1x'][()])[detected_full]
    s1y_det = np.array(mockDetections['injections']['spin1y'][()])[detected_full]
    s1z_det = np.array(mockDetections['injections']['spin1z'][()])[detected_full]
    s2x_det = np.array(mockDetections['injections']['spin2x'][()])[detected_full]
    s2y_det = np.array(mockDetections['injections']['spin2y'][()])[detected_full]
    s2z_det = np.array(mockDetections['injections']['spin2z'][()])[detected_full]
    z_det = np.array(mockDetections['injections']['redshift'][()])[detected_full]

    # This is dP_draw/(dm1*dm2*dz*ds1x*ds1y*ds1z*ds2x*ds2y*ds2z)
    precomputed_p_m1m2z_spin = np.array(mockDetections['injections']['sampling_pdf'][()])[detected_full]

    # In general, we'll want either dP_draw/(dm1*dm2*dz*da1*da2*dcost1*dcost2) or dP_draw/(dm1*dm2*dz*dchi_eff*dchi_p).
    # In preparation for computing these quantities, divide out by the component draw probabilities dP_draw/(ds1x*ds1y*ds1z*ds2x*ds2y*ds2z)
    # Note that injections are uniform in spin magnitude (up to a_max = 0.998) and isotropic, giving the following:
    dP_ds1x_ds1y_ds1z = (1./(4.*np.pi))*(1./0.998)/(s1x_det**2+s1y_det**2+s1z_det**2)
    dP_ds2x_ds2y_ds2z = (1./(4.*np.pi))*(1./0.998)/(s2x_det**2+s2y_det**2+s2z_det**2)
    precomputed_p_m1m2z = precomputed_p_m1m2z_spin/dP_ds1x_ds1y_ds1z/dP_ds2x_ds2y_ds2z

    return m1_det,m2_det,s1x_det,s1y_det,s1z_det,s2x_det,s2y_det,s2z_det,z_det,precomputed_p_m1m2z,nTrials

def genInjectionFile(ifar_threshold,snr_threshold,filename):

    # Load
    m1_det,m2_det,s1x_det,s1y_det,s1z_det,s2x_det,s2y_det,s2z_det,z_det,p_draw_m1m2z,nTrials = loadInjections(ifar_threshold,snr_threshold)

    # Derived parameters
    q_det = m2_det/m1_det
    Xeff_det = (m1_det*s1z_det + m2_det*s2z_det)/(m1_det+m2_det)
    Xp_det = np.maximum(np.sqrt(s1x_det**2+s1y_det**2),(3.+4.*q_det)/(4.+3.*q_det)*q_det*np.sqrt(s2x_det**2+s2y_det**2))
    a1_det = np.sqrt(s1x_det**2 + s1y_det**2 + s1z_det**2)
    a2_det = np.sqrt(s2x_det**2 + s2y_det**2 + s2z_det**2)
    cost1_det = s1z_det/a1_det
    cost2_det = s2z_det/a2_det

    # Compute marginal draw probabilities for chi_effective and joint chi_effective vs. chi_p probabilities
    p_draw_xeff = np.zeros(Xeff_det.size)
    p_draw_xeff_xp = np.zeros(Xeff_det.size)
    for i in range(p_draw_xeff.size):
        if i%500==0:
            print(i)
        p_draw_xeff[i] = chi_effective_prior_from_isotropic_spins(q_det[i],1.,Xeff_det[i])
        p_draw_xeff_xp[i] = joint_prior_from_isotropic_spins(q_det[i],1.,Xeff_det[i],Xp_det[i],ndraws=10000)

    # Draw probabilities for component spin magnitudes and tilts
    p_draw_a1a2cost1cost2 = (1./(4.*np.pi))**2*(1./0.998)**2*np.ones(a1_det.size)

    # Combine
    pop_reweight = 1./(p_draw_m1m2z*p_draw_xeff_xp)
    pop_reweight_XeffOnly = 1./(p_draw_m1m2z*p_draw_xeff)
    pop_reweight_noSpin = 1./p_draw_m1m2z

    # Also compute factors of dVdz that we will need to reweight these samples during inference later on
    dVdz = 4.*np.pi*Planck15.differential_comoving_volume(z_det).to(u.Gpc**3*u.sr**(-1)).value

    # Store and save
    injectionDict = {
            'm1':m1_det,
            'm2':m2_det,
            'Xeff':Xeff_det,
            'Xp':Xp_det,
            'z':z_det,
            's1z':s1z_det,
            's2z':s2z_det,
            'a1':a1_det,
            'a2':a2_det,
            'cost1':cost1_det,
            'cost2':cost2_det,
            'dVdz':dVdz,
            'p_draw_m1m2z':p_draw_m1m2z,
            'p_draw_chiEff':p_draw_xeff,
            'p_draw_chiEff_chiP':p_draw_xeff_xp,
            'p_draw_a1a2cost1cost2':p_draw_a1a2cost1cost2,
            'nTrials':nTrials
            }

    for key in injectionDict.keys():
        if key!='nTrials':
            injectionDict[key] = injectionDict[key].tolist()

    with open(filename,"w") as f:
        json.dump(injectionDict,f)

if __name__=="__main__":

    genInjectionFile(1,10.,'./injectionDict_FAR_1_in_1.pickle')
