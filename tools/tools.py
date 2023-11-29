import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import astropy.constants as const

def RA_DEC_fq(fx):
    
    # RA, DEC in deg.
    RA = np.zeros(fx['PRIMARY'].data.shape[1])
    RA_del = fx[0].header['CDELT1']
    for i in range(len(RA)):
        RA[i] = RA_del*(i+1-fx[0].header['CRPIX1']) + fx[0].header['CRVAL1']
    print("RA range      : {:.2f} to {:.2f}".format(min(RA), max(RA)))

    DEC = np.zeros(fx['PRIMARY'].data.shape[2])
    DEC_del = fx[0].header['CDELT2']
    for i in range(len(DEC)):
        DEC[i] = DEC_del*(i+1-fx[0].header['CRPIX2']) + fx[0].header['CRVAL2']
    print("DEC range     : {:.2f} to {:.2f}".format(min(DEC), max(DEC)))

    fq = np.zeros(fx['PRIMARY'].data.shape[0])
    fq_del = fx[0].header['CDELT3']
    for i in range(len(fq)):
        fq[i] = fq_del*(i+1-fx[0].header['CRPIX3']) + fx[0].header['CRVAL3']
    fq /= 10**6 # Hz into MHz

    # f= 1420/(1+z) MHz
    f_HI = 1420
    z = f_HI/fq -1
    print("redshift range: {:.2f} to {:.2f}".format(min(z), max(z)))
    
    return RA, DEC, fq


def pk_3D_to_2D(kper,T, RA, fq):
    
    # T: (RA,DEC) for each frequency
    # fq: an array for frequencies
    
    print(kper.shape)
    print(T.shape)
    print(RA.shape)
    print(fq.shape)

    k_num=len(kper)
    d_k = np.diff(kper)[0]
    k_min=kper[0] - d_k/2.
    k_max=kper[-1] + d_k/2. 
    
    p_k_f = np.zeros((k_num,k_num))
    N_k_f = np.zeros((k_num,k_num))

    n1=T.shape[1]
    n2=T.shape[0] # freq.
    
    f_c = fq[int(n2/2)] # central freq.

    f_HI = 1420
    z = f_HI/f_c -1
    cosmo = FlatLambdaCDM(H0=100.0, Om0=0.30964, Tcmb0=2.7255, Neff=3.046, Ob0=0.0490)
    comoving_distance = cosmo.comoving_distance(z).value #*cosmo.h # in Mpc
    
    xx=np.pi/180.*RA * comoving_distance # in Mpc
    rLbox1 = max(xx)-min(xx)
    print("  ", rLbox1, " Mpc")
    
    z = f_HI/fq -1
    comoving_distance = cosmo.comoving_distance(z).value #*cosmo.h # in Mpc
    zz=comoving_distance # in Mpc
    rLbox3 = max(zz)-min(zz)
    print("  ", rLbox3, " Mpc")

    k_F1=2.*np.pi/rLbox1 #RA
    k_F3=2.*np.pi/rLbox3 #freq
    
    for m in range(n1): # y
        for n in range(n1): # x
            for l in range(n2): # z
                
                cn = abs(n-int(n1/2))
                cm = abs(m-int(n1/2))
                cl = abs(l-int(n2/2))
        
                kpara = k_F3*cl
                kperp = k_F1*np.sqrt(cn**2 + cm**2)

                if((kpara>=k_min and kpara < k_max) and (kperp>=k_min and kperp < k_max)):
                    i=int((kperp-k_min)/d_k)
                    j=int((kpara-k_min)/d_k) 

                    p_k_f[j,i] = p_k_f[j,i] + abs(T[l,n,m])**2
                    N_k_f[j,i] = N_k_f[j,i] + 1 
            
    p_k_f = (p_k_f/N_k_f)*(rLbox1**2*rLbox3/(n1**4*n2**2))
    return p_k_f, N_k_f

def DensityFlux_to_Temp(fx):
    # This part comes from: https://gitlab.com/flomertens/ps_eor/-/blob/master/ps_eor/datacube.py#L1068
    freq_start = fx[0].header['CRVAL3']
    df = fx[0].header['CDELT3']
    nf = fx[0].header['NAXIS3']
    freq_end = freq_start + (nf - 1) * df
    fits_freqs = np.linspace(freq_start, freq_end, nf)
    lamb = const.c.value / fits_freqs

    cart_map = np.squeeze(fx['PRIMARY'].data)
    _, nx, ny = cart_map.shape
    res = abs(np.radians(fx[0].header['CDELT1']))

    bmaj = fx[0].header['BMAJ']
    bmin = fx[0].header['BMIN']
    omega = np.radians(bmaj) * np.radians(bmin) * np.pi / (4 * np.log(2))
    imager_scale_factor = omega / (res ** 2)

    jy2k = ((1e-26 * lamb ** 2) / (2 * const.k_B.value))
    jypsf2K = jy2k / omega # From Jy/PSF to K
    cart_map = cart_map * jypsf2K[:, None, None]
    return cart_map
