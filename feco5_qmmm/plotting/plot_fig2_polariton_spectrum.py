import numpy as np
from scipy import signal
from scipy import fftpack
import glob
import columnplots as clp
import MDAnalysis as mda

'''
Note that IR calculations should be run at NVE ensemble!!!

Outside the cavity, a very long simulation of a single Fe(CO)5 is presented

Dipole autocorrelaction function will be calculated!

Inside the cavity, a 1 ps simulation of 16 Fe(CO)5 will be given

Photon autocorrelation function will be calculated!
'''

'''
Universal parameters
'''

au2fs = 0.02418884254
au2eV = 27.211399
eV2cminv = 8065.540106923572
fsinv2eV = 4.135668 # 1fs-1 to 4.13 eV
fsinv2cminv = fsinv2eV * eV2cminv

'''
Parameters for controlling the Pade approximation
'''
sigma = 5e5
w_step = 1e-5
e_cutoff_ev = 0.5
e_cutoff_au = e_cutoff_ev / au2eV
e_start_ev = 0.01
e_start_au = e_start_ev / au2eV


def fft_pade(time,signal,sigma=sigma,max_len=None,w_min=e_start_au,w_max=e_cutoff_au,w_step=w_step,read_freq=None):
    """ Routine to take the Fourier transform of a time signal using the method
          of Pade approximants.
        Inputs:
          time:      (list or Numpy NDArray) signal sampling times
          signal:    (list or Numpy NDArray)
        Optional Inputs:
          sigma:     (float) signal damp factor, yields peaks with
                       FWHM of 2/sigma
          max_len:   (int) maximum number of points to use in Fourier transform
          w_min:     (float) lower returned frequency bound
          w_max:     (float) upper returned frequency bound
          w_step:    (float) returned frequency bin width
        Returns:
          fsignal:   (complex NDArray) transformed signal
          frequency: (NDArray) transformed signal frequencies
        From: Bruner, Adam, Daniel LaMaster, and Kenneth Lopata. "Accelerated
          broadband spectra using transition signal decomposition and Pade
          approximants." Journal of chemical theory and computation 12.8
          (2016): 3741-3750.
    """

    # center signal about zero
    signal = np.asarray(signal) - signal[0]

    stepsize = time[1] - time[0]

    # Damp the signal with an exponential decay.
    damp = np.exp(-(stepsize*np.arange(len(signal)))/float(sigma))
    signal *= damp

    M = len(signal)
    N = int(np.floor(M / 2))

    # Check signal length, and truncate if too long
    if max_len:
        if M > max_len:
            N = int(np.floor(max_len / 2))

    # G and d are (N-1) x (N-1)
    # d[k] = -signal[N+k] for k in range(1,N)
    d = -signal[N+1:2*N]

    try:
        from scipy.linalg import toeplitz, solve_toeplitz
        # Instead, form G = (c,r) as toeplitz
        #c = signal[N:2*N-1]
        #r = np.hstack((signal[1],signal[N-1:1:-1]))
        b = solve_toeplitz((signal[N:2*N-1],\
            np.hstack((signal[1],signal[N-1:1:-1]))),d,check_finite=False)
    except (ImportError,np.linalg.linalg.LinAlgError) as e:
        # OLD CODE: sometimes more stable
        # G[k,m] = signal[N - m + k] for m,k in range(1,N)
        G = signal[N + np.arange(1,N)[:,None] - np.arange(1,N)]
        b = np.linalg.solve(G,d)

    # Now make b Nx1 where b0 = 1
    b = np.hstack((1,b))

    # b[m]*signal[k-m] for k in range(0,N), for m in range(k)
    a = np.dot(np.tril(toeplitz(signal[0:N])),b)
    p = np.poly1d(np.flip(a))
    q = np.poly1d(np.flip(b))

    if read_freq is None:
        # choose frequencies to evaluate over
        frequency = np.arange(w_min,w_max,w_step)
    else:
        frequency = read_freq

    W = np.exp(-1j*frequency*stepsize)

    fsignal = p(W)/q(W)

    return frequency, fsignal

def fft(x, dtfs):
    # Adding zeros to the end of x
    #N = 1500
    #x = np.pad(x, (0, N), 'constant')
    lineshape = fftpack.dct(x, type=1)
    freq_au = np.linspace(0, 0.5/dtfs * 1e15, len(x))
    # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
    freq_cminverse = freq_au / (100.0 * 299792458.0)
    # Calculate spectra
    #field_description =  freq_au**2
    field_description =  freq_au**2
    spectra = lineshape * field_description
    return freq_cminverse, spectra
    #return freq_cminverse[0:spectra.size//2], spectra[0:spectra.size//2]

def smooth(x,window_len=1,window='hamming'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if window_len<3:
        return x

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len//2-1:-window_len//2]

def auto_correlation_function_simple(x):
    n = x.size
    if n % 2 == 0:
        x_shifted = np.zeros(n*2)
    else:
        x_shifted = np.zeros(n*2-1)
    x_shifted[n//2 : n//2+n] = x
    # Convolute the shifted array with the flipped array, which is equivalent to performing a correlation
    autocorr_full = (signal.fftconvolve(x_shifted, x[::-1], mode='same')[-n:]/ np.arange(n, 0, -1))
    # Truncate the autocorrelation array
    autocorr = autocorr_full[0:n//2]
    return autocorr

def get_dipole_spectrum(filename, method="fft", dt_fs=0.5, nmax=20000):
    data = np.loadtxt(filename)
    nsize = data[:,0].size
    if nsize < nmax:
        nmax = nsize
    nstart = int(nmax * 0.)
    print("outside cavity run for %.2f ps" %(nmax * 0.5 * 1e-3))
    dipole_x, dipole_y, dipole_z = data[nstart:nmax,0], data[nstart:nmax,1], data[nstart:nmax,2]
    dacf_x = auto_correlation_function_simple(dipole_x)
    dacf_y = auto_correlation_function_simple(dipole_y)
    dacf_z = auto_correlation_function_simple(dipole_z)
    if method == "fft":
        dacf_x_freq, dacf_x_sp = fft(dacf_x, dt_fs)
        dacf_y_freq, dacf_y_sp = fft(dacf_y, dt_fs)
        dacf_z_freq, dacf_z_sp = fft(dacf_z, dt_fs)
    elif method == "pade":
        t = np.linspace(0.0, dt_fs*(dipole_x.size -1), dipole_x.size) / au2fs
        dacf_x_freq, dacf_x_sp = fft_pade(t, dacf_x)
        dacf_y_freq, dacf_y_sp = fft_pade(t, dacf_y)
        dacf_z_freq, dacf_z_sp = fft_pade(t, dacf_z)
        dacf_x_freq *= au2eV * eV2cminv
        dacf_y_freq *= au2eV * eV2cminv
        dacf_z_freq *= au2eV * eV2cminv
        dacf_x_sp = np.abs(dacf_x_sp) * dacf_x_freq**2
        dacf_y_sp = np.abs(dacf_y_sp) * dacf_x_freq**2
        dacf_z_sp = np.abs(dacf_z_sp) * dacf_x_freq**2

    sp_tot = smooth(dacf_x_sp + dacf_y_sp + dacf_z_sp)
    #sp_tot = dacf_x_sp + dacf_y_sp + dacf_z_sp
    sp_tot /= np.max(sp_tot)
    return dacf_x_freq, sp_tot

def get_polariton_traj_from_xyz(xyz_filename, nmax=450):
    traj = mda.Universe(xyz_filename)
    frames = traj.trajectory
    nframes = len(frames)
    print("In %s nframes = %d" %(xyz_filename, nframes))
    qx_lst, qy_lst = [], []
    for idx, ts in enumerate(frames):
        current_coord = ts._pos.copy()
        qx = current_coord[176, 0]
        qy = current_coord[177, 1]
        qx_lst.append(qx)
        qy_lst.append(qy)
    qx_lst = np.array(qx_lst)
    qy_lst = np.array(qy_lst)
    return qx_lst[:nmax], qy_lst[:nmax]

def get_polariton_traj_from_out(filename, nmax=450):
    data = np.loadtxt(filename)
    qx, qy = data[:,6], data[:,10]
    print("Size of the system is", qx.size)
    return qx[:nmax], qy[:nmax]

def get_ph_spectrum(qx, qy, method="fft", dt_fs=2):
    dacf_x = auto_correlation_function_simple(qx)
    dacf_y = auto_correlation_function_simple(qy)
    if method == "fft":
        dacf_x_freq, dacf_x_sp = fft(dacf_x, dt_fs)
        dacf_y_freq, dacf_y_sp = fft(dacf_y, dt_fs)
    elif method == "pade":
        t = np.linspace(0.0, dt_fs*(qx.size -1), qx.size) / au2fs
        dacf_x_freq, dacf_x_sp = fft_pade(t, dacf_x)
        dacf_y_freq, dacf_y_sp = fft_pade(t, dacf_y)
        dacf_x_freq *= au2eV * eV2cminv
        dacf_y_freq *= au2eV * eV2cminv
        dacf_x_sp = np.abs(dacf_x_sp) * dacf_x_freq**2
        dacf_y_sp = np.abs(dacf_y_sp) * dacf_x_freq**2

    sp_tot = smooth(dacf_x_sp + dacf_y_sp)
    #sp_tot = dacf_x_sp + dacf_y_sp
    sp_tot /= np.max(sp_tot)
    return dacf_x_freq, sp_tot

fig, axes = clp.initialize(2, 2, width=7, fontsize=12, return_fig_args=True, LaTeX=True,
                            labelthem=True, labelthemPosition=[0.12, 0.94],
                            commonY=[-0.015, 0.5, "Normalized IR intensity [arb. units]"] )

# Prepare data outside the cavity
dipole_traj_outcav_thermostat = "../qm_solvent_free/outcav_eq_weak_bath_longsimulation/dipole_traj"
dipole_traj_outcav_qmmm = "../qmmm/outcav_eq_weak_bath_longsimulation_rerun/dipole_traj"
freq_d_thermostat, sp_d_thermostat = get_dipole_spectrum(dipole_traj_outcav_thermostat, method="fft")
freq_d_qmmm, sp_d_qmmm = get_dipole_spectrum(dipole_traj_outcav_qmmm, method="fft")
freq_d_thermostat_pade, sp_d_thermostat_pade = get_dipole_spectrum(dipole_traj_outcav_thermostat, method="pade")
freq_d_qmmm_pade, sp_d_qmmm_pade = get_dipole_spectrum(dipole_traj_outcav_qmmm, method="pade")

# plot IR spectrum of Fe(CO)5 outside the cavity
xs = [freq_d_thermostat_pade, freq_d_thermostat]
ys = [sp_d_thermostat_pade, sp_d_thermostat]
colors = ["k", "c"]
labels = [r"Pad\'e", "FFT"]
clp.plotone(xs, ys, axes[0,0], colors=colors, labels=labels, lw_lst=[1, 1, 0.5], xlim=[1905, 2120],
            alphaspacing=0.1,
            ylim=[0, 1.2])


xs = [freq_d_qmmm_pade, freq_d_qmmm]
ys = [sp_d_qmmm_pade, sp_d_qmmm]
clp.plotone(xs, ys, axes[1,0], colors=colors, labels=labels, lw_lst=[1, 1, 0.5], xlim=[1905, 2120],
            showlegend=False, alphaspacing=0.1,
            xlabel="frequency [cm$^{-1}$]",
            ylim=[0, 1.2])

# Prepare the data inside the cavity
photon_traj_incav_thermostat_xyz = "../qm_solvent_free/incav_eq_weak_bath/E0_1e-4_rerun/simu_1.xc.xyz"
photon_traj_incav_qmmm_out = "../qmmm/incav_eq_weak_bath/E0_1e-4_rerun/simu_1.out"

qx_thermostat, qy_thermostat = get_polariton_traj_from_xyz(photon_traj_incav_thermostat_xyz)
freq_q_thermostat, sp_q_thermostat = get_ph_spectrum(qx_thermostat, qy_thermostat, method="fft")
freq_q_thermostat_pade, sp_q_thermostat_pade = get_ph_spectrum(qx_thermostat, qy_thermostat, method="pade")

qx_qmmm, qy_qmmm = get_polariton_traj_from_out(photon_traj_incav_qmmm_out)
freq_q_qmmm, sp_q_qmmm = get_ph_spectrum(qx_qmmm, qy_qmmm, method="fft")
freq_q_qmmm_pade, sp_q_qmmm_pade = get_ph_spectrum(qx_qmmm, qy_qmmm, method="pade")

xs = [freq_q_thermostat_pade, freq_q_thermostat]
ys = [sp_q_thermostat_pade, sp_q_thermostat]

clp.plotone(xs, ys, axes[0,1], colors=colors, labels=labels, lw_lst=[1, 1, 0.5], xlim=[1705, 2320],
            alphaspacing=0.1,
            ylim=[0, 1.2], showlegend=False)

xs = [freq_q_qmmm_pade, freq_q_qmmm]
ys = [sp_q_qmmm_pade, sp_q_qmmm]

clp.plotone(xs, ys, axes[1,1], colors=colors, labels=labels, lw_lst=[1, 1, 0.5], xlim=[1705, 2320],
            showlegend=False, alphaspacing=0.1, xlabel="frequency [cm$^{-1}$]",
            ylim=[0, 1.2])


# add linear response CO frequency
side_to_axis_ratio = (978 + 971) / (1138)
axes[0,0].axvline(x=2010, c='b', linestyle="--", lw=1, ymax=1.0/1.2)
axes[0,0].axvline(x=2025, c='b', linestyle="--", lw=1, ymax=1.0/1.2/ side_to_axis_ratio)
axes[1,0].axvline(x=2010, c='b', linestyle="--", lw=1, ymax=1.0/1.2)
axes[1,0].axvline(x=2025, c='b', linestyle="--", lw=1, ymax=1.0/1.2/ side_to_axis_ratio)
# add cavity mode frequency
axes[0,1].axvline(x=2018, c=clp.red, linestyle="--", lw=1)
axes[1,1].axvline(x=2018, c=clp.red, linestyle="--", lw=1)

# add text showing out or inside the cavity
axes[0,0].text(0.15, 0.8, "vacuum\ncavity off", fontsize=12, transform=axes[0,0].transAxes)
axes[0,1].text(0.15, 0.8, "vacuum\ncavity on", fontsize=12, transform=axes[0,1].transAxes)
axes[1,0].text(0.15, 0.8, "$n$-dodecane\ncavity off", fontsize=12, transform=axes[1,0].transAxes)
axes[1,1].text(0.15, 0.8, "$n$-dodecane\ncavity on", fontsize=12, transform=axes[1,1].transAxes)


# finally, we insert some images
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

data_img = mpimg.imread('feco5_vib_2.tiff')
imagebox = OffsetImage(data_img, zoom=0.08)
ab = AnnotationBbox(imagebox, (2050.5, 0.5), frameon=False)
ab.set(zorder=-1)
axes[0,0].add_artist(ab)

data_img = mpimg.imread('feco5_vib_1.tiff')
imagebox = OffsetImage(data_img, zoom=0.08)
ab = AnnotationBbox(imagebox, (1970.5, 0.5), frameon=False)
ab.set(zorder=-1)
axes[0,0].add_artist(ab)

clp.adjust(tight_layout=True, savefile="spectrum.pdf")
