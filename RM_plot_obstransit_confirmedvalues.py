import numpy as np
# import emcee
# from tqdm import tqdm
# import corner
import starry

from otherscripts.ylm_rot import get_ylm_coeffs
print("Using `starry` version %s." % starry.__version__)
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams["savefig.dpi"] = 100
rcParams["figure.dpi"] = 100
matplotlib.rcParams['axes.formatter.useoffset'] = False

# # this reads the data from the example file for HD189733 (provided by Rod Loger)
# time, rv, err = np.loadtxt("../data/test_for_expectcalc/HD189733_rvs.txt", unpack=True, skiprows=1, delimiter=',')

## # the following line is slightly modified to read the weighted RV files produced by The Great Reduction or MARMOT
## frameid, time, rv, err = np.loadtxt("../data/51_Peg/RVs_ID2864_newweighted_marmot.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_63433_c/RVs_ID3185_newweighted_marmot.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HAT-P-2/RVs_ID3204_newweighted_marmot.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HAT-P-2_nocomb/RVs_ID3204_newweighted_marmot.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HAT-P-2_nocomb_fresh/RVs_ID3204_newweighted_marmot_ccf.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_149026/RVs_ID3207_newweighted_marmot.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_149026_fresh/RVs_ID3207_newweighted_marmot.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_149026_nobadords/RVs_ID3207_newweighted_marmot.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_149026_nocomb/RVs_ID3207_newweighted_marmot_ccf_Mar.txt", unpack=True, delimiter=' ')
frameid, time, rv, err = np.loadtxt("../data/HD_149026_nocomb_fresh/RVs_ID3207_newweighted_marmot_ccf.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_149026_newccf_nocomb/RVs_ID3207_newweighted_marmot_ccfopt_5kms.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_149026_popultest/RVs_ID3207_newweighted_marmot_ccfopt.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_189733_nocomb/RVs_ID3155_newweighted_marmot.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_189733_nocomb_fresh/RVs_ID3155_newweighted_marmot_ccf.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_189733_newccf/RVs_ID3155_weighted_marmot_ccffit.txt", unpack=True, delimiter=' ')
idx = np.argsort(time)
time = time[idx]
time_finesampl = np.linspace(time[0], time[-1], num=10000)
## time used as zero point on the plot axis
## tref = 2459162.25
tref = 2454456.78760  # value for HD149026 from Albrecht et al. 2012
Per = 2.87588874  # value for HD149026 from Zhang et al. 2018
## tref = 2455288.84910  # value for HAT-P-2 from Bonomo et al. 2017
## tref = 2456021.256  # value for 51 Peg from Martins et al. 2015
## tref = 2454279.436714  # value for HD189733 from Casasayas et al.

# frameid, time, rv, err = np.loadtxt("../data/HD189733/RVs_ID3155_weighted.txt", unpack=True, delimiter=' ')
# idx = np.argsort(time)
# time = time[idx]
# time_finesampl = np.linspace(time[0], time[-1], num=50)
# # time used as zero point on the plot axis
# # tref = 2459162.25
# tref = 2454279.0  # example value for HD189377 by Rod Luger

rv = rv[idx]
err = err[idx]


class Normal(object):
    """A normal distribution."""

    def __init__(self, mean, sigma):
        self.mean = mean
        self.sigma = sigma

    def sample(self):
        return self.mean + self.sigma * np.random.randn()

    def evaluate(self, x):
        return -0.5 * (x - self.mean) ** 2 / self.sigma ** 2


class Uniform(object):
    """A uniform distribution."""

    def __init__(self, low, high):
        self.low = low
        self.high = high

    def sample(self):
        return np.random.uniform(self.low, self.high)

    def evaluate(self, x):
        if (x < self.low) or (x > self.high):
            return -np.inf
        else:
            return 0


class Sine(object):
    """
    A sine distribution.

    This is an uninformative distribution for the
    inclination of the stellar rotation axis.
    """

    def __init__(self):
        pass

    def sample(self):
        x = np.random.random()
        y = np.arccos(1 - x) * 180 / np.pi
        z = np.random.random()
        if z < 0.5:
            return 180 - y
        else:
            return y

    def evaluate(self, x):
        if x < 0 or x > 180:
            return -np.inf
        else:
            return np.log10(np.sin(x * np.pi / 180))


class Prior:
    """A class containing all the information on our priors."""

    def __init__(self):
        self.params = ['veq',  # equatorial rotational velocity of the star
                       'obl',  # obliquity angle lambda
                       'inc',  # inclination angle of the stellar rotation axis with respect to the line of sight
                       'alpha',  # chemical composition parameter??
                       'q1',  # limb darkening parameter??
                       'q2',  # limb darkening parameter??
                       'baseline',  # RV of the baseline, probably the regular RV at mid-transit time or so
                       'K',  # radial velocity semi-amplitude
                       'b_t0',  # time of periastron
                       'b_per',  # period in days
                       'b_inc',  # inclination of the system in degree
                       'b_r',  # radius of the planet in stellar radii
                       'b_a',  # semi-major axis of the orbit in stellar radii
                       'ln_err']

        self.names = [r'$v_\mathrm{eq}$',
                      r'$\lambda$',
                      r'$i$',
                      r'$\alpha$',
                      r'$q_1$',
                      r'$q_2$',
                      r'$v_\mathrm{0}$',
                      r'$K$',
                      r'$t_\mathrm{0,p}$',
                      r'$P_\mathrm{p}$',
                      r'$i_\mathrm{p}$',
                      r'$r_\mathrm{p}$',
                      r'$a$',
                      r'$\ln\sigma$']

        # # # example value for HD189377 by Rod Luger:
        # # self.baseline_guess = Normal(-9950.0, 3.0)
        # self.baseline_guess = Normal(-4.0, 4.0)
        # self.K_guess = Normal(200.56, 0.88)
        # self.b_t0_guess = Normal(2454279.436714, 0.000015)
        # self.b_per_guess = Normal(2.21857567, 0.00000015)
        # self.b_inc_guess = Normal(85.710, 0.024)
        # self.b_r_guess = Normal(0.15667, 0.00012)
        # self.b_a_guess = Normal(8.863, 0.020)

        # # The actual priors we use.
        # # Some of these are taken from TESS data (ExoFOP/SG2 table)
        # # Others are uniform / uninformative priors
        # # self.obl = Uniform(-90, 90)
        # self.inc = Normal(92, 4)
        # self.q1 = Uniform(0, 1)
        # self.q2 = Uniform(0, 1)
        # self.baseline = Normal(-1, 0.2)  # actual data used (broad range of RV results)
        # self.ln_err = Uniform(-3, 3)

        self.q1 = Uniform(0, 1)
        self.q2 = Uniform(0, 1)
        self.ln_err = Uniform(-3, 3)

        # # # values for 51 Peg from Martins et al. 2015 and Butler et al. 2006:
        # self.baseline = Normal(-50, 2)  # actual data used (broad range of RV results)
        # self.inc = Normal(80.0, 10.0)
        # self.obl = Normal(0.2, 12.5)
        # self.veq = Normal(2600, 500)
        # self.alpha = Normal(0.0, 0.5)
        # self.K = Normal(55.9, 0.7)  # actual data used
        # self.b_t0 = Normal(2456021.936, 0.0005)  # actual data used
        # self.b_per = Normal(4.230785, 0.000036)  # actual data used
        # self.b_inc = Normal(80.0, 10.0)  # actual data used
        # self.b_r = Normal(0.06891, 0.0009)  # actual data used
        # self.b_a = Normal(10.28, 0.19)  # actual data used

        # # # values for HAT-P-2 from Loeillet et al.2008 (L08) or Albrecht et al. 2012 (A12) or Ment2018 (M18) or Bonomo2017 (B17):
        # self.baseline = Normal(-12700, 30)  # L08: -19855.1 5.8; FOC data used (broad range of RV results)
        # self.inc = Normal(90.0, 12.5)  # unknown, P_rot missing for estimate, value = inc of planet, error = error of lambda; not optimal...
        # self.obl = Normal(0.2, 12.5)  # L08, upper error only 12.2 but made symmetric
        # self.veq = Normal(22900, 1200)  # L08, upper error only 1100 but made symmetric
        # self.alpha = Normal(0.0, 0.5)  # Fe/H 0.11 0.10 (L08)
        # self.K = Normal(953.3, 3.6)  # M18
        # self.b_t0 = Normal(2455288.84910, 0.00037)  # B17
        # self.b_per = Normal(5.6334754, 0.0000026)  # B17
        # self.b_inc = Normal(90.0, 0.93)  # L08, upper error only 0.85 but made symmetric
        # self.b_r = Normal(0.06891, 0.0009)  # L08 used, A12: 0.0723 0.0006
        # self.b_a = Normal(10.28, 0.19)  # L08, upper error only 0.12 but made symmetric

        # # values for HD149026 from Albrecht et al. 2012 (A12), Carter et al. 2009 (C09) or Wolf et al. 2007 (W07) or Zhang2018 (Z18) or Ment2018 (M18):
        self.baseline = Normal(-14514, 5)  # actual data used (broad range of RV results) FOC
        self.inc = Normal(90, 0.1)  # unknown, P_rot missing for estimate; not optimal...
        self.obl = Normal(12, 7)  # A12
        self.veq = Normal(7700, 800)  # A12
        self.alpha = Normal(0.0, 0.5)  # Fe/H 0.36 0.08; used before: 0.0 0.5
        self.K = Normal(39.22, 0.68)  # M18
        self.b_t0 = Normal(2454456.78760, 0.00016)  # Z18
        self.b_per = Normal(2.87588874, 0.00000059)  # Z18
        self.b_inc = Normal(84.55, 0.81)  # C09, upper error only 0.35 but made symmetric
        self.b_r = Normal(0.0507, 0.0009)  # A12
        self.b_a = Normal(6.01, 0.23)  # C09, upper error only 0.17 but made symmetric

        # # # example value for HD189733 from Casasayas-Barris et al. 2017 (CB17):
        # self.baseline = Normal(22275, 5)  # actual data used (broad range of RV results) FOC
        # self.inc = Normal(92, 12)  # CB17; stellar rotation axis (lower error -4, but bigger one adopted to be symmetric)
        # self.obl = Normal(-0.31, 0.17)  # CB17
        # self.veq = Normal(3500, 1000)  # CB17
        # self.alpha = Normal(0.0, 0.5)  # some estimate by me
        # self.K = Normal(205.0, 6.0)  # CB17
        # self.b_t0 = Normal(2454279.436714, 0.000015)  # CB17
        # self.b_per = Normal(2.21857567, 0.00000015)  # CB17
        # self.b_inc = Normal(85.7100, 0.0023)  # CB17
        # self.b_r = Normal(0.1513, 0.0072)  # CB17; but computed from R_P and R_S by me
        # self.b_a = Normal(8.84, 0.27)  # CB17

        # Distributions for the initial MCMC step
        self.veq_guess = Normal(4500, 100)
        self.obl_guess = Normal(-0.4, 0.3)
        self.inc_guess = Normal(100, 2.0)
        self.alpha_guess = Uniform(-0.5, 0.5)
        self.q1_guess = Uniform(0, 1)
        self.q2_guess = Uniform(0, 1)
        self.baseline_guess = Normal(-50.0, 30.0)
        self.K_guess = Normal(415, 13)
        self.b_t0_guess = Normal(2458268.455, 0.002)
        self.b_per_guess = Normal(5.55149, 0.00001)
        self.b_inc_guess = Normal(87.6, 1.0)
        self.b_r_guess = Normal(0.092, 0.003)
        self.b_a_guess = Normal(9.5, 0.2)

        self.ln_err_guess = Normal(0, 0.1)

    def evaluate(self, x):
        """Evaluate the log prior."""
        return np.sum([getattr(self, p).evaluate(x[i]) for i, p in enumerate(self.params)])

    def sample(self):
        """Sample from the prior distribution."""
        return [getattr(self, p).sample() for i, p in enumerate(self.params)]

    def guess(self):
        """Sample from the `guess' distribution."""
        return [getattr(self, p + "_guess").sample() for i, p in enumerate(self.params)]


# Instantiate the prior
prior = Prior()

# Instantiate our `starry` system
star = starry.kepler.Primary(5)
planet = starry.kepler.Secondary(0)
system = starry.kepler.System(star, planet)
map_unif = starry.Map(2)

def compute(x):
    """Compute the RV model given a parameter vector `x`."""
    # Get our params
    veq, obl, inc, alpha, q1, q2, baseline, K, b_t0, b_per, b_inc, b_r, b_a, ln_err = x

    # Planet params
    planet.tref = b_t0
    planet.porb = b_per
    planet.inc = b_inc
    planet.r = b_r
    planet.a = b_a

    # Stellar brightness-weighted velocity profile
    star.reset()
    star[:3, :] = get_ylm_coeffs(veq=veq, inc=inc, obl=obl, alpha=np.abs(alpha))
    star[0, 0] = 1
    sqrtq1 = np.sqrt(q1)
    u1 = 2 * sqrtq1 * q2
    u2 = sqrtq1 * (1 - 2 * q2)
    star[1] = u1
    star[2] = u2

    # Compute the integral of the brightness-weighted velocity field.
    # As we explain in `DifferentialRotationWithSphericalHarmonics.ipynb`,
    # what we're actually computing here is the integral of (Iv + I)
    system.compute(time_finesampl)
    intIv_plus_I = star.lightcurve

    # Compute the integral of the brightness, I; this is the
    # RM effect normalization.
    map_unif[0, 0] = 1
    map_unif[1] = u1
    map_unif[2] = u2
    intI = map_unif.flux(xo=planet.X, yo=planet.Y, ro=planet.r)

    # The RM effect is the integral of Iv divided by the integral of I
    # Note that we must subtract out the I term from the numerator
    model = (intIv_plus_I - intI) / intI

    # Add a baseline
    model += baseline - K * np.sin(2 * np.pi / b_per * (time_finesampl - b_t0))

    return model

# calculate some quality check numbers
error_median = np.median(err)

fig = plt.figure(figsize=(14, 6))
plt.errorbar(time - tref, rv, yerr=err, fmt='o', label='median error bar: {}'.format(error_median))
plt.xlabel("Time [BJD - %d]" % tref, fontsize=16)
plt.ylabel("RV [m / s]", fontsize=16)
plt.legend()

for i in range(100):
    model = compute(prior.sample())
    plt.plot(time_finesampl - tref, model, 'C1', alpha=0.1)

############################################################
## ## the following part can be used to do phase-folding and binning of the data
# time_fold = time%Per
# time_finesampl_fold = time_finesampl%Per
# nbin = 40

# n, _ = np.histogram(time_fold - tref, bins=nbin)
# sy, _ = np.histogram(time_fold - tref, bins=nbin, weights=rv)
# sy2, _ = np.histogram(time_fold - tref, bins=nbin, weights=rv*rv)
# mean = sy / n
# std = np.sqrt(sy2/n - mean*mean)

# fig = plt.figure(figsize=(14, 6))
# plt.errorbar(time_fold - tref, rv, yerr=err, fmt='o', label='median error bar: {}'.format(error_median))
# plt.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, fmt='b-')
# plt.xlabel("Time [BJD - %d]" % tref, fontsize=16)
# plt.ylabel("RV [m / s]", fontsize=16)
# plt.legend()

# for i in range(100):
#     model = compute(prior.sample())
#     plt.plot(time_finesampl_fold - tref, model, 'C1', alpha=0.1)
############################################################

plt.show()


# for indxind, toi in enumerate(TOI_ID):
#     print(indxind, toi, V_rot[indxind])
#
#     # compute the right plotting window
#     dur = float(duration[indxind]) / 24  # convert hours to days
#     tref = float(T0[indxind])
#     time_finesampl = np.linspace(tref - 1.5 * dur, tref + 1.5 * dur, num=50)
#
#     # prepare all values to be used for the model
#     v_rot_val = float(V_rot[indxind])
#     v_rot_val_err = 0.1 * v_rot_val

# ##### OLD SNIPPETS #####


# # # example value for HD189377 by Rod Luger:
# # self.baseline = Uniform(-10000, -9900)  # actual data used (broad range of RV results)
# self.baseline = Normal(-1, 20)  # actual data used (broad range of RV results)
# self.K = Uniform(180, 220)  # actual data used
# self.b_t0 = Normal(2454279.436714, 0.000015)  # actual data used
# self.b_per = Normal(2.21857567, 0.00000015)  # actual data used
# self.b_inc = Normal(85.710, 0.024)  # actual data used
# self.b_r = Normal(0.15667, 0.00012)  # actual data used
# self.b_a = Normal(8.863, 0.020)  # actual data used

# # # example value for HD189377 by Rod Luger:
# # self.baseline_guess = Normal(-9950.0, 3.0)
# self.baseline_guess = Normal(-4.0, 4.0)
# self.K_guess = Normal(200.56, 0.88)
# self.b_t0_guess = Normal(2454279.436714, 0.000015)
# self.b_per_guess = Normal(2.21857567, 0.00000015)
# self.b_inc_guess = Normal(85.710, 0.024)
# self.b_r_guess = Normal(0.15667, 0.00012)
# self.b_a_guess = Normal(8.863, 0.020)
