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

# # for TOI1148 / TIC349827430 / KELT-24 use this transit times:
# duration = 4.315966 / 24  # convert hours to days
# tref = 2458684.82189

# # for TOI1153 / TIC154840461 / HD106018 use this transit times:
# duration = 0.978 / 24 * 2  # convert hours to days
# tref = 2458685.579346

# # for TOI1271 / TIC286923464 / HD118203 use this transit times:
# duration = 5.6457 / 24  # convert hours to days
# tref = 2458712.662354

# # for TOI1415 / TIC148782377 / HIP71409 use this transit times:
# duration = 5.163571 / 24  # convert hours to days
# tref = 2458933.177794

# # for TOI1447 / HD117173 / TIC298073824 use this transit times:
# duration = 4.474106 / 24  # convert hours to days
# tref = 2458684.296509

# # for TOI1665 / TIC354006740 / HD39315 use this transit times:
# duration = 3.086 / 24  # convert hours to days
# tref = 2458819.00417

# # for TOI1718 / HD58727 / TIC257241363 use this transit times:
# duration = 2.931427 / 24  # convert hours to days
# tref = 2458848.027133

# # for TOI1719 / HD76854 / TIC293617835 use this transit times:
# duration = 1.321287 / 24  # convert hours to days
# tref = 2458844.193788

# # for TOI1726 / TIC130181866 / HD63433 use this transit times:
# duration = 3.227344 / 24  # convert hours to days
# tref = 2458845.372818

# # for TOI1773 / 55Cnc / TIC332064670 use this transit times:
# duration = 1.459 / 24  # convert hours to days
# tref = 2458872.16257

# # for TOI1776 / HD95072 / TIC21535395 use this transit times:
# duration = 1.49854 / 24  # convert hours to days
# tref = 2458871.488566

# # for TOI1778 / TIC39699648 / HD77946 use this transit times:
# duration = 2.533431 / 24  # convert hours to days
# tref = 2458876.025023

# # for TOI1796 / TIC138819293 / GJ436 use this transit times:
# duration = 0.954893 / 24  # convert hours to days
# tref = 2458899.671523

# # for TOI1799 / TIC8967242 / HD96735 use this transit times:
# duration = 2.510763 / 24  # convert hours to days
# tref = 2458904.815773

# # for TOI1821 / TIC82308728 / HD97658 use this transit times:
# duration = 2.268 / 24  # convert hours to days
# tref = 2458904.92627

# # for TOI1831 / TIC27194429 / HD116206 use this transit times:
# duration = 0.848036 / 24  # convert hours to days
# tref = 2458928.117932

####### Eike's objects (START) #######
# # for TOI1181 / TIC229510866 use this transit times:
# toi_id = 'TOI1181'
# duration = 4.093308 / 24  # convert hours to days
# tref = 2458957.821395

# # for TOI1516 / TIC376637093 use this transit times:
# toi_id = 'TOI1516'
# duration = 2.829 / 24  # convert hours to days
# tref = 2458765.325886056

# # for TOI2046 / TIC468574941 use this transit times:
# toi_id = 'TOI2046'
# duration = 2.806 / 24  # convert hours to days
# tref = 2458792.3958005020

# for TOI1408 / TIC364186197 use this transit times:
toi_id = 'TOI1408'
duration = 1.566191 / 24  # convert hours to days
tref = 2458740.860909
####### Eike's objects (END) #######

time_finesampl = np.linspace(tref - 1.5 * duration, tref + 1.5 * duration, num=50)


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

        # The actual priors we use.
        # Some of these are taken from TESS data (ExoFOP/SG2 table)
        # Others are uniform / uninformative priors
        self.obl = Uniform(-1.0, 1.0)
        self.inc = Normal(90.0, 2.0)  # Sine()
        self.q1 = Uniform(0, 1)
        self.q2 = Uniform(0, 1)
        self.baseline = Normal(0.0, 0.1)  # actual data used (broad range of RV results)
        self.ln_err = Uniform(-3, 3)

        ########### These are the objects for Eike (START) ##############
        # # This is for TOI1181 / TIC229510866
        # self.veq = Normal(12000, 2000)
        # self.alpha = Normal(0.0, 0.1)
        # self.K = Normal(143.20, 10.0)  # actual data used
        # self.b_t0 = Normal(2458957.821395, 0.000222)  # actual data used
        # self.b_per = Normal(2.103193, 0.000011)  # actual data used
        # self.b_inc = Normal(90.0, 2.0)  # actual data used
        # self.b_r = Normal(0.07603077606858845, 0.006829311044773664)  # actual data used
        # self.b_a = Normal(3.77, 0.3)  # actual data used

        # # This is for TOI1516 / TIC376637093
        # self.veq = Normal(10380, 1000)
        # self.alpha = Normal(-0.06, 0.06)
        # self.K = Normal(132.12, 10.0)  # actual data used
        # self.b_t0 = Normal(2458765.32531, 0.00019)  # actual data used
        # self.b_per = Normal(2.05603, 0.00001)  # actual data used
        # self.b_inc = Normal(90.0, 1.0)  # actual data used
        # self.b_r = Normal(0.12107871816865484, 0.010425178804747654)  # actual data used
        # self.b_a = Normal(6.02, 0.6)  # actual data used

        # # This is for TOI2046 / TIC468574941
        # self.veq = Normal(7440, 1000)
        # self.alpha = Normal(0.16, 0.06)
        # self.K = Normal(162.34, 15.0)  # actual data used
        # self.b_t0 = Normal(2458792.39519, 0.00038)  # actual data used
        # self.b_per = Normal(1.4972, 0.00001)  # actual data used
        # self.b_inc = Normal(90.0, 1.0)  # actual data used
        # self.b_r = Normal(0.13534838407376606, 0.017880088859501816)  # actual data used
        # self.b_a = Normal(4.63, 0.45)  # actual data used

        # This is for TOI1408 / TIC364186197
        self.veq = Normal(12000, 1000)
        self.alpha = Normal(0.0, 0.01)
        self.K = Normal(29.38, 3.0)  # actual data used
        self.b_t0 = Normal(2458740.860909, 0.000354)  # actual data used
        self.b_per = Normal(4.424575, 0.00001)  # actual data used
        self.b_inc = Normal(90.0, 1.0)  # actual data used
        self.b_r = Normal(0.0724, 0.0141)  # 0.0281 actual data used 0.0719173387667067 0.02944764070789319
        self.b_a = Normal(8.86, 0.90)  # actual data used
        ########### These are the objects for Eike (END) ##############

        # # This is for TOI1148 / TIC349827430 / KELT-24
        # self.veq = Normal(19490, 1000)
        # self.alpha = Normal(0.0, 0.5)
        # self.K = Normal(462.0, 20.0)  # actual data used
        # self.b_t0 = Normal(2458684.82189, 0.00032)  # actual data used
        # self.b_per = Normal(5.55111, 0.00014)  # actual data used
        # self.b_inc = Normal(89.17, 0.75)  # actual data used
        # self.b_r = Normal(0.0866615898349809, 0.010195081567143037)  # actual data used
        # self.b_a = Normal(9.95, 0.18)  # actual data used

        # # This is for TOI1153 / TIC154840461 / HD106018
        # self.veq = Normal(10000, 500)
        # self.alpha = Normal(0.0, 0.5)
        # self.K = Normal(13.0, 1.3)  # actual data used
        # self.b_t0 = Normal(2458685.579346, 0.00119)  # actual data used
        # self.b_per = Normal(6.03614, 0.0006)  # actual data used
        # self.b_inc = Normal(84.34275272438427, 1.5)  # actual data used
        # self.b_r = Normal(0.04519303184484871, 0.00546114521038452)  # actual data used
        # self.b_a = Normal(9.84, 1.0)  # actual data used

        # # This is for TOI1271 / TIC286923464 / HD118203
        # self.veq = Normal(4700, 500)
        # self.alpha = Normal(0.211367, 0.0109545)
        # self.K = Normal(37.42, 2.0)  # actual data used
        # self.b_t0 = Normal(2458712.662354, 0.00093)  # actual data used
        # self.b_per = Normal(6.134842, 0.00056)  # actual data used
        # self.b_inc = Normal(88.75, 1.0)  # actual data used
        # self.b_r = Normal(0.054883422602936396, 0.003249664695925899)  # actual data used
        # self.b_a = Normal(7.24, 0.3)  # actual data used

        # # This is for TOI1415 / TIC148782377 / HIP71409
        # self.veq = Normal(22000, 1000)
        # self.alpha = Normal(0.0, 0.5)
        # self.K = Normal(4.59, 0.5)  # actual data used
        # self.b_t0 = Normal(2458933.177794, 0.001553)  # actual data used
        # self.b_per = Normal(14.419358, 0.002257)  # actual data used
        # self.b_inc = Normal(89.9721865137071, 1.5)  # actual data used
        # self.b_r = Normal(0.031092504234019224, 0.004602984713770016)  # actual data used
        # self.b_a = Normal(20.6, 3.5)  # actual data used

        # # This is for TOI1447 / HD117173 / TIC298073824
        # self.veq = Normal(85000, 1000)
        # self.alpha = Normal(0.0, 0.5)
        # self.K = Normal(8.21, 0.8)  # actual data used
        # self.b_t0 = Normal(2458684.296509, 0.003088)  # actual data used
        # self.b_per = Normal(1.229303, 0.000136)  # actual data used
        # self.b_inc = Normal(70.83084428107044, 3.541542214053522)  # actual data used
        # self.b_r = Normal(0.009861864156675873, 0.002757576103228098)  # actual data used
        # self.b_a = Normal(1.34, 0.4)  # actual data used

        # # This is for TOI1665 / TIC354006740 / HD39315
        # self.veq = Normal(8000, 800)
        # self.alpha = Normal(0.0, 0.5)
        # self.K = Normal(47.07, 3.0)  # actual data used
        # self.b_t0 = Normal(2458819.00417, 0.0022)  # actual data used
        # self.b_per = Normal(1.76406, 0.00029)  # actual data used
        # self.b_inc = Normal(88.0, 13.2)  # actual data used
        # self.b_r = Normal(0.055893409423545495, 0.021776639249061504)  # actual data used
        # self.b_a = Normal(3.26, 0.48899999999999993)  # actual data used

        # # This is for TOI1718 / HD58727 / TIC257241363
        # self.veq = Normal(4000, 400)
        # self.alpha = Normal(0.2, 0.1)
        # self.K = Normal(5.84, 0.5)  # actual data used
        # self.b_t0 = Normal(2458848.027133, 0.000695)  # actual data used
        # self.b_per = Normal(5.586947, 0.000372)  # actual data used
        # self.b_inc = Normal(89.642, 1.0)  # actual data used
        # self.b_r = Normal(0.0394, 0.0070)  # actual data used
        # self.b_a = Normal(14.4, 2.0)  # actual data used

        # # This is for TOI1719 / HD76854 / TIC293617835
        # self.veq = Normal(8000, 1000)
        # self.alpha = Normal(0.0, 0.5)
        # self.K = Normal(19.84, 2.0)  # actual data used
        # self.b_t0 = Normal(2458844.193788, 0.001695)  # actual data used
        # self.b_per = Normal(2.756792, 0.000348)  # actual data used
        # self.b_inc = Normal(89.4843310226913, 1.5)  # actual data used
        # self.b_r = Normal(0.017235364868821376, 0.003)  # actual data used
        # self.b_a = Normal(10.0, 2.5)  # actual data used

        # # This is for TOI1726 / TIC130181866 / HD63433
        # self.veq = Normal(7500, 800)
        # self.alpha = Normal(0.0169559, 0.0168932)
        # self.K = Normal(1.75, 0.5)  # actual data used
        # self.b_t0 = Normal(2458845.372818, 0.001305)  # actual data used
        # self.b_per = Normal(7.108153, 0.000651)  # actual data used
        # self.b_inc = Normal(89.38, 0.43)  # actual data used
        # self.b_r = Normal(0.02191711724293351, 0.004)  # actual data used
        # self.b_a = Normal(17.15, 0.3)  # actual data used

        # # This is for TOI1773 / 55Cnc / TIC332064670
        # self.veq = Normal(300, 50)
        # self.alpha = Normal(0.35, 0.1)
        # self.K = Normal(2.98, 0.5)  # actual data used
        # self.b_t0 = Normal(2458872.16257, 0.00186)  # actual data used
        # self.b_per = Normal(0.73665, 0.00008)  # actual data used
        # self.b_inc = Normal(83.3, 0.9)  # actual data used
        # self.b_r = Normal(0.017172410069594173, 0.002521391758569257)  # actual data used
        # self.b_a = Normal(3.05, 0.3)  # actual data used

        # # This is for TOI1776 / HD95072 / TIC21535395
        # self.veq = Normal(1000, 500)
        # self.alpha = Normal(0.0, 0.5)
        # self.K = Normal(1.14, 0.5)  # actual data used
        # self.b_t0 = Normal(2458871.488566, 0.003068)  # actual data used
        # self.b_per = Normal(2.800991, 0.00059)  # actual data used
        # self.b_inc = Normal(87.39, 4.366944108369561)  # actual data used
        # self.b_r = Normal(0.0136, 0.003)  # actual data used
        # self.b_a = Normal(11.2, 2.5)  # actual data used

        # # This is for TOI1778 / TIC39699648 / HD77946
        # self.veq = Normal(4000, 500)
        # self.alpha = Normal(0.0, 0.5)
        # self.K = Normal(2.67, 0.5)  # actual data used
        # self.b_t0 = Normal(2458876.025023, 0.002952)  # actual data used
        # self.b_per = Normal(6.515951, 0.001553)  # actual data used
        # self.b_inc = Normal(88.93041666122694, 1.5)  # actual data used
        # self.b_r = Normal(0.019515629487932274, 0.002)  # actual data used
        # self.b_a = Normal(15.0, 4.5)  # actual data used

        # # This is for TOI1796 / TIC138819293 / GJ436
        # self.veq = Normal(300, 50)
        # self.alpha = Normal(0.0, 0.5)
        # self.K = Normal(11.27, 1.0)  # actual data used
        # self.b_t0 = Normal(2458899.671523, 0.000169)  # actual data used
        # self.b_per = Normal(2.643973, 0.000027)  # actual data used
        # self.b_inc = Normal(86.44, 0.17)  # actual data used
        # self.b_r = Normal(0.08111400034755328, 0.004933055565378086)  # actual data used
        # self.b_a = Normal(13.9, 0.5)  # actual data used

        # # This is for TOI1799 / TIC8967242 / HD96735
        # self.veq = Normal(2000, 200)
        # self.alpha = Normal(0.0, 0.5)
        # self.K = Normal(1.09, 0.3)  # actual data used
        # self.b_t0 = Normal(2458904.815773, 0.003937)  # actual data used
        # self.b_per = Normal(7.093853, 0.002086)  # actual data used
        # self.b_inc = Normal(89.57781675164001, 1.5)  # actual data used
        # self.b_r = Normal(0.015594179328583412, 0.002)  # actual data used
        # self.b_a = Normal(19.0, 3.0)  # actual data used

        # # This is for TOI1821 / TIC82308728 / HD97658
        # self.veq = Normal(500, 50)
        # self.alpha = Normal(-0.24, 0.1)
        # self.K = Normal(2.9, 0.3)  # actual data used
        # self.b_t0 = Normal(2458904.92627, 0.008)  # actual data used
        # self.b_per = Normal(9.497, 0.00275)  # actual data used
        # self.b_inc = Normal(89.45, 0.4)  # actual data used
        # self.b_r = Normal(0.025579599431759848, 0.0038776857592502613)  # actual data used
        # self.b_a = Normal(23.9, 1.33)  # actual data used

        # # This is for TOI1831 / TIC27194429 / HD116206
        # self.veq = Normal(35500, 1000)
        # self.alpha = Normal(0.0, 0.5)
        # self.K = Normal(11.87, 1.0)  # actual data used
        # self.b_t0 = Normal(2458928.117932, 0.000601)  # actual data used
        # self.b_per = Normal(0.555293, 0.000023)  # actual data used
        # self.b_inc = Normal(88.0, 1.5)  # actual data used
        # self.b_r = Normal(0.018974231259244047, 0.003)  # actual data used
        # self.b_a = Normal(1.43, 0.1)  # actual data used

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


fig = plt.figure(figsize=(14, 6))
# plt.errorbar(time - tref, rv, yerr=err, fmt='o')
plt.xlabel("Time [BJD - %d]" % tref, fontsize=16)
plt.ylabel("RV [m / s]", fontsize=16)
plt.title('{}'.format(toi_id))

for i in range(100):
    model = compute(prior.sample())
    plt.plot(time_finesampl - tref, model, 'C1', alpha=0.1)

plt.savefig('../plots/{}_RM_expectation.png'.format(toi_id))
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

# this reads the data from the example file for HD189733 (provided by Rod Loger)
# time, rv, err = np.loadtxt("../data/test_for_expectcalc/HD189733_rvs.txt", unpack=True, skiprows=1, delimiter=',')
# this line is slightly modified to read the weighted RV files produced by The Great Reduction
# frameid, time, rv, err = np.loadtxt("../data/HD189733/RVs_ID3155_weighted.txt", unpack=True, delimiter=' ')
# idx = np.argsort(time)
# time = time[idx]
# time_finesampl = np.linspace(time[0], time[-1], num=50)
# # time used as zero point on the plot axis
# # tref = 2459162.25
# tref = 2454279.0  # example value for HD189377 by Rod Luger
#
# rv = rv[idx]
# err = err[idx]

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
