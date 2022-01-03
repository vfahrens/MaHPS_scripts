import numpy as np
# import emcee
# from tqdm import tqdm
# import corner
import starry
from astropy.time import Time, TimeDelta

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
## frameid, time, rv, err = np.loadtxt("../data/HD_149026_nocomb_fresh/RVs_ID3207_newweighted_marmot_ccf.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_149026_newccf_nocomb/RVs_ID3207_newweighted_marmot_ccfopt_5kms.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_149026_popultest/RVs_ID3207_newweighted_marmot_ccfopt.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_189733_nocomb/RVs_ID3155_newweighted_marmot.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_189733_nocomb_fresh/RVs_ID3155_newweighted_marmot_ccf.txt", unpack=True, delimiter=' ')
## frameid, time, rv, err = np.loadtxt("../data/HD_189733_newccf/RVs_ID3155_weighted_marmot_ccffit.txt", unpack=True, delimiter=' ')
## idx = np.argsort(time)
## time = time[idx]
## time_finesampl = np.linspace(time[0], time[-1], num=10000)
## time used as zero point on the plot axis
## tref = 2459162.25

## tref = 2456021.256  # value for 51 Peg from Martins et al. 2015


#####################################################################################
## this is new and more compact, for plotting the individual nights
#####################################################################################

## ## HD189733 ## ##
obsdate = '20201106'  ## dates: 20201106 (transit); 20201109 (templates)
path_for_date = '../data/HD_189733_nocomb_fresh/RVs_ID3155_newweighted_marmot_ccf_' + obsdate + '.txt'
outfile_plot = '../plots/HD189733_confparam/HD189733_RV_binnedmodel_{}.pdf'.format(obsdate)

frameid, time, rv, err = np.loadtxt(path_for_date, unpack=True, delimiter=' ')
plt_title = 'HD189733 ({}-{}-{})'.format(obsdate[:4], obsdate[4:6], obsdate[6:]) 

tref = 2454279.436714  # value for HD189733 from Casasayas et al.
tdur = 1.84  # [hours] value for HD189733 from Addison et al. 2019
Per = 2.21857567  # value for HD189733 from Casasayas et al.
exptime_frames = 900 / 84600  # exposure time, converted to days
idx = np.argsort(time)
time = time[idx]

obs_phase_start = ((time[0]-tref)%Per)
obs_phase_end = ((time[-1]-tref)%Per)
N_orbit = 1 + ( (time[0] - obs_phase_start - tref)/Per )

plot_time_window = ( tdur / 24 ) * 1.75
plot_phase_start = tref + (N_orbit*Per) - plot_time_window
plot_phase_end = tref + (N_orbit*Per) + plot_time_window

if obsdate == '20201109':
    obs_time_med = np.median(time - tref)
    offset_phase = obs_time_med%Per
    plot_phase_start = tref + (N_orbit*Per) - plot_time_window - (Per - offset_phase)
    plot_phase_end = tref + (N_orbit*Per) + plot_time_window - (Per - offset_phase)

res_time_finesampl = 10000
time_finesampl = np.linspace(plot_phase_start, plot_phase_end, num=res_time_finesampl)

template_bary_mean = -22119.436
rv = rv + template_bary_mean  #np.median(rv)
baseline_use = 0
baseline_err_use = 3.7

baselines_individ = {'20201106':22280, '20201109':22136}
#baseline_use = baselines_individ[obsdate]
# ## ## HD189733 END ## ##


# ## ## HD149026 ## ##
# obsdate = '20210426'  ## dates: 20210305 (transit); 20210306 (templates); 20210308 (transit); 20210423 (transit); 20210426 (transit)
# path_for_date = '../data/HD_149026_nocomb_fresh/RVs_ID3207_newweighted_marmot_ccf_' + obsdate + '.txt'
# outfile_plot = '../plots/HD149026_confparam/HD149026_RV_binnedmodel_{}.pdf'.format(obsdate)

# frameid, time, rv, err = np.loadtxt(path_for_date, unpack=True, delimiter=' ')
# plt_title = 'HD149026 ({}-{}-{})'.format(obsdate[:4], obsdate[4:6], obsdate[6:]) 

# tref = 2454456.78760  # value for HD149026 from Albrecht et al. 2012
# tdur = 3.230  # [hours] value for HD149026 from Albrecht et al. 2012
# Per = 2.87588874  # value for HD149026 from Zhang et al. 2018
# exptime_frames = 1200 / 84600  # exposure time, converted to days
# idx = np.argsort(time)
# time = time[idx]

# obs_phase_start = ((time[0]-tref)%Per)
# obs_phase_end = ((time[-1]-tref)%Per)
# N_orbit = 1 + ( (time[0] - obs_phase_start - tref)/Per )

# plot_time_window = ( tdur / 24 ) * 1.5
# plot_phase_start = tref + (N_orbit*Per) - plot_time_window
# plot_phase_end = tref + (N_orbit*Per) + plot_time_window

# if obsdate == '20210306':
#     obs_time_med = np.median(time - tref)
#     offset_phase = obs_time_med%Per
#     plot_phase_start = tref + (N_orbit*Per) - plot_time_window - (Per - offset_phase)
#     plot_phase_end = tref + (N_orbit*Per) + plot_time_window - (Per - offset_phase)

# res_time_finesampl = 10000
# time_finesampl = np.linspace(plot_phase_start, plot_phase_end, num=res_time_finesampl)

# template_bary_mean = 14498.537
# rv = rv + template_bary_mean  #np.median(rv)
# baseline_use = 0
# baseline_err_dates = {'20210305':7.0, '20210306':7.6, '20210308':8.9, '20210423':6.3, '20210426':6.6}
# baseline_err_use = baseline_err_dates[obsdate]

# baselines_individ = {'20210305':-14485, '20210306':-14503, '20210308':-14504, '20210423':-14516, '20210426':-14494}
# #baseline_use = baselines_individ[obsdate]
# ## ## HD149026 END ## ##


# ## HAT-P-2 ## ##
# obsdate = '20210225'  ## dates: 20210224 (out of transit); 20210225 (transit); 20210228 (templates); 20210324 (templates)
# path_for_date = '../data/HAT-P-2_nocomb_fresh/RVs_ID3204_newweighted_marmot_ccf_' + obsdate + '.txt'
# outfile_plot = '../plots/HAT-P-2_confparam/HAT-P-2_RV_binnedmodel_{}.pdf'.format(obsdate)

# frameid, time, rv, err = np.loadtxt(path_for_date, unpack=True, delimiter=' ')
# plt_title = 'HAT-P-2 ({}-{}-{})'.format(obsdate[:4], obsdate[4:6], obsdate[6:]) 

# tref = 2455288.84910  # value for HAT-P-2 from Bonomo et al. 2017
# tdur = 4.2888  # [hours] value for HAT-P-2 from Pal et al. 2010
# Per = 5.6334754  # value for HAT-P-2 from Bonomo et al. 2017
# exptime_frames = 1800 / 84600  # exposure time, converted to days
# idx = np.argsort(time)
# time = time[idx]

# obs_phase_start = ((time[0]-tref)%Per)
# obs_phase_end = ((time[-1]-tref)%Per)
# N_orbit = 1 + ( (time[0] - obs_phase_start - tref)/Per )

# plot_time_window = ( tdur / 24 ) * 1.5
# plot_phase_start = tref + (N_orbit*Per) - plot_time_window
# plot_phase_end = tref + (N_orbit*Per) + plot_time_window

# if obsdate == '20210224':
#     obs_time_med = np.median(time - tref)
#     offset_phase = obs_time_med%Per
#     plot_phase_start = tref + (N_orbit*Per) - plot_time_window - (Per - offset_phase)
#     plot_phase_end = tref + (N_orbit*Per) + plot_time_window - (Per - offset_phase)

# res_time_finesampl = 10000
# time_finesampl = np.linspace(plot_phase_start, plot_phase_end, num=res_time_finesampl)

# template_bary_mean = 11792.095
# rv = rv + template_bary_mean  #np.median(rv)
# baseline_use = 0
# baseline_err_dates = {'20210224':28.9, '20210225':25.6}
# baseline_err_use = baseline_err_dates[obsdate]

# baselines_individ = {'20210224':-11859, '20210225':-12584, '20210228':-11704, '20210324':-11851}
# #baseline_use = baselines_individ[obsdate]
# ## HAT-P-2 END ## ##






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
        # self.baseline = Normal(baseline_use, baseline_err_use)  # L08: -19855.1 5.8; FOC data used (broad range of RV results)
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

        # # # values for HD149026 from Albrecht et al. 2012 (A12), Carter et al. 2009 (C09) or Wolf et al. 2007 (W07) or Zhang2018 (Z18) or Ment2018 (M18):
        # self.baseline = Normal(baseline_use, baseline_err_use)  # actual data used (broad range of RV results) FOC
        # self.inc = Normal(90, 0.1)  # unknown, P_rot missing for estimate; not optimal...
        # self.obl = Normal(12, 7)  # A12
        # self.veq = Normal(7700, 800)  # A12
        # self.alpha = Normal(0.0, 0.5)  # Fe/H 0.36 0.08; used before: 0.0 0.5
        # self.K = Normal(39.22, 0.68)  # M18
        # self.b_t0 = Normal(2454456.78760, 0.00016)  # Z18
        # self.b_per = Normal(2.87588874, 0.00000059)  # Z18
        # self.b_inc = Normal(84.55, 0.81)  # C09, upper error only 0.35 but made symmetric
        # self.b_r = Normal(0.0507, 0.0009)  # A12
        # self.b_a = Normal(6.01, 0.23)  # C09, upper error only 0.17 but made symmetric

        # # example value for HD189733 from Casasayas-Barris et al. 2017 (CB17):
        self.baseline = Normal(baseline_use, baseline_err_use)  # actual data used (broad range of RV results) FOC
        self.inc = Normal(92, 12)  # CB17; stellar rotation axis (lower error -4, but bigger one adopted to be symmetric)
        self.obl = Normal(-0.31, 0.17)  # CB17
        self.veq = Normal(3500, 1000)  # CB17
        self.alpha = Normal(0.0, 0.5)  # some estimate by me
        self.K = Normal(205.0, 6.0)  # CB17
        self.b_t0 = Normal(2454279.436714, 0.000015)  # CB17
        self.b_per = Normal(2.21857567, 0.00000015)  # CB17
        self.b_inc = Normal(85.7100, 0.0023)  # CB17
        self.b_r = Normal(0.1513, 0.0072)  # CB17; but computed from R_P and R_S by me
        self.b_a = Normal(8.84, 0.27)  # CB17

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

########### original plotting routine
# fig = plt.figure(figsize=(14, 6))
# ax = fig.add_subplot(111)
# plt.errorbar(time - tref, rv, yerr=err, fmt='o', label='median error bar: {:.1f} m/s'.format(error_median))
# plt.xlabel("Time [BJD - %d]" % tref, fontsize=18)
# plt.ylabel("RV [m / s]", fontsize=18)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# plt.legend(fontsize=18)
# ax.set_title(plt_title, fontsize=18)
# plt.tight_layout()

# for i in range(100):
#     model = compute(prior.sample())
#     plt.plot(time_finesampl - tref, model, 'C1', alpha=0.1)
########### original ploting routine END


########################################
### binning of the models
readtime = 87.5 / 84600  ## readout time of the CCD in days

## make the models
num_of_models = 100
num_iter_rvoffset = 100  ## make 100x100 models, needed to evaluate the RV zeropoint offset
num_of_models_large = num_of_models * num_iter_rvoffset
all_models_large = np.empty(shape=(num_of_models_large,res_time_finesampl))
for iijii in range(num_of_models_large):
    model = compute(prior.sample())
    all_models_large[iijii] = model
    if iijii%200 == 0:
        print('Made {} models.'.format(iijii))

all_models = all_models_large[:num_of_models]  ## only use the first 100 models for plotting and calculating stuff

## calculate the median model value for each time_finesampl coordinate
model_rvs_median = []
#model_rvs_error = []
for timcor in range(res_time_finesampl):
    mod_rvs_onetimcor = all_models[:, timcor]
    model_rvs_median.append(np.median(mod_rvs_onetimcor))
    #model_rvs_error.append(np.std(mod_rvs_onetimcor))

model_rvs_median = np.array(model_rvs_median)

## calculate the binning of the model to reflect that of the data
num_bins_before_firstdata = np.floor( (time[0] - plot_phase_start - 0.5*exptime_frames) / (exptime_frames + readtime) )
first_bin_start = time[0] - 0.5*exptime_frames - num_bins_before_firstdata*(exptime_frames + readtime)
interval_covered = plot_phase_end - first_bin_start
num_of_bins = np.floor( interval_covered / (exptime_frames + readtime) )
num_of_bins = int(num_of_bins)

bins = []  ## list of bin edges: for each bin the left edge and for last bin additional the right edge
#times_bins = []
for incm in range(num_of_bins):
    bin_left_edge = time[0] - 0.5*exptime_frames - (num_bins_before_firstdata - incm)*(exptime_frames + readtime) - tref
    #bin_expmidtime = time[0] - (num_bins_before_firstdata - incm)*(exptime_frames + readtime)
    bins.append(bin_left_edge)
    #times_bins.append(bin_expmidtime)
    
bin_right_edge_last = time[0] + 0.5*exptime_frames + readtime - (num_bins_before_firstdata - num_of_bins + 1)*(exptime_frames + readtime) - tref
bins.append(bin_right_edge_last)
#times_bins = np.array(times_bins)

## create the binned model
n, _ = np.histogram(time_finesampl - tref, bins=bins)
#sy, _ = np.histogram(time_finesampl - tref, bins=bins, weights=model_rvs_median)
#median_rv_per_bin = np.median(sy)
#mean = sy / n


## get all the model RVs in a given bin
time_finesampl_arr = np.array(time_finesampl -tref)
indics_per_bin = []
mean_rv_per_bin = []
for bininin in range(len(bins)-1):
    part_in_bin = np.logical_and(bins[bininin] <= time_finesampl_arr, time_finesampl_arr < bins[bininin+1])
    indics_in_bin = np.where(part_in_bin)[0]
    indics_per_bin.append(indics_in_bin)
## get the stddev of all RV values in that bin, without doing the median for each time coordinate
errors_using_all_rvs = []
for bonon in range(len(bins)-1):
    rvs_only_one_bin = all_models[:, indics_per_bin[bonon]]
    rvs_only_one_bin_flat = rvs_only_one_bin.reshape(-1)
    std_one_bin = np.std(rvs_only_one_bin_flat)
    errors_using_all_rvs.append(std_one_bin)

    ## also get the RVs per bin by computing the mean(!, because exposure is time integration) for all times in the bin
    rvs_in_bin_model_med = model_rvs_median[indics_per_bin[bonon]]
    rv_mean_one_bin = np.mean(rvs_in_bin_model_med)
    mean_rv_per_bin.append(rv_mean_one_bin)

    
# ## calculate the offset from zero
bin_mid_time = [ (bins[lotr] + bins[lotr+1])/2 for lotr in range(len(bins)-1) ]
time_expstart = time - tref - 0.5*exptime_frames
pos_of_bins_with_obsdata = np.logical_and( bins > (time_expstart[0] - 0.25*exptime_frames), bins < (time_expstart[-1] + 0.25*exptime_frames) )
indics_of_bins_with_obsdata = np.where(pos_of_bins_with_obsdata)[0]

print(indics_of_bins_with_obsdata)
if obsdate == '20201106':
    indics_of_bins_with_obsdata = indics_of_bins_with_obsdata[:-8]
if obsdate == '20210308':
    indics_of_bins_with_obsdata = indics_of_bins_with_obsdata[:-7]
if obsdate == '20210423':
    indics_of_bins_with_obsdata = indics_of_bins_with_obsdata[:-5]
print(indics_of_bins_with_obsdata)

timecord_of_bins_with_obsdata = bin_mid_time[indics_of_bins_with_obsdata[0]:indics_of_bins_with_obsdata[-1]+1]
modelrvs_of_bins_with_obsdata = mean_rv_per_bin[indics_of_bins_with_obsdata[0]:indics_of_bins_with_obsdata[-1]+1]

## calculate the zeropoint offset several times
median_model_rv_alliters = []
for coat in range(num_iter_rvoffset):
    modpart = all_models_large[coat*num_of_models:(coat+1)*num_of_models]

    # ## calculate the median model value for each time_finesampl coordinate
    model_rvs_med_all = []
    for cortim in range(res_time_finesampl):
        mod_rvs_cortim = modpart[:, cortim]
        model_rvs_med_all.append(np.median(mod_rvs_cortim))
    # ## calculate the RV per bin, but only for those with data
    bin_rvs_subset = []
    for huhu in indics_of_bins_with_obsdata:
        model_rvs_med_in_bin = model_rvs_med_all[indics_per_bin[huhu][0]:indics_per_bin[huhu][-1]]
        bin_rv = np.mean(model_rvs_med_in_bin)
        bin_rvs_subset.append(bin_rv)

    # ## calculate the median in the used datarange for this subset
    median_model_rv_alliters.append(np.median(bin_rvs_subset))

## get the median offset of the RV model for several iterations
ref_rv_model = np.median(median_model_rv_alliters)
ref_rv_err_model = np.std(median_model_rv_alliters)
print(ref_rv_model, ref_rv_err_model)

## calculate the shift of the night
#diff_of_medians_rvs_mod = np.median(rv) - np.median(modelrvs_of_bins_with_obsdata)
diff_of_medians_rvs_mod = np.median(rv) - ref_rv_model
if obsdate == '20201106':
    diff_of_medians_rvs_mod = np.median(rv[:-8]) - ref_rv_model
if obsdate == '20210308':
    diff_of_medians_rvs_mod = np.median(rv[:-7]) - ref_rv_model
if obsdate == '20210423':
    diff_of_medians_rvs_mod = np.median(rv[:-5]) - ref_rv_model
print(diff_of_medians_rvs_mod, ref_rv_err_model)
rv_nightshift = rv - diff_of_medians_rvs_mod

model_errbar_median = np.median(errors_using_all_rvs)
## make the plot
fig = plt.figure(figsize=(14, 6))
ax = fig.add_subplot(111)
plt.errorbar(time - tref, rv_nightshift, yerr=err, fmt='o', label='MaHPS data \nmedian error bar: {:.1f} m/s'.format(error_median))
plt.errorbar((_[1:] + _[:-1])/2, mean_rv_per_bin, yerr=errors_using_all_rvs, fmt='r*', label='model RVs (binned) \nmedian error bar: {:.1f} m/s'.format(model_errbar_median))  ##  mean, yerr=std
#plt.scatter(timecord_of_bins_with_obsdata, modelrvs_of_bins_with_obsdata, color='g')
plt.plot([], [], ' ', label='nightly RV shift: {:.1f} ± {:.1f} m/s'.format(diff_of_medians_rvs_mod, ref_rv_err_model))

plt.xlabel("Time [BJD - %d]" % tref, fontsize=18)
plt.ylabel("RV [m / s]", fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=15)
ax.set_title(plt_title, fontsize=18)
plt.tight_layout()

for jjijj in range(num_of_models):
    model_single = all_models[jjijj]
    
    plt.plot(time_finesampl - tref, model_single, 'C1', alpha=0.1)
### binning of the models END
#########################################


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

#plt.show()
plt.savefig(outfile_plot)


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
