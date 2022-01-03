import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams["savefig.dpi"] = 100
rcParams["figure.dpi"] = 100
matplotlib.rcParams['axes.formatter.useoffset'] = False


frameid, time, rv, err, ord = np.loadtxt("../data/HD_149026_nocomb/RVs_ID3207_singleord_marmot_ccf_Mar.txt", unpack=True, delimiter=' ')
frameid2, time2, rv2, err2, ord2 = np.loadtxt("../data/HD_149026_fresh/RVs_ID3207_singleord_marmot.txt", unpack=True, delimiter=' ')

fig = plt.figure(figsize=(14, 6))


plt.errorbar(time - tref, rv, yerr=err, fmt='o')
plt.xlabel("Time [BJD - %d]" % tref, fontsize=16)
plt.ylabel("RV [m / s]", fontsize=16)

for i in range(100):
    model = compute(prior.sample())
    plt.plot(time_finesampl - tref, model, 'C1', alpha=0.1)

plt.show()
