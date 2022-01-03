# make a Keplerian fit configuration file

# Required packages for setup
# myradvel is available via Ubuntu, but not in the Windows anaconda environments
# import myradvel
import paths_and_files as pf


def make_radvel_conffile(redmine_id, n_cand, inst_list):
    with open(pf.radvel_config.format(redmine_id), 'w') as conffile:
        conffile.write("import myradvel\n")
        conffile.write("import paths_and_files as pf\n")
        conffile.write("import pandas as pd\n")
        conffile.write("import numpy as np\n")
        conffile.write("\n")
        conffile.write("redmine_id = {}\n".format(redmine_id))
        # Define global planetary system and dataset parameters
        conffile.write("starname = 'ID{}'\n".format(redmine_id))
        # number of planets in the system
        conffile.write("nplanets = {}\n".format(int(n_cand)))
        # list of instrument names. Can be whatever you like (no spaces)
        # but should match 'tel' column in the input file.
        instnames = list(inst_list)
        conffile.write("instnames = {}\n".format(list(inst_list)))
        # number of instruments with unique velocity zero-points
        conffile.write("ntels = {}\n".format(len(instnames)))
        # Fitting basis, see myradvel.basis.BASIS_NAMES for available basis names
        conffile.write("fitting_basis = 'per tp e w k'\n")
        # reference epoch for RV timestamps (i.e. this number has been subtracted off your timestamps)
        conffile.write("bjd0 = 0\n")
        # map the numbers in the Parameters keys to planet letters (for plotting and tables)
        conffile.write("planet_letters = {}\n")
        conffile.write("for n in range(nplanets):\n")
        conffile.write("    planet_letters.update({n + 1: chr(98 + n)})\n")
        # chr(98) = 'b', count characters up from there

        # Define prior centers (initial guesses) in a basis of your choice (need not be in the fitting basis)
        # initialize Parameters object
        conffile.write("anybasis_params = myradvel.Parameters(nplanets, "
                       "basis='per tp e w k', planet_letters=planet_letters)\n")

        conffile.write("lit_params = []\n")
        conffile.write("with open(pf.lit_planet_params.format(redmine_id)) as paramfile:\n")
        conffile.write("    for line in paramfile:\n")
        conffile.write("        line = line.strip()\n")
        conffile.write("        line = line.split()\n")
        conffile.write("        for entry in range(len(line)):\n")
        conffile.write("            line[entry] = float(line[entry])\n")
        conffile.write("        lit_params.append(line)\n")
        conffile.write("for n in range(nplanets):\n")
        # period of n-th planet
        conffile.write("    anybasis_params['per{}'.format(str(n+1))] = "
                       "myradvel.Parameter(value=float(lit_params[n][0]))\n")
        # time of inferior conjunction of n-th planet
        conffile.write("    anybasis_params['tp{}'.format(str(n+1))] = myradvel.Parameter(value=lit_params[n][1])\n")
        # eccentricity of n-th planet
        conffile.write("    anybasis_params['e{}'.format(str(n+1))] = myradvel.Parameter(value=lit_params[n][2])\n")
        # argument of periastron of the star's orbit for n-th planet
        conffile.write("    anybasis_params['w{}'.format(str(n+1))] = "
                       "myradvel.Parameter(value=(2.0 * np.pi * lit_params[n][3] / 360.0))\n")
        # velocity semi-amplitude for n-th planet
        conffile.write("    anybasis_params['k{}'.format(str(n+1))] = myradvel.Parameter(value=lit_params[n][4])\n")

        # abscissa for slope and curvature terms (should be near mid-point of time baseline)
        conffile.write("time_base = 2456778\n")
        # conffile.write("time_base = 2459097\n")

        # # add a part for 51 Peg fitting test like Hanna does it
        # conffile.write("mod = myradvel.RVModel(anybasis_params, time_base=time_base)\n")

        # slope: (If rv is m/s and time is days then [dvdt] is m/s/day)
        conffile.write("anybasis_params['dvdt'] = myradvel.Parameter(value=0.0)\n")
        # curvature: (If rv is m/s and time is days then [curv] is m/s/day^2)
        conffile.write("anybasis_params['curv'] = myradvel.Parameter(value=0.0)\n")

        # analytically calculate gamma if vary=False and linear=True
        conffile.write("for inst in instnames:\n")
        # velocity zero-point for each instrument
        conffile.write("    anybasis_params['gamma_{}'.format(inst)] = myradvel.Parameter(value=0.0)\n")
        # jitter for each instrument
        conffile.write("    anybasis_params['jit_{}'.format(inst)] = myradvel.Parameter(value=2.6)\n")

        # Convert input orbital parameters into the fitting basis
        conffile.write("params = anybasis_params.basis.to_any_basis(anybasis_params, fitting_basis)\n")

        # Set the 'vary' attributes of each of the parameters in the fitting basis.
        # A parameter's 'vary' attribute should be set to False if you wish to hold it fixed
        # during the fitting process. By default, all 'vary' parameters are set to True.
        conffile.write("params['dvdt'].vary = False\n")
        conffile.write("params['curv'].vary = False\n")
        # conffile.write("params['gamma'].vary = False\n")
        conffile.write("for n in range(nplanets):\n")
        conffile.write("    params['per{}'.format(str(n+1))].vary = False\n")
        conffile.write("    params['tp{}'.format(str(n+1))].vary = True\n")
        conffile.write("    params['e{}'.format(str(n+1))].vary = False\n")
        conffile.write("    params['w{}'.format(str(n+1))].vary = True\n")
        conffile.write("    params['k{}'.format(str(n+1))].vary = False\n")

        # Load radial velocity data, the data is contained in an ASCII file,
        # must have 'time', 'mnvel', 'errvel', and 'tel' keys
        # the velocities are expected to be in m/s
        conffile.write("data = pd.read_csv(pf.out_RVs_compare.format(redmine_id), sep=' ')\n")

        # Define prior shapes and widths here.
        conffile.write("priors = [\n")
        # Keeps eccentricity < 1
        conffile.write("    myradvel.prior.EccentricityPrior(nplanets),\n")
        # Gaussian prior on tc1 with center at tc1 and width 300 days
        conffile.write("    myradvel.prior.Gaussian('tp1', params['tp1'].value, 300.0),\n")
        # # # additional priors for 51 Peg fit
        # conffile.write("    myradvel.prior.HardBounds('e1', 0., 0.1),\n")
        # conffile.write("    myradvel.prior.HardBounds('k1', 50.0, 60.0),\n")
        # conffile.write("    myradvel.prior.HardBounds('per1', 4.2, 4.3)\n")
        conffile.write("]\n")

        conffile.write("for inst in instnames:\n")
        conffile.write("    priors.append(myradvel.prior.HardBounds('jit_{}'.format(inst), 0.0, 10.0))\n")

        # print(priors)

        # optional argument that can contain stellar mass in solar units (mstar) and
        # uncertainty (mstar_err). If not set, mstar will be set to nan.
        conffile.write("stellar = dict(mstar=1.12, mstar_err=0.06)\n")

    return pf.radvel_config.format(redmine_id)


# make_radvel_conffile(2864, 3, ['foces', 'lick', 'elodie'])
