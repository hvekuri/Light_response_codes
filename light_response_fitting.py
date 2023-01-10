import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import os 
import datetime as dt
from chamber_tools import fGPP
from lmfit import Model, Parameters
import config


class LightResponseFitting:
    def __init__(self, data_path, alpha_min, alpha_max, GPmax_min, GPmax_max, plotting):
        self.data_path = data_path
        self.alpha_min = alpha_min
        self.alpha_max = alpha_max
        self.GPmax_min = GPmax_min
        self.GPmax_max = GPmax_max
        self.plotting = plotting

    def fit_LR(self, GPP, par):
        # Define the fit model
        GPPmodel = Model(fGPP)

        # Use Parameter class for model params
        params = Parameters()
        params.add_many(("alpha", -0.001, True, self.alpha_min, self.alpha_max), ("GPmax", -1,
                                                                                              True, self.GPmax_min, self.GPmax_max))

        # Fit
        result = GPPmodel.fit(GPP, PAR=par,
                              alpha=params["alpha"], GPmax=params["GPmax"], method="leastsq")

        alpha_fit = result.params['alpha'].value
        GPmax_fit = result.params['GPmax'].value
        alpha_se = result.params['alpha'].stderr
        GPmax_se = result.params['GPmax'].stderr

        GP1200 = fGPP(1200, alpha_fit, GPmax_fit)
        if(result.covar is not None):
            GP1200_unc = result.eval_uncertainty(sigma=1.96, PAR=1200)[0]
        else:
            GP1200_unc = np.nan

        return alpha_fit, alpha_se, GPmax_fit, GPmax_se, GP1200, GP1200_unc, result

    def make_plot(self, par, GPP, result, alpha, GPmax, collar, date):
        fig, ax = plt.subplots()
        collar = str(collar)
        # Calculate the fitted function and 2-sigma uncertainty at arbitary x
        xp = np.linspace(0, 2000)
        GPP_fit = fGPP(xp, alpha, GPmax)
        if(result.covar is not None):
            GPP_unc = result.eval_uncertainty(sigma=1.96, PAR=xp)
            ax.fill_between(xp, GPP_fit-GPP_unc, GPP_fit +
                            GPP_unc, alpha=0.5, color="grey")
        ax.scatter(par, GPP, s=6, c='k')
        ax.plot(xp, GPP_fit, c='k', label="collar: "+str(collar))
        ax.set_title("Collar " + collar)
        ax.set_xlabel("PAR [$\mu$mol m$^-$$^2$ s$^-$$^1$]")
        ax.set_ylabel("GPP [mg CO$_2$ m$^-$$^2$ s$^-$$^1$]")
        ax.grid()
        fig.tight_layout()
        fig.savefig('Results/LR_'+str(date)+"_"+str(collar)+".png")


    def run(self):
        results = pd.DataFrame(columns=["Date", "Collar", "Alpha", "Alpha_se", "GPmax", "GPmax_se", "NEE1200", "GP1200", "GP1200unc", "Reco"])
        idx = 0

        data = pd.read_csv(self.data_path)
        data.Date = pd.to_datetime(data.Date).dt.date
        for date in data.Date.unique():
            for collar in data[data.Date==date].Collar.unique():
                cur = data[(data.Date==date) & (data.Collar==collar) & (data.PAR.notnull())]
                if len(cur) >= 3 and len(cur[cur.PAR<10])>0:

                    Reco = cur[cur.PAR<10].Flux.values[0]

                    par = np.array(cur.PAR)
                    GPP = np.array(cur.Flux-Reco)

                    alpha_fit, alpha_se, GPmax_fit, GPmax_se, GP1200, GP1200_unc,  result = self.fit_LR(
                        GPP, par)

                    results.loc[idx] = date, collar, alpha_fit, alpha_se, GPmax_fit, GPmax_se, GP1200 + \
                        Reco, GP1200, GP1200_unc, Reco

                    if self.plotting:
                        self.make_plot(par, GPP, result, alpha_fit, GPmax_fit, collar, date)

                    idx += 1

        results.to_csv('Results/light_responses.csv', index = False)

if __name__ == "__main__":
    lr_fitter = LightResponseFitting(config.flux_data_path, config.alpha_bounds[0], config.alpha_bounds[1], config.GPmax_bounds[0], config.GPmax_bounds[1], config.plotting)
    lr_fitter.run()