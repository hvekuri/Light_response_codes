import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from chamber_tools import fGPP, fGPP_VI, fRespiration_VI
from lmfit import Model, Parameters
import config


class LightResponseFitting:
    def __init__(self, data_path, param_bounds, GPP_model):
        self.data_path = data_path
        self.param_bounds = param_bounds
        self.GPP_model = GPP_model

    def fit_LR(self, GPP, par, rvi):
        # Define the fit model
        GPPmodel = Model(fGPP_VI, independent_vars = ['PAR', 'VI'])

        # Use Parameter class for model params
        params = Parameters()
        params.add_many(("alpha", -0.001, True, self.param_bounds[0], self.param_bounds[1]), ("GPmax", -1,
                                                                                              True, self.param_bounds[2], self.param_bounds[3]))

        # Fit
        result = GPPmodel.fit(GPP, PAR=par, VI=rvi,
                          alpha=params["alpha"], GPmax=params["GPmax"], method="leastsq")

        return result

    def fit_TR(self, Reco, tair, tsoil, rvi):
        # Define the fit model
        Rmodel = Model(fRespiration_VI, independent_vars = ['Tair', 'Tsoil', 'VI'])

        # Use Parameter class for model params
        params = Parameters()
        
        min_Rd0 = config.param_bounds[4]
        max_Rd0 = config.param_bounds[5]
        params.add('Rd0', 0.1, True, min_Rd0, max_Rd0)

        min_Rs0 = config.param_bounds[6]
        max_Rs0 = config.param_bounds[7]
        params.add('Rs0', value=0.1, vary=True, min=min_Rs0, max=max_Rs0)

        params.add('Es', value=config.Es, vary=False)

        params.add('bd', value=config.bd, vary=False)

        # Fit
        result = Rmodel.fit(Reco, Tair=tair, Tsoil = tsoil, VI = rvi,
                          Rs0=params["Rs0"], Rd0=params["Rd0"], Es = params['Es'], bd = params['bd'], method="leastsq")


        return result


    def run(self):
        results = pd.DataFrame()
        idx = 0

        data = pd.read_csv(self.data_path, delimiter=";")

        # Drop empty rows
        data = data[data.PAR.notnull()]
        data.Date = pd.to_datetime(data.Date, format="%d.%m.%Y").dt.date

        # TEMPERATURE RESPONSES
        resp_params = pd.DataFrame()

        for id in data.ID.unique():
            cur = data[data.ID==id]
            resp_data = cur[cur.PAR<10].copy()
            tair = resp_data['T_air']
            tsoil = resp_data['T_soil']
            reco = resp_data['Flux']
            rvi = resp_data['VI']

            result = self.fit_TR(reco, tair, tsoil, rvi)

            for p in ['Rd0', 'Rs0', 'Es', 'bd']:
                resp_params.loc[id, p] = result.params[p].value
                resp_params.loc[id, p+'_sd'] = result.params[p].stderr


        # LIGHT RESPONSES
        for date in data.Date.unique():
            for id in data[data.Date==date].ID.unique():
                cur = data[(data.Date==date) & (data.ID==id) & (data.VI>0)]
                if len(cur) >= 3 and len(cur[cur.PAR<10])>0:

                    Reco = cur[(cur.PAR<10)].Flux.values.mean()

                    par = np.array(cur.PAR)
                    GPP = np.array(cur.Flux-Reco)

                    if self.GPP_model == 'GPP_VI':
                        rvi = np.array(cur.VI)
                    else:
                        rvi = np.ones(len(par))

                   
                    if len(cur[(cur.PAR>10) & (cur.VI>0)]) >= 2:

                        result = self.fit_LR(
                            GPP, par, rvi)

                        results.loc[idx, 'Date'] = date
                        results.loc[idx, 'ID'] = id

                        for p in ['alpha', 'GPmax']:
                            results.loc[idx, p] = result.params[p].value
                            results.loc[idx, p+'_sd'] = result.params[p].stderr

                        results.loc[idx, 'Reco'] = Reco
     
                        idx += 1

        for id in resp_params.index.unique():
            for var in resp_params.columns:
                results.loc[results.ID==id, var] = resp_params.loc[id, var]

        
        results.to_csv('Results/fit_params.csv', index = False)
        return results


if __name__ == "__main__":
    lr_fitter = LightResponseFitting(config.flux_data_path, config.param_bounds, config.GPP_model)
    lr_fitter.run()