import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import os 
import datetime as dt
from chamber_tools import linear_fit, exponential_fit, calc_flux

import config


class ChamberFluxCalculator:
    def __init__(self, meta_data_path, data_path, fit_type, par_sd_limit, nrmse_limit, remove_from_start, min_length, plotting):
        self.meta_data_path = meta_data_path

        if fit_type == 'linear':
            self.fit_type = 'linear'
            self.fit = linear_fit
        else:
            self.fit_type = 'exponential'
            self.fit = exponential_fit

        self.data_path = data_path
        self.par_sd_limit = par_sd_limit
        self.nrmse_limit = nrmse_limit
        self.remove_from_start = remove_from_start
        self.min_length = min_length 
        self.plotting = plotting
        if 'Results' not in os.listdir():
            os.mkdir('Results/')

    def read_meta_data(self):
        meta_data = pd.read_csv(self.meta_data_path, delimiter=";")
        meta_data.Date = pd.to_datetime(meta_data.Date, format="%d.%m.%Y").dt.date
        for col in ['Start', 'End']:
            meta_data[col] = pd.to_datetime(meta_data[col], format="%H:%M:%S").dt.time
        return meta_data

    def calc_flux_and_stats(self, cur, mean_temp, pres, system_height):
        co2 = np.asarray(cur.co2)
        secs = np.asarray((cur.index-cur.index[0]).seconds)
        par = np.asarray(cur.par)

        slope, nrmse, co2_hat = self.fit(secs, co2)                   

        if np.isnan(mean_temp):
            mean_temp = np.mean(cur.temp)

        mean_par = np.mean(par)
        sd_par = np.std(par)

        flux = calc_flux(slope, pres, mean_temp, system_height)

        len_closure = (cur.index[-1]-cur.index[0]).seconds 

        return flux, nrmse, mean_par, sd_par, mean_temp, len_closure, co2_hat, slope

    # i, orig_data, cur_data, co2_hat, plot, date, seq, meas_ok, nrmse, par, tang
    def make_plot(self, fig, collar, data, acc_meas, i, meas_ok, co2_hat, slope, nrmse, sd_par):
        ax = fig.add_subplot(3, 2, i)
        secs = np.asarray((data.index-data.index[0]).seconds)
        acc_secs = (acc_meas.index-data.index[0]).seconds
        start_in_secs = (acc_meas.index[0]-data.index[0]).seconds
        end_in_secs = (acc_meas.index[-1]-data.index[0]).seconds
        ax.scatter(secs, data.co2, c='r', s=2, label="CO2")
        ax.scatter(acc_secs, acc_meas.co2, c='k', s=2, label="CO2")
        col = 'k'
        if meas_ok == False:
            col = 'r'
        ax.plot(acc_secs, co2_hat, c='k', label=self.fit_type+' fit')
        if self.fit_type == 'exponential':
            tang = co2_hat[0]+slope*np.arange(len(acc_secs))
            ax.plot(acc_secs, tang, c='g', label='tangent')
        ax.axvline(x=start_in_secs, c=col)
        ax.axvline(x=end_in_secs, c=col)
        ax.legend(loc='upper left', prop={'size': 6})
        ax.set_ylabel("CO$_2$ [ppm]", fontsize=6)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.set_title("NRMSE: "+str(round(nrmse, 2)) +
                    " PAR_sd: "+str(round(sd_par, 2)), fontsize=6)
        axx = ax.twinx()
        axx.plot(secs, data.par, c='b', label="PAR")
        axx.set_ylabel("PAR [$\mu$mol m$^-$$^2$ s$^-$$^1$]", fontsize=6)
        axx.set_ylim(0, 2000)
        axx.tick_params(axis='both', which='major', labelsize=6)
        axx.legend(loc='upper right', prop={'size': 6})
        date = data.index[0].date()
        plt.suptitle("Date: "+str(date)+" Collar: " + collar)
        
        return fig


    def calc_fluxes(self, meta_data):
        def check_nrmse(flux, nrmse):
            small_flux = 0.1
            if flux < small_flux:
                nrmse_not_ok = False
            else:
                if nrmse > self.nrmse_limit:
                    nrmse_not_ok = True
                else:
                    nrmse_not_ok = False
            return nrmse_not_ok

        results = pd.DataFrame(columns=['Date', 'Collar', 'Closure', 'Flux', 'NRMSE', 'PAR'])
        c = 0
        # Browse through dates
        for date in meta_data.Date.unique():
            data = pd.read_csv(self.data_path+'chamber_data_'+str(date)+'.csv', index_col=0)

            data.index = pd.to_datetime(data.index, format="%Y-%m-%d %H:%M:%S")
            
            # Browse through collars
            for collar in meta_data[meta_data.Date==date].Collar.unique():
                if self.plotting:
                    fig = plt.figure()

                cl = 1
                # Browse through closures
                for id in meta_data[(meta_data.Date==date) & (meta_data.Collar==collar)].index:

                    pres = meta_data.loc[id, "Air_pres (hPa)"]
                    system_height = meta_data.loc[id, "System_height (m)"]
                    mean_temp = meta_data.loc[id, "Temperature (K)"]

                    start_time = meta_data.iloc[id].Start
                    end_time = meta_data.iloc[id].End
                    orig_start = dt.datetime.combine(date, start_time)
                    orig_end = dt.datetime.combine(date, end_time)
                    orig_meas = data[orig_start:orig_end]
                    start = dt.datetime.combine(date, start_time)+dt.timedelta(seconds=self.remove_from_start)
                    end = dt.datetime.combine(date, end_time)
                    cur = data[start:end]

                    if len(cur) < 60:
                        print('Measurement too short: '+str(date)+', Collar: '+ str(collar) + ', Closure: '+str(cl))

                    else:
                        flux, nrmse, mean_par, sd_par, mean_temp, len_closure, co2_hat, slope = self.calc_flux_and_stats(cur, mean_temp, pres, system_height)

                        nrmse_not_ok = check_nrmse(flux, nrmse)

                        while len_closure >= self.min_length and (nrmse_not_ok or sd_par>self.par_sd_limit):
                            end -= dt.timedelta(seconds=10)
                            cur = data[start:end]

                            flux, nrmse, mean_par, sd_par, mean_temp, len_closure, co2_hat, slope = self.calc_flux_and_stats(cur, mean_temp, pres, system_height)
                            nrmse_not_ok = check_nrmse(flux, nrmse)

                        if np.isnan(mean_par):
                            print('PAR missing: '+str(date)+', Collar: '+ str(collar) + ', Closure: '+str(cl))

                        if np.isnan(mean_temp):
                            print('Temperature missing: '+str(date)+', Collar: '+ str(collar) + ', Closure: '+str(cl))

                        if len_closure > self.min_length and ~np.isnan(mean_temp):
                            results.loc[c] = date, collar, cl, flux, nrmse, mean_par
                            meas_ok = True
                            if self.plotting:
                                self.make_plot(fig, collar, orig_meas, cur, cl, meas_ok, co2_hat, slope, nrmse, sd_par)
                        else:
                            meas_ok = False

                        c+=1
                        cl+=1

                if self.plotting:
                    fig.tight_layout()
                    fig.savefig('Results/'+str(date)+'_'+str(collar)+'.png')
                    plt.close()

        return results 


    def run(self):
        meta_data = self.read_meta_data()
        results = self.calc_fluxes(meta_data)        
        results.to_csv('Results/flux_results.csv', index=False)


if __name__ == "__main__":
    calc = ChamberFluxCalculator(config.meta_data_path, config.data_path, config.fit_type, config.par_sd_limit, config.nrmse_limit, config.remove_from_start, config.min_length, config.plotting)
    calc.run()