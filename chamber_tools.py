import numpy as np
import datetime as dt
import math
from scipy import stats
from lmfit import Model, Parameters

def calc_flux(slope, pres, temp, system_height):
    # pres in [hPa]
    # temp in [K]
    # system_height in [m2]
    # -> returns flux in [mg m-2 s-1]
    if np.isnan(pres):
        flux = slope*(44.01/0.0224)*(273.15/(temp))*(system_height)
    else:
        flux = slope * 44.01*pres*100 * system_height/(8.31446*temp)*0.001
    return flux


def calc_NRMSE(y, y_hat):
    return math.sqrt(1/(len(y))*((y_hat-y)**2).sum())/(max(y)-min(y))



# Does linear fitting, returns slope, nrmse, rmse and predicted co2
def linear_fit(secs_lic, co2):
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        secs_lic, co2)
    co2_hat = intercept+slope*secs_lic
    nrmse = calc_NRMSE(np.asarray(co2), np.asarray(co2_hat))
    return slope, nrmse, co2_hat


# Does exponential fitting first with 17th order taylor series to get initial guesses
# Returns tangent's slope at s=0, nrmse, rmse and predicted co2
def exponential_fit(secs_lic, co2):
    co2_ = np.asarray(co2).reshape(-1, 1)
    t = np.asarray(secs_lic).reshape(-1, 1)

    # Fit Taylor series
    fluxmodel_taylor = Model(feksp_taylor)
    params_T = Parameters()
    params_T.add_many(("a_T", 300), ("b_T", 300), ("c_T", -0.000002))
    result_T = fluxmodel_taylor.fit(
        co2_, secs=t, a=params_T["a_T"], b=params_T["b_T"], c=params_T["c_T"])

    # Fit exponential model
    fluxmodel = Model(feksp)
    params = Parameters()
    params.add_many(("a", result_T.params["a"].value),
                    ("b", result_T.params["b"].value), ("c", result_T.params["c"].value))
    result = fluxmodel.fit(co2_, secs=t, a=params["a"],
                           b=params["b"], c=params["c"])

    a_fit = result.params["a"].value
    b_fit = result.params["b"].value
    c_fit = result.params["c"].value
    slope = (b_fit-a_fit)*c_fit
    co2_hat = feksp(secs_lic, a_fit, b_fit, c_fit)
    nrmse = calc_NRMSE(np.asarray(co2), np.asarray(co2_hat))

    return slope, nrmse, co2_hat


def feksp(secs, a, b, c):
    return a+(b-a)*np.exp(c*secs)


def feksp_taylor(secs, a, b, c):
    def ftaylorseries(c, x):
        series = 0
        for i in range(17):
            series += (c**i*x**i)/math.factorial(i)
        return series
    return a+(b-a)*ftaylorseries(c, secs)


def fGPP(PAR, alpha, GPmax):
    return alpha*GPmax*PAR/(alpha*PAR + GPmax)


