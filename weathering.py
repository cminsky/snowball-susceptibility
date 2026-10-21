import numpy as np

### WHAK weathering (Walker et al. 1981) in each grid cell,
###   W = (q/q_ref)^mu exp[gamma (T - T_ref)] (pCO2/pCO2_ref)^beta
### with q_ref = 300 mm/yr, T_ref = 288 K and pCO2_ref = 280 ppm, which cancel in W/W_ref

def whak(ts, runoff, co2, mu=0.75, gamma=0.095, beta=0.0):
    # for surface temperature ts [K], runoff [m/s] and CO2 [ppm]; zero where runoff <= 0. Default parameters are
    # West et al. (2005)
    q = np.maximum(runoff, 0.0) * 1e3 * 86400.0 * 365.0 # m/s -> mm/yr
    W = (q / 300.0) ** mu * np.exp(gamma * (ts - 288.0)) * (co2 / 280.0) ** beta
    return np.where(q > 0, W, 0.0)

### MAC weathering (Maher & Chamberlain 2014) in each grid cell, in the form of Graham & Pierrehumbert (2020), Eq. 24:
###   w = alpha / (1/k_eff + m A t_s + alpha/(q C_eq))
###   k_eff = k0 exp[(T - T_ref)/T_e] (pCO2/pCO2_ref)^beta,  C_eq = Lambda pCO2^n

# Graham & Pierrehumbert (2020) Table 1
MAC_PARAMS = {
    "alpha": 3.39e5,      # [-]
    "t_s": 1.0e5,         # soil age [yr]
    "m": 0.27,            # mineral molar mass [kg/mol]
    "A": 100.0,           # specific surface area [m2/kg]
    "Lambda": 1.4,        # [mol/m3] (1.4e-3 mol/L)
    "n": 0.316,           # thermodynamic pCO2 exponent
    "k0": 8.7e-6,         # reference rate constant [mol/m2/yr]
    "beta": 0.2,          # kinetic pCO2 exponent
    "T_e": 11.1,          # kinetic e-folding temperature [K]
    "T_ref": 288.0,       # [K]
    "pco2_ref": 280e-6,   # [bar]
}

def mac(ts, runoff, co2, p=MAC_PARAMS):
    # weathering rate [mol/m2/yr] for surface temperature ts [K], runoff [m/s] and CO2 [ppm]; zero where runoff <= 0
    q = np.maximum(runoff, 0.0) * 86400.0 * 365.0 # m/s -> m/yr
    pco2 = co2 * 1e-6 # ppm -> bar
    k_eff = p["k0"] * np.exp((ts - p["T_ref"]) / p["T_e"]) * (pco2 / p["pco2_ref"]) ** p["beta"]
    qC = q * (p["Lambda"] * pco2 ** p["n"]) # q C_eq
    w = p["alpha"] / (1 / k_eff + p["m"] * p["A"] * p["t_s"] + p["alpha"] / np.where(qC > 0, qC, 1.0))
    return np.where(qC > 0, w, 0.0)
