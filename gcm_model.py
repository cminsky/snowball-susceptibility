import numpy as np
from helpers import background_co2_via_temperature, reference_T

### GCM climate lookup at the nearest simulated solar constant: albedo and log(weathering) interpolated in log(pCO2)

def climate_curves(gcm_data):
    # albedo and log(normalized weathering) vs. log(pCO2) at each solar constant, from a GCM grid (helpers.load_grid)
    curves = {}
    for gsol, entry in gcm_data.items():
        co2 = np.asarray(entry["co2s"], dtype=float)
        order = np.argsort(co2)
        curves[float(gsol)] = {"log_co2": np.log(co2[order]), "albedo": np.asarray(entry["albedo"], dtype=float)[order],
                               "log_wsum": np.log(np.asarray(entry["wsum"], dtype=float)[order])}
    return curves

def band_curves(band_data):
    # log(weathering) vs. log(pCO2) in each latitude band at each solar constant, from helpers.band_means
    curves = {}
    for gsol, entry in band_data.items():
        curves[float(gsol)] = {}
        for band, b in entry.items():
            co2 = np.asarray(b["co2"], dtype=float)
            order = np.argsort(co2)
            curves[float(gsol)][band] = {"log_co2": np.log(co2[order]),
                                         "log_wsum": np.log(np.asarray(b["wsum"], dtype=float)[order])}
    return curves

def nearest_gsol(curves, gsol):
    gsols = np.array(sorted(curves))
    return float(gsols[np.argmin(np.abs(gsols - float(gsol)))])

def interp(curve, pCO2, field):
    return float(np.interp(np.log(float(pCO2)), curve["log_co2"], curve[field]))

### carbon cycle model coupled to the GCM climate

def run_model(
                gcm_data,         # GCM grid (helpers.load_grid), no OHT
                gsol,             # solar constant [W/m2]

                # model timing
                dt = 1e3/1e6,     # time step [Myr]
                t_max = 80,       # model run time [Myr]

                # carbon reference state
                pCO2_ref = 280,   # reference pCO2 [ppm]
                N_ref = 2.83,     # reference surficial carbon inventory [Emol]

                # volcanic degassing and weathering
                V_C = 8.5,        # volcanic degassing at wsum=1 [Emol/Myr], rescaled to balance weathering at the start
                W_sea = 1.6,      # seafloor weathering [Emol/Myr]

                # perturbation
                C_imb = 0,        # carbon removal [Emol/Myr]
                theta = 0,        # decay rate of forcing, F(t) = C_imb*exp(-theta*t) [1/Myr]
                snowball_albedo_threshold = 0.6,

                # latitude-dependent LIP weathering (off when lip_band is None)
                lip_band = None,  # "equatorial", "mid_latitude", or "polar"
                band_data = None, # helpers.band_means
            ):

    # --- setup ---
    t = 0
    curves = climate_curves(gcm_data)
    climate = curves[nearest_gsol(curves, gsol)]
    if lip_band is not None:
        bands = band_curves(band_data)
        band = bands[nearest_gsol(bands, gsol)][lip_band]

    # --- initial conditions: background pCO2, where global mean T equals that of the reference run ---
    pCO20 = background_co2_via_temperature(gcm_data, nearest_gsol(curves, gsol), reference_T())
    a = interp(climate, pCO20, "albedo")
    N = N0 = N_ref * np.sqrt(pCO20 / pCO2_ref)

    # degassing balances weathering at the initial state
    W_sil0 = V_C - W_sea
    V_C = W_sil0 * np.exp(interp(climate, pCO20, "log_wsum")) + W_sea

    # --- run model ---
    t_arr, N_arr = [t], [N]
    while t < t_max:

        # once in a Snowball, stay there until t_max
        if a > snowball_albedo_threshold:
            for _ in range(round((t_max - t) / dt)):
                t += dt
                t_arr.append(t); N_arr.append(N)
            break

        # --- update sources and sinks ---
        pCO2 = pCO2_ref * (N / N_ref) ** 2
        W_sil = W_sil0 * np.exp(interp(climate, pCO2, "log_wsum"))
        C = C_imb * np.exp(-theta*t)
        if lip_band is not None:
            # LIP weathering scales with the local weathering in its latitude band
            C *= float(np.exp(interp(band, pCO2, "log_wsum") - interp(band, pCO20, "log_wsum")))

        # --- update reservoir and albedo ---
        N = max(0, N + (V_C - W_sil - W_sea - C) * dt)
        a = interp(climate, pCO2_ref * (N / N_ref) ** 2, "albedo")

        t += dt
        t_arr.append(t); N_arr.append(N)

    return {"t": np.array(t_arr), "N": np.array(N_arr), "N0": N0}
