import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from weathering import whak, mac

### loading GCM output (data/)

REF_GSOL, REF_CO2 = 1361, 280 # reference run (0 Ma, 1361 W/m2, 280 ppm)

def load_runs(variant, topo):
    # every equilibrated run in data/{variant}_{topo}.nc, in file order (by gsol, then CO2)
    with Dataset(f"data/{variant}_{topo}.nc") as nc:
        nc.set_auto_mask(False)
        area = nc["cell_area"][:]
        gsol, co2, age = nc["gsol"][:], nc["co2"][:], nc["age"][:]
        albedo, temp = nc["albedo"][:], nc["temp"][:]
        ts, runoff = nc["ts"][:], nc["runoff"][:]
    runs = []
    for i in range(len(gsol)):
        runs.append({"gsol": int(gsol[i]), "co2": int(co2[i]), "age": float(age[i]), "albedo": float(albedo[i]),
                     "temp": float(temp[i]), "ts": ts[i], "runoff": runoff[i], "area": area})
    return runs

def reference_run(variant="default"):
    for r in load_runs(variant, "0Ma"):
        if r["gsol"] == REF_GSOL and r["co2"] == REF_CO2:
            return r

def global_weathering(r, law="whak", **params):
    # WHAK or MAC summed over all cells; both are zero where there is no runoff, including the ocean
    if law == "mac":
        return float(np.sum(mac(r["ts"], r["runoff"], r["co2"]) * r["area"]))
    return float(np.sum(whak(r["ts"], r["runoff"], r["co2"], **params) * r["area"]))

def load_grid(variant, topo, law="whak", **params):
    # {gsol: {"co2s", "wsum", "albedo", "temp"}}, with wsum normalized to the reference run with the
    # same OHT; law and params choose the weathering law (default: WHAK with the baseline parameters)
    w_ref = global_weathering(reference_run(variant), law, **params)
    grid = {}
    for r in load_runs(variant, topo):
        if r["gsol"] not in grid:
            grid[r["gsol"]] = {"co2s": [], "wsum": [], "albedo": [], "temp": []}
        g = grid[r["gsol"]]
        g["co2s"].append(r["co2"])
        g["wsum"].append(global_weathering(r, law, **params) / w_ref)
        for k in ("albedo", "temp"):
            g[k].append(r[k])
    return grid

def load_ages(variant, topo):
    # {gsol: time [Gyr before present]} of the simulated solar constants
    return {r["gsol"]: r["age"] for r in load_runs(variant, topo)}

def band_means(variant, topo):
    # {gsol: {band: {"co2", "wsum"}}}, with wsum the baseline WHAK sum over the equatorial (< 30),
    # mid-latitude (30-60), or polar (> 60) band
    with Dataset(f"data/{variant}_{topo}.nc") as nc:
        lats = nc["lat"][:]
    bands_rows = {"equatorial": np.abs(lats) < 30, "mid_latitude": (np.abs(lats) >= 30) & (np.abs(lats) < 60),
                  "polar": np.abs(lats) >= 60}
    out = {}
    for r in load_runs(variant, topo):
        if r["gsol"] not in out:
            out[r["gsol"]] = {name: {"co2": [], "wsum": []} for name in bands_rows}
        bands = out[r["gsol"]]
        W = whak(r["ts"], r["runoff"], r["co2"])
        for name, rows in bands_rows.items():
            a = r["area"][rows]
            bands[name]["co2"].append(r["co2"])
            bands[name]["wsum"].append(float(np.sum(W[rows] * a)))
    return out

### solar evolution (Gough 1981)

def F_solar(t, F0=1361.0, t0=4.54):
    # solar constant [W/m2] at t Gyr before present
    return F0 / (1 + 2/5 * (1 - (t0 - t) / t0))

def t_gyr_from_gsol(gsol, F0=1361.0, t0=4.54):
    # time [Gyr before present] for a given solar constant [W/m2]
    return t0 * (F0 / gsol - 1.0) / 0.4

### plotting

# figure sizes for AGU (up to 170 mm wide, text >= 8 pt)
FULL_WIDTH = 170 / 25.4   # [in]
MEDIUM_WIDTH = 120 / 25.4 # single-panel figures [in]
ROW_HEIGHT = 2.4          # height of one row of panels, including tick and axis labels [in]
MEDIUM_HEIGHT = 3.2       # single-panel figures [in]

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8.5, "axes.labelsize": 8.5, "axes.titlesize": 8.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8.5, "legend.title_fontsize": 8.5,
    "figure.titlesize": 8.5, "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "mathtext.cal": "Arial", "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.5, "ytick.minor.width": 0.5, "lines.linewidth": 1.0, "lines.markersize": 3,
    "patch.linewidth": 0.5, "legend.borderaxespad": 0.3, "figure.dpi": 150, "savefig.dpi": 300,
    "savefig.pad_inches": 0.02, "pdf.fonttype": 42,
})

# (key, label, color) for each paleogeography
TOPOS_SERIES = [
    ("0Ma",         "0 Ma",           "tab:green"),
    ("0Ma_flat",    "0 Ma (flat)",    "lightgreen"),
    ("300Ma",       "300 Ma",         "tab:blue"),
    ("300Ma_flat",  "300 Ma (flat)",  "tab:cyan"),
    ("720Ma_Mer21", "720 Ma (Mer21)", "tab:red"),
    ("720Ma_Li23",  "720 Ma (Li23)",  "tab:orange"),
]
TOPOS = [key for key, _, _ in TOPOS_SERIES]

def setup_time_ax(ax):
    ax.set_xlabel("Time (Gya)")
    ax.invert_xaxis()
    ax.set_xlim(3, -0.1)

def place_legends_below(fig, legends, gap=0.03, y=0.0):
    # lays a row of figure legends (created with loc="upper center") out beneath the axes, centered as a group
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_fig = fig.transFigure.inverted()
    widths = [leg.get_window_extent(renderer).transformed(to_fig).width for leg in legends]
    x = 0.5 - (sum(widths) + gap * (len(widths) - 1)) / 2
    for leg, w in zip(legends, widths):
        leg.set_bbox_to_anchor((x + w / 2, y), transform=fig.transFigure)
        x += w + gap
    fig.canvas.draw()

### carbon inventory

N_REF = 2.83   # pre-industrial surficial carbon inventory [Emol]
CO2_REF = 280  # pre-industrial pCO2 [ppm]

def carbon_inventory(co2_ppm):
    # N = N_ref * sqrt(CO2/CO2_ref) [Emol]
    return N_REF * np.sqrt(np.asarray(co2_ppm, dtype=float) / CO2_REF)

### background state and Snowball threshold

def reference_T(variant="default"):
    # global mean T [K] of the reference run (0 Ma, 1361 W/m2, 280 ppm) with the same OHT
    return reference_run(variant)["temp"]

def background_co2_via_temperature(grid, gsol, T_target):
    # CO2 [ppm] where global mean surface temperature equals T_target, interpolating in log(CO2); None if the
    # simulations don't reach T_target
    co2s = np.array(grid[gsol]["co2s"], dtype=float)
    temps = np.array(grid[gsol]["temp"], dtype=float)
    if temps.min() > T_target or temps.max() < T_target:
        return None
    order = np.argsort(temps)
    return float(np.exp(np.interp(T_target, temps[order], np.log(co2s)[order])))

def wsum_at_co2(grid, gsol, co2):
    # normalized weathering at a given CO2 [ppm], interpolating log(wsum) in log(CO2)
    co2s = np.array(grid[gsol]["co2s"], dtype=float)
    idx = np.argsort(co2s)
    vals = np.array(grid[gsol]["wsum"], dtype=float)[idx]
    return float(np.exp(np.interp(np.log(co2), np.log(co2s[idx]), np.log(vals))))

def snowball_threshold_co2(grid, gsol):
    # highest simulated CO2 [ppm] with global mean albedo above 0.6; NaN if there is no Snowball at this gsol
    co2s = np.array(grid[gsol]["co2s"], dtype=float)
    snowball = np.array(grid[gsol]["albedo"], dtype=float) > 0.6
    if not np.any(snowball):
        return np.nan
    return float(co2s[snowball].max())

def compute_stats(variant, topo):
    # background CO2, the weathering there, and Snowball threshold CO2 at each gsol
    grid = load_grid(variant, topo)
    ages = load_ages(variant, topo)
    T_target = reference_T(variant)
    gsols = np.array(sorted(grid))
    bg_co2, bg_wsum, sb_co2 = [], [], []
    for g in gsols:
        co2 = background_co2_via_temperature(grid, g, T_target)
        if co2 is None:
            bg_co2.append(np.nan)
            bg_wsum.append(np.nan)
        else:
            bg_co2.append(co2)
            bg_wsum.append(wsum_at_co2(grid, g, co2))
        sb_co2.append(snowball_threshold_co2(grid, g))
    return {"gsol": gsols, "age": np.array([ages[g] for g in gsols]),
            "background": {"co2": np.array(bg_co2), "wsum": np.array(bg_wsum)}, "snowball": {"co2": np.array(sb_co2)}}
