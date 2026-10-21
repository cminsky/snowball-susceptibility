import numpy as np

### idealized carbon cycle box model: dN/dt = V - W + f(t), W = W_b * y^n, y = (N/N_b)^2; f < 0 removes carbon

def run_model(
                # model timing
                dt = 1e3/1e6,     # time step [Myr]
                t_max = 80,       # model run time [Myr]

                # initial state
                N0 = 2.83,        # background surficial carbon inventory [Emol]
                V_C = 8.5,        # volcanic degassing [Emol/Myr]
                n = 0.5,          # weathering feedback strength

                # perturbation
                f_fn = lambda t: 0, # forcing [Emol/Myr]
            ):

    t = 0
    N = N0
    W0 = V_C # background is a steady state

    t_arr, N_arr = [t], [N]
    while t < t_max:
        W = W0 * ((N / N0) ** 2) ** n
        N = max(0, N + (V_C - W + f_fn(t)) * dt)
        t += dt
        t_arr.append(t)
        N_arr.append(N)

    return {"t": np.array(t_arr), "N": np.array(N_arr)}

### analytic maximum excursion |nu'_max| for f = F0*exp(-theta*t), linearized about the background

def nu_max_exp(k, theta, F0, alpha, V):
    # k = alpha*V/N_b
    return F0 / (alpha * V) * np.exp(theta * np.log(theta / k) / (k - theta))

def numerical_nu_max_exp(Nb, V, theta, F0, n):
    # |nu'_max| from the box model with f = -F0*exp(-theta*t)
    k = 2 * n * V / Nb
    t_run = min(200.0, 5 / min(k, theta)) # set run length time
    N = run_model(dt=0.01, t_max=t_run, N0=Nb, V_C=V, n=n, f_fn=lambda t: -F0 * np.exp(-theta * t))["N"]
    return float(np.max(np.abs(N - Nb)) / Nb), t_run
