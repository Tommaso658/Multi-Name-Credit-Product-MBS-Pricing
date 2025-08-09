# NOTE: The code below contains numerous pricing functions for credit tranches
# using different copula models (Gaussian, t-Student, Vasicek), and calibration
# routines. All function docstrings and inline comments have been translated to English
# without altering the original function logic.

# The code contains exact pricing, approximations (KL, HP, LHP), and calibration
# routines using optimization, Monte Carlo simulations, and analytical formulas.

# Please refer to the function docstrings for parameter descriptions and return values.

# Due to the large number of functions in this script, each one is independently documented.

# Main categories include:
# - Exact/HP pricing under Vasicek or t-Student copula
# - KL approximation pricing
# - LHP (Large Homogeneous Portfolio) approximation
# - Monte Carlo pricing using Gaussian or t-Student copula
# - Calibration of degrees of freedom ν via minimization of MSE
# - Plotting utilities for results
#
# The code supports analysis of tranche pricing across varying portfolio sizes (I),
# calibration of copula parameters, and comparison between models (Gaussian vs t-Student).

# ---- Libraries needed ----

import numpy as np
from scipy.stats import norm
from scipy.special import comb
from scipy.integrate import quad

from scipy.optimize import brentq
from scipy.optimize import fsolve
import numpy as np
from scipy.stats import t
from scipy.special import gammaln
from scipy.integrate import quad


from scipy.special import comb
from matplotlib import pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from scipy.integrate import fixed_quad

# ---- Functions ----
# -------------------






######## ---- Double-t-student prices and Model ------ #########


def price_cumul_tranche_HP_double_t_student_optimized(Kd, Ku, avg_recovery, rho, p, discount, I, nu):
    """
    Optimized pricing of a cumulative tranche [Kd, Ku] under the double t-student copula (HP approximation).

    Parameters:
        Kd         : float – lower attachment (e.g. 0.0)
        Ku         : float – upper attachment (e.g. 0.06)
        avg_recovery : float – average recovery rate
        rho        : float – asset correlation
        p          : float – marginal default probability
        discount   : float – discount factor
        I          : int   – number of loans
        nu         : float – degrees of freedom of the t-distribution

    Returns:
        price      : float – discounted cumulative price of the tranche
    """

    # Step 1: Normalize attachment points
    d = Kd / (1 - avg_recovery)
    u = Ku / (1 - avg_recovery)

    # Step 2: m from 0 to I
    m_vec = np.arange(I + 1)  # (I+1,)

    # Step 3: log binomial coefficients
    log_binom_coeff = gammaln(I + 1) - gammaln(m_vec + 1) - gammaln(I - m_vec + 1)

    # Step 4: Loss profile
    loss_tranche_vec = np.clip((m_vec / I - d) / (u - d), 0, 1)

    # Step 5: Pre-computed constants
    sqrt_rho = np.sqrt(rho)
    sqrt_1_minus_rho = np.sqrt(1 - rho)
    tinv_p = t.ppf(p, nu)

    # Step 6: Define integrand
    def integrand(y):
        py = t.cdf((tinv_p - sqrt_rho * y) / sqrt_1_minus_rho, nu)
        py = np.clip(py, 1e-15, 1 - 1e-15)  # avoid log(0)
        log_probs = (log_binom_coeff +
                     m_vec * np.log(py) +
                     (I - m_vec) * np.log(1 - py))
        weighted_losses = np.exp(log_probs) * loss_tranche_vec
        return np.sum(weighted_losses) * t.pdf(y, nu)

    # Step 7: Integration over y
    expected_loss, _ = quad(integrand, -6, 6, epsrel=1e-10, epsabs=1e-12)

    # Step 8: Discounted price
    price = discount * (1 - expected_loss)

    return price


def price_cumul_tranche_KL_double_t_optimized(kd, ku, recovery, rho, p, discount, I, nu):
    """
    Optimized KL-based pricing of a tranche [kd, ku] under a double t-Student copula.

    Parameters:
        kd (float): lower attachment point (e.g. 0.0)
        ku (float): upper attachment point (e.g. 0.06)
        recovery (float): recovery rate
        rho (float): asset correlation
        p (float): marginal default probability
        discount (float): discount factor
        I (int): portfolio size (number of loans)
        nu (int): degrees of freedom of the t-distribution

    Returns:
        price (float): discounted cumulative tranche price
    """

    # Normalize attachment points in z-space
    d = kd / (1 - recovery)
    u = ku / (1 - recovery)

    # Loss function
    def loss_tr(z):
        return np.clip((z - d) / (u - d), 0, 1)

    # Pre-computed constants
    sqrt_rho = np.sqrt(rho)
    sqrt_1_minus_rho = np.sqrt(1 - rho)
    tinv_p = t.ppf(p, nu)

    # Conditional default probability
    def p_y(y):
        return t.cdf((tinv_p - sqrt_rho * y) / sqrt_1_minus_rho, df=nu)

    # KL divergence and normalization factor
    def KL(z, y, py):
        epsilon = 1e-12
        z = np.clip(z, epsilon, 1 - epsilon)
        py = np.clip(py, epsilon, 1 - epsilon)
        return z * np.log(z / py) + (1 - z) * np.log((1 - z) / (1 - py))

    def C1(z):
        return np.sqrt(I / (2 * np.pi * z * (1 - z)))

    # Integrand for fixed y
    def integrand_z(z, y):
        return C1(z) * np.exp(-I * KL(z, y)) * loss_tr(z)

    # Midpoint integration over y
    dy = 0.1
    y_grid = np.arange(-6, 6 + dy, dy)
    y_mid = 0.5 * (y_grid[:-1] + y_grid[1:])
    dy_vec = np.diff(y_grid)

    integrals = []
    for y in y_mid:
        result, _ = quad(lambda z: integrand_z(z, y), 0, 1, epsrel=1e-8)
        integrals.append(result * t.pdf(y, nu))  # weight by density of y

    # Expected loss from smooth part
    exp_loss = np.sum(dy_vec * np.array(integrals))


    # Final discounted price
    price = discount * (1 - exp_loss)
    return price


def KL(z, y, py):
    epsilon = 1e-12
    z = np.clip(z, epsilon, 1 - epsilon)
    py = np.clip(py, epsilon, 1 - epsilon)
    return z * np.log(z / py) + (1 - z) * np.log((1 - z) / (1 - py))

def price_cumul_tranche_KL_double_t_optimized_2(kd, ku, recovery, rho, p, discount, I, nu):
    d = kd / (1 - recovery)
    u = ku / (1 - recovery)

    def loss_tr(z):
        return np.clip((z - d) / (u - d), 0, 1)

    sqrt_rho = np.sqrt(rho)
    sqrt_1_minus_rho = np.sqrt(1 - rho)
    tinv_p = t.ppf(p, nu)

    # Grid sur y (moins fin, mais suffisant)
    y_grid = np.linspace(-6, 6, 60)
    dy = y_grid[1] - y_grid[0]
    exp_loss = 0.0

    for y in y_grid:
        py = t.cdf((tinv_p - sqrt_rho * y) / sqrt_1_minus_rho, df=nu)

        def integrand(z):
            z = np.clip(z, 1e-12, 1 - 1e-12)
            kl = z * np.log(z / py) + (1 - z) * np.log((1 - z) / (1 - py))
            C = np.sqrt(I / (2 * np.pi * z * (1 - z)))
            return C * np.exp(-I * kl) * loss_tr(z)

        # fixed_quad = quadrature fixe rapide sur z ∈ [0,1]
        integral, _ = fixed_quad(integrand, 0, 1, n=30)
        exp_loss += integral * t.pdf(y, df=nu) * dy

    return discount * (1 - exp_loss)

def price_LHP_t_student_optimized(discounts, recovery, Ku_list, K_d, corr, p, nu):
    """
    Computes LHP tranche prices under a t-Student copula using vectorized integration.

    Parameters:
        discounts : float – discount factor
        recovery  : float – recovery rate
        Ku_list   : array-like – upper detachment points
        K_d       : float – lower detachment point (usually 0)
        corr      : float – correlation (rho)
        p         : float – marginal probability of default
        nu        : int   – degrees of freedom of the t-distribution

    Returns:
        LHP       : np.ndarray – prices of each tranche [K_d, Ku]
    """

    Ku_list = np.array(Ku_list).flatten()
    K = t.ppf(p, nu)
    u_vec = Ku_list / (1 - recovery)
    d = K_d / (1 - recovery)

    sqrt_1mc = np.sqrt(1 - corr)
    sqrt_c = np.sqrt(corr)

    # Define vectorized integrand for all tranches
    def vectorized_integrand(z):
        return compute_all_tranches_contribution(
            z, u_vec, d, K_d, Ku_list, sqrt_1mc, sqrt_c, K, nu, recovery
        )

    # Integrate once for all tranches
    expected_losses = np.array([
        quad(lambda z: vectorized_integrand(z)[i], 0, 1, epsrel=1e-10, epsabs=1e-12)[0]
        for i in range(len(Ku_list))
    ])

    price = discounts * (1 - expected_losses)

    return price[0]


def compute_all_tranches_contribution(z, u_vec, d, K_d, Ku_list, sqrt_1mc, sqrt_c, K, nu, recovery):
    """
    Computes the vector of integrand values for all tranches at a given z.
    """

    if z <= 1e-15 or z >= 1 - 1e-15:
        return np.zeros_like(u_vec)

    try:
        tinv_z = t.ppf(z, nu)
        y_star = (sqrt_1mc * tinv_z - K) / sqrt_c

        tpdf_y_neg = t.pdf(-y_star, nu)
        tpdf_tinv_z = t.pdf(tinv_z, nu)

        if tpdf_tinv_z < 1e-15:
            prob_z = 0
        else:
            prob_z = (sqrt_1mc / sqrt_c) * tpdf_y_neg / tpdf_tinv_z

    except Exception:
        return np.zeros_like(u_vec)

    # Tranche loss: min(max(((1 - R)*z - Kd)/(Ku - Kd), 0), 1)
    numerator = (1 - recovery) * z - K_d
    denominators = Ku_list - K_d

    tranche_losses = np.clip(numerator / denominators, 0, 1)
    return tranche_losses * prob_z


def tranche_prices_double_t_optimized(nu_optimal, I_list, Kd_calibration, Kd_list, Ku_list,
                                      recovery, corr_model, p, discount, flag):
    """
    Compute tranche prices under the double t-copula using different methods (Exact, KL, LHP).

    Parameters:
        nu_optimal : float – optimal degrees of freedom
        I_list     : list or array – portfolio sizes (granularity levels)
        Kd_calibration : float – lower detachment point (often 0)
        Kd_list    : list – not used here, included for signature compatibility
        Ku_list    : list – upper detachment points
        recovery   : float – recovery rate
        corr_model : list – correlation per tranche
        p          : float – marginal default probability
        discount   : float – discount factor
        flag       : int – 0 (Exact/MC), 1 (KL), 2 (LHP)

    Returns:
        tranche_prices : 2D numpy array – shape (len(I_list), len(Ku_list))
    """

    I_list_rounded = np.round(I_list).astype(int)
    Ku_array = np.array(Ku_list)
    n_I = len(I_list)
    n_Ku = len(Ku_array)
    tranche_prices = np.zeros((n_I, n_Ku))

    if flag == 0:
        # Exact pricing (Monte Carlo or HP double t-student)
        for j, I_val in enumerate(I_list_rounded):
            tranche_prices[j, :] = np.array([
                price_cumul_tranche_HP_double_t_student_optimized(
                    Kd_calibration, ku, recovery, corr_model, p, discount, I_val, nu_optimal
                )
                for i, ku in enumerate(Ku_array)
            ])

    elif flag == 1:
        # KL approximation
        for j, I_val in enumerate(I_list_rounded):
            tranche_prices[j, :] = np.array([
                price_cumul_tranche_KL_double_t_optimized_2(
                    Kd_calibration, ku, recovery, corr_model, p, discount, I_val, nu_optimal
                )
                for i, ku in enumerate(Ku_array)
            ])

    elif flag == 2:
        # LHP: vectorized once and replicated for all I
        lhp_prices = [price_LHP_t_student_optimized(
            discount, recovery, ku, Kd_calibration, corr_model, p, nu_optimal
        ) for ku in Ku_array]
        tranche_prices = np.tile(lhp_prices, (n_I, 1))

    return tranche_prices

import numpy as np
from scipy.stats import t, chi2

def t_copula_cumulative(discount, Kd, Ku_list, Nsim, rho, recovery, pd, I, nu):
    """
    Monte Carlo pricing of cumulative tranches under a t-Student copula.

    Parameters:
        discount   : float – discount factor
        Kd         : float – lower attachment point (shared)
        Ku_list    : array-like – upper detachment points for tranches
        Nsim       : int – number of simulations
        rho        : float – asset correlation
        recovery   : float – recovery rate
        pd         : float – marginal default probability
        I          : int – portfolio size (number of loans)
        nu         : int – degrees of freedom for the t-distribution

    Returns:
        price_cum      : array – expected price for each tranche (1D)
        price_cum_up   : array – price + 1σ
        price_cum_down : array – price - 1σ
        var_exp_loss   : array – variance of loss estimator
    """

    Ku = np.array(Ku_list).reshape(1, -1)  # 1 × m
    m = Ku.shape[1]

    # Step 1 – Default threshold
    K = t.ppf(pd, nu)

    # Step 2 – Systemic and idiosyncratic factors (t distributed)
    y = np.random.randn(Nsim, 1) / np.sqrt(chi2.rvs(df=nu, size=(Nsim, 1)) / nu)      # Nsim × 1
    eps = np.random.randn(Nsim, I) / np.sqrt(chi2.rvs(df=nu, size=(Nsim, I)) / nu)    # Nsim × I

    # Step 3 – Latent variables
    X = np.sqrt(rho) * y + np.sqrt(1 - rho) * eps                                    # Nsim × I
    defaults = X < K                                                                 # Nsim × I
    L = ((1 - recovery) / I) * np.sum(defaults, axis=1)                              # Nsim × 1

    # Step 4 – Tranche losses
    L = L.reshape(-1, 1)                                                             # Nsim × 1
    numer = np.maximum(L - Kd, 0) - np.maximum(L - Ku, 0)                            # Nsim × m
    denom = Ku - Kd                                                                 # 1 × m
    payoff = numer / denom                                                          # Nsim × m

    # Step 5 – Statistics
    exp_loss = np.mean(payoff, axis=0)                                              # 1 × m
    var_exp_loss = np.var(payoff, axis=0, ddof=0) / Nsim                            # 1 × m
    std_exp_loss = 1.96*np.sqrt(var_exp_loss)

    # Step 6 – Prices and bounds
    price_cum = discount * (1 - exp_loss)
    price_cum_up = discount * (1 - (exp_loss - std_exp_loss))
    price_cum_down = discount * (1 - (exp_loss + std_exp_loss))

    return price_cum, price_cum_up, price_cum_down, var_exp_loss





###### Vasicek Prices and Model ######


def HP_price_vasicek(discount, I, avg_recovery, Kd, Ku, p, rho):
    """
    Calcule le prix d'une tranche de crédit synthétique.

    Args:
        discount (float): facteur d'actualisation.
        I (int): nombre d'entités dans le portefeuille.
        avg_recovery (float): taux de recouvrement moyen.
        Kd (float): attachement (début de la tranche, e.g. 0.03).
        Ku (float): détachement (fin de la tranche, e.g. 0.07).
        p (float): probabilité de défaut marginale.
        rho (float): corrélation du facteur systémique.

    Returns:
        float: prix de la tranche.
    """

    d = Kd / (1 - avg_recovery)
    u = Ku / (1 - avg_recovery)

    def loss_tranche(m):
        return np.clip((m / I - d) / (u - d), 0.0, 1.0)


    def p_y(y):
        return norm.cdf((norm.ppf(p) - np.sqrt(rho) * y) / np.sqrt(1 - rho))


    def p_m_given_y(m, y):
        return comb(I, m) * p_y(y) ** m * (1 - p_y(y)) ** (I - m)


    expected_loss = 0.0
    for m in range(I + 1):
        integrand = lambda y: p_m_given_y(m, y) * norm.pdf(y)
        integral_result, _ = quad(integrand, -6, 6, limit=10)
        expected_loss += integral_result * loss_tranche(m)


    price = discount * (1 - expected_loss)
    return price


def HP_price_vasicek_fast(discount, I, avg_recovery, Kd, Ku, p, rho, n_y=61):
    """
    Version optimisée du calcul de prix de tranche sous modèle Vasicek.
    """

    d = Kd / (1 - avg_recovery)
    u = Ku / (1 - avg_recovery)

    def loss_tranche(m):
        return np.clip((m / I - d) / (u - d), 0.0, 1.0)


    y_grid = np.linspace(-6, 6, n_y)
    dy = y_grid[1] - y_grid[0]

    p_y_vec = norm.cdf((norm.ppf(p) - np.sqrt(rho) * y_grid) / np.sqrt(1 - rho))
    phi_y = norm.pdf(y_grid)


    comb_table = comb(I, np.arange(I + 1))

    expected_loss = 0.0

    for m in range(I + 1):
        py_m = p_y_vec ** m * (1 - p_y_vec) ** (I - m)
        integrand = comb_table[m] * py_m * phi_y
        integral = np.sum(integrand) * dy
        expected_loss += integral * loss_tranche(m)

    return discount * (1 - expected_loss)

import numpy as np
from scipy.stats import norm
from scipy.integrate import quad


def price_cumul_tranche_KL_vasicek_optimized(kd, ku, recovery, rho, p, disc, I):
    """
    KL approximation of the price of a credit tranche using Vasicek copula.

    Parameters:
        kd (float): lower detachment (e.g., 0.03)
        ku (float): upper detachment (e.g., 0.07)
        recovery (float): recovery rate (e.g., 0.4)
        rho (float): asset correlation (e.g., 0.2)
        p (float): unconditional default probability (e.g., 0.01)
        disc (float): discount factor (e.g., exp(-rT))
        I (int): number of obligors

    Returns:
        float: tranche price
    """

    # Normalize detachment points to loss space
    d = kd / (1 - recovery)
    u = ku / (1 - recovery)

    # Tranche loss function
    def loss_tr(z):
        return np.clip((z - d) / (u - d), 0.0, 1.0)

    # Precompute constants
    sqrt_rho = np.sqrt(rho)
    sqrt_1mrho = np.sqrt(1 - rho)
    norminv_p = norm.ppf(p)

    def p_y(y):
        return norm.cdf((norminv_p - sqrt_rho * y) / sqrt_1mrho)

    def KL(z, y):
        z = np.clip(z, 1e-10, 1 - 1e-10)
        py = np.clip(p_y(y), 1e-10, 1 - 1e-10)
        return z * np.log(z / py) + (1 - z) * np.log((1 - z) / (1 - py))

    def C1(z):
        z = np.clip(z, 1e-10, 1 - 1e-10)
        return np.sqrt(I / (2 * np.pi * z * (1 - z)))

    # Integrand for fixed y
    def integrand_z(z, y):
        return C1(z) * np.exp(-I * KL(z, y)) * loss_tr(z)

    # Outer integral over y ∈ [-6, 6] with normal density
    def integrand_y(y):
        val, _ = quad(integrand_z, 1e-4, 1 - 1e-4, args=(y,), epsrel=1e-6)
        return val * norm.pdf(y)

    # Integrate over y
    EL, _ = quad(integrand_y, -6, 6, epsrel=1e-6)

    # Final price
    price = disc * (1 - EL)
    return price

def price_cumul_tranche_KL_vasicek_fast(kd, ku, recovery, rho, p, disc, I, n_y=61, n_z=100):
    """
    Fast approximation of tranche price using Vasicek copula and KL method.
    """


    d = kd / (1 - recovery)
    u = ku / (1 - recovery)

    y_grid = np.linspace(-6, 6, n_y)
    dy = y_grid[1] - y_grid[0]
    norminv_p = norm.ppf(p)


    p_y = norm.cdf((norminv_p - np.sqrt(rho) * y_grid) / np.sqrt(1 - rho))
    phi_y = norm.pdf(y_grid)


    z_grid = np.linspace(1e-4, 1 - 1e-4, n_z)
    dz = z_grid[1] - z_grid[0]


    z_clipped = np.clip(z_grid, 1e-10, 1 - 1e-10)
    loss_tr = np.clip((z_clipped - d) / (u - d), 0.0, 1.0)
    C1 = np.sqrt(I / (2 * np.pi * z_clipped * (1 - z_clipped)))


    EL = 0.0
    for i, (py, w_y) in enumerate(zip(p_y, phi_y)):
        py = np.clip(py, 1e-10, 1 - 1e-10)
        KL = z_clipped * np.log(z_clipped / py) + (1 - z_clipped) * np.log((1 - z_clipped) / (1 - py))
        integrand_z = C1 * np.exp(-I * KL) * loss_tr
        integral_z = np.sum(integrand_z) * dz
        EL += integral_z * w_y * dy

    return disc * (1 - EL)

def price_LHP_vasicek(discounts, recovery, K_u, K_d, corr, p):
    """
    Computes the price of a tranche [K_d, K_u] using the Vasicek model
    under the Large Homogeneous Portfolio (LHP) assumption.

    Parameters:
        discounts (float): discount factor (e.g., 0.95)
        recovery (float): recovery rate (e.g., 0.4)
        K_u (float): upper attachment point (e.g., 0.06)
        K_d (float): lower attachment point (e.g., 0.03)
        corr (float): asset correlation (rho)
        p (float): marginal probability of default

    Returns:
        LHP (float): price of the tranche under the Vasicek LHP model
    """

    # 1. Default threshold in standard normal space
    K = norm.ppf(p)

    # 2. Normalize attachment points in loss space
    u = K_u / (1 - recovery)
    d = K_d / (1 - recovery)

    # 3. Compute y*(z) transformation
    def y_star(z):
        return ((np.sqrt(1 - corr)) * norm.ppf(z) - K) / np.sqrt(corr)

    # 4. Change-of-variable density: dz/dy
    def prob(z):
        return (np.sqrt(1 - corr) / np.sqrt(corr)) * \
               norm.pdf(-y_star(z)) / norm.pdf(norm.ppf(z))

    # 5. Tranche loss profile function
    def l(z):
        loss = ((1 - recovery) * z - K_d) / (K_u - K_d)
        return np.clip(loss, 0, 1)

    # 6. Integrand: loss(z) * density(z)
    def integrand(z):
        return l(z) * prob(z)

    # 7. Numerical integration over z ∈ [0, 1]
    ELHP, _ = quad(integrand, 0, 1, epsrel=1e-8)

    # 8. Price = discounted notional * (1 - expected loss)
    LHP = discounts * (1 - ELHP)

    return LHP

def price_LHP_vasicek_fast(discounts, recovery, K_u, K_d, corr, p, n_z=1000):
    """
    Fast approximation of LHP Vasicek tranche price using numerical quadrature.
    """

    # Default threshold in Gaussian space
    K = norm.ppf(p)

    # Normalize tranche boundaries
    u = K_u / (1 - recovery)
    d = K_d / (1 - recovery)

    # Grid over z ∈ [ε, 1 - ε] to avoid issues at 0 and 1
    eps = 1e-6
    z_grid = np.linspace(eps, 1 - eps, n_z)
    dz = z_grid[1] - z_grid[0]

    # Compute y*(z)
    y_star = ((np.sqrt(1 - corr)) * norm.ppf(z_grid) - K) / np.sqrt(corr)

    # Compute dz/dy
    numerator = (np.sqrt(1 - corr) / np.sqrt(corr)) * norm.pdf(-y_star)
    denominator = norm.pdf(norm.ppf(z_grid))
    prob_density = numerator / denominator

    # Compute tranche loss function l(z)
    tranche_loss = ((1 - recovery) * z_grid - K_d) / (K_u - K_d)
    tranche_loss = np.clip(tranche_loss, 0, 1)

    # Element-wise multiplication and integration
    integrand = tranche_loss * prob_density
    ELHP = np.sum(integrand) * dz

    # Discounted tranche price
    LHP = discounts * (1 - ELHP)
    return LHP

def gaussian_copula(discount, Kd, KuList, Nsim, rho, recovery, pd, I):
    """
    Gaussian one-factor copula – Monte Carlo pricing of cumulative tranches.

    Parameters:
        discount     : float – discount factor (e.g., 0.95)
        Kd           : float – lower attachment of the tranche (usually 0)
        KuList       : list or np.ndarray – upper detachments [0.03, 0.06, 0.09]
        Nsim         : int – number of Monte Carlo simulations
        rho          : float – asset correlation ρ (in (0,1))
        recovery     : float – recovery rate R
        pd           : float – default probability Q over horizon
        I            : int – number of homogeneous loans in the portfolio

    Returns:
        priceCum     : np.ndarray – cumulative tranche prices [0–Ku[i]]
        priceCumUp   : np.ndarray – upper bound (+1σ) of MC estimate
        priceCumDown : np.ndarray – lower bound (−1σ) of MC estimate
        varExpLoss   : np.ndarray – variance of MC loss estimator
    """
    KuList = np.array(KuList).flatten()  # Ensure 1D array
    m = len(KuList)

    # 1. Generate latent variables
    K = norm.ppf(pd)                           # default threshold
    Y = np.random.randn(Nsim, 1)               # systemic factor (N×1)
    Eps = np.random.randn(Nsim, I)             # idiosyncratic shocks (N×I)
    X = np.sqrt(rho) * Y + np.sqrt(1 - rho) * Eps  # latent variables (N×I)

    # 2. Defaults and portfolio loss
    Default = X < K                            # boolean (N×I)
    L = (1 - recovery) / I * np.sum(Default, axis=1)  # loss per scenario (N,)

    # 3. Compute tranche payoffs
    Ku = KuList.reshape(1, -1)                 # shape (1, m)
    L_mat = L.reshape(-1, 1)                   # shape (N, 1)
    numer = np.maximum(L_mat - Kd, 0) - np.maximum(L_mat - Ku, 0)  # (N, m)
    denom = Ku - Kd                            # (1, m)
    payoff = numer / denom                    # (N, m)

    # 4. Monte Carlo estimates
    expLoss = np.mean(payoff, axis=0)          # (m,)
    varExpLoss = np.var(payoff, axis=0, ddof=0) / Nsim  # variance of estimator
    stdExpLoss = 1.96*np.sqrt(varExpLoss)

    # 5. Discounted prices and confidence bounds
    priceCum = discount * (1 - expLoss)
    priceCumUp = discount * (1 - (expLoss - stdExpLoss))  # +1σ
    priceCumDown = discount * (1 - (expLoss + stdExpLoss))  # −1σ

    return priceCum, priceCumUp, priceCumDown, varExpLoss

def tranche_prices_vasicek(I_list, Kd_calibration, Kd_list, Ku_list,
                           recovery, corr_model, p, discount, flag):
    """
    Computes tranche prices [Kd, Ku] under Vasicek copula model.

    Parameters:
        I_list         : array-like – list of exposure sizes
        Kd_calibration : float – lower bound (fixed for all)
        Kd_list        : unused (for compatibility)
        Ku_list        : list of upper bounds for tranches
        recovery       : float – recovery rate
        corr_model     : float – correlation (rho)
        p              : float – marginal default probability
        discount       : float – discount factor
        flag           : int – method: 0=Exact, 1=KL, 2=LHP

    Returns:
        tranche_prices : ndarray [len(I_list), len(Ku_list)]
    """

    I_rounded = np.round(I_list).astype(int)
    Ku_array = np.array(Ku_list)
    n_I = len(I_list)
    n_Ku = len(Ku_array)
    tranche_prices = np.zeros((n_I, n_Ku))

    if flag == 0:
        # Exact pricing
        for j, I_val in enumerate(I_rounded):
            tranche_prices[j, :] = [
                HP_price_vasicek_fast(discount, I_val, recovery, Kd_calibration, ku, p, corr_model) for i, ku in enumerate(Ku_array)
            ]

    elif flag == 1:
        # KL approximation
        for j, I_val in enumerate(I_rounded):
            tranche_prices[j, :] = [
                price_cumul_tranche_KL_vasicek_fast(Kd_calibration, ku, recovery, corr_model, p, discount, I_val) for i, ku in enumerate(Ku_array)
            ]

    elif flag == 2:
        # LHP approximation
        for j in range(n_I):
            tranche_prices[j, :] = [
                price_LHP_vasicek_fast(
                    discount, recovery, ku, Kd_calibration, corr_model, p
                ) for i, ku in enumerate(Ku_array)
            ]

    return np.array(tranche_prices)




###### Calibration Functions ######


def calibration_freedom_t_student_LHP(discounts, corr_mkt, Kd_calibration, Ku_list,
                                       recovery, p, display_flag='off'):
    nu_list = np.linspace(2,25,24)
    MSE_list = np.full_like(nu_list, np.inf, dtype=float)
    MSE_min = np.inf
    idx_min = 0
    rho_model = 0.3

    price_equity_mkt = price_LHP_vasicek(discounts, recovery, Ku_list[0], Kd_calibration, corr_mkt[0], p)
    for i, nu in enumerate(nu_list):
        def equation(rho):
            return price_LHP_t_student_optimized(discounts, recovery, Ku_list, Kd_calibration, rho, p, nu) - price_equity_mkt

        try:
            rho_star = brentq(equation, 1e-4, 0.999, xtol=1e-4)
        except Exception:
            continue

        mse = 0
        for j in range(1, len(Ku_list)):
            price_t = price_LHP_t_student_optimized(discounts, recovery, Ku_list[j], Kd_calibration, rho_star, p, nu)
            price_v = price_LHP_vasicek(discounts, recovery, Ku_list[j],
                                        Kd_calibration, corr_mkt[j], p)
            mse += (price_t - price_v)**2

        MSE_list[i] = mse
        if mse < MSE_min:
            MSE_min = mse
            idx_min = i
            rho_model = rho_star

    return MSE_min, rho_model, idx_min, nu_list, MSE_list

def compute_for_nu_LHP(args):
    nu, discounts, recovery, Ku_list, Kd_calibration, p, price_equity_mkt, vasicek_prices = args
    try:
        def equation(rho):
            return price_LHP_t_student_optimized(
                discounts, recovery, Ku_list[0], Kd_calibration, rho, p, nu
            ) - price_equity_mkt

        rho_star = brentq(equation, 1e-4, 0.999, xtol=1e-4)

        mse = sum(
            (price_LHP_t_student_optimized(
                discounts, recovery, Ku_list[j], Kd_calibration, rho_star, p, nu
            ) - vasicek_prices[j]) ** 2
            for j in range(1, len(Ku_list))
        )

        return mse, rho_star, nu
    except Exception:
        return np.inf, None, nu


def calibration_freedom_t_student_LHP_optimized(discounts, corr_mkt, Kd_calibration, Ku_list,
                                                recovery, p, display_flag='off'):
    nu_list = np.linspace(2, 25,24)


    vasicek_prices = [
        price_LHP_vasicek(discounts, recovery, Ku, Kd_calibration, rho, p)
        for Ku, rho in zip(Ku_list, corr_mkt)
    ]
    price_equity_mkt = vasicek_prices[0]


    args_list = [
        (nu, discounts, recovery, Ku_list, Kd_calibration, p, price_equity_mkt, vasicek_prices)
        for nu in nu_list
    ]

    with ProcessPoolExecutor() as executor:
        results = list(tqdm(executor.map(compute_for_nu_LHP, args_list), total=len(args_list)))

    MSE_list = np.array([r[0] for r in results])
    rho_list = [r[1] for r in results]

    idx_min = int(np.argmin(MSE_list))

    return MSE_list[idx_min], rho_list[idx_min], idx_min, nu_list, MSE_list


def calibration_freedom_t_student_KL(discounts, corr_mkt,I, Kd_calibration, Ku_list,
                                       recovery, p, display_flag='off'):
    nu_list = np.linspace(2,25,24)
    MSE_list = np.full_like(nu_list, np.inf, dtype=float)
    MSE_min = np.inf
    idx_min = 0
    rho_model = 0.3


    price_equity_mkt = price_LHP_vasicek(discounts, recovery, Ku_list[0], Kd_calibration, corr_mkt[0], p)

    print(price_equity_mkt)

    for i, nu in enumerate(nu_list):
        def equation(rho):
            return  price_cumul_tranche_KL_double_t_optimized(Kd_calibration, Ku_list[0],
                                                              recovery, rho, p,
                                                              discounts, I, nu) - price_equity_mkt

        try:
            rho_star = brentq(equation, 1e-4, 0.999, xtol=1e-4)
            print(rho_star)
        except Exception:
            continue

        mse = 0
        for j in range(1, len(Ku_list)):

            price_t = price_cumul_tranche_KL_double_t_optimized(Kd_calibration, Ku_list[j],
                                                                     recovery, rho_star, p,
                                                                     discounts, I, nu)

            price_v = price_LHP_vasicek(discounts, recovery, Ku_list[j],
                                        Kd_calibration, corr_mkt[j], p)
            mse += (price_t - price_v)**2

        MSE_list[i] = mse
        if mse < MSE_min:
            MSE_min = mse
            idx_min = i
            rho_model = rho_star

    return MSE_min, rho_model, idx_min, nu_list, MSE_list


def compute_for_nu_KL(args):
    nu, Kd_calibration, Ku_list, recovery, p, discounts, I, price_equity_mkt, vasicek_prices = args
    try:
        def equation(rho):
            return price_cumul_tranche_KL_double_t_optimized_2(
                Kd_calibration, Ku_list[0], recovery, rho, p, discounts, I, nu
            ) - price_equity_mkt

        rho_star = brentq(equation, 1e-4, 0.999, xtol=1e-4)

        mse = sum(
            (price_cumul_tranche_KL_double_t_optimized_2(
                Kd_calibration, Ku_list[j], recovery, rho_star, p, discounts, I, nu
            ) - vasicek_prices[j]) ** 2
            for j in range(1, len(Ku_list))
        )

        return mse, rho_star, nu
    except Exception as e:
        print(f"[ERROR] ν = {nu}: {e}")
        return np.inf, None, nu

def calibration_freedom_t_student_KL_optimized(discounts, corr_mkt, I, Kd_calibration, Ku_list,
                                               recovery, p, display_flag='off'):
    nu_list = np.linspace(2, 25, 24)


    vasicek_prices = [
        price_LHP_vasicek(discounts, recovery, Ku, Kd_calibration, rho, p)
        for Ku, rho in zip(Ku_list, corr_mkt)
    ]
    price_equity_mkt = vasicek_prices[0]


    args_list = [
        (nu, Kd_calibration, Ku_list, recovery, p, discounts, I, price_equity_mkt, vasicek_prices)
        for nu in nu_list
    ]


    with ProcessPoolExecutor(max_workers=1) as executor:
        results = list(tqdm(executor.map(compute_for_nu_KL, args_list), total=len(args_list)))

    MSE_list = np.array([r[0] for r in results])
    rho_list = [r[1] for r in results]
    idx_min = int(np.argmin(MSE_list))

    return MSE_list[idx_min], rho_list[idx_min], idx_min, nu_list, MSE_list


###### Plot functions ######


def plot_mse(nu_values, mse_list, title):
    """
    Plot the Mean Squared Error (MSE) as a function of ν values, with minimum highlighted.

    Parameters:
        nu_values (array-like): List or array of ν (degrees of freedom) tested.
        mse_list (array-like): Corresponding list of MSE values.
        title (str): Plot title (optional).
    """
    nu_values = np.array(nu_values)
    mse_list = np.array(mse_list)

    # Identify minimum
    idx_min = np.argmin(mse_list)
    nu_star = nu_values[idx_min]
    mse_min = mse_list[idx_min]

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(nu_values, mse_list, marker='o', linestyle='-', linewidth=2,
             label='MSE(ν)', color='steelblue')

    # Highlight minimum
    plt.plot(nu_star, mse_min, 'ro', markersize=8,
             label=f'Minimum: ν = {nu_star}, MSE = {mse_min:.2e}')
    plt.axvline(x=nu_star, linestyle='--', color='red', alpha=0.6)

    # Labels and grid
    plt.title(title, fontsize=13)
    plt.xlabel('Degrees of Freedom ν', fontsize=11)
    plt.ylabel('Mean Squared Error (MSE)', fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()



def plot_tranche_prices_vs_I(I_list, HP_prices, KL_prices, LHP_prices, Ku_list):
    """
    Create subplots of tranche prices (HP, KL, LHP) vs number of exposures I for 3 tranches.

    Parameters:
        I_list (array): List of number of exposures.
        HP_prices (array): Matrix of HP prices [n_I × n_tranches].
        KL_prices (array): Matrix of KL prices [n_I × n_tranches].
        LHP_prices (array): Matrix of LHP prices [1 × n_tranches] or [n_set × n_tranches].
        Ku_list (list): List of upper attachment points for each tranche.
    """
    tranche_names = ['Equity', 'Junior Mezzanine', 'Senior Mezzanine']
    colors = ['magenta', 'green', 'navy']

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

    for k in range(3):
        ax = axes[k]

        # Plot HP and KL curves
        ax.plot(I_list, HP_prices[:, k], 'o-', label='HP (t-Student)',
                linewidth=2, color=colors[0])
        ax.plot(I_list, KL_prices[:, k], 's--', label='KL (t-Student)',
                linewidth=2, color=colors[1])

        # Horizontal LHP line
        ax.axhline(y=LHP_prices[0, k] if LHP_prices.ndim == 2 else LHP_prices[k],
                   linestyle='--', color=colors[2],
                   linewidth=2, label='LHP (Vasicek)')

        # Axis formatting
        ax.set_xscale('log')
        ax.set_xlabel('Number of Exposures I (log scale)', fontsize=11)
        ax.set_title(f'{tranche_names[k]} Tranche\n[0–{int(100 * Ku_list[k])}%]', fontsize=12)
        ax.grid(True, which='both', linestyle='--', alpha=0.5)
        ax.legend(loc='best', fontsize=9)

    axes[0].set_ylabel('Price (% of Notional)', fontsize=11)
    plt.suptitle('CDO Tranche Prices vs Number of Exposures (Double t-Student Model)', fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


def price_tranches_over_I(
    I_list,                    # list or array of I values
    copula_type,               # 0 = Gaussian, 1 = t-Student
    discount,
    Kd,
    Ku_list,
    Nsim,
    rho,
    recovery,
    pd,
    nu=None                    # required only for t-Student
):
    """
    Computes cumulative tranche prices for a list of I values using Gaussian or t-Student copula.

    Parameters:
        I_list (list or array): List of exposure counts to test (I values).
        copula_type (int): 0 for Gaussian copula, 1 for t-Student copula.
        discount (float): Discount factor.
        Kd (float): Lower attachment point.
        Ku_list (list): Upper detachment points.
        Nsim (int): Number of Monte Carlo scenarios.
        rho (float): Asset correlation.
        recovery (float): Recovery rate.
        pd (float): Default probability.
        nu (float, optional): Degrees of freedom (only for t-Student copula).

    Returns:
        Tuple of 3 arrays [n_I × n_tranches]: price, price_up, price_down
    """

    n_I = len(I_list)
    n_tranches = len(Ku_list)

    price_matrix = np.zeros((n_I, n_tranches))
    price_up_matrix = np.zeros_like(price_matrix)
    price_down_matrix = np.zeros_like(price_matrix)

    for idx, I in enumerate(I_list):
        if copula_type == 0:
            # Gaussian copula
            price, price_up, price_down, _ = gaussian_copula(
                discount, Kd, Ku_list, Nsim, rho, recovery, pd, I
            )
        elif copula_type == 1:
            if nu is None:
                raise ValueError("Degrees of freedom 'nu' must be provided for t-Student copula.")
            price, price_up, price_down, _ = t_copula_cumulative(
                discount, Kd, Ku_list, Nsim, rho, recovery, pd, I, nu
            )
        else:
            raise ValueError("copula_type must be 0 (Gaussian) or 1 (t-Student)")

        price_matrix[idx, :] = price
        price_up_matrix[idx, :] = price_up
        price_down_matrix[idx, :] = price_down

    return price_matrix, price_up_matrix, price_down_matrix


def subplot_tranche_comparison(
    I_list,
    Ku_list,
    hp_prices,
    copula_prices,
    copula_prices_up,
    copula_prices_down,
    hp_flag,       # 0 = Gaussian HP, 1 = t-Student HP
    copula_flag,   # 0 = Gaussian MC, 1 = t-Student MC
    lhp_prices=None,  # New: vector of LHP prices for each tranche
    lhp_flag=None     # New: 0 = Vasicek, 1 = t-Student
):
    """
    Plots HP vs Copula prices (with CI bands) across I for each tranche,
    and optionally adds the LHP asymptotic limit (Vasicek or t-Student).
    """

    colors = {
        'hp': 'magenta' if hp_flag else 'blue',
        'mc': 'darkgreen' if copula_flag else 'black',
        'band': (0.8, 1.0, 0.8) if copula_flag else (0.85, 0.85, 0.85),
        'lhp': 'red' if lhp_flag == 1 else 'orange'
    }

    n_tranches = len(Ku_list)
    fig, axs = plt.subplots(1, n_tranches, figsize=(5 * n_tranches, 4), sharey=True)

    if n_tranches == 1:
        axs = [axs]

    for k in range(n_tranches):
        ax = axs[k]
        x = np.array(I_list)

        # Plot confidence band
        x_fill = np.concatenate([x, x[::-1]])
        y_fill = np.concatenate([copula_prices_up[:, k], copula_prices_down[::-1, k]])
        ax.fill_between(x_fill, y_fill, facecolor=colors['band'], alpha=0.5, label='Monte Carlo CI')

        # Plot Monte Carlo mean
        ax.plot(x, copula_prices[:, k], 's--', linewidth=2, color=colors['mc'], label='Monte Carlo mean')

        # Plot HP prices
        ax.plot(x, hp_prices[:, k], 'o-', linewidth=2, color=colors['hp'], label='HP price')

        # Plot LHP asymptotic line (optional)
        if lhp_prices is not None:
            ax.axhline(lhp_prices[1,k], color=colors['lhp'], linestyle=':', linewidth=2,
                       label='LHP (∞, ' + ('t-Student' if lhp_flag else 'Vasicek') + ')')

        # Format
        ax.set_xscale('log')
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.set_xlabel('Number of Exposures $I$ (log scale)', fontsize=11)
        if k == 0:
            ax.set_ylabel('Tranche Price (% of Notional)', fontsize=11)
        ax.set_title(f'Tranche $[0\\%, {100 * Ku_list[k]:.0f}\\%]$', fontsize=12)
        ax.legend()

    fig.suptitle('Tranche Prices vs Portfolio Size – HP vs Monte Carlo Copula', fontsize=14)
    fig.tight_layout()
    plt.show()



def score_function(hp_prices, copula_prices_down, copula_prices_up):
    """
    Calcule le taux d'inclusion des prix HP dans les intervalles de confiance Monte Carlo.

    Parameters:
        hp_prices: np.ndarray [n_I x n_tranches]
        copula_prices_down: np.ndarray [n_I x n_tranches]
        copula_prices_up: np.ndarray [n_I x n_tranches]
        verbose (bool): affiche les erreurs tranche par tranche

    Returns:
        coverage_count: int, nombre de cas couverts
        total: int, nombre total de cas
        coverage_ratio: float, pourcentage de couverture
    """

    inside = (hp_prices >= copula_prices_down) & (hp_prices <= copula_prices_up)
    coverage_count = np.sum(inside,axis=1)
    total = len(hp_prices)
    coverage_ratio = coverage_count / total

    return coverage_count, total, coverage_ratio