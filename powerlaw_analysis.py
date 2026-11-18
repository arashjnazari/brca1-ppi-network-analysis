"""Discrete tail inference for a supplied network degree sequence.

The bootstrap and likelihood comparisons use an independent-observation
reference model. They do not account for dependence between graph degrees
or for selection of a protein-centered subnetwork.
"""

import json
import platform
from importlib.metadata import version
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import powerlaw
from scipy.optimize import minimize_scalar
from scipy.special import log_ndtr, zeta


ALPHA_BOUNDS = (1.000001, 20.0)
PRIMARY_MIN_TAIL = 100
SENSITIVITY_MIN_TAIL = 50
BOOTSTRAPS = 1000
SEED = 2026


def validate_degrees(values):
    """Require a nonempty sample of finite positive integer degrees."""
    values = np.asarray(values, dtype=float)
    if (values.ndim != 1 or values.size == 0 or not np.isfinite(values).all()
            or np.any(values < 1) or np.any(values != np.floor(values))):
        raise ValueError("Degrees must be a nonempty one-dimensional array of positive integers.")
    if np.any(values >= np.iinfo(np.int64).max):
        raise ValueError("Degree exceeds the supported integer range.")
    return values.astype(np.int64)


def estimate_alpha(tail, xmin):
    """Numerically maximize the exact Hurwitz-zeta discrete likelihood."""
    tail = np.asarray(tail)
    if len(tail) < 2 or np.any(tail < xmin) or np.all(tail == xmin):
        raise ValueError("A finite tail MLE needs observations above xmin.")
    mean_log = np.log(tail).mean()
    result = minimize_scalar(
        lambda alpha: alpha * mean_log + np.log(zeta(alpha, xmin)),
        bounds=ALPHA_BOUNDS, method="bounded", options={"xatol": 1e-8},
    )
    if not result.success or not np.isfinite(result.fun):
        raise RuntimeError("Discrete power-law optimization failed.")
    alpha = float(result.x)
    if min(alpha - ALPHA_BOUNDS[0], ALPHA_BOUNDS[1] - alpha) < 1e-4:
        raise RuntimeError("The fitted alpha reached an optimization boundary; expand the bounds.")
    return alpha


def ks_distance(tail, xmin, alpha):
    """Exact discrete KS supremum, including unobserved integers in gaps.

    Between sample jumps, the empirical CDF is constant and the model CDF
    increases. Evaluating both k and k-1 at every observed k therefore finds
    the full supremum without allocating an array up to the largest degree.
    """
    values, counts = np.unique(tail, return_counts=True)
    cumulative = np.cumsum(counts) / len(tail)
    before = (np.cumsum(counts) - counts) / len(tail)
    normalizer = zeta(alpha, xmin)
    model_at = 1 - zeta(alpha, values.astype(float) + 1) / normalizer
    model_before = 1 - zeta(alpha, values.astype(float)) / normalizer
    return float(max(np.max(np.abs(cumulative - model_at)),
                     np.max(np.abs(before - model_before))))


def scan_xmin(degrees, min_tail):
    """Fit every observed eligible cutoff; retain at least min_tail values."""
    degrees = validate_degrees(degrees)
    if min_tail < 2:
        raise ValueError("min_tail must be at least 2.")
    rows = []
    for xmin in np.unique(degrees)[:-1]:
        tail = degrees[degrees >= xmin]
        if len(tail) < min_tail:
            break
        alpha = estimate_alpha(tail, int(xmin))
        rows.append({"xmin": int(xmin), "alpha": alpha, "n_tail": len(tail),
                     "tail_fraction": len(tail) / len(degrees),
                     "ks": ks_distance(tail, int(xmin), alpha)})
    if not rows:
        raise ValueError(f"No nonconstant tail contains at least {min_tail} observations.")
    return rows


def fit_powerlaw(degrees, min_tail=PRIMARY_MIN_TAIL):
    """Select xmin by minimum KS, resolving numerical ties by smaller xmin."""
    return min(scan_xmin(degrees, min_tail), key=lambda row: (row["ks"], row["xmin"]))


def sample_powerlaw(alpha, xmin, size, rng):
    """Exact integer inverse-CDF sampling for the unbounded zeta tail.

    Compare survival probabilities directly, using integer bracketing and
    binary search. Draws are not clipped at the observed maximum or N-1.
    """
    if not alpha > 1 or xmin < 1 or size < 0:
        raise ValueError("Require alpha > 1, xmin >= 1 and a nonnegative sample size.")
    if size == 0:
        return np.empty(0, dtype=np.int64)
    target = 1.0 - rng.random(size)
    normalizer = zeta(alpha, xmin)
    lo = np.full(size, xmin - 1, dtype=np.int64)
    hi = np.full(size, 2 * xmin, dtype=np.int64)
    needs_more = zeta(alpha, hi.astype(float) + 1) / normalizer > target
    while np.any(needs_more):
        if np.any(hi[needs_more] > (2**52) // 2):
            raise OverflowError("An exact draw exceeds the supported floating-point integer range.")
        lo[needs_more] = hi[needs_more]
        hi[needs_more] *= 2
        needs_more = zeta(alpha, hi.astype(float) + 1) / normalizer > target
    while np.any(hi - lo > 1):
        active = hi - lo > 1
        midpoint = lo[active] + (hi[active] - lo[active]) // 2
        below = zeta(alpha, midpoint.astype(float) + 1) / normalizer > target[active]
        indices = np.flatnonzero(active)
        lo[indices[below]] = midpoint[below]
        hi[indices[~below]] = midpoint[~below]
    return hi


def wilson_interval(successes, trials):
    """95% Wilson interval for Monte Carlo exceedance probability."""
    z = 1.959963984540054
    fraction = successes / trials
    denominator = 1 + z * z / trials
    center = (fraction + z * z / (2 * trials)) / denominator
    half = z * np.sqrt(fraction * (1 - fraction) / trials + z * z / (4 * trials**2)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def goodness_of_fit(degrees, fitted, min_tail, bootstraps, seed, progress=None):
    """CSN semiparametric bootstrap, reselecting alpha and xmin each time."""
    if bootstraps < 1:
        raise ValueError("bootstraps must be positive.")
    degrees = validate_degrees(degrees)
    rng = np.random.default_rng(seed)
    body = degrees[degrees < fitted["xmin"]]
    tail_probability = fitted["n_tail"] / len(degrees)
    distances = np.empty(bootstraps)
    for iteration in range(bootstraps):
        n_tail = int(rng.binomial(len(degrees), tail_probability))
        n_body = len(degrees) - n_tail
        synthetic = np.concatenate((
            rng.choice(body, size=n_body, replace=True) if n_body else np.empty(0, dtype=int),
            sample_powerlaw(fitted["alpha"], fitted["xmin"], n_tail, rng),
        ))
        distances[iteration] = fit_powerlaw(synthetic, min_tail)["ks"]
        if progress and ((iteration + 1) % 100 == 0 or iteration + 1 == bootstraps):
            progress(iteration + 1, bootstraps)
    exceedances = int(np.count_nonzero(distances >= fitted["ks"]))
    # A finite simulation must not be presented as an exactly zero probability.
    p_value = (exceedances + 1) / (bootstraps + 1)
    low, high = wilson_interval(exceedances, bootstraps)
    return {"p_bootstrap": p_value, "bootstrap_exceedances": exceedances,
            "bootstraps": bootstraps, "seed": seed,
            "mc_se": float(np.sqrt(p_value * (1 - p_value) / (bootstraps + 1))),
            "mc_ci_low": low, "mc_ci_high": high}, distances


def compare_models(degrees, fitted, analysis):
    """Compare three library models on exactly the same fitted tail.

    The pinned library uses an exact discrete power law and exponential,
    and integrates a continuous lognormal over half-integer bins. R and p_lr
    assess relative likelihood, not absolute goodness of fit.
    """
    if version("powerlaw") != "1.5":
        raise RuntimeError("Install powerlaw==1.5 to reproduce the validated model-comparison API.")
    tail = degrees[degrees >= fitted["xmin"]]
    package = powerlaw.Fit(tail, xmin=fitted["xmin"], discrete=True,
                           estimate_discrete=False, discrete_approximation="round", verbose=False)
    alpha_package = float(package.power_law.alpha)
    alpha_difference = abs(alpha_package - fitted["alpha"])
    if alpha_difference > 0.002:
        raise RuntimeError("Independent discrete MLE implementations disagree.")
    rows = []
    for model, parameters in (("power_law", 1), ("exponential", 1), ("lognormal", 2)):
        distribution = getattr(package, model)
        log_probabilities = np.asarray(distribution.loglikelihoods(tail), dtype=float)
        if distribution.noise_flag or not np.isfinite(log_probabilities).all():
            raise RuntimeError(f"Invalid {model} likelihood fit.")
        log_likelihood = float(log_probabilities.sum())
        if model == "power_law":
            ratio, comparison_p = np.nan, np.nan
        else:
            ratio, comparison_p = package.distribution_compare(
                "power_law", model, normalized_ratio=False, nested=False)
        rows.append({
            "analysis": analysis, "model": model, "xmin": fitted["xmin"],
            "n_tail": fitted["n_tail"], "log_likelihood": log_likelihood,
            "n_parameters": parameters, "AIC": 2 * parameters - 2 * log_likelihood,
            "R_powerlaw_vs_model": ratio, "p_lr": comparison_p,
            "alpha": getattr(distribution, "alpha", np.nan),
            "lambda": getattr(distribution, "Lambda", np.nan),
            "mu": getattr(distribution, "mu", np.nan),
            "sigma": getattr(distribution, "sigma", np.nan) if model == "lognormal" else np.nan,
        })
    table = pd.DataFrame(rows)
    table["delta_AIC"] = table["AIC"] - table["AIC"].min()
    return table, alpha_package, alpha_difference


def plot_tail_fits(degrees, fits, models, label, path):
    """Draw conditional empirical and model CCDFs for the two declared rules."""
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2))
    for ax, fitted in zip(axes, fits):
        tail = degrees[degrees >= fitted["xmin"]]
        values, counts = np.unique(tail, return_counts=True)
        empirical = np.cumsum(counts[::-1])[::-1] / len(tail)
        ax.loglog(values, empirical, "o", ms=4, color="#263b50", label="Empirical tail")
        k = np.arange(fitted["xmin"], int(tail.max()) + 1)
        ax.loglog(k, zeta(fitted["alpha"], k) / zeta(fitted["alpha"], fitted["xmin"]),
                  color="#2166ac", lw=1.8, label="Discrete power law")
        selected = models[models["analysis"] == fitted["analysis"]].set_index("model")
        rate = selected.loc["exponential", "lambda"]
        ax.loglog(k, np.exp(-rate * (k - fitted["xmin"])), color="#bc5145", lw=1.7, label="Discrete exponential")
        mu, sigma = selected.loc["lognormal", ["mu", "sigma"]]
        survival = np.exp(log_ndtr((mu - np.log(k - 0.5)) / sigma)
                          - log_ndtr((mu - np.log(fitted["xmin"] - 0.5)) / sigma))
        ax.loglog(k, survival, color="#28816d", lw=1.7, linestyle="--", label="Rounded discrete lognormal")
        ax.set(xlabel="Degree k", ylabel="P(K >= k | K >= xmin)",
               title=f"{fitted['analysis'].replace('_', ' ').title()}: minimum tail {fitted['min_tail']}\n"
                     f"xmin = {fitted['xmin']}, alpha = {fitted['alpha']:.3f}, n_tail = {fitted['n_tail']}")
        ax.text(0.03, 0.04,
                f"Tail fraction: {fitted['tail_fraction']:.1%}\n"
                f"Bootstrap GOF p = {fitted['p_bootstrap']:.3f} (i.i.d. reference)",
                transform=ax.transAxes, fontsize=9,
                bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "#dddddd"})
        ax.grid(alpha=0.2, which="major")
        ax.legend(fontsize=8, loc="upper right")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"{label}: degree-tail fits in a selected STRING association network", fontsize=13)
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_xmin_sensitivity(scan, fits, label, path):
    """Show the fit-distance/sample-size tradeoff without additional testing."""
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.5), sharex=True)
    for ax, column, ylabel in zip(axes, ("alpha", "ks", "n_tail"),
                                ("MLE alpha", "Discrete KS distance", "Tail observations")):
        ax.plot(scan["xmin"], scan[column], "o-", ms=4, color="#263b50", lw=1.3)
        for fitted, color, style in zip(fits, ("#2166ac", "#bc5145"), ("--", ":")):
            ax.axvline(fitted["xmin"], color=color, linestyle=style, lw=1.5,
                       label=f"{fitted['analysis']}: xmin={fitted['xmin']}")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(fontsize=9)
    for value in (PRIMARY_MIN_TAIL, SENSITIVITY_MIN_TAIL):
        axes[2].axhline(value, color="#8c8c8c", lw=0.8, linestyle="--")
    axes[2].set_xlabel("Lower cutoff xmin")
    fig.suptitle(f"{label}: cutoff sensitivity (candidate tails with at least 50 observations)")
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def run_analysis(degrees, output_dir, label, input_sha256, bootstraps=BOOTSTRAPS, seed=SEED):
    """Export the declared primary analysis and minimum-tail sensitivity check."""
    degrees = validate_degrees(degrees)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    scan = pd.DataFrame(scan_xmin(degrees, SENSITIVITY_MIN_TAIL))
    scan["eligible_primary"] = scan["n_tail"] >= PRIMARY_MIN_TAIL
    scan["eligible_sensitivity"] = scan["n_tail"] >= SENSITIVITY_MIN_TAIL
    scan.to_csv(output_dir / "xmin_sensitivity.csv", index=False)
    fitted_rows, comparison_tables, bootstrap_tables = [], [], []
    for index, (analysis, minimum) in enumerate((("primary", PRIMARY_MIN_TAIL),
                                                ("sensitivity", SENSITIVITY_MIN_TAIL))):
        fitted = fit_powerlaw(degrees, minimum)
        run_seed = seed + index
        print(f"{label} {analysis}: xmin={fitted['xmin']}, alpha={fitted['alpha']:.4f}, "
              f"n_tail={fitted['n_tail']}/{len(degrees)}, KS={fitted['ks']:.4f}", flush=True)
        gof, distances = goodness_of_fit(
            degrees, fitted, minimum, bootstraps, run_seed,
            progress=lambda completed, total: print(f"  bootstrap {completed}/{total}", flush=True),
        )
        comparison, package_alpha, difference = compare_models(degrees, fitted, analysis)
        fitted_rows.append({"analysis": analysis, "min_tail": minimum,
                            "n_observations": len(degrees), **fitted, **gof,
                            "alpha_powerlaw_package": package_alpha,
                            "alpha_abs_difference": difference})
        comparison_tables.append(comparison)
        bootstrap_tables.append(pd.DataFrame({"analysis": analysis,
                                              "replication": np.arange(1, bootstraps + 1),
                                              "ks": distances}))
        print(f"  GOF p={gof['p_bootstrap']:.4f}; exceedances={gof['bootstrap_exceedances']}/{bootstraps}", flush=True)
    summary = pd.DataFrame(fitted_rows)
    models = pd.concat(comparison_tables, ignore_index=True)
    summary.to_csv(output_dir / "fit_summary.csv", index=False)
    models.to_csv(output_dir / "model_comparison.csv", index=False)
    pd.concat(bootstrap_tables, ignore_index=True).to_csv(output_dir / "bootstrap_ks.csv", index=False)
    plot_tail_fits(degrees, fitted_rows, models, label, output_dir / "tail_fits.png")
    plot_xmin_sensitivity(scan, fitted_rows, label, output_dir / "xmin_sensitivity.png")
    metadata = {
        "network": label, "input_sha256": input_sha256, "python": platform.python_version(),
        "packages": {name: version(name) for name in
                     ("numpy", "scipy", "pandas", "matplotlib", "networkx", "powerlaw")},
        "primary_min_tail": PRIMARY_MIN_TAIL, "sensitivity_min_tail": SENSITIVITY_MIN_TAIL,
        "bootstraps_per_analysis": bootstraps, "primary_seed": seed, "sensitivity_seed": seed + 1,
        "alpha_bounds": list(ALPHA_BOUNDS), "powerlaw_support": "integers xmin,...,infinity",
        "ks_convention": "supremum over all integer CDF points; includes observed jumps and preceding integers",
        "bootstrap_model": "independent semiparametric observations; empirical body and exact discrete power-law tail",
        "bootstrap_refits": ["xmin", "alpha"], "p_bootstrap_formula": "(exceedances + 1) / (B + 1)",
        "mc_ci": "95% Wilson interval for the bootstrap exceedance probability, not a confidence interval for alpha",
        "alternative_models": ["exact discrete exponential", "rounded discrete lognormal (half-integer bins)"],
        "interpretation": "Exploratory i.i.d.-reference diagnostics, not a graph-level or disease-mechanism test.",
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Saved inference results to {output_dir}")
    return summary, models
