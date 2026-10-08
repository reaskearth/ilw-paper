"""Reproduce the landfall likelihood-ratio figure in the Supporting Information.

Comola et al., "Climate-Conditioned Catastrophe Modeling for Dynamic Risk
Assessment" (Research Square rs-9124834/v2). The figure asks whether the observed record
of continental US (CONUS) tropical cyclone landfalls is more likely under seasonally
forced event sets (SEAS5, initialized in mid-June or mid-April) than under static event
sets pooled over many seasons (ERA5, or the same SEAS5 initialization).

Panels:

(a) total log10 likelihood ratio over 1981-2024, for each forcing x static baseline,
    at three aggregation levels (CONUS, 7 regions, 72 gates) and three intensity groups
    (tropical storm and stronger, Cat 1+, Cat 3+);
(b) the one-sided Diebold-Mariano p-value for equal expected log score per season, with
    a Newey-West variance (3-season lag) and the Harvey-Leybourne-Newbold small-sample
    correction (t distribution, n - 1 degrees of freedom);
(c) one cell of (a) season by season (gate level, Cat 1+, SEAS5 June init vs static
    ERA5), with its running total and the observed Cat 1+ landfalls per season;
(d) the landfall gates and regions.

Method. Likelihoods are empirical: the probability of an observed count in a unit and
season is the fraction of simulated seasons that produced it, floored at 1/(N+1). Each
static pool leaves out the season being scored. Each storm is counted once, at its
maximum-intensity landfall among the CONUS gates. For the June initialization, which
tracks storms forming from July 1, the event sets and the observations alike only count
storms forming on or after July 1 (day of year 181, Jan 1 = 0).

Inputs. The paths at the top of the script mirror ``adjust-ylt.py``: the repository's
``DATA/`` folder, with the restricted files under ``DATA/RESTRICTED/``.

- ``DATA/RESTRICTED/Reask_landfall-data.parquet``: per-simulated-season landfall counts
  by gate and Saffir-Simpson category for ERA5 (2,500 simulations per season) and SEAS5
  (1,000 per ensemble member and season; 25 members in 1981-2016, 51 from 2017), binned
  at 18, 33, 43, 50, 58 and 70 m/s;
- ``DATA/Reask_climate-indices.parquet``: the year-ensembles each event set covers;
- ``DATA/Reask_gates.parquet``: the landfall gates and their regions;
- ``DATA/RESTRICTED/Reask_ibtracs-landfall-data.parquet``: the observed IBTrACS (USA
  agency) landfalls at those gates, binned here with the same breaks.

Usage::

    python likelihood-test-figure.py

Writes ``OUTPUT/likelihood-test-figure.png`` (300 dpi) and ``.pdf``, the values in every
cell of (a) and (b) as ``OUTPUT/likelihood-test-figure_data.csv``, and panel (c)'s
series as ``OUTPUT/likelihood-test-figure_panel-c-data.csv``. Cartopy downloads the
Natural Earth coastlines and state boundaries for panel (d) on first use.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import polars as pl

# change this to the correct path to your data
proj_dir = Path("../")

data_dir = proj_dir / "DATA"
restricted_dir = data_dir / "RESTRICTED"

COUNTS_PATH = restricted_dir / "Reask_landfall-data.parquet"
OBS_PATH = restricted_dir / "Reask_ibtracs-landfall-data.parquet"
METRICS_PATH = data_dir / "Reask_climate-indices.parquet"
GATES_PATH = data_dir / "Reask_gates.parquet"

# where the figure and its data tables are written
OUT_DIR = proj_dir / "OUTPUT"
STEM = "likelihood-test-figure"

# the seasons the package holds for both SEAS5 initializations: all are scored, and the
# static pools are built from all of them (less the season being scored)
YEARS = list(range(1981, 2025))
# the published Saffir-Simpson breaks of the event-set counts, in m/s
PACKAGE_BREAKS = [18, 33, 43, 50, 58, 70]

# Cumulative Saffir-Simpson intensity groups, keyed by minimum category code (Cat0 is a
# tropical storm)
INTENSITY_GROUPS = {"Cat0+": 0, "Cat1+": 1, "Cat3+": 3}
AGG_LEVELS = ("CONUS", "region", "gate")

# (key, forcing label, static baseline label) for each comparison, top to bottom
ROWS = [
    ("SEAS5June__ERA5static", "SEAS5\nJune init", "ERA5"),
    ("SEAS5June__SEAS5Junestatic", "SEAS5\nJune init", "SEAS5 June"),
    ("SEAS5April__ERA5static", "SEAS5\nApril init", "ERA5"),
    ("SEAS5April__SEAS5Aprilstatic", "SEAS5\nApril init", "SEAS5 April"),
]
# the cell broken down season by season in panel (c)
FOCUS = ("SEAS5June__ERA5static", "Cat1+", "gate")
# smallest p the colour scale of panel (b) resolves
P_FLOOR = 1e-4


# --------------------------------------------------------------------------- #
# Loading, with CBRA v0.30.11
# --------------------------------------------------------------------------- #
def import_cbra():
    """Import CBRA v0.30.11 (``pylt``), checking it bins as the counts were binned."""
    from pylt import io, settings

    breaks = settings.CATEGORIES["Vm"].to("m/s").magnitude
    if not np.allclose(breaks, PACKAGE_BREAKS):
        raise RuntimeError(
            f"pylt's default Saffir-Simpson breaks are {breaks.tolist()} m/s, not the "
            f"{PACKAGE_BREAKS} m/s the package was binned with; use CBRA v0.30.11"
        )
    return io, settings


def load_inputs(io, settings):
    """Gates, event sets and observations, as the scoring below expects them.

    Returns (region_map, gates, datasets, obs): region_map is (gate, region) for the 72
    CONUS gates; gates the same as a GeoDataFrame for the map; datasets maps each event
    set to (per-simulated-season counts, N per season); obs maps the genesis window
    ("Jun" or "full") to (observed counts, 1 per season).
    """
    kw = dict(
        ldf_data_path=COUNTS_PATH, metrics_path=METRICS_PATH, gates_path=GATES_PATH
    )
    gates = io.load_gates(kw["gates_path"], collapse_reversible_gates=True).filter(
        pl.col.basin == "CONUS"
    )
    region_map = gates.select(
        gate=pl.col.gate.cast(pl.String), region=pl.col.region.cast(pl.String)
    )
    conus = region_map["gate"].to_list()
    # three South Florida gates have reversed-direction twins, folded into them here
    reversed_gates = io.get_reversible_gates(kw["gates_path"])
    jun = settings.DOG_LB_SEAS5["MIDJUNE"]

    def tidy(counts):
        return (
            counts.with_columns(
                gate=pl.col.gate.cast(pl.String).replace(reversed_gates)
            )
            .filter(pl.col.gate.is_in(conus))
            .group_by("season", "ensemble", "sample", "gate", "category")
            .agg(pl.col("count").sum())
        )

    def event_set(source, timing, cutoff):
        counts, vye = io.load_and_process_ldfs(
            source,
            YEARS,
            timing=timing,
            day_of_year_cutoff=cutoff,
            cutoff_type="genesis",
            gates_subset=conus + list(reversed_gates),
            return_valid_year_ensembles=True,
            extra_filters=[],  # a fresh list: v0.30.11's default leaks between calls
            **kw,
        )
        spe = settings.SAMPLES_PER_ENSEMBLEYEAR[source]
        n = {
            int(s): int(k) * spe
            for s, k in vye.group_by("season")
            .agg(pl.col.ensemble.n_unique())
            .iter_rows()
        }
        if sorted(n) != YEARS:
            raise ValueError(f"{source} {timing}: seasons {sorted(n)}")
        return tidy(counts), n

    def observed(cutoff):
        counts, _ = io.load_and_process_ldfs(
            "IBTrACS_USA",
            YEARS,
            trackset_vers="IBTRACS_LINFILLED",
            day_of_year_cutoff=cutoff,
            cutoff_type="genesis",
            gates_subset=conus,
            collapse_reversible_gates=True,
            ldf_data_path=OBS_PATH,
            gates_path=kw["gates_path"],
            metrics_path=kw["metrics_path"],
            extra_filters=[],
            return_valid_year_ensembles=True,
        )
        return tidy(counts), {y: 1 for y in YEARS}

    datasets = {
        "ERA5_Jun": event_set("ERA5", None, jun),
        "ERA5_full": event_set("ERA5", None, None),
        "SEAS5_MIDJUNE": event_set("SEAS5", "MIDJUNE", jun),
        "SEAS5_MIDAPRIL": event_set("SEAS5", "MIDAPRIL", None),
    }
    obs = {"Jun": observed(jun), "full": observed(None)}
    gates_gdf = gates.with_columns(
        region=pl.col.region.cast(pl.String)
    ).st.to_geopandas()
    return region_map, gates_gdf, datasets, obs


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #
def member_unit_counts(counts, level, region_map, min_cat) -> pl.DataFrame:
    """Per-member landfall count in each unit at the given intensity floor.

    Only members with a non-zero count in a unit appear (k >= 1); zero-count
    members are recovered later against the event-set size N.
    Returns (season, ensemble, sample, unit, k).
    """
    df = counts.filter(pl.col.category >= min_cat)
    if level == "CONUS":
        df = df.with_columns(unit=pl.lit("CONUS"))
    elif level == "region":
        df = df.join(region_map, on="gate", how="left").rename({"region": "unit"})
    elif level == "gate":
        df = df.with_columns(unit=pl.col.gate)
    else:
        raise ValueError(level)
    return df.group_by("season", "ensemble", "sample", "unit").agg(
        k=pl.col.count.sum().cast(pl.Int64)
    )


def units_for_level(level, region_map) -> list[str]:
    """List the aggregation units (gates/regions/'CONUS') for a level."""
    if level == "CONUS":
        return ["CONUS"]
    if level == "region":
        return sorted(region_map["region"].unique().to_list())
    return sorted(region_map["gate"].unique().to_list())


def _prob(obs, nz, tot_nz, by_season):
    """Empirical P(count == o) for each obs row, with a 1/(N+1) floor.

    Parameters
    ----------
    obs : pl.DataFrame with (unit, season, o) [+ N already joined as 'N']
    nz : pl.DataFrame of exact-count member tallies
        alt: (unit, season, k, m); null: (unit, k, m)
    tot_nz : pl.DataFrame of members with k >= 1 per group
        alt: (unit, season, tot_nz); null: (unit, tot_nz)
    by_season : bool
        True for the (per-year) alternative, False for the pooled null.

    Returns obs with added columns: p (probability) and floored (bool).
    """
    join_grp = ["unit", "season"] if by_season else ["unit"]

    out = (
        obs.join(
            nz.rename({"k": "o"}),
            on=join_grp + ["o"],
            how="left",
        )
        .join(tot_nz, on=join_grp, how="left")
        .with_columns(
            m=pl.col.m.fill_null(0),
            tot_nz=pl.col.tot_nz.fill_null(0),
        )
    )
    # number of members with exactly o:
    #   o == 0  -> N - (members with k >= 1)
    #   o >= 1  -> m (tally of exact matches)
    n_match = pl.when(pl.col.o == 0).then(pl.col.N - pl.col.tot_nz).otherwise(pl.col.m)
    # floor any count the event set never produced (incl. o == 0 when every
    # member had a landfall) at 1/(N+1) so the log-likelihood stays finite
    floored = n_match == 0
    return out.with_columns(
        floored=floored,
        p=pl.when(floored).then(1.0 / (pl.col.N + 1.0)).otherwise(n_match / pl.col.N),
    ).drop("m", "tot_nz")


def alt_probabilities(muc_alt, n_by_year, obs):
    """Likelihood of each observed (unit, year) count under the per-year (alt) dist."""
    nz = muc_alt.group_by("unit", "season", "k").agg(m=pl.len())
    tot_nz = muc_alt.group_by("unit", "season").agg(tot_nz=pl.len())
    # Every evaluated season must have an event-set size; a null N from the
    # left-join would silently poison p_alt/loglr, so fail fast instead.
    missing = sorted(set(obs["season"].to_list()) - set(n_by_year))
    if missing:
        raise ValueError(f"n_by_year is missing event-set sizes for seasons: {missing}")
    n_df = pl.DataFrame(
        {"season": list(n_by_year), "N": list(n_by_year.values())},
        schema={"season": obs.schema["season"], "N": pl.Int64},
    )
    obs = obs.join(n_df, on="season", how="left")
    return _prob(obs, nz, tot_nz, by_season=True).rename(
        {"p": "p_alt", "floored": "floored_alt", "N": "N_alt"}
    )


def null_probabilities(muc_null, n_null, obs, null_n_by_year=None):
    """Likelihood of each observed (unit, year) count under the pooled (null) dist.

    ``muc_null`` pools all members across the baseline years; ``n_null`` is the
    summed event-set size over those years.

    With ``null_n_by_year`` (event-set size per baseline season), the null is
    leave-one-season-out: season Y is scored against the pool of every baseline season
    except Y, so the static climatology never contains the season it is scoring.
    Seasons outside the baseline are scored against the full pool.
    """
    # _prob divides by N; a zero pooled event-set size (e.g. no baseline years
    # present in null_n_by_year) would produce inf/nan probabilities.
    if n_null <= 0:
        raise ValueError(f"pooled null event-set size must be positive, got {n_null}")
    if null_n_by_year is not None:
        return _null_probabilities_loo(muc_null, n_null, obs, null_n_by_year)
    nz = muc_null.group_by("unit", "k").agg(m=pl.len())
    tot_nz = muc_null.group_by("unit").agg(tot_nz=pl.len())
    obs = obs.with_columns(N=pl.lit(n_null, dtype=pl.Int64))
    return _prob(obs, nz, tot_nz, by_season=False).rename(
        {"p": "p_null", "floored": "floored_null", "N": "N_null"}
    )


def _null_probabilities_loo(muc_null, n_null, obs, null_n_by_year):
    """Leave-one-season-out version of ``null_probabilities``.

    The pooled tallies minus season Y's own give the pool without Y, exactly, so no
    per-season re-pooling is needed.
    """
    season_dtype = obs.schema["season"]
    muc_null = muc_null.with_columns(pl.col.season.cast(season_dtype))
    seasons = obs.select("season").unique()

    # members with exactly k landfalls, and with any, in the pool without season Y
    nz = (
        muc_null.group_by("unit", "k")
        .agg(m_all=pl.len())
        .join(seasons, how="cross")
        .join(
            muc_null.group_by("unit", "season", "k").agg(m_own=pl.len()),
            on=["unit", "season", "k"],
            how="left",
        )
        .select("unit", "season", "k", m=pl.col.m_all - pl.col.m_own.fill_null(0))
    )
    tot_nz = (
        muc_null.group_by("unit")
        .agg(t_all=pl.len())
        .join(seasons, how="cross")
        .join(
            muc_null.group_by("unit", "season").agg(t_own=pl.len()),
            on=["unit", "season"],
            how="left",
        )
        .select("unit", "season", tot_nz=pl.col.t_all - pl.col.t_own.fill_null(0))
    )
    n_loo = {int(y): n_null - null_n_by_year.get(int(y), 0) for y in seasons["season"]}
    if min(n_loo.values()) <= 0:
        raise ValueError("leave-one-season-out pool is empty for some season")
    n_df = pl.DataFrame(
        {"season": list(n_loo), "N": list(n_loo.values())},
        schema={"season": season_dtype, "N": pl.Int64},
    )
    obs = obs.join(n_df, on="season", how="left")
    return _prob(obs, nz, tot_nz, by_season=True).rename(
        {"p": "p_null", "floored": "floored_null", "N": "N_null"}
    )


def observed_counts(obs_counts, level, region_map, min_cat, units, eval_years):
    """Observed count for every (unit, year) in the eval grid, zero-filled.

    Returns (unit, season, o).
    """
    muc = member_unit_counts(obs_counts, level, region_map, min_cat)
    obs = muc.select("unit", "season", o=pl.col.k)
    grid = pl.DataFrame(
        [(u, y) for u in units for y in eval_years],
        schema={"unit": pl.String, "season": obs.schema["season"]},
        orient="row",
    )
    return grid.join(obs, on=["unit", "season"], how="left").with_columns(
        o=pl.col.o.fill_null(0)
    )


def likelihood_ratio(
    obs_probs,
):
    """Combine per-(unit, year) alt/null probabilities into a per-row logLR.

    logLR > 0 favours the (seasonally-forced) alternative.
    """
    return obs_probs.with_columns(
        loglr=(pl.col.p_alt.log() - pl.col.p_null.log()),
    )


def summarize(lr_df, meta):
    """One headline row for a comparison: summed/mean logLR and floor diagnostics."""
    agg = lr_df.select(
        n_obs=pl.len(),
        sum_loglr=pl.col.loglr.sum(),
        mean_loglr=pl.col.loglr.mean(),
        n_floor_alt=pl.col.floored_alt.sum(),
        n_floor_null=pl.col.floored_null.sum(),
    ).to_dicts()[0]
    # log10 likelihood ratio is easier to read as "orders of magnitude"
    agg["sum_log10lr"] = agg["sum_loglr"] / np.log(10)
    return {**meta, **agg}


def _loo_sizes(spec, null_n):
    """Per-season null event-set sizes if ``spec`` asks for a leave-one-out null."""
    if not spec.get("leave_one_out", False):
        return None
    return {y: null_n[y] for y in spec["baseline_years"] if y in null_n}


def run_all(specs, datasets, obs_datasets, region_map):
    """Run every (spec x intensity x level) comparison.

    ``datasets`` maps dataset key -> (counts_df, n_by_year) for the model event
    sets; ``obs_datasets`` maps window -> (counts_df, _) for the observations.

    Returns (summary_df, lr_frames) where ``lr_frames`` keys
    (spec_key, intensity, level) -> the per-(unit, year) logLR frame.
    """
    muc_cache = {}

    def muc(dataset_key, level, mincat):
        ck = (dataset_key, level, mincat)
        if ck not in muc_cache:
            cts = datasets[dataset_key][0]
            muc_cache[ck] = member_unit_counts(cts, level, region_map, mincat)
        return muc_cache[ck]

    rows, lr_frames = [], {}
    for spec in specs:
        alt_n = datasets[spec["alt_dataset"]][1]
        null_n = datasets[spec["null_dataset"]][1]
        obs_cts = obs_datasets[spec["obs_window"]][0]
        for intensity, mincat in INTENSITY_GROUPS.items():
            for level in AGG_LEVELS:
                units = units_for_level(level, region_map)

                alt_muc = muc(spec["alt_dataset"], level, mincat).filter(
                    pl.col.season.is_in(spec["eval_years"])
                )
                null_muc = muc(spec["null_dataset"], level, mincat).filter(
                    pl.col.season.is_in(spec["baseline_years"])
                )
                n_alt = {y: alt_n[y] for y in spec["eval_years"] if y in alt_n}
                n_null = sum(null_n[y] for y in spec["baseline_years"] if y in null_n)

                obs = observed_counts(
                    obs_cts, level, region_map, mincat, units, spec["eval_years"]
                )
                pa = alt_probabilities(alt_muc, n_alt, obs)
                pn = null_probabilities(
                    null_muc,
                    n_null,
                    obs.select("unit", "season", "o"),
                    null_n_by_year=_loo_sizes(spec, null_n),
                )
                lr = likelihood_ratio(
                    pa.join(pn, on=["unit", "season", "o"], how="left")
                )
                lr_frames[(spec["key"], intensity, level)] = lr
                rows.append(
                    summarize(
                        lr,
                        dict(
                            key=spec["key"],
                            alt=spec["alt_label"],
                            null=spec["null_label"],
                            baseline=spec["baseline"],
                            intensity=intensity,
                            level=level,
                            eval_start=min(spec["eval_years"]),
                            eval_end=max(spec["eval_years"]),
                            N_null=n_null,
                        ),
                    )
                )
    return pl.DataFrame(rows), lr_frames


def year_scores(lr_frame) -> np.ndarray:
    """Yearly total log-LR series D_y (sum of per-unit logLR within each season).

    Aggregating over units within a season is what makes the seasons the
    independent replicates for the significance tests: it absorbs the
    cross-sectional dependence between gates/regions in the same year (one storm
    is assigned to one gate; a busy season lifts them together).
    """
    return (
        lr_frame.group_by("season")
        .agg(D=pl.col.loglr.sum())
        .sort("season")["D"]
        .to_numpy()
    )


def _newey_west_lag(n):
    """Automatic Bartlett-kernel truncation lag (Newey-West rule of thumb)."""
    return int(np.floor(4 * (n / 100) ** (2 / 9)))


def dm_test(lr_frame, hac_lag=None):
    """Diebold-Mariano test for equal log-score, on the yearly log-LR series.

    H0: the static and seasonally-forced forecasts have equal expected log score
    (mean D_y = 0). Returns a signed z (positive => the seasonally-forced model
    scores better), its two-sided normal p-value, and the HAC (Newey-West)
    long-run-variance lag used to allow for serial correlation across seasons.
    """
    from scipy import stats

    d = year_scores(lr_frame)
    n = len(d)
    dbar = float(d.mean())
    dm = d - dbar
    g0 = float(dm @ dm) / n
    if hac_lag is None:
        hac_lag = _newey_west_lag(n)
    lrv = g0
    for k in range(1, min(hac_lag, n - 1) + 1):
        gk = float(dm[k:] @ dm[:-k]) / n
        lrv += 2 * (1 - k / (hac_lag + 1)) * gk
    lrv = max(lrv, 0.0)
    se = np.sqrt(lrv / n) if lrv > 0 else np.nan
    z = dbar / se if se and np.isfinite(se) and se > 0 else np.nan
    p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
    return dict(
        z_dm=float(z),
        p_dm=p,
        mean_D=dbar,
        sum_logLR=float(d.sum()),
        n_years=int(n),
        hac_lag=int(hac_lag),
    )


def specs():
    """Return the four comparisons, each scored on every season against a leave-one-out pool."""

    def spec(key, alt, null, window):
        return dict(
            key=key,
            alt_label=key.split("__")[0],
            null_label=key.split("__")[1],
            alt_dataset=alt,
            null_dataset=null,
            obs_window=window,
            baseline=f"{YEARS[0]}-{YEARS[-1]}",
            baseline_years=YEARS,
            eval_years=YEARS,
            leave_one_out=True,
        )

    return [
        spec("SEAS5June__ERA5static", "SEAS5_MIDJUNE", "ERA5_Jun", "Jun"),
        spec("SEAS5June__SEAS5Junestatic", "SEAS5_MIDJUNE", "SEAS5_MIDJUNE", "Jun"),
        spec("SEAS5April__ERA5static", "SEAS5_MIDAPRIL", "ERA5_full", "full"),
        spec(
            "SEAS5April__SEAS5Aprilstatic", "SEAS5_MIDAPRIL", "SEAS5_MIDAPRIL", "full"
        ),
    ]


def score(datasets, obs, region_map):
    """Score every cell; return (results per cell, panel c's season-by-season series).

    The p-value is the one-sided Diebold-Mariano test with the Harvey-Leybourne-Newbold
    correction: the statistic scaled by sqrt((n - 1) / n) and referred to a t
    distribution with n - 1 degrees of freedom.
    """
    from scipy import stats

    sp = specs()
    summary, lr = run_all(sp, datasets, obs, region_map)
    dm = pl.DataFrame(
        [
            dict(
                key=s["key"],
                intensity=it,
                level=lvl,
                **dm_test(lr[(s["key"], it, lvl)]),
            )
            for s in sp
            for it in INTENSITY_GROUPS
            for lvl in AGG_LEVELS
        ]
    )
    n = len(YEARS)
    res = summary.join(
        dm.select("key", "intensity", "level", "z_dm", "hac_lag"),
        on=["key", "intensity", "level"],
    )
    res = res.with_columns(
        p_dm_hln=pl.Series(
            stats.t.sf(res["z_dm"].to_numpy() * np.sqrt((n - 1) / n), n - 1)
        )
    )
    ser = (
        lr[FOCUS]
        .group_by("season")
        .agg(loglr=pl.col.loglr.sum(), o=pl.col.o.sum())
        .sort("season")
        .with_columns(log10lr=pl.col.loglr / np.log(10))
    )
    cell = res.filter(
        pl.col.key == FOCUS[0], pl.col.intensity == FOCUS[1], pl.col.level == FOCUS[2]
    )["sum_log10lr"].item()
    assert abs(ser["log10lr"].sum() - cell) < 1e-9
    return res, ser


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #
LEVELS = list(AGG_LEVELS)
LEVEL_LABELS = {"CONUS": "CONUS", "region": "Region", "gate": "Gate"}
INTENS = list(INTENSITY_GROUPS)
INTENS_LABELS = {"Cat0+": "TS+", "Cat1+": "Cat 1+", "Cat3+": "Cat 3+"}

INK = "#0b0b0b"
INK_2 = "#52514e"
SURFACE = "#ffffff"
# diverging blue <-> red with a neutral gray midpoint (reference palette); blue = the
# seasonally-forced model favoured, red = the static model favoured
DIVERGING = ["#8f2424", "#e34948", "#f0efec", "#3987e5", "#104281"]
# single-hue sequential for significance; orange, the second sequential hue, so blue
# is never asked to mean both "favoured" and "significant"
SEQUENTIAL = ["#fbeee8", "#f5b89b", "#eb6834", "#b54315", "#6e2708"]


def _cmap(colors, name):
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(name, colors)


def _text_color(rgba):
    r, g, b_ = rgba[:3]
    return SURFACE if 0.2126 * r + 0.7152 * g + 0.0722 * b_ < 0.45 else INK


def _pow10(x):
    """1e-4 -> '10^-4', 2e-4 -> '2x10^-4', as mathtext."""
    e = int(np.floor(np.log10(x)))
    m = round(x / 10**e)
    if m == 10:
        m, e = 1, e + 1
    return f"10$^{{{e}}}$" if m == 1 else f"{m}×10$^{{{e}}}$"


def _fmt_p(p, floor):
    # cell labels stay short: below 0.001 only the order of magnitude is printed (the
    # colour carries the gradation; the data table has the values)
    if p <= floor * 1.0001:
        return f"<{_pow10(floor)}"
    if p < 0.001:
        return "<10$^{-3}$"
    # one significant figure below 0.1, unless that rounds up to 0.1 itself, which is
    # then printed as 0.10 like the values above it
    if p < 0.1 and float(f"{p:.1g}") < 0.1:
        return f"{p:.1g}"
    return f"{p:.2f}"


POS, NEG = DIVERGING[3], DIVERGING[1]


def _grid(res, col, row_keys):
    lut = {
        (r["key"], r["intensity"], r["level"]): r[col]
        for r in res.iter_rows(named=True)
    }
    cols = [(lvl, it) for lvl in LEVELS for it in INTENS]
    return np.array([[lut[(k, it, lvl)] for lvl, it in cols] for k in row_keys])


def _heatmap(ax, mat, cmap, norm, labels, focus_ij, rows, show_rows):
    """One heatmap panel, styled as make_lrt_figure's, with the focus cell outlined."""
    from matplotlib.patches import Rectangle
    from matplotlib.transforms import blended_transform_factory

    nr, nc = mat.shape
    for i in range(nr):
        for j in range(nc):
            rgba = cmap(norm(mat[i, j]))
            ax.add_patch(
                Rectangle(
                    (j + 0.03, i + 0.04), 0.94, 0.92, facecolor=rgba, edgecolor="none"
                )
            )
            ax.text(
                j + 0.5,
                i + 0.52,
                labels[i][j],
                ha="center",
                va="center",
                fontsize=6.5,
                color=_text_color(rgba),
            )
    fi, fj = focus_ij
    ax.add_patch(
        Rectangle(
            (fj + 0.01, fi + 0.02),
            0.98,
            0.96,
            facecolor="none",
            edgecolor=INK,
            lw=1.3,
            zorder=5,
        )
    )
    ax.set_xlim(0, nc)
    ax.set_ylim(nr, 0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0, pad=2)

    ax.set_xticks(np.arange(nc) + 0.5)
    ax.set_xticklabels(
        [INTENS_LABELS[it] for _ in LEVELS for it in INTENS], fontsize=6.5, color=INK_2
    )
    ax.xaxis.tick_top()
    for g, lvl in enumerate(LEVELS):
        x0 = g * len(INTENS)
        ax.text(
            x0 + len(INTENS) / 2,
            -0.78,
            LEVEL_LABELS[lvl],
            ha="center",
            va="bottom",
            fontsize=7,
            color=INK,
            fontweight="bold",
        )
        ax.plot(
            [x0 + 0.1, x0 + len(INTENS) - 0.1],
            [-0.68, -0.68],
            color=INK_2,
            lw=0.6,
            clip_on=False,
        )
        if g:
            ax.axvline(x0, color=SURFACE, lw=3)

    groups = [group for _, group, _ in rows]
    ax.set_yticks(np.arange(nr) + 0.5)
    ax.set_yticklabels(
        [baseline for _, _, baseline in rows] if show_rows else [],
        fontsize=6.5,
        color=INK,
    )
    starts = [i for i in range(nr) if i == 0 or groups[i] != groups[i - 1]]
    tr = blended_transform_factory(ax.figure.transFigure, ax.transData)
    texts = []
    for k, s in enumerate(starts):
        e = starts[k + 1] if k + 1 < len(starts) else nr
        if s:
            ax.axhline(s, color=SURFACE, lw=3)
        if show_rows:
            texts.append(
                ax.text(
                    0.012,
                    (s + e) / 2,
                    groups[s],
                    ha="left",
                    va="center",
                    fontsize=7,
                    color=INK,
                    fontweight="bold",
                    linespacing=1.1,
                    transform=tr,
                    clip_on=False,
                )
            )
            ax.plot(
                [0.098, 0.098],
                [s + 0.12, e - 0.12],
                color=INK_2,
                lw=0.6,
                transform=tr,
                clip_on=False,
            )
    return texts


def _column_headers(fig, ax, group_texts):
    """Head the two row-label columns, centred over each, as the level headers are."""
    from matplotlib.transforms import blended_transform_factory

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    tr = blended_transform_factory(fig.transFigure, ax.transData)
    for header, artists in [
        ("Seasonal\nforcing", group_texts),
        ("Static\nbaseline", [t for t in ax.get_yticklabels() if t.get_text()]),
    ]:
        boxes = [inv.transform(a.get_window_extent(renderer)) for a in artists]
        x0, x1 = min(b[0][0] for b in boxes), max(b[1][0] for b in boxes)
        ax.text(
            (x0 + x1) / 2,
            -0.78,
            header,
            ha="center",
            va="bottom",
            fontsize=7,
            color=INK,
            fontweight="bold",
            linespacing=1.1,
            transform=tr,
            clip_on=False,
        )
        ax.plot(
            [x0, x1], [-0.68, -0.68], color=INK_2, lw=0.6, transform=tr, clip_on=False
        )


# coastline order, west to east, and where each region's name sits (lon, lat, ha)
REGIONS = {
    "WEST_GULF": ("West Gulf", -93.7, 26.5, "center"),
    "CENTRAL_GULF": ("Central Gulf", -89.4, 28.15, "center"),
    "WEST_FLORIDA": ("West\nFlorida", -84.7, 26.0, "center"),
    "SOUTH_FLORIDA": ("South Florida", -81.6, 23.75, "center"),
    "EAST_FLORIDA": ("East\nFlorida", -78.9, 28.2, "left"),
    "MID-ATLANTIC": ("Mid-Atlantic", -75.3, 33.6, "left"),
    "NORTH-EAST": ("Northeast", -70.2, 39.6, "left"),
}
REGION_TONES = (INK, "#85847f")  # alternate along the coast, so each region reads


def _gate_map(fig, rect, gates):
    """Map the 72 CONUS gates, each region in alternating tone and named in place.

    ``gates`` is a GeoDataFrame of the collapsed CONUS gates with a string ``region``
    column. No hue: blue and orange already mean "seasonal favoured" and "significant"
    in the panels above. Gates are separated by small gaps at their ends.
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    assert set(gates["region"]) == set(REGIONS), set(gates["region"])

    proj = ccrs.LambertConformal(central_longitude=-83, standard_parallels=(28, 42))
    pc = ccrs.PlateCarree()
    ax = fig.add_axes(rect, projection=proj)
    ax.set_extent([-99.5, -66.0, 23.0, 46.2], crs=pc)
    land = cfeature.NaturalEarthFeature("physical", "land", "50m")
    states = cfeature.NaturalEarthFeature(
        "cultural", "admin_1_states_provinces_lakes", "50m"
    )
    ax.add_feature(land, facecolor="#ecebe7", edgecolor="#c9c8c3", lw=0.3)
    ax.add_feature(states, facecolor="none", edgecolor=SURFACE, lw=0.4)
    ax.spines["geo"].set_visible(False)

    for k, (region, (name, lon, lat, ha)) in enumerate(REGIONS.items()):
        tone = REGION_TONES[k % 2]
        for geom in gates.loc[gates["region"] == region, "geometry"]:
            xs, ys = geom.xy
            ax.plot(
                xs,
                ys,
                color=tone,
                lw=2.2,
                solid_capstyle="butt",
                transform=pc,
                zorder=3,
            )
            # a surface-coloured dot at each end opens a gap between gates
            ax.plot(
                [xs[0], xs[-1]],
                [ys[0], ys[-1]],
                ls="none",
                marker="o",
                ms=1.0,
                mfc=SURFACE,
                mec="none",
                transform=pc,
                zorder=4,
            )
        ax.text(
            lon,
            lat,
            name,
            ha=ha,
            va="center",
            fontsize=6,
            color=INK,
            linespacing=1.0,
            transform=pc,
            zorder=5,
        )
    return ax


def draw(
    res,
    ser,
    gates,
    rows,
    focus,
    p_col,
    p_title,
    p_floor,
    out_stem,
    panel_c_title="Gate-level Cat 1+ by season, SEAS5 June init vs static ERA5",
    obs_label="Observed\nCat 1+",
):
    """Draw the four-panel figure and write it as ``out_stem``.pdf and a 300 dpi .png.

    Parameters
    ----------
    res : pl.DataFrame
        One row per (key, intensity, level), with ``sum_log10lr`` and ``p_col``.
    ser : pl.DataFrame
        The ``focus`` cell season by season: ``season``, ``log10lr`` and ``o`` (the
        observed count shown under panel c).
    gates : geopandas.GeoDataFrame
        The collapsed CONUS gates, with a string ``region`` column, for panel d.
    rows : list of (key, forcing label, static baseline label)
        The heatmap rows, top to bottom.
    focus : (key, intensity, level)
        The cell broken down in panel c and outlined in panels a and b.
    p_col, p_title, p_floor
        The p-value column for panel b, its title, and the smallest p the colour scale
        resolves (labelled "<floor").
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from matplotlib.transforms import ScaledTranslation, blended_transform_factory

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial",
                "Helvetica",
                "Liberation Sans",
                "Nimbus Sans",
                "DejaVu Sans",
            ],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Liberation Sans",
            "mathtext.it": "Liberation Sans:italic",
            "mathtext.bf": "Liberation Sans:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.edgecolor": INK_2,
            "axes.labelcolor": INK,
            "xtick.color": INK_2,
            "ytick.color": INK_2,
        }
    )
    floor = p_floor
    keys = [k for k, _, _ in rows]
    lr = _grid(res, "sum_log10lr", keys)
    p = _grid(res, p_col, keys)
    focus_ij = (
        keys.index(focus[0]),
        LEVELS.index(focus[2]) * len(INTENS) + INTENS.index(focus[1]),
    )

    lim = float(min(max(np.ceil(np.nanmax(np.abs(lr))), 2.0), 6.0))
    div, seq = _cmap(DIVERGING, "lrt_div"), _cmap(SEQUENTIAL, "lrt_seq")
    n_lr = Normalize(-lim, lim, clip=True)
    n_p = Normalize(0, -np.log10(floor), clip=True)

    # Layout in inches, top to bottom, so the spacing does not depend on the figure
    # height: each heatmap needs room above it for its column headers and the panel
    # title above those.
    W = 7.2
    hm_h, c_h = 1.2, 2.05
    title_gap = 0.16  # panel title baseline above the top of its headers
    hdr_a, hdr_b = 0.5, 0.38  # headers above (a) (two-line row headers), (b)
    stack = [
        ("top", 0.06),
        ("title_a", title_gap),
        ("hdr_a", hdr_a),
        ("a", hm_h),
        ("gap_ab", 0.2),
        ("title_b", title_gap),
        ("hdr_b", hdr_b),
        ("b", hm_h),
        ("gap_bc", 0.22),
        ("title_c", 0.3),
        ("c", c_h),
        ("bottom", 0.4),
    ]
    H = sum(h for _, h in stack)
    y, top_of = H, {}
    for name, h in stack:
        top_of[name] = y
        y -= h
    fig = plt.figure(figsize=(W, H), facecolor=SURFACE)

    hm_left, hm_right, cb_w = 0.235, 0.895, 0.012

    def axes_at(name, left, right):
        y1 = top_of[name] / H
        return fig.add_axes(
            [left, y1 - dict(stack)[name] / H, right - left, dict(stack)[name] / H]
        )

    ax_a = axes_at("a", hm_left, hm_right)
    ax_b = axes_at("b", hm_left, hm_right)
    texts = _heatmap(
        ax_a,
        lr,
        div,
        n_lr,
        [[f"{v:.2f}" for v in r] for r in lr],
        focus_ij,
        rows,
        True,
    )
    _heatmap(
        ax_b,
        -np.log10(p),
        seq,
        n_p,
        [[_fmt_p(v, floor) for v in r] for r in p],
        focus_ij,
        rows,
        True,
    )
    # the row-label column headers once, over (a); (b) has the same rows
    _column_headers(fig, ax_a, texts)

    # colourbars, vertical, to the right of each heatmap; strongest evidence at the
    # top of both (large LR, small p)
    for ax, cmap, norm, ticks, ticklabels, label in [
        (
            ax_a,
            div,
            n_lr,
            np.arange(-lim, lim + 1, 2 if lim > 3 else 1),
            None,
            "log$_{10}$ LR",
        ),
        (
            ax_b,
            seq,
            n_p,
            [-np.log10(t) for t in [1, 0.1, 0.01, 0.001, floor]],
            ["1", "0.1", "0.01", "0.001", _pow10(floor)],
            "p",
        ),
    ]:
        pos = ax.get_position()
        cax = fig.add_axes([hm_right + 0.02, pos.y0 + 0.004, cb_w, pos.height - 0.008])
        cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
        cb.set_ticks(ticks)
        if ticklabels is not None:
            cb.set_ticklabels(ticklabels)
        cb.ax.tick_params(labelsize=6, color=INK_2, labelcolor=INK_2, length=2)
        cb.outline.set_visible(False)
        cb.ax.set_title(label, fontsize=6.5, color=INK_2, pad=6, loc="left")

    # (c) season by season, across the width of the heatmap columns and row labels
    c_right, d_left = 0.585, 0.635
    ax_c = axes_at("c", 0.085, c_right)
    pdf = ser.to_pandas()
    seasons, v = pdf["season"].to_numpy(), pdf["log10lr"].to_numpy()
    ax_c.bar(
        seasons,
        v,
        width=0.72,
        color=[POS if x >= 0 else NEG for x in v],
        linewidth=0,
        zorder=2,
        label="Season",
    )
    cum = np.cumsum(v)
    ax_c.plot(
        seasons,
        cum,
        color=INK,
        lw=1.1,
        zorder=3,
        marker="o",
        ms=1.8,
        label="Running total",
    )
    ax_c.axhline(0, color=INK_2, lw=0.6, zorder=1)
    ax_c.annotate(
        f"{cum[-1]:.2f}",
        (seasons[-1], cum[-1]),
        xytext=(4, 0),
        textcoords="offset points",
        fontsize=6.5,
        color=INK,
        va="center",
        fontweight="bold",
    )
    ax_c.set_xlim(seasons[0] - 0.8, seasons[-1] + 2.2)
    ymin = min(v.min(), cum.min(), 0)
    ymax = max(v.max(), cum.max(), 0)
    pad = 0.08 * (ymax - ymin)
    ax_c.set_ylim(ymin - pad, ymax + pad)
    ax_c.set_ylabel("log$_{10}$ LR", fontsize=7)
    ax_c.set_xticks(np.arange(1985, seasons[-1] + 1, 5))
    ax_c.set_xticks(seasons, minor=True)
    # year labels pushed down to make room for the observed-count strip
    ax_c.tick_params(labelsize=6.5, length=2.5, width=0.6)
    ax_c.tick_params(axis="x", pad=11)
    ax_c.tick_params(which="minor", length=1.5, width=0.4)
    ax_c.grid(axis="y", color="#e6e5e2", lw=0.5, zorder=0)
    for s in ("top", "right"):
        ax_c.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax_c.spines[s].set_linewidth(0.6)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    ax_c.legend(
        handles=[
            Patch(color=POS, lw=0, label="Season, seasonal forcing favoured"),
            Patch(color=NEG, lw=0, label="Season, static climate favoured"),
            Line2D(
                [], [], color=INK, lw=1.1, marker="o", ms=1.8, label="Running total"
            ),
        ],
        loc="upper left",
        frameon=False,
        fontsize=6.5,
        handlelength=1.4,
        borderaxespad=0.2,
    )
    # observed Cat 1+ CONUS landfalls per season, as a strip just under the axis
    strip = blended_transform_factory(
        ax_c.transData, ax_c.transAxes
    ) + ScaledTranslation(0, -5 / 72, fig.dpi_scale_trans)
    for s, o in zip(seasons, pdf["o"].to_numpy()):
        ax_c.text(
            s,
            0,
            f"{int(o)}",
            ha="center",
            va="top",
            fontsize=5.3,
            color=INK_2,
            transform=strip,
        )
    ax_c.text(
        seasons[0] - 1.4,
        0,
        obs_label,
        ha="right",
        va="top",
        fontsize=5.5,
        color=INK_2,
        linespacing=1.0,
        transform=strip,
    )

    # (d) the gates and regions every aggregation level in (a)-(c) is built from
    y1 = top_of["c"] / H
    _gate_map(fig, [d_left, y1 - c_h / H, 0.985 - d_left, c_h / H], gates)

    for name, x, letter, title in [
        ("title_a", 0.012, "a", "Log$_{10}$ likelihood ratio, seasonal vs static"),
        (
            "title_b",
            0.012,
            "b",
            p_title,
        ),
        (
            "title_c",
            0.012,
            "c",
            panel_c_title,
        ),
        ("title_c", d_left - 0.005, "d", "Landfall gates and regions"),
    ]:
        yb = (top_of[name] - 0.13) / H
        fig.text(x, yb, letter, fontsize=9, fontweight="bold", color=INK)
        fig.text(x + 0.028, yb, title, fontsize=7.5, color=INK)

    out_stem = Path(out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext, kw in [("pdf", {}), ("png", {"dpi": 300})]:
        fig.savefig(f"{out_stem}.{ext}", facecolor=SURFACE, **kw)
    plt.close(fig)


def main():
    """Load, score, and write the figure and its data tables."""
    # CBRA v0.30.11 on current libraries raises deprecation warnings; none matter here
    warnings.filterwarnings("ignore")
    io, settings = import_cbra()
    region_map, gates, datasets, obs = load_inputs(io, settings)
    res, ser = score(datasets, obs, region_map)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / STEM
    res.select(
        "key",
        "intensity",
        "level",
        "n_obs",
        "sum_log10lr",
        "z_dm",
        "p_dm_hln",
        "n_floor_alt",
        "n_floor_null",
    ).write_csv(f"{out}_data.csv")
    ser.select("season", "log10lr", observed="o").write_csv(f"{out}_panel-c-data.csv")
    draw(
        res,
        ser,
        gates,
        ROWS,
        FOCUS,
        "p_dm_hln",
        "Diebold–Mariano p-value for equal skill (one-sided)",
        P_FLOOR,
        out,
    )
    print(f"wrote {out}.png, .pdf, _data.csv and _panel-c-data.csv")


if __name__ == "__main__":
    main()
