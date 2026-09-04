"""Run exploratory OS/PFS sensitivity analyses for the PRINCE Mito3 cohort."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2


def cox_fit(time: np.ndarray, event: np.ndarray, score: np.ndarray, arm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    def negative_partial_log_likelihood(beta: np.ndarray) -> float:
        linear_predictor = score * beta[0] + arm * beta[1]
        value = 0.0
        for event_time in np.unique(time[event == 1]):
            deaths = (time == event_time) & (event == 1)
            risk = time >= event_time
            value += linear_predictor[deaths].sum() - deaths.sum() * np.log(np.exp(linear_predictor[risk]).sum())
        return -value

    estimate = minimize(negative_partial_log_likelihood, np.zeros(2), method="BFGS").x
    step = 1e-4
    hessian = np.zeros((2, 2))
    baseline = negative_partial_log_likelihood(estimate)
    for first in range(2):
        first_step = np.zeros(2)
        first_step[first] = step
        hessian[first, first] = (
            negative_partial_log_likelihood(estimate + first_step)
            - 2 * baseline
            + negative_partial_log_likelihood(estimate - first_step)
        ) / step**2
        for second in range(first):
            second_step = np.zeros(2)
            second_step[second] = step
            hessian[first, second] = hessian[second, first] = (
                negative_partial_log_likelihood(estimate + first_step + second_step)
                - negative_partial_log_likelihood(estimate + first_step - second_step)
                - negative_partial_log_likelihood(estimate - first_step + second_step)
                + negative_partial_log_likelihood(estimate - first_step - second_step)
            ) / (4 * step**2)
    standard_error = np.sqrt(np.diag(np.linalg.inv(hessian)))
    return estimate, standard_error


def logrank_pvalue(time: np.ndarray, event: np.ndarray, group: np.ndarray) -> float:
    observed_minus_expected = 0.0
    variance = 0.0
    for event_time in np.unique(time[event == 1]):
        risk = time >= event_time
        deaths = (time == event_time) & (event == 1)
        total_at_risk = risk.sum()
        total_deaths = deaths.sum()
        group_at_risk = (risk & (group == 1)).sum()
        group_deaths = (deaths & (group == 1)).sum()
        if total_at_risk > 1:
            fraction = group_at_risk / total_at_risk
            observed_minus_expected += group_deaths - total_deaths * fraction
            variance += total_deaths * fraction * (1 - fraction) * (total_at_risk - total_deaths) / (total_at_risk - 1)
    return float(chi2.sf(observed_minus_expected**2 / variance, 1))


def km_curve(time: np.ndarray, event: np.ndarray, group: np.ndarray) -> pd.DataFrame:
    records: list[dict[str, float | int]] = []
    for group_value, group_name in ((0, "low"), (1, "high")):
        mask = group == group_value
        survival = 1.0
        records.append({"group": group_name, "time_days": 0.0, "survival": survival})
        for event_time in np.unique(time[mask & (event == 1)]):
            at_risk = int(np.sum(mask & (time >= event_time)))
            deaths = int(np.sum(mask & (time == event_time) & (event == 1)))
            survival *= 1 - deaths / at_risk
            records.append({"group": group_name, "time_days": float(event_time), "survival": survival})
    return pd.DataFrame(records)


def analyse_endpoint(frame: pd.DataFrame, time_column: str, event_column: str, endpoint: str) -> tuple[dict[str, float | int | str], pd.DataFrame]:
    time = frame[time_column].to_numpy(float)
    event = frame[event_column].to_numpy(int)
    score = frame["Mito3"].to_numpy(float)
    score_sd = score.std(ddof=1)
    score_per_sd = score / score_sd
    arm = (frame["arm"].eq("C2")).to_numpy(int)
    median = float(np.median(score))
    group = (score >= median).astype(int)
    estimate, standard_error = cox_fit(time, event, score_per_sd, arm)
    cox_p = float(chi2.sf((estimate[0] / standard_error[0]) ** 2, 1))
    result = {
        "endpoint": endpoint,
        "n": int(len(frame)),
        "events": int(event.sum()),
        "mito3_median_split": median,
        "low_group_n": int(np.sum(group == 0)),
        "high_group_n": int(np.sum(group == 1)),
        "cox_hr_per_score_sd_arm_adjusted": float(np.exp(estimate[0])),
        "cox_ci_low": float(np.exp(estimate[0] - 1.96 * standard_error[0])),
        "cox_ci_high": float(np.exp(estimate[0] + 1.96 * standard_error[0])),
        "cox_p": cox_p,
        "median_split_logrank_p": logrank_pvalue(time, event, group),
        "analysis_role": "post hoc exploratory survival sensitivity; no optimized cutoff",
    }
    return result, km_curve(time, event, group).assign(endpoint=endpoint)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    args = parser.parse_args()
    frame = pd.read_csv(args.scores, sep="\t")
    frame = frame.loc[frame["received_nivolumab"].eq("Y")].copy()
    required = {"Mito3", "arm", "os_time", "os_event", "clinical.observation.pfs", "clinical.observation.pfs.event"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    frame["os_event"] = frame["os_event"].astype(float).astype(int)
    frame["clinical.observation.pfs.event"] = frame["clinical.observation.pfs.event"].astype(str).str.lower().eq("true").astype(int)
    results = []
    curves = []
    for time_column, event_column, endpoint in (
        ("os_time", "os_event", "OS"),
        ("clinical.observation.pfs", "clinical.observation.pfs.event", "PFS"),
    ):
        result, curve = analyse_endpoint(frame, time_column, event_column, endpoint)
        results.append(result)
        curves.append(curve)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(args.output, sep="\t", index=False)
    curve_frame = pd.concat(curves, ignore_index=True)
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True)
    colors = {"low": "#0072B2", "high": "#D55E00"}
    for axis, endpoint in zip(axes, ("OS", "PFS")):
        subset = curve_frame.loc[curve_frame["endpoint"].eq(endpoint)]
        for group in ("low", "high"):
            group_data = subset.loc[subset["group"].eq(group)]
            axis.step(group_data["time_days"], group_data["survival"], where="post", label=f"Mito3 {group}", color=colors[group], linewidth=1.6)
        pvalue = next(row["median_split_logrank_p"] for row in results if row["endpoint"] == endpoint)
        axis.set_title(endpoint)
        axis.set_xlabel("Days")
        axis.text(0.04, 0.08, f"log-rank P={pvalue:.3f}", transform=axis.transAxes, fontsize=8)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Kaplan–Meier survival")
    axes[1].legend(frameon=False, fontsize=8, loc="lower left")
    figure.tight_layout()
    figure.savefig(args.figure, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
