#!/usr/bin/env python3

import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from pathlib import Path

# paths
INPUT_CSV = Path("/maps/jh2589/eping/project_rates.csv")
OUTPUT_DIR = Path("/maps/epr26/ex_ante_forecast_out")
OUTPUT_DIR.mkdir(exist_ok=True)

# data
df = pd.read_csv(INPUT_CSV).dropna()

# titles
axis_labels = {
    "s_exante_rate": "Time-Shifted Historical Match",
    "k_exante_rate": "Recent Project",
    "regional_exante_rate": "Recent Regional",
    "s_expost_rate": ""
}

def make_scatter(counterfactual_cols, filename, observed_col="k_rate", colors=None):
    n = len(counterfactual_cols)
    fig, axs = plt.subplots(1, n, figsize=(6*n, 5), sharey=True)

    all_x = df[counterfactual_cols].values.flatten()
    all_y = df[observed_col].values.flatten()
    all_vals = np.concatenate([all_x, all_y])
    min_val = np.floor(all_vals.min() / 5) * 5
    max_val = np.ceil(all_vals.max() / 5) * 5

    for i, x_col in enumerate(counterfactual_cols):
        ax = axs[i] if n > 1 else axs
        x = df[x_col]
        y = df[observed_col]

        # stats
        mae = np.mean(np.abs(x - y))
        mape = np.mean(np.abs(x - y) / y) * 100

        # Run OLS regression and extract slope and confidence interval
        model = sm.OLS(y, x).fit()
        slope_new = model.params.iloc[0]  # Since there's no intercept, the first parameter is the slope
        conf = model.conf_int(alpha=0.05)  # 95% confidence interval (pandas)
        lower, upper = conf.iloc[0]  # Extract confidence interval for slope

        slope = np.sum(x * y) / np.sum(x * x)
        prediction = (1 - np.sum((y - x)**2) / np.sum((y - np.mean(y))**2))

        # points
        ax.scatter(x, y, color=colors[i])
        mins = min(x.min(), y.min() + 0.5)
        maxs = max(x.max(), y.max() + 0.5)
        line = np.linspace(mins, maxs, 100)
        ax.plot(line, line, linestyle="--", color="gray")
        ax.plot(line, line * slope_new, linestyle=":", linewidth = 1, color=colors[i]) # Plot the estimated slope
        ax.plot(line, line * lower, linestyle=":", linewidth = 0.5, color=colors[i]) # Plot the lower bound of CI of slope
        ax.plot(line, line * upper, linestyle=":", linewidth = 0.5, color=colors[i]) # Plot the upper bound of CI of slope
        ax.fill_between(line, line * lower, line * upper, color=colors[i], alpha=0.25)
        ax.set_xlim(mins, maxs)
        ax.set_ylim(mins, maxs)
        ax.set_aspect("equal", adjustable="box")

        # text
        ax.set_title(axis_labels.get(x_col, x_col), fontsize=16)
        ax.set_xlabel("Counterfactual deforestation rate (%)", fontsize=14)
        if i == 0:
            ax.set_ylabel("Observed deforestation rate (%)", fontsize=14)

        ax.text(
            mins + 0.02*(maxs - mins),
            maxs - 0.1*(maxs - mins),
#           f"MAE: {mae:.3f}\nSlope: {slope_new:.3f} [{lower:.3f}, {upper:.3f}]\nGoodness-of-Fit: {prediction:.2f}",
           f"MAE: {mae:.3f}%\nSlope: {slope_new:.3f}\nGoodness-of-Fit: {prediction:.2f}",
            fontsize=14,
            va="top"
        )

    plt.subplots_adjust(wspace=0.02)
    plt.savefig(OUTPUT_DIR / filename, bbox_inches="tight")
    plt.close()
    print(f"✅ Saved {filename}")

# === plots ===

make_scatter(
    ["regional_exante_rate", "k_exante_rate", "s_exante_rate"],
    "out_figure4a_ex_ante.png",
    colors=["#006CD1", "#40B0A6", "#CDAC60"]
)

make_scatter(
    ["s_expost_rate"],
    "out_figure4b_ex_post.png",
    colors=["#C13C3C"]
)