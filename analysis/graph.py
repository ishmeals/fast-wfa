import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
import numpy as np

# Check for command-line argument
if len(sys.argv) != 2:
    print("Usage: python graph.py <path_to_exp_results.csv>")
    sys.exit(1)

csv_file_path = sys.argv[1]

# Set output directory to the current working directory
output_dir = os.getcwd()

# Load data
try:
    data = pd.read_csv(csv_file_path)
except FileNotFoundError:
    print(f"Error: {csv_file_path} not found. Ensure the experiment results are generated.")
    exit(1)

# Set seaborn aesthetics
sns.set_theme(style="whitegrid")

# Configure matplotlib font settings
plt.rcParams['font.family'] = 'DejaVu Serif'  # Use a font similar to Computer Modern
plt.rcParams['mathtext.fontset'] = 'cm'  # Use Computer Modern for math


def plot_line(data, x, y, hue, title, x_label, y_label, output_file, sample_count, seq_length):
    """Generate a line plot with dynamic title."""
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=data, x=x, y=y, hue=hue, marker="o")
    plt.title(f"{title}\n(Sample Count: {sample_count}, Sequence Length: {seq_length})")
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.yscale("log")
    plt.legend(title=hue)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def plot_heatmap_grid(data, x, y, z, title, x_label, y_label, output_file, sample_count, seq_length):
    """Generate a grid of heatmaps for each algorithm with dynamic title and numeric gradient description."""
    algorithms = data["Algorithm"].unique()
    num_algorithms = len(algorithms)
    cols = 3  # Number of columns in the grid
    rows = -(-num_algorithms // cols)  # Ceiling division for rows

    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5 * rows), constrained_layout=True)

    for i, algorithm in enumerate(algorithms):
        row, col = divmod(i, cols)
        ax = axes[row, col] if rows > 1 else axes[col]
        algo_data = data[data["Algorithm"] == algorithm]
        pivot_table = algo_data.pivot_table(index=y, columns=x, values=z)
        sns.heatmap(
            pivot_table,
            annot=True,
            fmt=".2f",  # Regular numeric format
            cmap="viridis",
            ax=ax,
            cbar_kws={"label": "Time (s)"},  # Add label to color bar
        )
        ax.set_title(algorithm)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

    # Remove empty subplots
    for i in range(len(algorithms), rows * cols):
        row, col = divmod(i, cols)
        ax = axes[row, col] if rows > 1 else axes[col]
        fig.delaxes(ax)

    fig.suptitle(f"{title}\n(Sample Count: {sample_count}, Sequence Length: {seq_length})", fontsize=16)
    plt.savefig(output_file)
    plt.close()


def plot_heatmap_single(data, x, y, z, title, x_label, y_label, output_file):
    """Generate a single heatmap with a gradient description."""
    plt.figure(figsize=(10, 6))
    pivot_table = data.pivot_table(index=y, columns=x, values=z)
    sns.heatmap(
        pivot_table,
        annot=True,
        fmt=".2f",  # Regular numeric format
        cmap="viridis",
        cbar_kws={"label": "Time (s)"},  # Add label to color bar
    )
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


# Define experiments with parameters
experiments = {
    "Error v Time": {
        "type": "line",
        "x": "Error Rate",
        "y": "Avg Time",
        "hue": "Algorithm",
        "title": "Error Rate vs Average Time",
        "x_label": "Error Rate",
        "y_label": "Average Time (s)",
    },
    "Sequence Length v Time": {
        "type": "line",
        "x": "Sequence Length",
        "y": "Avg Time",
        "hue": "Algorithm",
        "title": "Sequence Length vs Average Time",
        "x_label": "Sequence Length",
        "y_label": "Average Time (s)",
    },
    "Gap Opening v Time": {
        "type": "line",
        "x": "Gap Opening Cost",
        "y": "Avg Time",
        "hue": "Algorithm",
        "title": "Gap Opening Cost vs Average Time",
        "x_label": "Gap Opening Cost",
        "y_label": "Average Time (s)",
    },
    "Gap Extension v Time": {
        "type": "line",
        "x": "Gap Extension Cost",
        "y": "Avg Time",
        "hue": "Algorithm",
        "title": "Gap Extension Cost vs Average Time",
        "x_label": "Gap Extension Cost",
        "y_label": "Average Time (s)",
    },
    "Mismatch Penalty v Time": {
        "type": "line",
        "x": "Mismatch Penalty",
        "y": "Avg Time",
        "hue": "Algorithm",
        "title": "Mismatch Penalty vs Average Time",
        "x_label": "Mismatch Penalty",
        "y_label": "Average Time (s)",
    },
    "Joint Error & Length": {
        "type": "heatmap_grid",
        "x": "Error Rate",
        "y": "Sequence Length",
        "z": "Avg Time",
        "title": "Joint Impact of Error Rate and Sequence Length on Time",
        "x_label": "Error Rate",
        "y_label": "Sequence Length",
    },
    "Gap Costs Interaction": {
        "type": "heatmap_grid",
        "x": "Gap Opening Cost",
        "y": "Gap Extension Cost",
        "z": "Avg Time",
        "title": "Interaction of Gap Costs",
        "x_label": "Gap Opening Cost",
        "y_label": "Gap Extension Cost",
    },
    "Sensitivity Analysis": {
        "type": "heatmap_grid",
        "x": "Gap Opening Cost",
        "y": "Gap Extension Cost",
        "z": "Avg Time",
        "title": "Sensitivity Analysis",
        "x_label": "Gap Opening Cost",
        "y_label": "Gap Extension Cost",
    },
    "Error Rate & Complexity": {
        "type": "line",
        "x": "Error Rate",
        "y": "Avg Time",
        "hue": "Algorithm",
        "title": "Error Rate vs Time with Complexity",
        "x_label": "Error Rate",
        "y_label": "Average Time (s)",
    },
    "Length & Gap Penalties": {
        "type": "heatmap_grid",
        "x": "Sequence Length",
        "y": "Gap Opening Cost",
        "z": "Avg Time",
        "title": "Combination of Sequence Length and Gap Penalties",
        "x_label": "Sequence Length",
        "y_label": "Gap Opening Cost",
    },
}

# Generate plots for each experiment
for experiment, params in experiments.items():
    filtered_data = data[data["Experiment"] == experiment]
    if filtered_data.empty:
        continue
    output_file = os.path.join(output_dir, f"{experiment.replace(' ', '_')}.png")
    sample_count = filtered_data["Sample Count"].iloc[0] if "Sample Count" in filtered_data.columns else "N/A"
    seq_length = filtered_data["Sequence Length"].iloc[0] if "Sequence Length" in filtered_data.columns else "N/A"

    if params["type"] == "line":
        plot_line(
            filtered_data,
            x=params["x"],
            y=params["y"],
            hue=params["hue"],
            title=params["title"],
            x_label=params["x_label"],
            y_label=params["y_label"],
            output_file=output_file,
            sample_count=sample_count,
            seq_length=seq_length,
        )
    elif params["type"] == "heatmap_grid":
        plot_heatmap_grid(
            filtered_data,
            x=params["x"],
            y=params["y"],
            z=params["z"],
            title=params["title"],
            x_label=params["x_label"],
            y_label=params["y_label"],
            output_file=output_file,
            sample_count=sample_count,
            seq_length=seq_length,
        )

print(f"Graphs saved to {output_dir}")
