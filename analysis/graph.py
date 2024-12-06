import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
file_path = '/mnt/data/exp_results.csv'
data = pd.read_csv(file_path)

# Create output directory for graphs
output_dir = "/mnt/data/graphs"
import os
os.makedirs(output_dir, exist_ok=True)

# Set seaborn aesthetics
sns.set_theme(style="whitegrid")


def plot_line(data, x, y, hue, title, x_label, y_label, output_file):
    """Generate a line plot."""
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=data, x=x, y=y, hue=hue, marker="o")
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend(title=hue)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def plot_heatmap(data, x, y, z, title, x_label, y_label, output_file):
    """Generate a heatmap."""
    plt.figure(figsize=(10, 6))
    pivot_table = data.pivot_table(index=y, columns=x, values=z)
    sns.heatmap(pivot_table, annot=True, fmt=".2f", cmap="viridis")
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


# Filter and plot for each experiment
experiments = {
    "Error Rate vs Time": {
        "type": "line",
        "x": "Error Rate",
        "y": "Avg Time",
        "hue": "Algorithm",
        "title": "Error Rate vs Average Time",
        "x_label": "Error Rate",
        "y_label": "Average Time (s)",
    },
    "Sequence Length vs Time": {
        "type": "line",
        "x": "Sequence Length",
        "y": "Avg Time",
        "hue": "Algorithm",
        "title": "Sequence Length vs Average Time",
        "x_label": "Sequence Length",
        "y_label": "Average Time (s)",
    },
    "Gap Opening Cost vs Time": {
        "type": "line",
        "x": "Gap Opening Cost",
        "y": "Avg Time",
        "hue": "Algorithm",
        "title": "Gap Opening Cost vs Average Time",
        "x_label": "Gap Opening Cost",
        "y_label": "Average Time (s)",
    },
    "Gap Extension Cost vs Time": {
        "type": "line",
        "x": "Gap Extension Cost",
        "y": "Avg Time",
        "hue": "Algorithm",
        "title": "Gap Extension Cost vs Average Time",
        "x_label": "Gap Extension Cost",
        "y_label": "Average Time (s)",
    },
    "Mismatch Penalty vs Time": {
        "type": "line",
        "x": "Mismatch Penalty",
        "y": "Avg Time",
        "hue": "Algorithm",
        "title": "Mismatch Penalty vs Average Time",
        "x_label": "Mismatch Penalty",
        "y_label": "Average Time (s)",
    },
    "Joint Error & Length": {
        "type": "heatmap",
        "x": "Error Rate",
        "y": "Sequence Length",
        "z": "Avg Time",
        "title": "Joint Impact of Error Rate and Sequence Length on Time",
        "x_label": "Error Rate",
        "y_label": "Sequence Length",
    },
    "Gap Costs Interaction": {
        "type": "heatmap",
        "x": "Gap Opening Cost",
        "y": "Gap Extension Cost",
        "z": "Avg Time",
        "title": "Interaction of Gap Costs",
        "x_label": "Gap Opening Cost",
        "y_label": "Gap Extension Cost",
    },
    "Sensitivity Analysis": {
        "type": "heatmap",
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
        "type": "heatmap",
        "x": "Sequence Length",
        "y": "Gap Opening Cost",
        "z": "Avg Time",
        "title": "Combination of Sequence Length and Gap Penalties",
        "x_label": "Sequence Length",
        "y_label": "Gap Opening Cost",
    },
}

for experiment, params in experiments.items():
    filtered_data = data[data["Experiment"] == experiment]
    output_file = os.path.join(output_dir, f"{experiment.replace(' ', '_')}.png")
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
        )
    elif params["type"] == "heatmap":
        plot_heatmap(
            filtered_data,
            x=params["x"],
            y=params["y"],
            z=params["z"],
            title=params["title"],
            x_label=params["x_label"],
            y_label=params["y_label"],
            output_file=output_file,
        )

print(f"Graphs saved to {output_dir}")
