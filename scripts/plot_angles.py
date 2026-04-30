# plot angle distributions grouped by residue size class
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

df = pd.read_csv("final/angles.tsv", sep="\t")

# make sure angles are in [-180, 180]
df["angle"] = ((df["angle"] + 180) % 360) - 180

size_order = ["Tiny", "Small", "Intermediate", "Large", "Bulky"]
colors = {
    "Tiny":         "#c8c4b5",
    "Small":        "#e8b44c",
    "Intermediate": "#e07840",
    "Large":        "#c93d1d",
    "Bulky":        "#8b0a1a",
}

fig, ax = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor("#bdbdbd")
ax.set_facecolor("#bdbdbd")

xs = np.linspace(-180, 180, 500)

for sz in size_order:
    vals = df[df["size_class"] == sz]["angle"].values
    if len(vals) < 2:
        continue
    kde = gaussian_kde(vals, bw_method=0.15)
    ax.plot(xs, kde(xs), color=colors[sz], linewidth=1.5, label=sz)

# vertical grid lines
for x in range(-150, 200, 50):
    ax.axvline(x, color="blue", linestyle=":", linewidth=0.5, alpha=0.6)

ax.set_xlim(-180, 180)
ax.set_xticks(range(-150, 200, 50))
ax.set_xlabel(u"Angle between adjacent C-\u03b1 -> Centroid vectors [\u00b0]", fontsize=13)
ax.set_ylabel("Norm. Freq. [A.U.]", fontsize=13)
ax.set_title(f"Tripeptide (XRX) in Helix (n = {len(df)})", fontsize=14)
ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="black")
ax.tick_params(axis="y", left=False, labelleft=True)
ax.grid(False)

plt.tight_layout()
plt.savefig("final/angle_plot.png", dpi=300, facecolor="#bdbdbd")
print(f"saved plot to final/angle_plot.png (n={len(df)})")
