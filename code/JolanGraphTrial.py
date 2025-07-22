import seaborn as sborn
import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker
import pandas as panda
import numpy as npy
import scipy 
from scipy.stats import gaussian_kde

bacData = panda.read_csv("OutputBacteria.csv")
arcData = panda.read_csv("OutputArchaea.csv")
combinedData = panda.concat([bacData, arcData], ignore_index=True)

combinedData["LogGeoMeanVolume"] = npy.log10(combinedData["GeoMeanVolume"])
combinedData = combinedData.dropna(subset=["LogGeoMeanVolume"])


filtered = combinedData.groupby("Phylum").filter(lambda x: len(x)>30)

global_min = filtered["LogGeoMeanVolume"].min()
global_max = filtered["LogGeoMeanVolume"].max()
bin_edges = npy.linspace(global_min, global_max, 41)

plot = sborn.FacetGrid(filtered, col="Phylum", col_wrap=3, sharex = True, sharey= False)
plot.map_dataframe(sborn.histplot, x="LogGeoMeanVolume", stat="count", color="skyblue", bins = bin_edges)

for ax, phylum in zip(plot.axes.flat, plot.col_names):
    subdf = filtered[filtered["Phylum"] == phylum]
    data = subdf["LogGeoMeanVolume"].values    

    #hist_vals, _ = npy.histogram(data,bins=bin_edges)

    #kde = gaussian_kde(data, bw_method=0.5)
    #x_vals = npy.linspace(data.min(), data.max(), 200)
    #y_vals = kde(x_vals)

    #y_scaled = y_vals * (hist_vals.max() / y_vals.max())
    #ax.plot(x_vals, y_scaled, color = "Red", linewidth =2)

    mean_val = npy.mean(data)
    ax.axvline(mean_val, color = "green", linestyle="--", linewidth=2)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    ax.tick_params(axis="x", which="minor", length= 3, color = "gray")

for ax in plot.axes.flat:
    ax.tick_params(axis="x", labelbottom=True)

plot.fig.suptitle("Prokaryote Log₁₀ Normalised Volume by Phylum", fontsize = 16)
plot.set_axis_labels("Log₁₀ Volume", "Number of Species")
pyplot.xticks(rotation=0)
pyplot.tight_layout()
plot.savefig("Figure_6.pdf", format="pdf")

pyplot.figure(figsize=(10,6))

colour = sborn.color_palette("husl", len(filtered["Phylum"].unique()))

for (phylum, subdf), colour in zip(filtered.groupby("Phylum"), colour):
    data = subdf["LogGeoMeanVolume"].values    

    hist_vals, _ = npy.histogram(data,bins=bin_edges)

    kde = gaussian_kde(data, bw_method=0.5)
    x_vals = npy.linspace(global_min, global_max, 300)
    y_vals = kde(x_vals)
    y_vals /= y_vals.max()
    #y_scaled = y_vals * (hist_vals.max() / y_vals.max())
    pyplot.plot(x_vals, y_vals, label=phylum, linewidth = 2, color=colour)
    pyplot.fill_between(x_vals, y_vals, alpha = 0.4, color=colour)

    mean_val = npy.mean(data)
    pyplot.axvline(mean_val, color = colour, linestyle="--", linewidth=2)

pyplot.gca().xaxis.set_major_locator(ticker.MultipleLocator(1))
pyplot.gca().xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

pyplot.tick_params(axis="x", which="minor", length= 3, color = "gray")

pyplot.xlabel("Log₁₀ Volume")
pyplot.ylabel("Relative Density")
pyplot.title("Normalized KDE of Log₁₀ Normalised Volume by Phylum")
pyplot.grid(True)
pyplot.legend()
pyplot.ylim(0, 1.05)
pyplot.tight_layout()
pyplot.savefig("Figure_7.pdf", format="pdf")