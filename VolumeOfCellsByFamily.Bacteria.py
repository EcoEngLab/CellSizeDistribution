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


filtered = combinedData.groupby("Family").filter(lambda x: len(x)>30)

# Before the loop, define global min/max and bin edges
global_min = filtered["LogGeoMeanVolume"].min()
global_max = filtered["LogGeoMeanVolume"].max()
bin_edges = npy.linspace(global_min, global_max, 36)  # 30 bins, for example
print("Bin edges:", bin_edges)
print("Bin width:", (global_max - global_min) / 35)

plot = sborn.FacetGrid(filtered, col="Family", col_wrap=3, sharex = True, sharey= False)
plot.map_dataframe(sborn.histplot, x="LogGeoMeanVolume", stat="count", color="skyblue", bins = bin_edges)

for ax, family in zip(plot.axes.flat, plot.col_names):
    subdf = filtered[filtered["Family"] == family]
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

plot.fig.suptitle("Prokaryote Log₁₀ Normalised Volume by Family", fontsize = 16)
plot.set_axis_labels("Log₁₀ Volume", "Number of Species")
pyplot.xticks(rotation=0)
pyplot.tight_layout()
plot.fig.subplots_adjust(top=0.95)
plot.savefig("FamilyHist.pdf", format="pdf")

pyplot.figure(figsize=(10, 8))  # Taller figure

colour = sborn.color_palette("husl", len(filtered["Family"].unique()))

for (family, subdf), colour in zip(filtered.groupby("Family"), colour):
    data = subdf["LogGeoMeanVolume"].values    

    kde = gaussian_kde(data, bw_method=0.5)
    x_vals = npy.linspace(global_min, global_max, 300)
    y_vals = kde(x_vals)
    y_vals /= y_vals.max()
    pyplot.plot(x_vals, y_vals, label=family, linewidth=2, color=colour)
    pyplot.fill_between(x_vals, y_vals, alpha=0.4, color=colour)

    mean_val = npy.mean(data)
    pyplot.axvline(mean_val, color=colour, linestyle="--", linewidth=2)

pyplot.gca().xaxis.set_major_locator(ticker.MultipleLocator(1))
pyplot.gca().xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

pyplot.tick_params(axis="x", which="minor", length=3, color="gray")

pyplot.xlabel("Log₁₀ Volume")
pyplot.ylabel("Relative Density")
pyplot.title("Normalized KDE of Log₁₀ Normalised Volume by Family")
pyplot.grid(True)
# Make legend smaller and move it outside the plot
pyplot.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0.)
pyplot.ylim(0, 1.05)
x_margin = 0.2  # Adjust this value as needed
pyplot.xlim(global_min - x_margin, global_max + x_margin)
pyplot.tight_layout(rect=[0, 0, 0.85, 0.95])
pyplot.subplots_adjust(top=0.9)
pyplot.savefig("FamilyKDE.pdf", format="pdf", bbox_inches='tight')
