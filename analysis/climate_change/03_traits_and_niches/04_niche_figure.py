import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import matplotlib

matplotlib.use("MacOSX")
from pathlib import Path

data_path = Path("../modve_data/modve_output/regua")
filling_fname = "a5_niche_range_filling_cc_vs_no_cc_long_v3.csv"
dist_fname = "opt_max_distances_cc_vs_no_cc.csv"

niche_filling = pd.read_csv(data_path / filling_fname)
distances = pd.read_csv(data_path / dist_fname)
# Add ScenarioNames column for plotting
distances["ScenarioNames"] = distances["Scenario"].map({"CC": "Climate change", "No CC": "Baseline"})
# Sort scenarios
distances["ScenarioNames"] = pd.Categorical(
    distances.ScenarioNames, categories=["Baseline", "Climate change"], ordered=True
)
niche_filling["ScenarioNames"] = pd.Categorical(
    niche_filling.ScenarioNames, categories=["Baseline", "Climate change"], ordered=True
)

# =============================================================================
# 0.  CONFIG
# =============================================================================
#COLORS = {"Baseline": "#004aad", "Climate change": "#f7766e"}   # adjust as needed
COLORS = {"Baseline": "#7792CA", "Climate change": "#FAADA8"}   # adjust as needed
SCENARIO_MAP = {"No CC": "Baseline", "CC": "Climate change"}

# =============================================================================
# 1.  DATA PREPARATION  (mirrors the R wrangling)
# =============================================================================
# -- Load / assume niche_filling is already a DataFrame ----------------------
# niche_filling = pd.read_csv("your_data.csv")

niche_filling["ScenarioNames"] = niche_filling["Scenario"].map(SCENARIO_MAP)
scenario_order = ["Baseline", "Climate change"]

# -- FN / PN filling subset --------------------------------------------------
fn_pn_mask = niche_filling["Variable"].str.contains("FN filling|PN filling", regex=True)
fn_pn_data = niche_filling[fn_pn_mask].copy()
fn_pn_data["FillingType"] = fn_pn_data["Variable"].apply(
    lambda x: "FN" if "FN" in x else "PN"
)
fn_pn_data["Variable"] = fn_pn_data["Variable"].str.replace(
    r" FN filling| PN filling", "", regex=True
)
fn_pn_data["Group"] = fn_pn_data["Variable"] + "\n" + fn_pn_data["FillingType"]
fn_pn_data.loc[fn_pn_data["Group"] == "Relative humidity\nFN", "Group"] = "Relative\nhumidity\nFN"
fn_pn_data.loc[fn_pn_data["Group"] == "Relative humidity\nPN", "Group"] = "Relative\nhumidity\nPN"

GROUP_ORDER = [
    "Temperature\nFN",      "Temperature\nPN",
    "Relative\nhumidity\nFN", "Relative\nhumidity\nPN",
    "Light\nFN",            "Light\nPN",
]
fn_pn_data["Group"] = pd.Categorical(fn_pn_data["Group"], categories=GROUP_ORDER, ordered=True)

# -- Geo range filling subset ------------------------------------------------
GEO_VARS_ORIG = ["Vertical range filling", "Horizontal range filling", "3D range filling"]
geo_data = niche_filling[niche_filling["Variable"].isin(GEO_VARS_ORIG)].copy()
GEO_VARS = ["Vertical\nrange filling", "Horizontal\nrange filling", "3D\nrange filling"]
geo_data["Variable"] = geo_data["Variable"].map({"Vertical range filling": "Vertical\nrange filling",
                                                 "Horizontal range filling": "Horizontal\nrange filling",
                                                 "3D range filling": "3D\nrange filling"})
geo_data["Variable"] = pd.Categorical(geo_data["Variable"], categories=GEO_VARS, ordered=True)

# =============================================================================
# 2.  HELPER – unified split-violin function
# =============================================================================
def split_violin(
    ax, data, x, y,
    hue="ScenarioNames",
    order=None,
    palette=COLORS,
    hue_order=None,
    ylabel="",
    xlabel_rotation=0,
    tick_fontsize=14,
    label_fontsize=17,
):
    if hue_order is None:
        hue_order = scenario_order

    sns.violinplot(
        data=data,
        x=x, y=y,
        hue=hue,
        hue_order=hue_order,
        order=order,
        split=True,          # ← the key argument for a split / mirrored violin
        inner="box",         # shows median + IQR inside each half
        palette=palette,
        linewidth=0.8,
        ax=ax,
        legend=False,
        gap=.05,
        inner_kws=dict(box_width=8, color=".3")
    )

    ax.set_xlabel("", fontsize=label_fontsize)
    ax.set(xlabel=None)
    ax.set_ylabel(ylabel, fontsize=label_fontsize)
    ax.tick_params(axis="x", labelsize=label_fontsize, rotation=xlabel_rotation)
    ax.tick_params(axis="y", labelsize=tick_fontsize)
    ax.yaxis.grid(True, linestyle="-", alpha=0.2)
    ax.xaxis.grid(False)
    ax.set_facecolor("white")
    ax.spines[["top", "right"]].set_visible(False)

    # Remove per-axes legend (collected at figure level later)
    legend = ax.get_legend()
    if legend:
        legend.remove()


# =============================================================================
#  SINGLE-VARIABLE SPLIT VIOLIN  (e.g. p_max_T, p_pns_T, p_opt_T ...)
# =============================================================================
def split_violin_single(
        ax, data, y,
        hue="ScenarioNames",
        palette=COLORS,
        hue_order=None,
        ylabel="",
        xlabel="",
        tick_fontsize=14,
        label_fontsize=17,
        width=0.7,
):
    """
    A single split violin (one Baseline half + one CC half, merged at the
    centre, no x category, no x-axis label/ticks).
    """
    if hue_order is None:
        hue_order = scenario_order

    data = data.copy()
    data["_x"] = ""  # constant dummy x -> a single violin position

    ax.yaxis.grid(True, linestyle="-", alpha=0.2)
    ax.set_facecolor("white")
    ax.spines[["top", "right"]].set_visible(False)

    sns.violinplot(
        data=data,
        x="_x", y=y,
        hue=hue,
        hue_order=hue_order,
        split=True,  # halves are mirrored & touch at the centre line
        inner="box",
        palette=palette,
        width=width,  # controls how "fat" the merged violin is
        linewidth=0.8,
        ax=ax,
        legend=False,
        gap=.05,
        inner_kws=dict(box_width=8, color=".3")
    )

    ax.set_xlabel(xlabel, fontsize=label_fontsize)
    ax.set_xticks([])  # no x-axis label/ticks at all
    ax.set_ylabel(ylabel, fontsize=label_fontsize)
    ax.tick_params(axis="y", labelsize=tick_fontsize)

    legend = ax.get_legend()
    if legend:
        legend.remove()

def split_box(ax, data, x, y, hue="ScenarioNames", order=None, palette=COLORS,
              hue_order=None, ylabel="", xlabel_rotation=0,
              tick_fontsize=14, label_fontsize=17):
    ax.yaxis.grid(True, linestyle="-", alpha=0.2)
    ax.xaxis.grid(False)
    ax.set_facecolor("white")

    if hue_order is None:
        hue_order = scenario_order
    sns.boxplot(
        data=data, x=x, y=y, hue=hue, order=order, hue_order=hue_order,
        palette=palette,
        notch=True,
        flierprops=dict(marker="o", markersize=6,
                         markerfacecolor="lightgrey", markeredgecolor="black",
                         alpha=0.5),
        linewidth=1,
        legend=False,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontsize=label_fontsize)
    ax.tick_params(axis="x", labelsize=tick_fontsize)
    ax.tick_params(axis="y", labelsize=tick_fontsize)
    if xlabel_rotation:
        plt.setp(ax.get_xticklabels(), rotation=xlabel_rotation, ha="right", rotation_mode="anchor")
    ax.spines[["top", "right"]].set_visible(False)

# =============================================================================
# 3.  LAYOUT  –  3 rows × 4 cols, rows 0-1 span cols 0-2 for the wide plots
# =============================================================================
fig = plt.figure(figsize=(15*0.8, 10*0.8), layout="constrained")   # bigger canvas, constrained layout
gs = gridspec.GridSpec(3, 6, figure=fig)                   # no manual hspace/wspace here
fig.get_layout_engine().set(h_pad=0.08, w_pad=0.08, hspace=0.08, wspace=0.06)

#ax_geo_rf    = fig.add_subplot(gs[0, 0:4])
ax_v_rf    = fig.add_subplot(gs[0, 0:2])
ax_h_rf    = fig.add_subplot(gs[0, 2:4])
ax_3d_rf    = fig.add_subplot(gs[0, 4:6])
ax_fn_pn     = fig.add_subplot(gs[1, 0:6])   # was 0:7 -> only 6 cols exist

ax_pns_T     = fig.add_subplot(gs[2, 0])
ax_pns_RH    = fig.add_subplot(gs[2, 1])
ax_pns_Light = fig.add_subplot(gs[2, 2])
ax_opt_T     = fig.add_subplot(gs[2, 3])
ax_max_T     = fig.add_subplot(gs[2, 4])
ax_max_RH    = fig.add_subplot(gs[2, 5])

# =============================================================================
# 4.  PANELS DEFINED IN THE R SNIPPET
# =============================================================================

# --- p_geo_rf ----------------------------------------------------------------
# split_violin(
#     ax_geo_rf, geo_data,
#     x="Variable", y="Value",
#     order=GEO_VARS,
#     ylabel="Range filling (%)"
# )
# ax_geo_rf.set_yscale("log")
# split_box(
#     ax_geo_rf, geo_data,
#     x="Variable", y="Value",
#     order=GEO_VARS,
#     ylabel="Range filling (%)"
# )
#ax_h_rf.set_title("Geographical range filling", fontsize=20, pad=6)

split_violin_single(ax_v_rf, geo_data[geo_data["Variable"] == "Vertical\nrange filling"],
                    "Value", ylabel="Range filling (%)", xlabel="Vertical range filling", width=0.6)
split_violin_single(ax_h_rf, geo_data[geo_data["Variable"] == "Horizontal\nrange filling"],
                    "Value", xlabel="Horizontal range filling", width=0.6)
split_violin_single(ax_3d_rf, geo_data[geo_data["Variable"] == "3D\nrange filling"],
                    "Value", xlabel="3D range filling", width=0.6)

# --- p_fn_pn -----------------------------------------------------------------
split_violin(
    ax_fn_pn, fn_pn_data,
    x="Group", y="Value",
    order=GROUP_ORDER,
    ylabel="Niche filling (%)"
)
#ax_fn_pn.set_title("FN / PN per environmental axis", fontsize=20, pad=6)

# ── Uncomment and replace each block with e.g.:
#    split_violin(ax_max_T, your_data, x="Variable", y="Value", ylabel="°C")

split_violin_single(ax_max_T, distances[distances["Variable"] == "Temperature"], "DistanceToMax",
             ylabel="Distance to max.\ntemperature (°C)")
split_violin_single(ax_max_RH, distances[distances["Variable"] == "Relative humidity"], "DistanceToMax",
             ylabel="Distance to max.\nhumidity (%)")
split_violin_single(ax_pns_T, niche_filling[niche_filling["Variable"] == "Temperature PN size"], "Value",
             ylabel="Temperature\nPN size (°C)")
split_violin_single(ax_pns_RH, niche_filling[niche_filling["Variable"] == "Relative humidity PN size"],
             "Value", ylabel="Relative humidity\nPN size (%)")
split_violin_single(ax_pns_Light, niche_filling[niche_filling["Variable"] == "Light PN size"],
             "Value", ylabel="Light PN size (lux)")
split_violin_single(ax_opt_T, distances[distances["Variable"] == "Temperature"], "DistanceToOpt",
             ylabel="Distance to opt.\ntemperature (°C)")

# =============================================================================
# 6.  SHARED LEGEND
# =============================================================================

from matplotlib.patches import Patch

legend_handles = [Patch(facecolor=COLORS[s], edgecolor="black", label=s) for s in scenario_order]

fig.legend(
    handles=legend_handles,
    loc="lower center",
    ncol=2,
    fontsize=16,
    bbox_to_anchor=(0.5, 0.64),
    frameon=False,
    title=None,
)

plt.tight_layout()

# =============================================================================
# 7.  SAVE
# =============================================================================
output_fname = "niche_range_filling_v11.pdf"
output_path = Path("../../figs/a5_plots_test/cc_vs_no_cc")
fig.savefig(output_path / output_fname, bbox_inches="tight", dpi=150)
print(f"Saved → {output_path}")
plt.show()