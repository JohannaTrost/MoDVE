import xarray as xr
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("MacOSX")
from pathlib import Path

# ------ 1. Load and process CMIP6 data ------

out_dir = Path("../../data/mc_input/climate/cmip6_ceda")

daily_clim = None

for year in range(1981, 2101):

    print(f"Processing year: {year}")

    experim = "ssp245" if year >= 2015 else "historical"

    # -- tas
    ds_tas = xr.open_dataset(out_dir / f"baf_tas_day_{experim}_{year}.nc")

    # average over lat & lon (space), avg/sd models
    tas_space_mean = ds_tas['tas'].mean(dim=["lat", "lon"])
    tas_all_mean = tas_space_mean.mean(dim="model")
    tas_model_sd = tas_space_mean.std(dim="model")

    # -- hurs
    ds_hurs = xr.open_dataset(out_dir / f"baf_hurs_day_{experim}_{year}.nc")

    # average over lat & lon (space), avg/sd models
    hurs_space_mean = ds_hurs['hurs'].mean(dim=["lat", "lon"])
    hurs_all_mean = hurs_space_mean.mean(dim="model")
    hurs_model_sd = hurs_space_mean.std(dim="model")

    # Step 3: put both into a dataframe
    df = pd.DataFrame({
        "time": ds_tas["time"].values,
        "temp_mean": tas_all_mean.values,
        "temp_model_sd": tas_model_sd.values,
        "relhum_mean": hurs_all_mean.values,
        "relhum_model_sd": hurs_model_sd.values
    }).set_index("time")

    daily_clim = pd.concat([daily_clim, df]) if daily_clim is not None else df

# --- Vizulaize

# Set the style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

# Create figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), sharex=True)

# Subplot 1: Temperature
ax1.plot(daily_clim.index, daily_clim['temp_mean'],
         color='red', alpha=0.8, linewidth=1, label='Mean Temperature')

# Add fill_between for temperature standard deviation
ax1.fill_between(daily_clim.index,
                 daily_clim['temp_mean'] - daily_clim['temp_model_sd'],
                 daily_clim['temp_mean'] + daily_clim['temp_model_sd'],
                 alpha=0.3, color='red', label='±1 SD')

ax1.set_ylabel('Temperature (°C)', fontsize=12)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Subplot 2: Relative Humidity
ax2.plot(daily_clim.index, daily_clim['relhum_mean'],
         color='blue', alpha=0.8, linewidth=1, label='Mean Relative Humidity')

# Add fill_between for relative humidity standard deviation
ax2.fill_between(daily_clim.index,
                 daily_clim['relhum_mean'] - daily_clim['relhum_model_sd'],
                 daily_clim['relhum_mean'] + daily_clim['relhum_model_sd'],
                 alpha=0.3, color='blue', label='±1 SD')

ax2.set_ylabel('Relative Humidity (%)', fontsize=12)
ax2.set_xlabel('Year', fontsize=12)
ax2.legend()
ax2.grid(True, alpha=0.3)

# Adjust layout and display
plt.tight_layout()
plt.show()

plt.savefig('../../figs/mc_input/relhum_temp_daily_ts_cmip6.png', dpi=700, bbox_inches='tight')

# -- Annual aggregated data viz

# Aggregate data to annual means
annual_clim = daily_clim.groupby(daily_clim.index.year).agg({
    'temp_mean': 'mean',
    'temp_model_sd': 'mean',
    'relhum_mean': 'mean',
    'relhum_model_sd': 'mean'
})

# Create figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), sharex=True)

# Subplot 1: Temperature
ax1.plot(annual_clim.index, annual_clim['temp_mean'],
         color='red', alpha=0.8, linewidth=2, label='Annual Mean Temperature', marker='o', markersize=3)

# Add fill_between for temperature standard deviation
ax1.fill_between(annual_clim.index,
                 annual_clim['temp_mean'] - annual_clim['temp_model_sd'],
                 annual_clim['temp_mean'] + annual_clim['temp_model_sd'],
                 alpha=0.3, color='red', label='±1 SD')

ax1.set_ylabel('Temperature (°C)', fontsize=12)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Subplot 2: Relative Humidity
ax2.plot(annual_clim.index, annual_clim['relhum_mean'],
         color='blue', alpha=0.8, linewidth=2, label='Annual Mean Relative Humidity', marker='o', markersize=3)

# Add fill_between for relative humidity standard deviation
ax2.fill_between(annual_clim.index,
                 annual_clim['relhum_mean'] - annual_clim['relhum_model_sd'],
                 annual_clim['relhum_mean'] + annual_clim['relhum_model_sd'],
                 alpha=0.3, color='blue', label='±1 SD')

ax2.set_ylabel('Relative Humidity (%)', fontsize=12)
ax2.set_xlabel('Year', fontsize=12)
ax2.legend()
ax2.grid(True, alpha=0.3)

# Adjust layout and display
plt.tight_layout()
plt.show()

plt.savefig('../../figs/mc_input/relhum_temp_annual_ts_cmip6.png', dpi=700, bbox_inches='tight')

# ----- Year selection
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from datetime import datetime

# ----------------------------
# 2. Functions for VPD + yearly features
# ----------------------------
def calc_vpd(temp_c, rh):
    # Tetens formula for saturation vapor pressure (kPa)
    es = 0.6108 * np.exp((17.27*temp_c)/(temp_c+237.3))
    ea = es * (rh/100.0)
    return es - ea

def yearly_features(df):
    df = df.copy()
    df["year"] = df["time"].dt.year
    df["month"] = df["time"].dt.month
    df["vpd"] = calc_vpd(df["temp_mean"], df["relhum_mean"])

    feats = df.groupby(["year"]).apply(lambda g: pd.Series({
        "tmean": g["temp_mean"].mean(),
        "t95": np.percentile(g["temp_mean"], 95),
        "t05": np.percentile(g["temp_mean"], 5),
        "rhmean": g["relhum_mean"].mean(),
        "rh05": np.percentile(g["relhum_mean"], 5),
        "rh95": np.percentile(g["relhum_mean"], 95),
        "vpdmean": g["vpd"].mean(),
        "vpd95": np.percentile(g["vpd"], 95),
        "drydays": (g["relhum_mean"] < 60).sum()
    })).reset_index()

    # ensemble mean per year
    feats_ens = feats.groupby("year").mean().reset_index()
    return feats_ens

features_df = yearly_features(daily_clim.reset_index())

# get extreme years
print(features_df[features_df["drydays"] == features_df["drydays"].max()])
print(features_df[features_df["t95"] == features_df["t95"].max()])
print(features_df[features_df["vpd95"] == features_df["vpd95"].max()])

# ----------------------------
# 3. Constrained medoids selection
# ----------------------------
target = 8  # number of representative years you want
hist_end = 2024  # maximum year with ERA5 inputs available

# scale features
X = features_df.drop(columns=["year"]).values
scaler = StandardScaler()
Xz = scaler.fit_transform(X)

# cluster all years (past + future if present)
km = KMeans(n_clusters=target, n_init=50, random_state=0).fit(Xz)
centers = km.cluster_centers_
cluster_labels = km.labels_

# match cluster centers to nearest HISTORICAL year
hist_mask = features_df["year"] <= hist_end
Xz_hist = Xz[hist_mask.values]
years_hist = features_df.loc[hist_mask, "year"].values
years_all = features_df["year"].values

cluster_labels_hist = cluster_labels[hist_mask.values]
cl_with_hist = np.unique(cluster_labels_hist)
print("Clusters with historical years in them:", cl_with_hist)

selected_years = []
center_years = []
for i, c in enumerate(centers):
    print(i)
    dists_all = np.linalg.norm(Xz[cluster_labels == i] - c, axis=1)
    best_all = years_all[cluster_labels == i][np.argmin(dists_all)]
    center_years.append(int(best_all))

    if i in cl_with_hist:
        dists = np.linalg.norm(Xz_hist[cluster_labels_hist == i] - c, axis=1)
        best = years_hist[cluster_labels_hist == i][np.argmin(dists)]
        selected_years.append(int(best))

        print(f"Cluster {i}: Selected historical year {best} (center year {best_all})")
    else:
        selected_years.append(int(best_all))
        print(f"Cluster {i}: No historical year, center year {best_all}")

print("Selected historical years for microclimc:", selected_years)

print("Center years for microclimc:", center_years)

# ----------------------------
# 4. Coverage diagnostic: distance of each year to nearest selected hist year
# ----------------------------
sel_mask = features_df["year"].isin(selected_years)
Xz_sel = Xz[sel_mask.values]
years_all = features_df["year"].values

cover_dist = []
for i, y in enumerate(years_all):
    d = np.min(np.linalg.norm(Xz_sel - Xz[i], axis=1))
    cover_dist.append(d)

coverage_df = pd.DataFrame({"year": years_all, "dist": cover_dist})
print("\nWorst-covered years:")
print(coverage_df.sort_values("dist", ascending=False).head(10))

# ----------------------------
# 5. Visualize clusters and selected years
# ----------------------------


import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# X = yearly feature matrix (n_years x n_features)


# Reduce to 2D for plotting
pca = PCA(n_components=2)
X_2d = pca.fit_transform(X)

print("Explained variance ratio by PCA components:", pca.explained_variance_ratio_)
print("Total explained variance:", round(sum(pca.explained_variance_ratio_), 3))

plt.figure(figsize=(10, 7))
scatter = plt.scatter(X_2d[:, 0], X_2d[:, 1], c=cluster_labels, cmap='tab10', s=25)

# Annotate points with year
for i, year in enumerate(years_all):
    plt.annotate(year, (X_2d[i, 0], X_2d[i, 1]), fontsize=8,
                 color='gray', alpha=0.5)

for c in range(len(centers)):
    center_year = center_years[c]
    sel_year = selected_years[c]

    i_sel = np.where(years_all == sel_year)[0][0]
    i_center = np.where(years_all == center_year)[0][0]

    plt.annotate(center_year, (X_2d[i_center, 0], X_2d[i_center, 1]), fontsize=8, color="grey",
                 alpha=0.5, weight='bold')
    plt.annotate(sel_year, (X_2d[i_sel, 0], X_2d[i_sel, 1]), fontsize=9, color='black', weight='bold')


plt.xlabel('PCA 1')
plt.ylabel('PCA 2')
plt.title('Yearly climate clusters (PCA 2D projection)')
plt.colorbar(scatter, label='Cluster')
plt.tight_layout()
plt.show()

plt.savefig('../../figs/mc_input/year_selection_kmeans.png', dpi=700, bbox_inches='tight')

# Print features of selected years
print(features_df[features_df["year"].isin(selected_years)])

# 2006, 2008, 2017, 2023
# 2092, 2078, 2042, 2089
# -> era5 here: 2024, 2022, 2021, 2020