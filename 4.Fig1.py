"""
4. Plotting Figure 1.
"""

import os
import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
from obspy import Stream, read
from obspy.taup import TauPyModel

from helpers import Event, process_traces, trace_label, trace_times


def read_traces(ev_old, ev_new, freqmin, freqmax, tmark, tstart, tend):
    """
    Read in the traces for the old and new events.

    Parameters
    ----------
    ev_old : Event
        Old event.
    ev_new : Event
        New event.
    freqmin : float
        Minimum frequency for bandpass filter.
    freqmax : float
        Maximum frequency for bandpass filter.
    tmark : str
        Time marker for the phase arrival.
    tstart : float
        Start time relative to the time marker.
    tend : float
        End time relative to the time marker.

    Returns
    -------
    st : obspy.core.stream.Stream
        Stream object containing the traces.
    """
    AAKold = read(f"waveforms/{ev_old.id}/II.AAK.00.BHZ.SAC")[0]
    AAKnew = read(f"waveforms/{ev_new.id}/II.AAK.00.BHZ.SAC")[0]
    KZAnew = read(f"waveforms/{ev_new.id}/KN.KZA..BHZ.SAC")[0]
    HORSnew = read(f"waveforms/{ev_new.id}/XP.HORS.01.HHZ.SAC")[0]

    st = Stream([AAKold, AAKnew, KZAnew, HORSnew])
    process_traces(
        st, freqmin=freqmin, freqmax=freqmax, tmark=tmark, tstart=tstart, tend=tend
    )
    return st


def get_PKiKP_reflection_point(evla, evlo, evdp, stla, stlo):
    """
    Get the reflection point for PKiKP phase.

    Parameters
    ----------
    evla : float
        Event latitude.
    evlo : float
        Event longitude.
    evdp : float
        Event depth.
    stla : float
        Station latitude.
    stlo : float
        Station longitude.

    Returns
    -------
    rlat : float
        Reflection point latitude.
    rlon : float
        Reflection point longitude.
    """
    model = TauPyModel(model="iasp91")
    # PKiKP reflection point at ICB
    ar = model.get_pierce_points_geo(
        source_depth_in_km=evdp,
        source_latitude_in_deg=evla,
        source_longitude_in_deg=evlo,
        receiver_latitude_in_deg=stla,
        receiver_longitude_in_deg=stlo,
        phase_list=("PKiKP",),
    )[0]
    idx = np.argmax(ar.pierce["depth"])  # deepest point
    rlat, rlon = ar.pierce["lat"][idx], ar.pierce["lon"][idx]
    return rlat, rlon


# Setting matplotlib parameters
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.labelsize"] = 16
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["xtick.labelsize"] = 13
plt.rcParams["ytick.labelsize"] = 13

fig = plt.figure(figsize=(12, 10))
# four subplots (left for geographic maps, and right for waveforms)
ax1 = fig.add_subplot(2, 2, 1, projection=ccrs.Mercator(central_longitude=75.0))
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3, projection=ccrs.Mercator(central_longitude=-56))
ax4 = fig.add_subplot(2, 2, 4)

axs = [ax1, ax2, ax3, ax4]

# Parameters for panels A and B
ev_old = Event("1999-11-21T10:39:53.080", -55.726, -26.949, 33.0, 4.8)
ev_new = Event("2007-01-21T15:03:16.370", -55.764, -27.016, 47.7, 5.6)
# None means no filtering because the data are already filtered
freqmin, freqmax = None, None
# event location for the new event
evla, evlo, evdp = ev_new.latitude, ev_new.longitude, ev_new.depth
colors = ["blue", "cyan", "red"]  # colors for AAK, KZA, and HORS
linestyles = ["--", "-", "-."]  # line styles AAK, KZA, and HORS


######################################################################################
# Panel A: Map of stations
######################################################################################
ax = ax1
[AAKold, AAKnew, KZAnew, HORSnew] = read_traces(
    ev_old, ev_new, freqmin=freqmin, freqmax=freqmax, tmark="t0", tstart=-3, tend=2
)
st = [AAKnew, KZAnew, HORSnew]

# add stations
ax.scatter(
    [tr.stats.sac.stlo for tr in st],
    [tr.stats.sac.stla for tr in st],
    s=80,
    marker="^",
    c=colors,
)
# add station names
for tr in st:
    ax.text(
        tr.stats.sac.stlo + 0.07,
        tr.stats.sac.stla,
        f"{tr.stats.station}\ndist={tr.stats.sac.gcarc:.3f}°\naz={tr.stats.sac.az:.2f}°",
        fontsize="medium",
        ha="left",
        va="bottom" if tr.stats.station == "HORS" else "center",
    )
# add raypaths
for tr in st:
    ax.plot(
        [evlo, tr.stats.sac.stlo], [evla, tr.stats.sac.stla], color="gray", zorder=-1
    )
# frame settings
ax.set_xlim(73.5, 77.0)
ax.set_ylim(40.0, 43.5)
ax.set_xticks(np.arange(73.5, 77.1, 1.0))
ax.set_yticks(np.arange(40.0, 43.6, 0.5))
ax.xaxis.set_major_formatter("{x:.1f}°E")
ax.yaxis.set_major_formatter("{x:.1f}°N")
ax.xaxis.set_minor_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.25))

# Top inset map: Events and stations on global maps.
# Set a finer threshold to have better great circle path.
plateCr = ccrs.PlateCarree()
plateCr.threshold /= 50.0

axin = ax.inset_axes([0.655, 0.64, 0.35, 0.35], projection=plateCr)
axin.coastlines()
axin.set_extent([-50, 90, -60, 50])
for tr in [AAKnew, KZAnew, HORSnew]:  # raypaths
    axin.plot(
        [evlo, tr.stats.sac.stlo],
        [evla, tr.stats.sac.stla],
        color="gray",
        linewidth=0.5,
        transform=ccrs.Geodetic(),
        zorder=1,
    )
axin.scatter(evlo, evla, c="red", marker="*", s=80, transform=ccrs.Geodetic())  # event
axin.text(
    evlo, evla + 5, "SSI", fontweight="bold", fontsize=8, ha="center", va="bottom"
)
axin.scatter(  # stations
    tr.stats.sac.stlo,
    tr.stats.sac.stla,
    c="black",
    marker="^",
    transform=ccrs.Geodetic(),
)
rlat, rlon = get_PKiKP_reflection_point(
    evla, evlo, evdp, tr.stats.sac.stla, tr.stats.sac.stlo
)
axin.scatter(rlon, rlat, c="red", marker="o", s=20)  # PKiKP reflection points

# Bottom inset map for PKiKP Fresnel zone
axin = ax.inset_axes([0.632, 0.05, 0.35, 0.35], projection=plateCr)
# add the PKiKP reflection point
for tr, color in zip([AAKnew, KZAnew, HORSnew], colors):
    rlat, rlon = get_PKiKP_reflection_point(
        evla, evlo, evdp, tr.stats.sac.stla, tr.stats.sac.stlo
    )
    axin.scatter(rlon, rlat, marker="o", s=6, color=color)

# add the Fresnel zone ellipse
if not os.path.exists("fresnelzone.npz"):
    sys.exit("Run fresnelzone.py first to generate fresnelzone.npz")
npzfile = np.load("fresnelzone.npz")
x, y, dts = npzfile["x"], npzfile["y"], npzfile["dts"]
# draw the dt=0.5 s line as the Fresnel zone
axin.contourf(x, y, dts, levels=[0, 0.5], colors=["blue"], alpha=0.2, zorder=-5)

# Frame settings
axin.set_xlim(23, 43)
axin.set_ylim(-20, 0)
axin.set_xticks(np.arange(23, 43, 5.0))
axin.set_yticks(np.arange(-20, 0 + 0.1, 4.0))
axin.tick_params(axis="x", labelsize=8)
axin.tick_params(axis="y", labelsize=8)
axin.xaxis.set_major_formatter("{x:.1f}°")
axin.yaxis.set_major_formatter("{x:.1f}°")
axin.set_title("PKiKP Fresnel zone at ICB", fontsize=7, fontweight="bold")


######################################################################################
# Panel B: Comparisons of waveforms
######################################################################################
ax = ax2
[AAKold, AAKnew, KZAnew, HORSnew] = read_traces(
    ev_old, ev_new, freqmin=freqmin, freqmax=freqmax, tmark="t0", tstart=-3.0, tend=2.0
)
# comparison between new HORS and new KZA
yshift = 0.0
tr = HORSnew
(line3,) = ax.plot(
    trace_times(tr),
    tr.data - yshift,
    label=trace_label(tr),
    color=colors[2],
    linestyle=linestyles[2],
)
tr = KZAnew
(line2,) = ax.plot(
    trace_times(tr),
    tr.data - yshift,
    label=trace_label(tr),
    color=colors[1],
    linestyle=linestyles[1],
)
ax.text(-1.5, 0.1 - yshift, "DF", fontweight="bold", fontsize=16)
ax.text(0.5, 0.45 - yshift, "CD", fontweight="bold", fontsize=16)
ax.add_artist(ax.legend(handles=[line2, line3], loc=(0.05, 0.825)))

# comparison between old AAK and new HORS
yshift = 1.8
tr = AAKold
(line1,) = ax.plot(
    trace_times(tr),
    tr.data - yshift,
    label=trace_label(tr),
    color=colors[0],
    linestyle=linestyles[0],
)
tr = KZAnew
(line2,) = ax.plot(
    trace_times(tr),
    tr.data - yshift,
    label=trace_label(tr),
    color=colors[1],
    linestyle=linestyles[1],
)
ax.text(-1.5, 0.1 - yshift, "DF", fontweight="bold", fontsize=16)
ax.text(0.5, 0.45 - yshift, "CD", fontweight="bold", fontsize=16)
ax.add_artist(ax.legend(handles=[line1, line2], loc=(0.05, 0.44)))

# comparison between old AAK and new HORS
yshift = 3.0
tr = AAKold
(line1,) = ax.plot(
    trace_times(tr),
    tr.data - yshift,
    label=trace_label(tr),
    color=colors[0],
    linestyle=linestyles[0],
)
tr = HORSnew
(line3,) = ax.plot(
    trace_times(tr),
    tr.data - yshift,
    label=trace_label(tr),
    color=colors[2],
    linestyle=linestyles[2],
)
ax.text(-1.5, 0.1 - yshift, "DF", fontweight="bold", fontsize=16)
ax.text(0.5, 0.45 - yshift, "CD", fontweight="bold", fontsize=16)
ax.add_artist(ax.legend(handles=[line1, line3], loc=(0.05, 0.18)))

# frame settings
ax.set_yticks([])
ax.set_xlabel("Time (s)")
ax.set_xlim(-3, 2.0)
ax.set_ylim(top=0.9)


# Parameters for Panel C and D
events = [
    Event("2006-01-23T13:52:25.290", -55.943, -26.630, 35.9, 4.9),
    Event("2006-05-13T23:53:31.780", -56.046, -27.610, 101.8, 5.0),
    Event("2007-01-21T15:03:16.370", -55.764, -27.016, 47.7, 5.6),
]
evlos = [ev.longitude for ev in events]
evlas = [ev.latitude for ev in events]
evnames = ["20060123", "20060513", "Doublet 99-07"]
colors = ["black", "black", "red"]

# Panel C: Map of the events
ax = ax3
ax.scatter(evlos, evlas, marker="*", c=colors, s=80)
for i in range(len(evlos)):
    ax.text(
        evlos[i] + 0.04, evlas[i], evnames[i], color=colors[i], va="top", fontsize=12
    )

ax.set_xlim(-28, -26)
ax.set_ylim(-57, -55)
ax.set_xticks(np.arange(-28, -26 + 0.1, 0.5))
ax.set_yticks(np.arange(-57, -55 + 0.1, 0.5))
ax.xaxis.set_major_formatter("{x:.1f}°")
ax.yaxis.set_major_formatter("{x:.1f}°")

# Panel D: Waveform comparisons
ax = ax4
yshift = 0.0
for ev in events:
    evdir = f"waveforms/{ev.id}"
    sta1, sta2 = "KZA", "HORS"
    tr1 = read(f"{evdir}/*.{sta1}.*.SAC")[0]
    tr2 = read(f"{evdir}/*.{sta2}.*.SAC")[0]
    process_traces(
        [tr1, tr2], freqmin=None, freqmax=None, tmark="t0", tstart=-3, tend=2.5
    )
    (line1,) = ax.plot(
        trace_times(tr1), tr1.data - yshift, color="cyan", lw=1.2, label=sta1
    )
    (line2,) = ax.plot(
        trace_times(tr2), tr2.data - yshift, color="red", ls="-.", lw=1.2, label=sta2
    )
    yshift += 1.4
ax.add_artist(ax.legend(handles=[line1, line2], loc="upper right"))
ax.text(
    -2.4,
    0.1,
    "2006/01/23 13:52:25\nLat: -55.943°\nLon: -26.63°\nDepth: 35.9 km",
    ha="left",
    fontsize=12,
)
ax.text(
    -2.4,
    -1.2,
    "2006/05/13 23:53:32\nLat: -56.046°\nLon: -27.61°\nDepth: 101.8 km",
    ha="left",
    fontsize=12,
)
# frame settings
ax.set_xlabel("Time (s)")
ax.set_ylim([-1.9, 0.8])
ax.set_xlim(-2.5, 2.5)
ax.set_yticks([])

for ax, label in zip(axs, "abcd"):
    ax.text(
        0.005,
        0.995,
        f"({label})",
        transform=ax.transAxes,
        fontsize="xx-large",
        fontweight="bold",
        verticalalignment="top",
        fontfamily="serif",
    )

fig.tight_layout(w_pad=-5.0, h_pad=4.0)
fig.savefig("Fig1.pdf")
