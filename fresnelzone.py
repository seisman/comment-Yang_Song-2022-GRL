"""
Calculate the Frensel zone of PKiKP reflection at ICB.
"""
import numpy as np
from obspy.geodetics import locations2degrees
from obspy.taup import TauPyModel

model = TauPyModel(model="prem")

# event 2007
evla, evlo, evdp = -55.764, -27.016, 47.7
# station AAK
stla, stlo = 42.6375, 74.4942
# ICB depth
dICB = 5149.5

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
# reflection point and PKiKP traveltime
rlon, rlat = ar.pierce["lon"][idx], ar.pierce["lat"][idx]
PKiKP_time = ar.time


def scattering_time(evla, evlo, evdp, stla, stlo, dICB, slat, slon):
    """
    Calculate traveltime for a scattering wave at ICB.

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
    dICB : float
        ICB depth.
    slat : float
        Scattering point latitude.
    slon : float
        Scattering point longitude.

    Returns
    -------
    t : float
        Total traveltime for the scatter wave.
    """
    # from ICB to the source
    t1 = model.get_travel_times(
        source_depth_in_km=dICB,
        receiver_depth_in_km=evdp,  # depth of ICB
        distance_in_degree=locations2degrees(evla, evlo, slat, slon),
        phase_list=("p",),
    )[0].time
    # from ICB to the receiver
    t2 = model.get_travel_times(
        source_depth_in_km=dICB,
        receiver_depth_in_km=0.0,  # depth of ICB
        distance_in_degree=locations2degrees(stla, stlo, slat, slon),
        phase_list=("p",),
    )[0].time
    return t1 + t2


# Calculate the traveltimes for a grid of points around the reflection point
rlon = np.round(rlon, decimals=4)
rlat = np.round(rlat, decimals=4)
dts = []
for lon in np.arange(rlon - 8.0, rlon + 8.0, 0.5):
    for lat in np.arange(rlat - 8.0, rlat + 8.0, 0.5):
        dt = scattering_time(evla, evlo, evdp, stla, stlo, dICB, lat, lon) - PKiKP_time
        dts.append(dt)
x, y = np.meshgrid(
    np.arange(rlon - 8.0, rlon + 8.0, 0.5), np.arange(rlat - 8.0, rlat + 8.0, 0.5)
)
dts = np.array(dts).reshape(x.shape)

# Save the data to a file
np.savez("fresnelzone", x=x, y=y, dts=dts)


"""
# How to read and visualize the Fresnel zone.
import matplotlib.pyplot as plt

# Load the data from a file and plot
npzfile = np.load("fresnelzone.npz")
x, y, dts = npzfile["x"], npzfile["y"], npzfile["dts"]
plt.contour(x, y, dts, levels=[0.0, 0.5])
plt.show()
"""
