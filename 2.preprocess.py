"""
2. Preprocess data

Data processing:

1. Reading waveforms
2. Remove mean and linear trend
3. Tapering
4. Integrate from velocity to displacement
5. Apply the WWSP instrument response
6. Interpolate to a 1000 Hz sampling rate
7. Filter the data
"""
import glob
from pathlib import Path

from obspy import read

from helpers import Event

# Poles & Zeros for the WWSSN short-period sensor, i.e., wwsp in SAC.
# The poles and zeros information are from SAC source code `sac/src/icm/wwsp.c`.
paz_wwsp = {
    "poles": [
        -5.0136607 + 6.4615109j,
        -5.0136607 - 6.4615109j,
        -8.2981509 + 0.0j,
        -8.6940765 - 7.1968661j,
        -8.6940765 + 7.1968661j,
    ],
    "zeros": [0j, 0j, 0j],
    "gain": 397.54767,
    "sensitivity": 1.0,
}

# frequency band in filtering
freqmin, freqmax = 0.5, 1.0

# List of events
events = [
    Event("1999-11-21T10:39:53.080", -55.726, -26.949, 33.0, 4.8),
    Event("2007-01-21T15:03:16.370", -55.764, -27.016, 47.7, 5.6),
    Event("2006-01-23T13:52:25.290", -55.943, -26.630, 35.9, 4.9),
    Event("2006-05-13T23:53:31.780", -56.046, -27.610, 101.8, 5.0),
]

for ev in events:
    if not Path(f"waveforms/{ev.id}").exists():
        print(f"{ev.id}: directory for SAC files not found.")
        continue

    print(f"Working on event: {ev.id}")
    for sacfile in glob.glob(f"waveforms/{ev.id}/*.SAC"):
        tr = read(sacfile)[0]
        tr.detrend("demean")
        tr.detrend("linear")
        tr.taper(max_percentage=0.05, type="hann")
        tr.integrate()  # the integration is slightly different from SAC
        tr.simulate(paz_simulate=paz_wwsp)  # apply the WWSP sensor response
        tr.interpolate(sampling_rate=1000)  # interpolate to 1000 Hz
        tr.filter(
            "bandpass", freqmin=freqmin, freqmax=freqmax, corners=2, zerophase=False
        )
        tr.write(sacfile, format="SAC")
