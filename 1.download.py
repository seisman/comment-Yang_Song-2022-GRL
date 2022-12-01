"""
1. Download the waveform data.

This script downloads the waveform data of three neighboring stations
(II.AAK, KN.KZA, and XP.HORS) for the 99-07 doublet and another two nearby
earthquakes from the IRIS DMC, and save them as SAC files.
"""
from pathlib import Path

from obspy.clients.fdsn import Client

from helpers import Event, save_stream_as_sac

# save data in the waveforms directory
Path("waveforms").mkdir(exist_ok=True)

# List of events to download
events = [
    Event("1999-11-21T10:39:53.080", -55.726, -26.949, 33.0, 4.8),
    Event("2007-01-21T15:03:16.370", -55.764, -27.016, 47.7, 5.6),
    Event("2006-01-23T13:52:25.290", -55.943, -26.630, 35.9, 4.9),
    Event("2006-05-13T23:53:31.780", -56.046, -27.610, 101.8, 5.0),
]

# Initialize an FDSN client
client = Client("IRIS")
for ev in events:
    print(f"Downloading event {ev.id}")
    # PKiKP and PKIKP at 130 degrees arrive at ~1150 seconds after the origin time.
    starttime, endtime = ev.origin + 1000, ev.origin + 1300
    bulk = [
        ("II", "AAK", "00", "BHZ", starttime, endtime),
        ("KN", "KZA", "", "BHZ", starttime, endtime),
        ("XP", "HORS", "01", "HHZ", starttime, endtime),
    ]
    # Request data and save to SAC files
    st = client.get_waveforms_bulk(bulk)
    inv = client.get_stations_bulk(bulk, level="channel")
    save_stream_as_sac(st, ev, inv)
