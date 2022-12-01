"""
3. Picking phases using SAC.

Mark a time to align phases. The time marker should be save to SAC's "T0" header.
"""
import subprocess as sp

from helpers import Event

# List of events
events = [
    Event("1999-11-21T10:39:53.080", -55.726, -26.949, 33.0, 4.8),
    Event("2007-01-21T15:03:16.370", -55.764, -27.016, 47.7, 5.6),
    Event("2006-01-23T13:52:25.290", -55.943, -26.630, 35.9, 4.9),
    Event("2006-05-13T23:53:31.780", -56.046, -27.610, 101.8, 5.0),
]

for ev in events:
    print(f"Working on event: {ev.id}")
    s = "wild echo off \n"
    s += f"r waveforms/{ev.id}/*.SAC \n"
    s += "ppk \n"
    s += "wh \n"
    s += "q \n"
    sp.Popen(["sac"], stdin=sp.PIPE).communicate(s.encode())
