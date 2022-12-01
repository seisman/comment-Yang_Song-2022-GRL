"""
Helper functions.
"""
from pathlib import Path
from typing import NamedTuple

from obspy import UTCDateTime
from obspy.io.sac import SACTrace
from obspy.io.sac.util import get_sac_reftime


class Event:
    """
    Class for event information.
    """

    def __init__(self, origin, latitude, longitude, depth, magnitude):
        self.origin = UTCDateTime(origin)
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth
        self.magnitude = magnitude
        self.id = self.origin.strftime("%Y%m%d%H%M%S")


class Doublet(NamedTuple):
    """
    Doublet (old and new event) and the corresponding filter bands.
    """

    ev_old: Event
    ev_new: Event
    freqmin: float
    freqmax: float


def save_stream_as_sac(st, ev, inv):
    """
    Save a stream object as SAC files in the "waveforms" directory.

    Parameters
    ----------
    st : obspy.core.stream.Stream
        Stream object.
    ev : Event
        Event object.
    inv : obspy.core.inventory.inventory.Inventory
        Inventory object.
    """
    st.merge(fill_value="interpolate")
    Path(f"waveforms/{ev.id}").mkdir(parents=True, exist_ok=True)
    for tr in st:
        sac = SACTrace.from_obspy_trace(tr)

        # set event information
        sac.evla = ev.latitude
        sac.evlo = ev.longitude
        sac.evdp = ev.depth
        sac.mag = ev.magnitude
        sac.reftime = ev.origin
        sac.o = 0.0
        sac.iztype = "io"

        # set station information
        coord = inv.get_coordinates(tr.id, datetime=ev.origin)
        sac.stla = coord["latitude"]
        sac.stlo = coord["longitude"]
        sac.stel = coord["elevation"]
        sac.stdp = coord["local_depth"]

        # set SAC header
        sac.lcalda = True  # calculate distance, azimuth and back-azimuth in saving
        sac.write(f"waveforms/{ev.id}/{tr.id}.SAC")


def process_traces(st, freqmin, freqmax, tmark, tstart, tend):
    """
    Filter, trim, and normalize traces.

    Parameters
    ----------
    st : obspy.core.stream.Stream
        Stream object that should be processed.
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
    """
    # band pass filtering
    if freqmin is not None and freqmax is not None:
        st.filter(
            "bandpass", freqmin=freqmin, freqmax=freqmax, corners=2, zerophase=False
        )
    for tr in st:
        # save necessary information in the stats dictionary
        tr.stats.distance = tr.stats.sac.dist
        tr.stats.tmark = get_sac_reftime(tr.stats.sac) + tr.stats.sac[tmark]
        # trim the data to the desired time window and normalize
        tr.trim(tr.stats.tmark + tstart, tr.stats.tmark + tend)
        # normalize the data by the peak-to-peak amplitude
        tr.normalize(norm=tr.data.max() - tr.data.min())


def trace_label(tr):
    """
    Return the label of a trace for plotting.

    Parameters
    ----------
    st : obspy.core.trace.Trace
        Trace object.

    Returns
    -------
    label : str
        Label of the trace, e.g., "99 AAK".
    """
    return f"{str(tr.stats.starttime.year)[2:]} {tr.stats.station}"


def trace_times(tr):
    """
    Return the relative times of a trace for plotting.

    Parameters
    ----------
    st : obspy.core.trace.Trace
        Trace object.

    Returns
    -------
    times : numpy.ndarray
        Times relative to the time marker.
    """
    return tr.times(reftime=tr.stats.tmark)
