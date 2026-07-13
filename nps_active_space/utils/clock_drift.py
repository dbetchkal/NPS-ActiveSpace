import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import os
from scipy.signal import correlate, find_peaks
from scipy.ndimage import sobel
import pandas as pd
import numpy as np
from shapely.geometry import Point

from nps_active_space.utils.enums import TrackSource
from nps_active_space.utils.models import FAAReleasable, Nvspl, Tracks
from nps_active_space.utils.computation import expected_Lp, NMSIM_bbox_utm, audible_time_delay, ambience_from_nvspl, interpolate_spline
from nps_active_space.utils.helpers import load_studyarea, get_deployment

logger = logging.getLogger(__name__)


def logsum(df, axis=1):
    """Take a dataframe with columns representing SPL time series and logsum them together"""
    pressure_sum = np.sum(10 ** (df / 10), axis=axis)
    return 10 * np.log10(pressure_sum)


def _format_missing_time_ranges(missing_times: pd.DatetimeIndex, max_ranges: int = 5) -> str:
    """Summarize missing timestamps as contiguous ranges for log messages."""
    if len(missing_times) == 0:
        return "none"
    missing_times = missing_times.sort_values()
    ranges: list[tuple[pd.Timestamp, pd.Timestamp]] = []
    range_start = missing_times[0]
    range_end = missing_times[0]
    for timestamp in missing_times[1:]:
        if timestamp - range_end <= pd.Timedelta(seconds=1):
            range_end = timestamp
        else:
            ranges.append((range_start, range_end + pd.Timedelta(seconds=1)))
            range_start = timestamp
            range_end = timestamp
    ranges.append((range_start, range_end + pd.Timedelta(seconds=1)))

    parts = [f"{start} to {end}" for start, end in ranges[:max_ranges]]
    if len(ranges) > max_ranges:
        parts.append(f"... and {len(ranges) - max_ranges} more range(s)")
    return "; ".join(parts)


def _find_acoustic_peaks(sobel_signal: pd.Series, min_prominence_db: float, min_isolation_sec: int) -> pd.DataFrame:
    """Find prominent, well-isolated peaks in a Sobel-filtered NVSPL signal.

    Used by the "peak_match" clock drift method as an alternative to whole-day
    cross-correlation: rather than forcing a match every day (which can lock onto
    wind/ambient noise on quiet days), this only surfaces acoustic events that
    are unambiguously louder than their surroundings.

    Returns
    -------
    pd.DataFrame
        Columns "time" and "prominence", sorted by descending prominence so the
        most credible candidate event is first. Empty if no peak clears the
        prominence/isolation thresholds.
    """
    peak_positions, properties = find_peaks(
        sobel_signal.to_numpy(),
        prominence=min_prominence_db,
        distance=max(min_isolation_sec, 1),
    )
    if len(peak_positions) == 0:
        return pd.DataFrame(columns=["time", "prominence"])
    peaks = pd.DataFrame({
        "time": sobel_signal.index[peak_positions],
        "prominence": properties["prominences"],
    })
    return peaks.sort_values("prominence", ascending=False).reset_index(drop=True)


def _predicted_flight_peaks(pts_section: gpd.GeoDataFrame) -> pd.DataFrame:
    """For each flight, find its predicted closest-approach time (peak predicted Lp).

    Returns
    -------
    pd.DataFrame
        Columns "track_id", "time", "Lp_est" -- one row per flight in `pts_section`.
    """
    rows = []
    for track_id, group in pts_section.groupby("track_id"):
        peak_idx = group["Lp_est"].idxmax()
        rows.append({
            "track_id": track_id,
            "time": group.loc[peak_idx, "time_audible"],
            "Lp_est": group.loc[peak_idx, "Lp_est"],
        })
    return pd.DataFrame(rows, columns=["track_id", "time", "Lp_est"])


def match_single_event(
    sobel_signal: pd.Series,
    pts_section: gpd.GeoDataFrame,
    max_clock_drift: pd.Timedelta,
    min_prominence_db: float = 6.0,
    min_isolation_sec: int = 30,
) -> tuple[float | None, dict | None]:
    """Estimate clock drift by matching a single confident acoustic peak to the
    nearest predicted flight arrival, instead of cross-correlating an entire day.

    Whole-day correlation is forced to pick a "best" lag even on quiet or windy
    days with no real aircraft signal, which is a likely source of daily drift
    estimates that implausibly flip sign. This alternative only returns a drift
    estimate when there is one clear, isolated acoustic peak with an unambiguous
    predicted flight nearby; otherwise it returns None, the same way boundary-limited
    correlation matches are rejected today.

    Parameters
    ----------
    sobel_signal: pd.Series
        Sobel-filtered, transportation-band NVSPL signal indexed by time (as computed
        in `ClockDriftFixer.get_clock_drift`).
    pts_section: gpd.GeoDataFrame
        Track points for the period, with "track_id", "time_audible", and "Lp_est" columns.
    max_clock_drift: pd.Timedelta
        Only consider predicted flight peaks within this distance of the acoustic peak.
    min_prominence_db: float, default 6.0
        Minimum peak prominence (dB) for an acoustic peak to be considered credible.
    min_isolation_sec: int, default 30
        Minimum separation (seconds) enforced between acoustic peaks.

    Returns
    -------
    drift_sec: float or None
        Seconds to add to the track time to align with the acoustic time, or None
        if no confident single-event match was found.
        match_info: dict or None
        Details about the match (acoustic_time, predicted_time, track_id, prominence,
        confidence), useful for QC plotting and logging. None if no match was found.
        ``confidence`` is the acoustic peak prominence in dB.
    """
    acoustic_peaks = _find_acoustic_peaks(sobel_signal, min_prominence_db, min_isolation_sec)
    if acoustic_peaks.empty:
        return None, None

    flight_peaks = _predicted_flight_peaks(pts_section)
    if flight_peaks.empty:
        return None, None

    for _, acoustic_peak in acoustic_peaks.iterrows():
        window_start = acoustic_peak["time"] - max_clock_drift
        window_end = acoustic_peak["time"] + max_clock_drift
        candidates = flight_peaks[
            (flight_peaks["time"] >= window_start) & (flight_peaks["time"] <= window_end)
        ]
        if candidates.empty:
            continue
        # the loudest predicted flight in the window is our best guess for
        # which flight produced this acoustic peak
        best_candidate = candidates.loc[candidates["Lp_est"].idxmax()]
        drift_sec = (acoustic_peak["time"] - best_candidate["time"]).total_seconds()
        prominence = float(acoustic_peak["prominence"])
        match_info = {
            "acoustic_time": acoustic_peak["time"],
            "predicted_time": best_candidate["time"],
            "track_id": best_candidate["track_id"],
            "prominence": prominence,
            "confidence": prominence,
        }
        return drift_sec, match_info

    return None, None


class ClockDriftFixer():
    """A tool to fix clock drift between the geographic and acoustic records by using cross-correlation.
    
    Usage:
    (1) Load NVSPL and track data
    (2) Instantiate this class, e.g. `fixer = ClockDriftFixer(...)`
    (3) Call fixer.drift_time_series(...)
    (4) Examine plots to determine which clock drift estimations are credible.
    (5) Call fixer.fit_drift_lines(...) to use the credible points to fit linear models of clock drift.
        The clock drift predictions based on the models get saved to a clock drift file.
    (6) The ground truthing module will automatically find the clock drift file and apply the corrections.
    """
    
    def __init__(self, project_dir: str, unit: str, site: str, year: str,
                 pts: Tracks, nvspl: Nvspl, database_type: TrackSource,
                 faa_path: str, faa_corrections_path: str, plot_dir: str = None):
        """Constructor for the ClockDriftFixer class. Loads data and computes predicted audibility of causal data.

        Parameters
        ----------
        project_dir: str
            Path to directory containing site directories (e.g. containing DENATRLA/)
        unit: str
            NPS unit code, e.g. "DENA"
        site: str
            Deployment site code, e.g. "TRLA"
        year: str
            Deployment year, e.g. "2025"
        pts: Tracks
            A GeoDataFrame
        pts: Tracks
            A Tracks object containing causal track points and times, with a row for each point.
            The times in this GeoDataFrame have not been corrected for clock drift.
        nvspl: Nvspl
            An NVSPL object containing the acoustic record for the period of interest. Should overlap
            with pts as much as possible. Any partial days will be excluded from analysis.
        database_type: TrackSource
            Causal track feed (GPS, ADSB, or AIS). Used only for the default clock drift filename.
        faa_path: str
            Path to the FAA releasable aircraft database (MASTER.txt).
        faa_corrections_path: str
            Path to FAA aircraft type corrections JSON.
        plot_dir: str, Optional
            If provided, intermediate plots will be saved to this directory. This is very helpful for
            doing quality control, and is highly recommended.
        """
        self.deployment = unit + site + year
        self.default_clock_drift_file = os.path.join(
            project_dir, unit+site, f"{unit}{site}{year}_clock_drift_{database_type}.csv")
        self.plot_dir = plot_dir
        if plot_dir is not None:
            os.makedirs(plot_dir, exist_ok=True)

        logger.info("Loading site data")
        studyarea = load_studyarea(project_dir, unit, site, year)
        utm_zone = NMSIM_bbox_utm(studyarea)
        mic = get_deployment(project_dir, unit, site, year).to_crs(utm_zone)
        mic_pt = Point(mic.x, mic.y, mic.z)
        self.nvspl = nvspl

        logger.info("Filtering for only fixed-wing and helicopters")
        ids = pts["track_id"].str.split("_").str[0]  # ids will be either N-numbers or ICAO addresses
        uses_n_numbers = ids[0].startswith("N")
        if uses_n_numbers:
            ids = ids.str[1:]  # the track id lists N-numbers like: N12345, but FAA uses just 12345
        faa = FAAReleasable(faa_path,
                            faa_corrections_path,
                            n_numbers=ids.unique() if uses_n_numbers else [],
                            icao_addresses=ids.unique() if not uses_n_numbers else [],
                            warnings=False).data
        if faa.empty:
            raise Exception("No matching FAA data found")
        faa_filtered = faa[faa["TYPE AIRCRAFT"].isin(["Fixed-wing", "Helicopter"])]
        id_col = "N-NUMBER" if uses_n_numbers else "MODE S CODE HEX"
        ids_to_keep = faa_filtered[id_col].unique()
        pts = pts[ids.isin(ids_to_keep)]
        
        # convert to 3D points in the correct CRS
        pts = pts.to_crs(utm_zone)
        pts.geometry = gpd.points_from_xy(
            pts.geometry.x, pts.geometry.y, pts.z)
        
        # Interpolate points. This is essential for GPS, and seems to help smooth stuff out for ADSB too
        logger.info("Interpolating points")
        interp_list = []
        for track_id, group in pts.groupby("track_id"):
            if len(group) < 3:  # need at least 3 points to interpolate
                continue
            df = interpolate_spline(group)
            df["track_id"] = track_id
            interp_list.append(df)
        pts = pd.concat(interp_list, axis=0, ignore_index=True)

        # compute expected time and magnitude of audiblity
        logger.info("Computing audibility")
        pts = expected_Lp(pts, mic_pt)
        pts = audible_time_delay(pts, "point_dt", mic_pt)
        self.pts = pts

        logger.info("Initialization done")


    def get_clock_drift(self, start_dt, end_dt, max_clock_drift=pd.Timedelta(minutes=5),
                        method: str = "correlation", min_prominence_db: float = 6.0,
                        min_isolation_sec: int = 30):
        """
        Gets clock drift during a period of time by comparing the acoustic and causal records.
        
        Parameters
        ----------
        start_dt: pd.Timestamp or datetime
            Start time of the period in which to compute clock drift.
        end_dt: pd.Timestamp or datetime
            End time of the period in which to compute clock drift.
        max_clock_drift: pd.Timedelta
            The maximum expected magnitude of clock drift. Will only check for drifts less than this
            in magnitude. This helps to avoid situations where the most powerful predicted Lp peaks happen
            to line up well two hours in the future, but we know clock drift isn't two hours.
        method: {"correlation", "peak_match"}, default "correlation"
            "correlation" cross-correlates the whole day's Sobel-filtered NVSPL signal against the
            predicted SPL signal, picking the lag with the highest correlation. This can produce
            implausible day-to-day sign flips on quiet or windy days with no real aircraft signal
            to lock onto, since a "best" lag is always chosen even when it's just noise.
            "peak_match" instead finds the single most prominent, isolated acoustic peak that day
            and pairs it with the nearest predicted flight arrival (see `match_single_event`),
            returning None on days without one clear, unambiguous event rather than guessing.
        min_prominence_db: float, default 6.0
            Only used by method="peak_match". Minimum peak prominence (dB) for an acoustic peak
            to be considered credible.
        min_isolation_sec: int, default 30
            Only used by method="peak_match". Minimum separation (seconds) enforced between
            acoustic peaks.
        """
        # get the nvspl data corresponding to this time period
        full_index = pd.date_range(start_dt, end_dt, freq="s", inclusive="left")
        available_index = self.nvspl.index.intersection(full_index)
        if len(available_index) < len(full_index):
            missing_index = full_index.difference(available_index)
            logger.warning(
                "NVSPL doesn't cover full period %s to %s: missing %d of %d seconds (%s)",
                start_dt,
                end_dt,
                len(missing_index),
                len(full_index),
                _format_missing_time_ranges(missing_index),
            )
            return np.nan
        nvspl_section = self.nvspl.loc[full_index]

        # get the track points corresponding to this time period
        pts_section = self.pts[(self.pts["point_dt"] > start_dt) & (self.pts["point_dt"] < end_dt)]
        if pts_section.empty:
            logger.warning(
                "No track points for period %s to %s (exclusive end)",
                start_dt,
                end_dt,
            )
            return np.nan

        # Correlating with NVSPL broadband SPL has issues with wind noise obscuring the aircraft signal.
        # To fix this, first process the NVSPL spectrogram with a vertical sobel filter,
        # which calculates the derivative across bands. Wind gusts tend to stretch across most
        # bands with gradual changes, but tonal sounds like fixed-wing aircraft have sharply differing
        # power levels between nearby bands. We then logsum this derivative value, which emphasizes
        # aircraft peaks more than summing. We only logsum over the transportation band to reduce noise.
        # Empirically, this works very well.
        spectro = nvspl_section.loc[:,"12.5":"20000"]
        sobel_signal = np.abs(sobel(spectro, axis=1))
        sobel_signal = pd.DataFrame(sobel_signal, index=spectro.index, columns=spectro.columns)
        sobel_signal = logsum(sobel_signal.loc[:,"80":"1250"])  # transportation band

        if method == "peak_match":
            drift_sec, match_info = match_single_event(
                sobel_signal, pts_section, max_clock_drift,
                min_prominence_db=min_prominence_db, min_isolation_sec=min_isolation_sec,
            )
            if drift_sec is None:
                logger.warning(
                    "No confident single-event peak match found for period %s to %s",
                    start_dt,
                    end_dt,
                )
                return None
            logger.info(
                "Peak match for %s to %s: acoustic_time=%s, drift=%.1f s, "
                "track=%s, confidence=%.1f dB prominence",
                start_dt,
                end_dt,
                match_info["acoustic_time"],
                drift_sec,
                match_info["track_id"],
                match_info["confidence"],
            )
            if self.plot_dir is not None:
                self._plot_peak_match(start_dt, end_dt, sobel_signal, match_info, drift_sec)
            return drift_sec

        # Set up an index that will be used for the predicted signal.
        # This contains each second of the time period, minus a period at the beginning and end
        # corresponding to the max clock drift size, to allow shifting this signal over the full NVSPL day
        # while always fully overlapping the NVSPL data
        pred_index = pd.date_range(start_dt + max_clock_drift, end_dt - max_clock_drift, freq="s", inclusive="left")

        # initialize a list of predicted SPL time series that we will later logsum to get the overall predicted SPL time series
        # put ambience in it because we want to add in ambience
        ambience = logsum(ambience_from_nvspl(nvspl_section, 90), axis=0)
        Lp_series_list = [pd.Series(index=pred_index, data=ambience)]

        # add in each flight's Lp contribution to the Lp_series_list
        for _, group in pts_section.groupby("track_id"):
            # get a series of predicted Lp values at various times for this flight
            # we would like to interpolate this series to make it line up with pred_index,
            # which we need in order to cross-correlate with the NVSPL data
            orig_series = pd.Series(
                index=group["time_audible"], data=group["Lp_est"].values)

            # figure out which times we need to calculate interpolated Lp_est values
            # within the temporal extent of this flight
            sec_index = pd.date_range(group["time_audible"].min().floor(
                "s"), group["time_audible"].max().ceil("s"), freq="s")
            seconds_to_calc = pd.Series(
                index=sec_index.difference(orig_series.index), data=np.nan)
            # .difference() because we don't want to calculate interpolated Lp_est if we
            # already have Lp_est for that time

            # interpolate
            series = pd.concat((orig_series, seconds_to_calc)).sort_index()
            series = series.interpolate(method="time", limit=1, limit_direction="both")
            # just pick out the seconds that match with the pred_index
            series = series.loc[sec_index]
            # expand the index to the full pred_index
            series = series.reindex(pred_index)
            # set any na values to no sound
            series.loc[series.isna()] = -100  # -100 dB for "no sound"
            Lp_series_list.append(series)

        Lp_pred = logsum(pd.concat(Lp_series_list, axis=1, ignore_index=True))

        # correlate to determine max clock drift
        correlation = correlate(sobel_signal, Lp_pred, mode="valid")
        assert len(correlation) == 2 * max_clock_drift.total_seconds() + 1, \
            f"{len(correlation)} != {2 * max_clock_drift.total_seconds() + 1}"
        max_idx = np.argmax(correlation)
        if max_idx == 0 or max_idx == len(correlation)-1:
            logger.warning(
                "Max correlation at boundary for period %s to %s, not true maximum",
                start_dt,
                end_dt,
            )
            return None  # boundary-limited, not true peak
        drifts = np.arange(-max_clock_drift.total_seconds(), max_clock_drift.total_seconds() + 1)
        clock_drift = drifts[max_idx]

        # plot
        if self.plot_dir is not None:
            fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(18, 10),
                                gridspec_kw={'height_ratios': [1, 1, 2]})
            fig.suptitle(f"{self.deployment}: {start_dt} - {end_dt}")

            ax[0].set_title("NVSPL")
            ax[0].set_ylabel("Sobel Filter Signal")
            sobel_signal.plot(ax=ax[0], color="black", linewidth=0.5)

            ax[1].set_title("Predicted SPL")
            ax[1].set_ylabel("SPL (dB)")
            Lp_pred.plot(ax=ax[1], color="black", linewidth=0.5)
            ax[1].sharex(ax[0])

            ax[2].set_title("Correlation")
            ax[2].plot(drifts, correlation, color="black")
            ax[2].vlines(clock_drift, correlation.min(), correlation.max(), linestyle="--", color="black")
            ax[2].set_xlabel("Possible Clock Drifts (sec)")
            ax[2].set_ylabel("Correlation")

            fig.tight_layout()
            fig.savefig(os.path.join(self.plot_dir, f"{start_dt.date()}.png"))
            plt.close()
        
        return clock_drift


    def _plot_peak_match(self, start_dt, end_dt, sobel_signal: pd.Series, match_info: dict, drift_sec: float) -> None:
        """Save a QC plot for a single-event peak match (used when method="peak_match")."""
        fig, ax = plt.subplots(figsize=(18, 6))
        fig.suptitle(f"{self.deployment}: {start_dt} - {end_dt} (peak match, drift={drift_sec:.1f}s)")
        ax.set_title("NVSPL Sobel Filter Signal")
        ax.set_ylabel("Sobel Filter Signal")
        sobel_signal.plot(ax=ax, color="black", linewidth=0.5)
        ax.axvline(match_info["acoustic_time"], color="red", linestyle="--", label="Acoustic peak")
        ax.axvline(match_info["predicted_time"], color="blue", linestyle="--",
                  label=f"Predicted arrival ({match_info['track_id']})")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(self.plot_dir, f"{start_dt.date()}_peak_match.png"))
        plt.close(fig)
    

    def _init_time_series_plot(self):
        """Utility for setting up a time series plot of clock drifts"""
        plt.subplots(figsize=(14, 6))
        plt.title(f"Estimated Daily Clock Drift, {self.deployment}")
        plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=5))
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        plt.ylabel("Estimated Clock Drift (sec)")
    

    def drift_time_series(self, start_dt=None, end_dt=None, max_clock_drift=pd.Timedelta(minutes=5),
                          show_plots: bool = True, method: str = "correlation",
                          min_prominence_db: float = 6.0, min_isolation_sec: int = 30):
        """Compute clock drift each day over an extended period of time to determine patterns over time.

        See `ClockDriftFixer.get_clock_drift` for a description of `method`, `min_prominence_db`,
        and `min_isolation_sec`.
        """
        if start_dt is None:
            start_dt = self.nvspl.index.min()
        if end_dt is None:
            end_dt = self.nvspl.index.max()

        logger.info(
            "Estimating daily clock drift using method=%s from %s to %s "
            "(max_drift=%s, min_prominence_db=%.1f, min_isolation_sec=%d)",
            method,
            start_dt,
            end_dt,
            max_clock_drift,
            min_prominence_db,
            min_isolation_sec,
        )

        freq = "d"
        period_start = start_dt.ceil(freq)
        period_end = end_dt.floor(freq) + pd.Timedelta(days=1)
        if period_end <= period_start:
            period_end = period_start + pd.Timedelta(days=1)
        period_bounds = pd.date_range(period_start, period_end, freq=freq)
        drifts = []
        for i in range(len(period_bounds)-1):
            logger.debug(f"{period_bounds[i]} to {period_bounds[i+1]}")
            drifts.append(self.get_clock_drift(
                period_bounds[i], period_bounds[i+1], max_clock_drift,
                method=method, min_prominence_db=min_prominence_db, min_isolation_sec=min_isolation_sec,
            ))
        
        # use the center of the period as the time anchor for each period's drift
        self.times = period_bounds[:-1] + (period_bounds.diff()[1:] / 2)
        self.drifts = np.array(drifts).astype(float)

        # plot
        self._init_time_series_plot()
        plt.plot(self.times, self.drifts, marker="o", color="black")
        # add labels
        label = 0
        for i in range(len(self.drifts)):
            if not np.isnan(self.drifts[i]):
                plt.text(self.times[i], self.drifts[i], f"  {label}", fontsize=9)
                label += 1
        if self.plot_dir is not None:
            plt.savefig(os.path.join(self.plot_dir, "Clock Drift Time Series.png"))
        if show_plots:
            plt.show()
        else:
            plt.close()

        return self.times, self.drifts
    

    def _fit_drift_line(self, times, drifts):
        """Utility to fit a line, used by fit_drift_lines()"""
        times_ns = times.astype(np.int64)
        times_days = times_ns / (1e9 * 60 * 60 * 24)
        return np.polyfit(times_days, drifts, deg=1)
    

    def fit_drift_lines(self, indices_to_use, clock_drift_file=None,
                        maintenance_times=[], start_dt=None, end_dt=None,
                        show_plots: bool = True):
        """
        Clock drift tends to be linear (for ADSB particularly).
        This function fits lines for each time period between maintenance visits,
        and creates a clock drift file based on these lines.

        Parameters
        ----------
        indices_to_use: list of ints
            A list of credible clock drift estimation points to use for fitting the line.
            The indices correspond to those shown on the time series plot produced by self.drift_time_series()
        clock_drift_file: str, optional
            Filename for saving the clock drift file. If not provided, a default path is used.
        maintenance_times: list of pd.Timestamp or datetime, default []
            Times when maintenance was performed on the system. Linear models will be fit for the periods between
            maintenance visits, since resetting the clock during visits can cause discontinuities in the clock drift.
        start_dt: pd.Timestamp or datetime, optional
            Start time for estimating clock drift. Should precede all track data.
            Default is midnight before the earliest track point.
        end_dt: pd.Timestamp or datetime, optional
            End time for estimating clock drift. Should come after all track data.
            Default is midnight after the last track point.
        """
        assert hasattr(self, "times") and hasattr(self, "drifts")
        self._init_time_series_plot()

        # fill in default values for clock drift file name and time period
        if clock_drift_file is None:
            clock_drift_file = self.default_clock_drift_file
        if start_dt is None:
            start_dt = self.pts["point_dt"].min().floor("d")  # make sure to encompass track points
        if end_dt is None:
            end_dt = self.pts["point_dt"].max().ceil("d")

        # select points we will use to fit the lines
        valid = ~np.isnan(self.drifts)
        x = self.times[valid][indices_to_use]
        y = self.drifts[valid][indices_to_use]
        plt.scatter(self.times, self.drifts, c="lightgray")
        plt.scatter(x, y, color="black")
        # add labels
        for i, (xi, yi) in enumerate(zip(x, y)):
            plt.text(xi, yi, f"  {indices_to_use[i]}", fontsize=9)

        # fit a line for each period between maintenance times
        fits_list = []
        file_entries_list = []
        period_boundaries = [start_dt] + sorted(maintenance_times) + [end_dt]
        for i in range(len(period_boundaries)-1):
            period_start = period_boundaries[i]
            period_end = period_boundaries[i+1]

            period_mask = (x >= period_start) & (x < period_end)
            slope, intercept = self._fit_drift_line(x[period_mask], y[period_mask])
            row = {"start_dt": period_start, "end_dt": period_end, "slope": slope, "intercept": intercept}
            fits_list.append(row)

            # Get entries for the clock drift file.
            # Subtract a teensy bit from period end time to avoid getting a duplicate time
            # with the start of the next time period, while retaining discontinuity behavior.
            # Include times from all the days we had overlapping nvspl and tracks,
            # for better human readability of the file (though unnecessary, will get linearly
            # interpolated later)
            new_times = np.concat([
                [period_start],
                self.times[(self.times > period_start) & (self.times < period_end)],
                [period_end - pd.Timedelta(seconds=0.1)]  
            ]).astype('datetime64[ns]')
            x_new = new_times.astype(np.int64) / (1e9 * 60 * 60 * 24)
            y_new = slope * x_new + intercept
            file_entries_list.append(pd.DataFrame({"Time": new_times, "Seconds": y_new}))

            plt.plot(new_times, y_new, color="black")
            logger.info(f"{period_start} to {period_end}, clock drift {slope:.3f} sec per day")
        
        self.drift_fits = pd.DataFrame(fits_list)
        file_entries = pd.concat(file_entries_list, ignore_index=True)
        file_entries.to_csv(clock_drift_file, index=False)
        logger.debug(f"Clock drift file entries:\n{file_entries}")

        if self.plot_dir is not None:
            plt.savefig(os.path.join(self.plot_dir, "Fitted Lines.png"))
        if show_plots:
            plt.show()
        else:
            plt.close()

        return self.drift_fits


def correct_clock_drift(tracks: Tracks, clock_drift_file: str, inplace: bool=True):
    """Fix the point_dt field of a set of tracks by correcting for clock drift.
    
    Parameters
    ----------
    tracks: Tracks
        A GeoDataFrame containing a "point_dt" field, with a row for each point.
        The times in this GeoDataFrame have not been corrected for clock drift.
    clock_drift_file: str
        Path to the clock drift file. This is a CSV file with columns "Time" and "Seconds",
        representing the drift at various points in time. Drift is defined as what you would
        need to add to the track point time to arrive at the corresponding NVSPL time.
        This function linearly interpolates between drifts in the CSV file to find the drift
        at any moment in time. Note that the range of the CSV file times must encompass
        the range of the track point times for interpolation to work.
    inplace: bool
        Whether to modify tracks["point_dt"] in place.
    
    Returns
    -------
    tracks: Tracks
        A GeoDataFrame of track points with times shifted to account for clock drift.
    """

    # read file into a pd.Series
    drifts = pd.read_csv(clock_drift_file, index_col="Time")["Seconds"]
    drifts.index = pd.to_datetime(drifts.index)
    # make sure the clock drifts encompass the track point times, so we can interpolate
    if tracks["point_dt"].min() < drifts.index.min() or tracks["point_dt"].max() > drifts.index.max():
        raise Exception(f"Clock drift corrections must encompass the whole track period ({tracks["point_dt"].min()} - {tracks["point_dt"].max()})")
    
    # add in entries corresponding to the track point times in between existing clock drift entries 
    # then interpolate to fill them in, then extract those interpolated values to get time adjustments
    # make sure there are no duplicate times in drifts_augmented, so that when we query drifts_augmented to get the adjustments,
    # each point_dt value will only correspond to one entry in drifts_augmented
    times_to_interpolate = pd.Series(data=np.nan, index=tracks["point_dt"].unique())
    times_to_interpolate = times_to_interpolate[~times_to_interpolate.index.isin(drifts.index)]
    drifts_augmented = pd.concat([drifts, times_to_interpolate])
    drifts_augmented.sort_index(inplace=True)
    drifts_augmented.interpolate(method="time", inplace=True)
    adjustments = drifts_augmented[tracks["point_dt"]]
    adjustments = pd.to_timedelta(adjustments, unit="s")
    correct_point_dt = tracks["point_dt"] + adjustments.values

    if inplace:
        tracks["point_dt"] = correct_point_dt
        return tracks
    else:
        track_copy = tracks.copy()
        track_copy["point_dt"] = correct_point_dt
        return track_copy