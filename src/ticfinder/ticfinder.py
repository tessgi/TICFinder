"""Find TIC numbers from RA, Dec and Mag"""
from dataclasses import dataclass
from typing import List, Optional

import astropy.units as u
import numpy as np
import pandas as pd
import ticgen
from astropy.coordinates import Distance, SkyCoord
from astropy.time import Time
from astroquery.utils.tap.core import TapPlus
from tqdm import tqdm

from . import log


@dataclass
class TICFinder:
    """TICFinder class converts RA, Dec and optionally magnitude into TIC numbers.

    TICFinder queries TIC and finds closest matches. Matches are weighted by

    a) Their inverse square distance to the input RA/Dec
    b) Their relative flux compared to the input magnitude if magnitudes are supplied.

    **Note** You must input RA and Dec in J2000 epoch.

    Parameters
    ----------
    ra: List or npt.NDArray
        Array of RA values. MUST be in J2000 epoch
    dec: List or npt.NDArray
        Array of RA values. MUST be in J2000 epoch
    magnitude: Optional, List or npt.NDArray
        Array of magnitude values. Presumed to be close to TESS Magnitude and Gaia RP magnitude.
    """

    ra: List
    dec: List
    magnitude: Optional = None

    def __post_init__(self):
        self.tap = TapPlus(url="http://TAPVizieR.u-strasbg.fr/TAPVizieR/tap")

        if (self.ra is None) | (self.dec is None):
            raise ValueError("Please pass RA and Dec")
        if not hasattr(self.ra, "__iter__"):
            self.ra = [self.ra]
        if not hasattr(self.dec, "__iter__"):
            self.dec = [self.dec]

        if self.magnitude is None:
            log.warning("No magnitudes supplied")
        else:
            if not hasattr(self.magnitude, "__iter__"):
                self.magnitude = [self.magnitude]

        self.epoch_decimal = Time.now().decimalyear
        self.epoch = np.floor(self.epoch_decimal).astype(int)

    def __repr__(self):
        return f"TICFinder ({len(self.ra)} Targets)"

    @staticmethod
    def from_csv(fname):
        """Read directly from a csv file. Must include RA, Dec and optionally magnitude.

        Parameters
        ----------
        fname : str
            Filename
        """
        return TICFinder.from_pandas(pd.read_csv(fname))

    @staticmethod
    def from_pandas(df):
        """Read directly from a pandas dataframe. Must include RA, Dec and optionally magnitude.

        Parameters
        ----------
        df : pd.DataFrame
            pandas dataframe
        """

        def get_column(keys):
            hits = [key in df.columns for key in keys]
            if np.any(hits):
                return np.asarray(
                    df[np.asarray(keys)[[key in df.columns for key in keys]][0]]
                )
            else:
                return None

        return TICFinder(
            ra=get_column(["ra", "RA", "RAJ2000", "ra [ deg ]"]),
            dec=get_column(["dec", "Dec", "DEC", "DEJ2000", "dec [ deg ]"]),
            magnitude=get_column(
                ["tmag", "magnitude", "mag", "tess_mag [ mag ]", "tess_mag", "Tmag"]
            ),
        )

    def get_tics(self, non_null_motion_only=False, show_progress=True):
        """Populates the `tic` attribute. Run this to fetch the TIC numbers

        Parameters
        ----------
        non_null_motion_only: bool
            Whether to return targets that have non NULL values for pmRA, pmDE, and Plx.
            If True, will only return targets where all three are valued.
        show_progress: bool
            Whether to show a progress bar. You can silence the progress bar by setting this to False.
        """
        gdf = []
        for idx in tqdm(
            range(len(self.ra)), disable=~show_progress, desc="Cross matching targets"
        ):
            ra, dec, magnitude = (
                self.ra[idx],
                self.dec[idx],
                self.magnitude if self.magnitude is None else self.magnitude[idx],
            )
            for idx in np.arange(1, 4):
                radius = 50 * idx**2 / 3600.0
                query = get_query(
                    ra,
                    dec,
                    magnitude,
                    radius=radius,
                    non_null_motion_only=non_null_motion_only,
                    ntarg=1,
                )
                job = self.tap.launch_job(query)
                series = job.get_results().to_pandas()
                if len(series) != 0:
                    break
            if len(series) == 0:
                series = pd.Series(
                    index=[
                        "ra_input",
                        "dec_input",
                        "magnitude_input",
                        "TIC",
                        "RAJ2000",
                        "DEJ2000",
                        "pmRA",
                        "pmDE",
                        "Plx",
                        "e_Plx",
                        "RPmag",
                        "Bmag",
                        "Vmag",
                        "Hmag",
                        "Kmag",
                        "Jmag",
                        "pix_sep",
                        "weight",
                        "RAJ2022",
                        "DEJ2022",
                        "tmag",
                        "motion_from_2000_to_2022_in_pixels",
                    ],
                    dtype=np.float64,
                )
                series[["ra_input", "dec_input", "magnitude_input"]] = (
                    ra,
                    dec,
                    magnitude if magnitude is not None else np.nan,
                )
            gdf.append(series)
        gdf = pd.concat(gdf).reset_index(drop=True)
        self.gdf = self._clean_dataframe(gdf)
        self.tic = np.asarray(gdf.TIC)
        self.best_fit_ra = np.asarray(gdf["RAJ2000"])
        self.best_fit_dec = np.asarray(gdf["DEJ2000"])
        self.best_fit_magnitude = np.asarray(gdf.tmag)
        self.motion_in_pixels = np.asarray(
            gdf[f"motion_from_2000_to_{self.epoch}_in_pixels"]
        )

    def _clean_dataframe(self, df):
        df[["pmRA", "pmDE", "Plx", "e_Plx"]] = df[
            ["pmRA", "pmDE", "Plx", "e_Plx"]
        ].fillna(0)
        df[[f"RAJ{self.epoch}", f"DEJ{self.epoch}", "tmag"]] = np.nan
        df.loc[np.abs(df.Plx / df.e_Plx) < 2, "Plx"] = 1e-5
        df.loc[df.Plx == 0, "Plx"] = 1e-5

        def clean(x):
            return x if np.isfinite(x) else None

        for idx in range(len(df)):
            c = SkyCoord(
                df.loc[idx]["RAJ2000"] * u.deg,
                df.loc[idx]["DEJ2000"] * u.deg,
                frame="icrs",
                obstime="2020-01-01T00:00:00",
                pm_ra_cosdec=df.loc[idx]["pmRA"] * u.mas / u.yr,
                pm_dec=df.loc[idx]["pmDE"] * u.mas / u.yr,
                distance=Distance(
                    parallax=df.loc[idx]["Plx"] * u.mas, allow_negative=True
                ),
            ).apply_space_motion(dt=(self.epoch_decimal - 2000) * u.year)
            tmag = ticgen.Star(
                Bmag=clean(df.loc[idx, "Bmag"]),
                Vmag=clean(df.loc[idx, "Vmag"]),
                Jmag=clean(df.loc[idx, "Jmag"]),
                Hmag=clean(df.loc[idx, "Hmag"]),
                Ksmag=clean(df.loc[idx, "Kmag"]),
            ).Tmag
            df.loc[idx, [f"RAJ{self.epoch}", f"DEJ{self.epoch}", "tmag"]] = np.hstack(
                [c.ra.deg, c.dec.deg, tmag]
            )
        df[f"motion_from_2000_to_{self.epoch}_in_pixels"] = (
            np.asarray(
                np.hypot(
                    df.RAJ2000 - df[f"RAJ{self.epoch}"],
                    df.DEJ2000 - df[f"DEJ{self.epoch}"],
                )
            )
            * u.deg.to(u.arcsecond)
            / 21
        )
        return df

    def to_csv(self):
        """Converts object to a pandas dataframe for easy storing"""
        if not hasattr(self, "tic"):
            log.warning("Run `get_tics` first.")
            return
        df = pd.DataFrame(
            np.vstack(
                [self.tic, self.best_fit_ra, self.best_fit_dec, self.best_fit_magnitude]
            ).T,
            columns=["TIC", "RA", "Dec", "Tmag"],
        )
        if self.magnitude is not None:
            df["input_magnitude"] = self.magnitude
        df["pix_sep"] = self.gdf["pix_sep"]
        df[f"motion_from_2000_to_{self.epoch}_in_pixels"] = self.motion_in_pixels
        return df

    def get_target_matches(
        self, idx, radius=0.014, ntarg=None, non_null_motion_only=False
    ):
        """Returns an ordered dataframe of the matches to the input

        Parameters
        ----------
        idx: int
            The index of the target to calculate the matches for.
        radius: float
            Radius in degrees to query around
        ntarg: Optional, int
            Number of targets to return. If None, will return all matches.
        non_null_motion_only: bool
            Whether to return targets that have non NULL values for pmRA, pmDE, and Plx.
            If True, will only return targets where all three are valued.

        Returns
        -------
        df : pd.DataFrame
            pandas containing all the information on all matches around the `idx`
            index target.
        """
        ra, dec, mag = (
            self.ra[idx],
            self.dec[idx],
            self.magnitude[idx] if self.magnitude is not None else None,
        )
        log.info(f"RA:{ra}, Dec:{dec}, Mag:{mag}")
        query = get_query(
            ra,
            dec,
            mag,
            radius=radius,
            ntarg=ntarg,
            non_null_motion_only=non_null_motion_only,
        )
        df = self.tap.launch_job(query).get_results().to_pandas()
        return self._clean_dataframe(df)


def get_query(
    ra, dec, magnitude=None, ntarg=None, radius=0.014, non_null_motion_only=False
):
    """Builds a TAP query to get matches for a given RA, Dec and magnitude."""
    motion = [
        """\nAND pmRA IS NOT NULL
            AND pmDE IS NOT NULL
            AND Plx IS NOT NULL
            """
        if non_null_motion_only
        else ""
    ][0]
    ntarg_statement = f"""TOP {ntarg} """ if ntarg is not None else ""
    mag0 = f"""{magnitude} as magnitude_input,""" if magnitude is not None else ""
    mag1 = (
        f"""* 1/(2*POWER(10, ABS({magnitude} - RPMag) * -0.4))"""
        if magnitude is not None
        else ""
    )
    mag2 = (
        f"""AND RPmag >= {magnitude - 3} AND RPmag <= {magnitude + 3}"""
        if magnitude is not None
        else ""
    )
    query = f"""SELECT {ntarg_statement}
            {ra} as ra_input,
            {dec} as dec_input,{mag0} TIC,
            RAJ2000, DEJ2000, pmRA, pmDE, Plx, e_Plx, RPmag, Bmag, Vmag, Hmag, Kmag, Jmag,
            DISTANCE(
               POINT('ICRS', RAJ2000, DEJ2000),
               POINT('ICRS', {ra}, {dec})) * 3600/21 as pix_sep,
            1/POWER(DISTANCE(
               POINT('ICRS', RAJ2000, DEJ2000),
               POINT('ICRS', {ra}, {dec})) * 3600/21, 2) {mag1} as weight
            FROM "IV/38/tic"
            WHERE 1=CONTAINS(POINT('ICRS',RAJ2000,DEJ2000), CIRCLE('ICRS', {ra}, {dec}, {radius}))
            {mag2} {motion}
            ORDER BY weight DESC"""
    return query
