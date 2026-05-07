#!/usr/bin/env python3
"""Fetch small public DACE RV test data sets for Agatha.

The output tables use Agatha's plain-text convention:

    Time RV eRV BIS FWHM Contrast

Time is DACE RJD. RV, eRV, BIS, and FWHM are in m/s as returned by DACE.
RV is median-centered per target/instrument so the files are convenient demos.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from dace_query.spectroscopy import Spectroscopy


DATASETS = [
    ("HD40307", "HARPS03", "DRS-3.3.6-CCF", "HARPS"),
    ("GJ436", "HARPS03", "DRS-3.3.6-CCF", "HARPS"),
    ("HD189567", "HARPS03", "DRS-3.3.6-CCF", "HARPS"),
]


def clean_name(value: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in value).strip("_")


def fetch_one(target: str, instrument: str, drs: str, mode: str) -> np.ndarray:
    all_data = Spectroscopy.get_timeseries(
        target,
        limit=10000,
        sorted_by_instrument=True,
    )
    data = all_data[instrument][drs][mode]

    cols = [
        np.asarray(data["rjd"], dtype=float),
        np.asarray(data["rv"], dtype=float),
        np.asarray(data["rv_err"], dtype=float),
        np.asarray(data["ccf_bispan"], dtype=float),
        np.asarray(data["ccf_fwhm"], dtype=float),
        np.asarray(data["ccf_contrast"], dtype=float),
    ]
    table = np.column_stack(cols)
    keep = np.all(np.isfinite(table[:, :3]), axis=1)
    table = table[keep]
    table = table[np.argsort(table[:, 0])]
    table[:, 1] -= np.nanmedian(table[:, 1])
    return table


def main() -> None:
    outdir = Path("data")
    outdir.mkdir(exist_ok=True)

    for target, instrument, drs, mode in DATASETS:
        table = fetch_one(target, instrument, drs, mode)
        outfile = outdir / f"{target}_DACE_{clean_name(instrument)}.dat"
        header = "Time RV eRV BIS FWHM Contrast"
        np.savetxt(
            outfile,
            table,
            fmt=["%.8f", "%.8f", "%.8f", "%.8f", "%.8f", "%.8f"],
            header=header,
            comments="",
        )
        print(f"Wrote {outfile} ({len(table)} rows)")


if __name__ == "__main__":
    main()
