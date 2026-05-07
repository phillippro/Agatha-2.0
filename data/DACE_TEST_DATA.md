# DACE Test Data

These files were fetched from the public DACE spectroscopy radial-velocity service with `dace-query` on 2026-05-07:

- `HD40307_DACE_HARPS03.dat`
- `GJ436_DACE_HARPS03.dat`
- `HD189567_DACE_HARPS03.dat`

The source call is reproducible with:

```sh
python3 scripts/fetch_dace_test_data.py
```

Each file is a plain-text Agatha table with columns:

```text
Time RV eRV BIS FWHM Contrast
```

`Time` is DACE `rjd`. `RV`, `eRV`, `BIS`, and `FWHM` are in m/s as returned by DACE. `RV` is median-centered per target to make the files convenient for demos and tests. `Contrast` is the DACE CCF contrast column.
