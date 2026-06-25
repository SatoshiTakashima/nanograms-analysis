# NanoGRAMS Test-Pulse Gain Pipeline

This directory turns the manual `testpulse.ipynb` workflow into repeatable scripts.

## 1. Fit test-pulse peaks

For the interactive human-check workflow, use the GUI:

```bash
python3 review_testpulse_gain_gui.py testpulse_gain_config.yaml
```

For each test-pulse time, the GUI shows FEC0-FEC3 together. Turn on each FEC's
`Peak` switch only when a clear test-pulse peak is visible. The single `Fit`
button fits every FEC whose `Peak` switch is on. You can edit the Gaussian
center, height, width, and fit-window half width independently for FEC0-FEC3;
input fields are editable only for FECs with `Peak` turned on. A `Peak` switch
that is on means that FEC is intended to be used in the compact gain table,
provided the fit also passes the configured quality cuts.

The fit is performed with `iminuit` using an unbinned maximum-likelihood model
over the configured fit window. The reported errors are the symmetric HESSE
errors. The plotted red curve is drawn only over that fit range. The GUI shows
both the full ADU-CMN range and peak-centered `+/- 30 ADU` zoom panels for all
four FECs. Fit overlays and fitted `Mean` labels are drawn only in the zoom
panels; the full-range panels remain plain histograms. Orange frames mark the
zoom panels whose `Peak` switch is on. `Save summary` is blocked when any FEC in
the displayed test-pulse file has `Peak` turned on but has not yet been fitted,
or when that fit failed the quality cuts. When the summary is saved, the current
GUI screen is also saved as a PNG in `products/testpulse_gain/review_plots/`.
`Fit` changes color after fitting, and `Save summary` changes color after a
successful save or blocked save attempt. `Prev` and `Next` are greyed out at the
first and last configured test-pulse files.

The GUI reads each ROOT file once per displayed test-pulse time and caches the
FEC0-FEC3 channel values, so moving within the same time and refitting selected
FECs should not reopen the ROOT file.

The GUI writes:

- `products/testpulse_gain/manual_selection.csv`: click-by-click review state.
- `products/testpulse_gain/run10_testpulse_data.csv`: accepted fit centers only.
- `products/testpulse_gain/review_plots/*_gui.png`: saved GUI screenshots.

The non-interactive batch script is still available:

```bash
python3 analyze_testpulse_gain.py testpulse_gain_config.yaml
```

Outputs:

- `products/testpulse_gain/testpulse_fits.csv`: all FEC fit diagnostics.
- `products/testpulse_gain/plots/`: per-file histograms and `time_trend.png`.
- `products/testpulse_gain/review_plots/`: one 2x2 FEC review image per test-pulse file.
- `products/testpulse_gain/manual_selection.csv`: manual OK/NG table for each time/FEC.

By default, `review.use_manual_selection` is `true`. The first run creates
`manual_selection.csv` with an empty `use_fit` column. Inspect the corresponding
`review_plots/*.png` images, fill `use_fit` with `1` for usable FECs and `0` for
rejected FECs, then rerun the same command. Only after every row has a manual
choice does the script write:

- `products/testpulse_gain/run10_testpulse_data.csv`: compact `time_id,FEC0,FEC1,FEC2,FEC3` table.

The compact table stores `NaN` for FECs rejected in `manual_selection.csv`. The
automatic fit quality is still written as `suggested_use_fit`, but the manual
`use_fit` column is the final decision when manual selection is enabled.

## 2. Interpolate gains for gamma runs

```bash
python3 interpolate_gain_for_gamma.py interpolate_gain_config.yaml
```

The output CSV contains one row per gamma run and interpolated `FEC0` to `FEC3`
peak ADU values. These values are meant to be passed to `NanoGRAMSCalibration`
as the `gain_tp_hash` for that run.

## Notes

- The default input file name is `outfile00001_000.root`. If the VATA converter
  output name changes, edit `input.file_name`.
- The default test-pulse channel is `17`, and the config records `ccal: 8`.
- A FEC is treated as alive only when the peak is strong enough, the Gaussian
  fit converges, and the fitted peak/sigma/chi-square pass the configured cuts.
