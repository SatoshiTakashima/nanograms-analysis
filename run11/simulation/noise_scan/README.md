# NanoGRAMS Simulation Noise Scan

This directory scans `database/detector_parameters.xml`:

```xml
<noise_level param0="5.0" param1="0.036" param2="0.0" />
```

For each configured `param1`, `run_noise_param_scan.py` creates an isolated run
directory, writes a modified detector-parameter XML, generates a Ruby simulation
runner, and optionally executes it.

## Dry run / prepare files

```bash
python3 run_noise_param_scan.py noise_scan_config.yaml
```

With `simulation.execute: false`, this prepares:

```text
products/noise_param1_scan/runs/param1_*/detector_parameters.xml
products/noise_param1_scan/runs/param1_*/run_simulation.rb
```

If matching hittree ROOT files already exist in those run directories, the script
will compare them to the experimental data immediately.

## Run simulations

Set:

```yaml
simulation:
  execute: true
```

Then run the same command. If the raw simulation output needs conversion to a
ComptonSoft `hittree`, set `simulation.postprocess_command` in the YAML. The
available placeholders are:

- `{run_dir}`
- `{sim_root}`
- `{sim_hittree}`

## Outputs

- `scan_results.csv`: param1, reduced chi-square, normalization, input files.
- `chi2_scan.png`: reduced chi-square vs param1.
- `best_summary.txt`: best param1.
- `best_comparison.png`: data spectrum and best normalized simulation spectrum.

By default, spectra use `event-sum` mode and only events with exactly one hit are
included, matching the "1-hit event" comparison requested for the resolution scan.
