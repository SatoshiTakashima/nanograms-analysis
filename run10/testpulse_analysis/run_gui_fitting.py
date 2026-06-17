#!/usr/bin/env python3
"""Click-through GUI for NanoGRAMS test-pulse gain review.

The intended flow is deliberately human-in-the-loop:

  1. inspect the CMN-subtracted histogram for one time/FEC,
  2. mark whether a visible peak exists,
  3. run the Gaussian fit only for peak-like histograms,
  4. inspect the fit overlay,
  5. accept only fits that also pass the configured quality cuts.

Accepted fits are written to the same compact gain table used by the
interpolation script: time_id,datetime,FEC0,FEC1,FEC2,FEC3.
"""

from __future__ import annotations

import os
import sys
sys.path.append("../../mymodule")
import analyze_testpulse_gain as tp
from fit_testpulse_data_gui import TestPulseReviewGUI

if __name__ == "__main__":
    config_file_path = "../metadata/config_testpulse_fit.yaml"
    cfg = tp.load_config(config_file_path)

    cfg.outdir.mkdir(parents=True, exist_ok=True)

    mpl_cache_dir = cfg.outdir / ".matplotlib"
    mpl_cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache_dir))

    app = TestPulseReviewGUI(cfg)
    app.show()
