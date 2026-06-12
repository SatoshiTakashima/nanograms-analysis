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

import argparse
import csv
import math
import os
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import analyze_testpulse_gain as tp


@dataclass
class ReviewState:
    peak_ok: bool | None = None
    use_fit: bool | None = None
    result: tp.FitResult | None = None
    message: str = ""


@dataclass
class InitialGuess:
    mean: float     = math.nan
    height: float   = math.nan
    width: float    = math.nan
    fit_window: float = math.nan


class RoundedButton:
    def __init__(self, fig, rect, label: str, callback,
                 facecolor: str = "#eef2f6", hovercolor: str = "#dde8f2"):
        from matplotlib.patches import FancyBboxPatch

        self.fig = fig
        self.ax  = fig.add_axes(rect)
        self.ax.set_axis_off()
        self.ax.set_xlim(0.0, 1.0)
        self.ax.set_ylim(0.0, 1.0)
        self.callback  = callback
        self.facecolor = facecolor
        self.base_facecolor = facecolor
        self.hovercolor = hovercolor
        self.enabled = True
        self.patch = FancyBboxPatch(
            (0.0, 0.0),
            1.0,
            1.0,
            boxstyle="round,pad=0.015,rounding_size=0.18",
            transform=self.ax.transAxes,
            linewidth=1.2,
            edgecolor="#c7d0da",
            facecolor=facecolor,
        )
        self.ax.add_patch(self.patch)
        self.label = self.ax.text(0.5, 0.5, label, ha="center", va="center", fontsize=11)
        self.cid = fig.canvas.mpl_connect("button_press_event", self.on_click)

    def on_click(self, event) -> None:
        if event.inaxes is self.ax and self.enabled:
            self.callback(event)

    def set_label(self, label: str) -> None:
        self.label.set_text(label)

    def set_facecolor(self, color: str) -> None:
        self.patch.set_facecolor(color)

    def set_enabled(self, enabled: bool) -> None:
        self.enabled = enabled
        self.patch.set_facecolor(self.base_facecolor if enabled else "#d8dee5")
        self.patch.set_edgecolor("#c7d0da" if enabled else "#c4cbd3")
        self.label.set_color("#111111" if enabled else "#7b848e")


class ToggleSwitch:
    def __init__(self, fig, label_pos, switch_rect, label: str, callback):
        from matplotlib.patches import Circle, FancyBboxPatch

        self.fig = fig
        self.ax = fig.add_axes(switch_rect)
        self.ax.set_axis_off()
        self.ax.set_xlim(0.0, 2.0)
        self.ax.set_ylim(0.0, 1.0)
        self.ax.set_aspect("equal")
        self.label_text = label
        self.callback = callback
        self.state: bool | None = None
        self.text = fig.text(label_pos[0], label_pos[1], label, ha="center", va="center", fontsize=11)
        self.track = FancyBboxPatch(
            (0.04, 0.08),
            1.92,
            0.84,
            boxstyle="round,pad=0.0,rounding_size=0.42",
            linewidth=1.0,
            edgecolor="#b7c0c8",
            facecolor="#cfd6dd",
        )
        self.knob = Circle((0.50, 0.50), 0.36,
                           facecolor="white", edgecolor="#b7c0c8", linewidth=1.0)
        self.ax.add_patch(self.track)
        self.ax.add_patch(self.knob)
        self.cid = fig.canvas.mpl_connect("button_press_event", self.on_click)

    def on_click(self, event) -> None:
        if event.inaxes is self.ax:
            self.callback(event)

    def set_state(self, state: bool | None, label: str | None = None) -> None:
        self.state = state
        if label is not None:
            self.text.set_text(label)
        is_on = state is True
        self.track.set_facecolor("#42c77b" if is_on else "#cfd6dd")
        self.track.set_edgecolor("#30a965" if is_on else "#b7c0c8")
        self.knob.center = (1.50 if is_on else 0.50, 0.50)


    #def parse_args() -> argparse.Namespace:
    #    parser = argparse.ArgumentParser(description="GUI review for NanoGRAMS test-pulse fits.")
    #    parser.add_argument("config", help="YAML configuration file.")
    #    parser.add_argument(
    #        "--start",
    #        default=None,
    #        help="Optional start item, for example 20251001/0124_04:FEC2.",
    #    )
    #    return parser.parse_args()


def parse_optional_bool(value: str | None) -> bool | None:
    return tp.parse_manual_bool(value)


def parse_optional_float(value: str | None) -> float:
    if value is None or str(value).strip() == "":
        return math.nan
    return float(value)


def load_review_states(cfg: tp.FitConfig) -> dict[tuple[str, int], ReviewState]:
    states: dict[tuple[str, int], ReviewState] = {}
    if not cfg.review_csv.exists():
        return states

    with cfg.review_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            time_id = row["time_id"]
            fec = int(row["fec"])
            key = (time_id, fec)
            use_fit = parse_optional_bool(row.get("use_fit"))
            peak_ok = parse_optional_bool(row.get("peak_ok"))
            if peak_ok is None and use_fit is not None:
                peak_ok = use_fit

            file_path = Path(row.get("file_path", ""))
            result = tp.FitResult(time_id=time_id, file_path=file_path, fec=fec)
            fit_passed = parse_optional_bool(row.get("fit_passed"))
            if fit_passed is None:
                fit_passed = parse_optional_bool(row.get("suggested_use_fit"))
            result.alive = bool(fit_passed) if fit_passed is not None else row.get("reason", "") == "ok"
            result.peak_adu = parse_optional_float(row.get("peak_adu"))
            result.peak_adu_err = parse_optional_float(row.get("peak_adu_err"))
            result.sigma_adu = parse_optional_float(row.get("sigma_adu"))
            result.sigma_adu_err = parse_optional_float(row.get("sigma_adu_err"))
            result.amplitude = parse_optional_float(row.get("amplitude"))
            result.baseline = parse_optional_float(row.get("baseline"))
            result.chi2_ndf = parse_optional_float(row.get("chi2_ndf"))
            result.n_events = int(float(row.get("n_events") or 0))
            result.n_fit = int(float(row.get("n_fit") or 0))
            result.reason = row.get("reason", "")
            result.fit_min_adu = parse_optional_float(row.get("fit_min_adu"))
            result.fit_max_adu = parse_optional_float(row.get("fit_max_adu"))
            states[key] = ReviewState(peak_ok=peak_ok, use_fit=use_fit, result=result)
    return states


def bool_to_cell(value: bool | None) -> str:
    if value is None:
        return ""
    return "1" if value else "0"


def result_value(result: tp.FitResult | None, name: str) -> float | int | str:
    if result is None:
        return ""
    return getattr(result, name)


def save_review_states(cfg: tp.FitConfig,
                       items: list[tuple[str, int]],
                       states: dict[tuple[str, int], ReviewState]) -> None:
    cfg.review_csv.parent.mkdir(parents=True, exist_ok=True)
    with cfg.review_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "time_id", "datetime", "fec", "peak_ok", "use_fit", "fit_passed",
            "peak_adu", "peak_adu_err", "sigma_adu", "sigma_adu_err",
            "amplitude", "baseline", "n_events", "n_fit", "chi2_ndf",
            "fit_min_adu", "fit_max_adu", "reason", "file_path",
        ])
        for time_id, fec in items:
            key = (time_id, fec)
            state = states.get(key, ReviewState())
            result = state.result
            file_path = cfg.data_root / time_id / cfg.file_name
            writer.writerow([
                time_id,
                tp.parse_time_id(time_id).isoformat(sep=" "),
                fec,
                bool_to_cell(state.peak_ok),
                bool_to_cell(state.use_fit),
                bool_to_cell(result.alive) if result is not None else "",
                result_value(result, "peak_adu"),
                result_value(result, "peak_adu_err"),
                result_value(result, "sigma_adu"),
                result_value(result, "sigma_adu_err"),
                result_value(result, "amplitude"),
                result_value(result, "baseline"),
                result_value(result, "n_events"),
                result_value(result, "n_fit"),
                result_value(result, "chi2_ndf"),
                result_value(result, "fit_min_adu"),
                result_value(result, "fit_max_adu"),
                result_value(result, "reason"),
                file_path,
            ])


def write_summary_from_review(cfg: tp.FitConfig,
                              items: list[tuple[str, int]],
                              states: dict[tuple[str, int], ReviewState]) -> None:
    by_time: dict[str, dict[int, ReviewState]] = {}
    for time_id, fec in items:
        by_time.setdefault(time_id, {})[fec] = states.get((time_id, fec), ReviewState())

    cfg.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with cfg.summary_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_id", "datetime", *[f"FEC{fec}" for fec in cfg.fec_ids]])
        for time_id in sorted(by_time, key=tp.parse_time_id):
            row = [time_id, tp.parse_time_id(time_id).isoformat(sep=" ")]
            for fec in cfg.fec_ids:
                state = by_time[time_id].get(fec, ReviewState())
                if state.peak_ok is True and state.result is not None and state.result.alive:
                    row.append(state.result.peak_adu)
                else:
                    row.append(float("nan"))
            writer.writerow(row)


class TestPulseReviewGUI:
    #def __init__(self, cfg: tp.FitConfig, start: str | None = None):
    def __init__(self, cfg: tp.FitConfig):
        import matplotlib

        if "MPLBACKEND" not in os.environ:
            try:
                matplotlib.use("TkAgg")
            except Exception:
                pass

        import matplotlib.pyplot as plt
        from matplotlib.widgets import TextBox

        self.plt = plt
        self.TextBox = TextBox
        self.cfg = cfg
        self.time_ids = tp.discover_time_ids(cfg)
        if not self.time_ids:
            raise RuntimeError("No test-pulse input files were found.")
        self.all_items = [
            (time_id, fec)
            for time_id in self.time_ids
            for fec in cfg.fec_ids
        ]
        self.items = self.all_items

        self.states = load_review_states(cfg)
        self.values_cache: dict[tuple[str, int], np.ndarray] = {}
        self.initial_boxes: dict[int, dict[str, TextBox]] = {}
        self.initial_guesses: dict[tuple[str, int], InitialGuess] = {}
        #self.time_index, self.selected_fec = self.start_position(start)
        self.time_index, self.selected_fec = 0, self.cfg.fec_ids[0]


        self.fig = plt.figure(figsize=(14.0, 9.0))
        self.title_text = self.fig.text(
            0.5, 0.985, "", ha="center", va="top", fontsize=13, weight="bold",
        )
        self.path_text = self.fig.text(
            0.5, 0.958, "", ha="center", va="top", fontsize=8.5,
        )
        self.fig.text(0.265, 0.915, "Full range", ha="center", va="center", fontsize=12, weight="bold")
        self.fig.text(0.735, 0.915, "Peak +/- 30 ADU", ha="center", va="center", fontsize=12, weight="bold")
        self.ax_wide: dict[int, object] = {}
        self.ax_zoom: dict[int, object] = {}
        self.axis_to_fec: dict[object, int] = {}
        left_x = [0.055, 0.285]
        right_x = [0.535, 0.765]
        y_pos = [0.700, 0.455]
        for index, fec in enumerate(self.cfg.fec_ids):
            row = index // 2
            col = index % 2
            self.ax_wide[fec] = self.fig.add_axes([left_x[col], y_pos[row], 0.195, 0.178])
            self.ax_zoom[fec] = self.fig.add_axes([right_x[col], y_pos[row], 0.195, 0.178])
            self.axis_to_fec[self.ax_wide[fec]] = fec
            self.axis_to_fec[self.ax_zoom[fec]] = fec

        self.status_ax = self.fig.add_axes([0.055, 0.255, 0.885, 0.135])
        self.status_ax.set_xticks([])
        self.status_ax.set_yticks([])
        self.status_ax.patch.set_facecolor("#f7f7f7")
        self.status_ax.patch.set_edgecolor("#cccccc")
        self.status_text = self.status_ax.text(
            0.01, 0.94, "", ha="left", va="top", fontsize=8.5,
            family="monospace", transform=self.status_ax.transAxes,
        )
        self.status_text.set_clip_on(True)

        self.controls = []
        self.peak_switches: dict[int, ToggleSwitch] = {}
        self.prev_button = self.add_button([0.055, 0.075, 0.095, 0.070], "Prev", self.on_prev)
        peak_positions = {
            0: ((0.205, 0.178), [0.245, 0.155, 0.060, 0.052]),
            1: ((0.205, 0.132), [0.245, 0.109, 0.060, 0.052]),
            2: ((0.205, 0.086), [0.245, 0.063, 0.060, 0.052]),
            3: ((0.205, 0.040), [0.245, 0.017, 0.060, 0.052]),
        }
        for fec in self.cfg.fec_ids:
            label_pos, switch_rect = peak_positions.get(
                fec,
                ((0.205, 0.178 - 0.046 * fec),
                 [0.245, 0.155 - 0.046 * fec, 0.060, 0.052]),
            )
            self.peak_switches[fec] = self.add_switch(
                label_pos,
                switch_rect,
                f"FEC{fec}",
                lambda event, fec=fec: self.on_toggle_peak(event, fec),
            )

        self.fig.text(0.350, 0.218, "Mean init", ha="left", va="bottom", fontsize=10)
        self.fig.text(0.435, 0.218, "Height init", ha="left", va="bottom", fontsize=10)
        self.fig.text(0.525, 0.218, "Width init", ha="left", va="bottom", fontsize=10)
        self.fig.text(0.615, 0.218, "Half Fit Window", ha="left", va="bottom", fontsize=10)
        input_y = {0: 0.170, 1: 0.124, 2: 0.078, 3: 0.032}
        for fec in self.cfg.fec_ids:
            y = input_y.get(fec, 0.170 - 0.046 * fec)
            self.initial_boxes[fec] = {
                "mean": self.add_textbox([0.350, y, 0.070, 0.032]),
                "height": self.add_textbox([0.435, y, 0.070, 0.032]),
                "width": self.add_textbox([0.525, y, 0.070, 0.032]),
                "fit_window": self.add_textbox([0.615, y, 0.070, 0.032]),
            }

        self.fit_button = self.add_button([0.735, 0.150, 0.075, 0.065], "Fit", self.on_fit, facecolor="#e8f1ff")
        self.save_button = self.add_button([0.825, 0.150, 0.125, 0.065], "Save summary", self.on_save_summary)
        self.next_button = self.add_button([0.910, 0.075, 0.075, 0.070], "Next", self.on_next)
        self.fig.canvas.mpl_connect("button_press_event", self.on_panel_click)
        self.set_guess_boxes_for_current_time(force=True)
        self.draw()

    def add_button(self, rect, label, callback,
                   facecolor: str = "#eef2f6") -> RoundedButton:
        button = RoundedButton(self.fig, rect, label, callback, facecolor=facecolor)
        self.controls.append(button)
        return button

    def add_switch(self, label_pos, switch_rect, label, callback) -> ToggleSwitch:
        switch = ToggleSwitch(self.fig, label_pos, switch_rect, label, callback)
        self.controls.append(switch)
        return switch

    def add_textbox(self, rect, label: str | None = None) -> TextBox:
        ax_box = self.fig.add_axes(rect)
        ax_box.set_facecolor("#f7f9fc")
        if label is not None:
            self.fig.text(rect[0], rect[1] + rect[3] + 0.014, label,
                          ha="left", va="bottom", fontsize=10)
        textbox = self.TextBox(ax_box, "", initial="")
        return textbox

        #def start_position(self, start: str | None) -> tuple[int, int]:
        #    if start is None:
        #        return 0, self.cfg.fec_ids[0]
        #    if ":FEC" not in start:
        #        raise ValueError("--start must look like 20251001/0124_04:FEC2")
        #    time_id, fec_text = start.split(":FEC", 1)
        #    fec = int(fec_text)
        #    item = (time_id, fec)
        #    if item not in self.all_items:
        #        raise ValueError(f"Start item is not in the configured input list: {start}")
        #    return self.time_ids.index(time_id), fec

    def current_key(self) -> tuple[str, int]:
        return self.time_ids[self.time_index], self.selected_fec

    def current_time_id(self) -> str:
        return self.time_ids[self.time_index]

    def current_file_path(self) -> Path:
        return self.cfg.data_root / self.current_time_id() / self.cfg.file_name

    def current_item_range(self) -> tuple[int, int]:
        start = self.time_index * len(self.cfg.fec_ids) + 1
        end = start + len(self.cfg.fec_ids) - 1
        return start, end

    def current_screen_image_path(self) -> Path:
        safe_time = self.current_time_id().replace("/", "_")
        return self.cfg.review_plot_dir / f"{safe_time}_gui.png"

    def save_current_screen_image(self) -> Path:
        image_path = self.current_screen_image_path()
        image_path.parent.mkdir(parents=True, exist_ok=True)
        self.fig.canvas.draw()
        self.fig.savefig(image_path, dpi=150)
        return image_path

    def current_state(self) -> ReviewState:
        key = self.current_key()
        self.states.setdefault(key, ReviewState())
        return self.states[key]

    def state_for(self, time_id: str, fec: int) -> ReviewState:
        key = (time_id, fec)
        self.states.setdefault(key, ReviewState())
        return self.states[key]

    def values_for(self, time_id: str, fec: int) -> np.ndarray:
        key = (time_id, fec)
        if key not in self.values_cache:
            file_path = self.cfg.data_root / time_id / self.cfg.file_name
            for loaded_fec, values in tp.read_fec_channel_values_by_fec(file_path, self.cfg).items():
                self.values_cache[(time_id, loaded_fec)] = values
        return self.values_cache[key]

    def current_values(self) -> np.ndarray:
        time_id, fec = self.current_key()
        return self.values_for(time_id, fec)

    def histogram_guess_for(self, time_id: str, fec: int) -> InitialGuess:
        values = self.values_for(time_id, fec)
        edges = np.arange(self.cfg.histogram_min_adu,
                          self.cfg.histogram_max_adu + self.cfg.histogram_bin_width_adu,
                          self.cfg.histogram_bin_width_adu)
        counts, edges = np.histogram(values, bins=edges)
        if counts.size == 0 or np.max(counts) <= 0:
            return InitialGuess()
        centers = 0.5 * (edges[:-1] + edges[1:])
        peak_index = int(np.argmax(counts))
        mean = float(centers[peak_index])
        height = float(counts[peak_index])
        near_peak = values[np.abs(values - mean) <= self.cfg.fit_window_adu]
        width = float(np.std(near_peak)) if near_peak.size > 1 else math.nan
        if not np.isfinite(width) or width <= 0.0:
            width = max(1.0, self.cfg.sigma_min_adu)
        width = float(np.clip(width, self.cfg.sigma_min_adu, self.cfg.sigma_max_adu))
        return InitialGuess(mean=mean, height=height, width=width, fit_window=self.cfg.fit_window_adu)

    def default_initial_guess_for(self, time_id: str, fec: int) -> InitialGuess:
        state = self.state_for(time_id, fec)
        if state.result is not None and np.isfinite(state.result.peak_adu):
            fit_window = self.cfg.fit_window_adu
            if np.isfinite(state.result.fit_min_adu) and np.isfinite(state.result.fit_max_adu):
                fit_window = 0.5 * (state.result.fit_max_adu - state.result.fit_min_adu)
            return InitialGuess(
                mean=state.result.peak_adu,
                height=state.result.amplitude,
                width=state.result.sigma_adu,
                fit_window=fit_window,
            )
        return self.histogram_guess_for(time_id, fec)

    def set_guess_boxes_for_current_time(self, force: bool = False) -> None:
        time_id = self.current_time_id()
        for fec in self.cfg.fec_ids:
            key = (time_id, fec)
            if force or key not in self.initial_guesses:
                self.initial_guesses[key] = self.default_initial_guess_for(time_id, fec)
            self.set_guess_boxes_for_fec(fec, self.initial_guesses[key])
        self.update_input_enabled_states()

    def set_guess_boxes_for_fec(self, fec: int, guess: InitialGuess) -> None:
        boxes = self.initial_boxes[fec]
        boxes["mean"].set_val("" if not np.isfinite(guess.mean) else f"{guess.mean:.3f}")
        boxes["height"].set_val("" if not np.isfinite(guess.height) else f"{guess.height:.3f}")
        boxes["width"].set_val("" if not np.isfinite(guess.width) else f"{guess.width:.3f}")
        boxes["fit_window"].set_val("" if not np.isfinite(guess.fit_window) else f"{guess.fit_window:.3f}")

    def initial_guess_from_boxes(self, fec: int) -> InitialGuess:
        boxes = self.initial_boxes[fec]
        return InitialGuess(
            mean=parse_optional_float(getattr(boxes["mean"], "text", "")),
            height=parse_optional_float(getattr(boxes["height"], "text", "")),
            width=parse_optional_float(getattr(boxes["width"], "text", "")),
            fit_window=parse_optional_float(getattr(boxes["fit_window"], "text", "")),
        )

    def sync_initial_guess_from_boxes(self, fec: int) -> InitialGuess:
        guess = self.initial_guess_from_boxes(fec)
        self.initial_guesses[(self.current_time_id(), fec)] = guess
        return guess

    def sync_current_time_initial_guesses_from_boxes(self) -> None:
        time_id = self.current_time_id()
        for fec in self.cfg.fec_ids:
            self.initial_guesses[(time_id, fec)] = self.initial_guess_from_boxes(fec)

    def set_textbox_enabled(self, textbox: TextBox, enabled: bool) -> None:
        textbox.set_active(enabled)
        textbox.ax.set_facecolor("#ffffff" if enabled else "#e9ecef")
        try:
            textbox.text_disp.set_color("#111111" if enabled else "#888888")
        except AttributeError:
            pass

    def update_input_enabled_states(self) -> None:
        time_id = self.current_time_id()
        for fec in self.cfg.fec_ids:
            enabled = self.state_for(time_id, fec).peak_ok is True
            for textbox in self.initial_boxes[fec].values():
                self.set_textbox_enabled(textbox, enabled)

    def update_button_labels(self) -> None:
        time_id = self.current_time_id()
        for fec, switch in self.peak_switches.items():
            switch.set_state(self.state_for(time_id, fec).peak_ok, f"FEC{fec}")
        self.update_input_enabled_states()
        self.update_navigation_buttons()

    def update_navigation_buttons(self) -> None:
        self.prev_button.set_enabled(self.time_index > 0)
        self.next_button.set_enabled(self.time_index < len(self.time_ids) - 1)

    def set_button_feedback(self, button: RoundedButton, color: str) -> None:
        button.set_facecolor(color)
        self.fig.canvas.draw_idle()
        try:
            self.fig.canvas.flush_events()
        except Exception:
            pass

    def reset_action_button_colors(self) -> None:
        self.fit_button.set_facecolor(self.fit_button.base_facecolor)
        self.save_button.set_facecolor(self.save_button.base_facecolor)

    def on_panel_click(self, event) -> None:
        if event.inaxes not in self.axis_to_fec:
            return
        self.selected_fec = self.axis_to_fec[event.inaxes]
        self.draw()

    def plot_histogram_panel(self,
                             ax,
                             counts: np.ndarray,
                             edges: np.ndarray,
                             title: str,
                             xlim: tuple[float, float],
                             result: tp.FitResult | None,
                             selected: bool,
                             draw_fit: bool,
                             show_xlabel: bool,
                             show_ylabel: bool,
                             show_legend: bool = False) -> None:
        ax.clear()
        ax.stairs(counts, edges, label="histogram")
        if draw_fit and result is not None and np.isfinite(result.peak_adu) and np.isfinite(result.sigma_adu):
            fit_min = result.fit_min_adu if np.isfinite(result.fit_min_adu) else result.peak_adu - self.cfg.fit_window_adu
            fit_max = result.fit_max_adu if np.isfinite(result.fit_max_adu) else result.peak_adu + self.cfg.fit_window_adu
            x = np.linspace(fit_min, fit_max, 400)
            if np.isfinite(result.amplitude) and np.isfinite(result.baseline):
                ax.plot(
                    x,
                    tp.gaussian_with_const(x,
                                           result.amplitude,
                                           result.peak_adu,
                                           result.sigma_adu,
                                           result.baseline),
                    color="tab:red",
                    lw=2.0,
                    label="fit function",
                )
            ax.axvline(result.peak_adu, color="tab:red", ls="--", alpha=0.75)
            ax.axvspan(fit_min, fit_max, color="tab:red", alpha=0.05, label="fit range")
            ax.text(
                0.03,
                0.91,
                f"Mean {result.peak_adu:.2f} ADU",
                ha="left",
                va="top",
                color="tab:red",
                fontsize=8,
                transform=ax.transAxes,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.70, "pad": 1.5},
            )

        ax.set_title(title, fontsize=9, pad=2)
        ax.set_xlabel("ADU - CMN" if show_xlabel else "", fontsize=7)
        ax.set_ylabel("Counts" if show_ylabel else "", fontsize=7)
        ax.set_xlim(*xlim)
        ax.grid(alpha=0.22, linewidth=0.6)
        ax.tick_params(labelsize=7, labelbottom=True, labelleft=True)
        for spine in ax.spines.values():
            spine.set_linewidth(2.2 if selected else 1.0)
            spine.set_edgecolor("#f08c00" if selected else "black")
        if show_legend:
            ax.legend(loc="upper right", fontsize=6, framealpha=0.85)

    def draw(self) -> None:
        time_id = self.current_time_id()
        file_path = self.current_file_path()
        self.title_text.set_text(
            f"Test-pulse ch{self.cfg.test_pulse_channel}  |  {time_id}  "
            f"({self.time_index + 1}/{len(self.time_ids)} files)"
        )
        self.path_text.set_text(self.format_display_path(file_path))
        try:
            self.fig.canvas.manager.set_window_title(
                f"Test-pulse ch{self.cfg.test_pulse_channel}: {time_id}"
            )
        except Exception:
            pass

        edges = np.arange(self.cfg.histogram_min_adu,
                          self.cfg.histogram_max_adu + self.cfg.histogram_bin_width_adu,
                          self.cfg.histogram_bin_width_adu)

        for panel_index, fec in enumerate(self.cfg.fec_ids):
            row = panel_index // 2
            col = panel_index % 2
            show_xlabel = row == 1
            show_ylabel = col == 0
            state = self.state_for(time_id, fec)
            values = self.values_for(time_id, fec)
            counts, hist_edges = np.histogram(values, bins=edges)
            centers = 0.5 * (hist_edges[:-1] + hist_edges[1:])
            if state.result is not None and np.isfinite(state.result.peak_adu):
                zoom_center = state.result.peak_adu
            elif counts.size > 0 and np.max(counts) > 0:
                zoom_center = float(centers[int(np.argmax(counts))])
            else:
                zoom_center = 0.0

            is_peak_on = state.peak_ok is True
            self.plot_histogram_panel(
                self.ax_wide[fec],
                counts,
                hist_edges,
                f"FEC{fec}",
                (self.cfg.histogram_min_adu, self.cfg.histogram_max_adu),
                state.result,
                selected=False,
                draw_fit=False,
                show_xlabel=show_xlabel,
                show_ylabel=show_ylabel,
                show_legend=False,
            )
            self.plot_histogram_panel(
                self.ax_zoom[fec],
                counts,
                hist_edges,
                f"FEC{fec}",
                (zoom_center - 30.0, zoom_center + 30.0),
                state.result,
                selected=is_peak_on,
                draw_fit=True,
                show_xlabel=show_xlabel,
                show_ylabel=show_ylabel,
                show_legend=True,
            )

        self.status_text.set_text(self.status_message())
        self.update_button_labels()
        self.fig.canvas.draw_idle()

    def status_message(self) -> str:
        state = self.current_state()
        time_id = self.current_time_id()
        item_start, item_end = self.current_item_range()
        active_fecs = [
            fec
            for fec in self.cfg.fec_ids
            if self.state_for(time_id, fec).peak_ok is True
        ]
        fitted_active = [
            fec
            for fec in active_fecs
            if self.state_for(time_id, fec).result is not None
        ]
        accepted = sum(
            1
            for item in self.all_items
            if self.states.get(item, ReviewState()).peak_ok is True
            and self.states.get(item, ReviewState()).result is not None
            and self.states.get(item, ReviewState()).result.alive
        )
        decided = sum(1 for item in self.all_items if self.states.get(item, ReviewState()).peak_ok is not None)
        result = state.result
        guess = self.initial_guess_from_boxes(self.selected_fec)
        active_text = ", ".join(f"FEC{fec}" for fec in active_fecs) if active_fecs else "none"
        lines = [
            f"file={self.time_index + 1}/{len(self.time_ids)}, FEC items={item_start}-{item_end}/{len(self.all_items)}, selected=FEC{self.selected_fec}",
            f"Peak ON: {active_text}; fitted active FECs={len(fitted_active)}/{len(active_fecs)}, decided={decided}/{len(self.all_items)}, accepted={accepted}",
            f"selected initial: mean={guess.mean:.3g}, height={guess.height:.3g}, width={guess.width:.3g}, fit_window={guess.fit_window:.3g}",
        ]
        if result is not None:
            lines.append(
                f"fit: peak={result.peak_adu:.3g} +/- {result.peak_adu_err:.2g} ADU, "
                f"sigma={result.sigma_adu:.3g} +/- {result.sigma_adu_err:.2g}, "
                f"chi2/ndf={result.chi2_ndf:.3g}, reason={result.reason}"
            )
        if state.message:
            lines.append(state.message)
        return self.wrap_status_lines(lines)

    def format_display_path(self, file_path: Path) -> str:
        text = str(file_path)
        if len(text) <= 130:
            return text
        split_at = text.rfind("/", 0, 100)
        if split_at < 50:
            split_at = 100
        return text[:split_at + 1] + "\n" + text[split_at + 1:]

    def wrap_status_lines(self, lines: list[str]) -> str:
        wrapped: list[str] = []
        for text in lines:
            for line in str(text).splitlines():
                if len(line) <= 120:
                    wrapped.append(line)
                    continue
                wrapped.extend(textwrap.wrap(line, width=120, subsequent_indent="  "))
        return "\n".join(wrapped[:6])

    def unfitted_peak_fecs_for_current_time(self) -> list[int]:
        time_id = self.current_time_id()
        missing: list[int] = []
        for fec in self.cfg.fec_ids:
            state = self.state_for(time_id, fec)
            if state.peak_ok is not True:
                continue
            if state.result is None or not np.isfinite(state.result.peak_adu):
                missing.append(fec)
        return missing

    def failed_peak_fecs_for_current_time(self) -> list[int]:
        time_id = self.current_time_id()
        failed: list[int] = []
        for fec in self.cfg.fec_ids:
            state = self.state_for(time_id, fec)
            if state.peak_ok is not True or state.result is None:
                continue
            if not state.result.alive:
                failed.append(fec)
        return failed

    def persist(self) -> None:
        save_review_states(self.cfg, self.all_items, self.states)

    def move(self, step: int) -> None:
        self.sync_current_time_initial_guesses_from_boxes()
        self.time_index = max(0, min(len(self.time_ids) - 1, self.time_index + step))
        self.reset_action_button_colors()
        self.set_guess_boxes_for_current_time(force=False)
        self.draw()

    def on_prev(self, _event) -> None:
        self.move(-1)

    def on_next(self, _event) -> None:
        self.move(1)

    def on_toggle_peak(self, _event, fec: int | None = None) -> None:
        if fec is not None:
            self.selected_fec = fec
        state = self.current_state()
        if state.peak_ok is not True:
            state.peak_ok = True
            state.use_fit = None
            state.message = f"FEC{self.selected_fec} peak marked OK. Press Fit."
        else:
            state.peak_ok = False
            state.use_fit = False
            state.message = f"FEC{self.selected_fec} peak marked NG. This FEC is rejected before fitting."
        self.update_input_enabled_states()
        self.persist()
        self.draw()

    def on_fit(self, _event) -> None:
        self.set_button_feedback(self.fit_button, "#ffe8a3")
        time_id = self.current_time_id()
        active_fecs = [
            fec
            for fec in self.cfg.fec_ids
            if self.state_for(time_id, fec).peak_ok is True
        ]
        if not active_fecs:
            self.current_state().message = "Turn on at least one Peak switch before fitting."
            self.set_button_feedback(self.fit_button, "#ffd6d6")
            self.draw()
            return

        file_path = self.cfg.data_root / time_id / self.cfg.file_name
        messages: list[str] = []
        for fec in active_fecs:
            state = self.state_for(time_id, fec)
            guess = self.sync_initial_guess_from_boxes(fec)
            state.result = tp.fit_fec(
                time_id,
                file_path,
                self.cfg,
                fec,
                initial_mean_adu=guess.mean,
                initial_height=guess.height,
                initial_width_adu=guess.width,
                fit_window_adu=guess.fit_window,
                values=self.values_for(time_id, fec),
            )
            if state.result is not None and np.isfinite(state.result.peak_adu):
                self.initial_guesses[(time_id, fec)] = self.default_initial_guess_for(time_id, fec)
            state.use_fit = bool(state.result.alive)
            if state.result.reason == "fit_failed:iminuit_missing":
                messages.append(
                    f"FEC{fec}: iminuit is missing in {sys.executable}."
                )
            else:
                messages.append(f"FEC{fec}: {state.result.reason}")
        if active_fecs:
            self.selected_fec = active_fecs[-1]
        if any("iminuit is missing" in message for message in messages):
            self.current_state().message = (
                "Fit failed because this Python cannot import iminuit. "
                "Install iminuit there, or run this GUI with the Python that has iminuit."
            )
            self.set_button_feedback(self.fit_button, "#ffd6d6")
        else:
            self.current_state().message = "Fit done for Peak ON FECs: " + "; ".join(messages)
            self.set_button_feedback(self.fit_button, "#cfe4ff")
        self.set_guess_boxes_for_current_time(force=False)
        self.persist()
        self.draw()

    def on_save_summary(self, _event) -> None:
        self.set_button_feedback(self.save_button, "#ffe8a3")
        missing = self.unfitted_peak_fecs_for_current_time()
        if missing:
            state = self.current_state()
            missing_text = ", ".join(f"FEC{fec}" for fec in missing)
            state.message = (
                f"Save blocked: Peak is ON but fit is missing for {missing_text}. "
                "Run Fit for them, or turn Peak OFF before saving."
            )
            self.set_button_feedback(self.save_button, "#ffd6d6")
            self.draw()
            return

        failed = self.failed_peak_fecs_for_current_time()
        if failed:
            state = self.current_state()
            failed_text = ", ".join(
                f"FEC{fec}({self.state_for(self.current_time_id(), fec).result.reason})"
                for fec in failed
            )
            state.message = (
                f"Save blocked: Peak is ON but fit did not pass quality cuts for {failed_text}. "
                "Turn Peak OFF or adjust/refit before saving."
            )
            self.set_button_feedback(self.save_button, "#ffd6d6")
            self.draw()
            return

        self.persist()
        write_summary_from_review(self.cfg, self.all_items, self.states)
        state = self.current_state()
        image_path = self.current_screen_image_path()
        state.message = (
            f"Saved summary CSV for all FECs: {self.cfg.summary_csv}\n"
            f"Saved GUI screen image: {image_path}"
        )
        self.set_button_feedback(self.save_button, "#d7f5dd")
        self.draw()
        self.save_current_screen_image()

    def show(self) -> None:
        self.plt.show()


    #def main() -> None:
    #args = parse_args()
    #cfg = tp.load_config(args.config)
    #cfg.outdir.mkdir(parents=True, exist_ok=True)
    #mpl_cache_dir = cfg.outdir / ".matplotlib"
    #mpl_cache_dir.mkdir(parents=True, exist_ok=True)
    #os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache_dir))

    #app = TestPulseReviewGUI(cfg, args.start)
    #app.show()


if __name__ == "__main__":
    #main()
    config_file_path = "config.yaml"
    cfg = tp.load_config(config_file_path)

    cfg.outdir.mkdir(parents=True, exist_ok=True)

    mpl_cache_dir = cfg.outdir / ".matplotlib"
    mpl_cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache_dir))

    #app = TestPulseReviewGUI(cfg, args.start)
    app = TestPulseReviewGUI(cfg)
    app.show()
