#!/usr/bin/env python3
"""Lightweight customtkinter GUI for NanoGRAMS test-pulse gain review."""

from __future__ import annotations

import math
import tkinter as tk
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

import customtkinter as ctk
import numpy as np

import analyze_testpulse_gain as tp
from fit_testpulse_data_gui import (
    InitialGuess,
    ReviewState,
    load_review_states,
    parse_optional_float,
    save_review_states,
    write_summary_from_review,
)


@dataclass
class HistogramData:
    counts: np.ndarray
    edges: np.ndarray
    centers: np.ndarray
    peak_center: float = math.nan
    peak_count: float = 0.0
    peak_width: float = math.nan


class HistogramPanel:
    def __init__(self, parent, title: str, width: int = 310, height: int = 170):
        self.frame = ctk.CTkFrame(parent, fg_color="#f8fafc", corner_radius=8)
        self.label = ctk.CTkLabel(
            self.frame,
            text=title,
            anchor="center",
            font=ctk.CTkFont(size=18, weight="bold"),
            text_color="#172033",
        )
        self.canvas = tk.Canvas(
            self.frame,
            width=width,
            height=height,
            bg="#ffffff",
            highlightthickness=1,
            highlightbackground="#d6dee8",
        )
        self.label.pack(fill="x", padx=8, pady=(6, 0))
        self.canvas.pack(fill="both", expand=True, padx=8, pady=(2, 6))
        self.width = width
        self.height = height
        self._click_callback = None
        self._last_draw_args = None
        self.canvas.bind("<Configure>", self.on_resize)

    def bind_click(self, callback):
        self._click_callback = callback
        self.canvas.bind("<Button-1>", lambda _event: callback())

    def on_resize(self, _event) -> None:
        if self._last_draw_args is not None:
            self.draw(*self._last_draw_args)

    @staticmethod
    def nice_step(span: float, target_ticks: int = 5) -> float:
        if span <= 0.0 or not np.isfinite(span):
            return 1.0
        raw_step = span / max(target_ticks, 1)
        exponent = math.floor(math.log10(raw_step))
        fraction = raw_step / (10.0 ** exponent)
        if fraction <= 1.0:
            nice_fraction = 1.0
        elif fraction <= 2.0:
            nice_fraction = 2.0
        elif fraction <= 5.0:
            nice_fraction = 5.0
        else:
            nice_fraction = 10.0
        return nice_fraction * (10.0 ** exponent)

    @classmethod
    def nice_ticks(cls, min_value: float, max_value: float, target_ticks: int = 5) -> list[float]:
        step = cls.nice_step(max_value - min_value, target_ticks)
        start = math.ceil(min_value / step) * step
        stop = math.floor(max_value / step) * step
        ticks = []
        value = start
        while value <= stop + 0.5 * step:
            ticks.append(0.0 if abs(value) < 1.0e-12 else value)
            value += step
        if not ticks:
            ticks = [min_value, max_value]
        return ticks

    def draw(self, histogram: HistogramData, xlim: tuple[float, float],
             result: tp.FitResult | None = None, selected: bool = False,
             draw_fit: bool = False) -> None:
        self._last_draw_args = (histogram, xlim, result, selected, draw_fit)
        c = self.canvas
        c.delete("all")
        w = max(c.winfo_width(), self.width)
        h = max(c.winfo_height(), self.height)
        left, right, top, bottom = 58, w - 18, 18, h - 38
        plot_w = max(1, right - left)
        plot_h = max(1, bottom - top)

        x_min, x_max = xlim
        centers = histogram.centers
        counts = histogram.counts
        mask = (centers >= x_min) & (centers <= x_max)
        ymax_data = float(np.max(counts[mask])) if np.any(mask) else float(np.max(counts) if counts.size else 1.0)
        ymax = 1.1 * ymax_data
        ymax = max(ymax, 1.0)

        border = "#f59e0b" if selected else "#64748b"
        c.create_rectangle(left, top, right, bottom, outline=border, width=3 if selected else 1)

        def sx(x):
            return left + (x - x_min) / (x_max - x_min) * plot_w

        def sy(y):
            return bottom - y / ymax * plot_h

        for i in range(5):
            y = top + plot_h * i / 4.0
            c.create_line(left, y, right, y, fill="#e8eef5")
        for value in self.nice_ticks(0.0, ymax, target_ticks=5):
            y = bottom - value / ymax * plot_h
            c.create_line(left - 4, y, left, y, fill="#475569")
            c.create_text(left - 8, y, text=f"{value:.0f}", anchor="e", fill="#475569", font=("TkDefaultFont", 8))
        for value in self.nice_ticks(x_min, x_max, target_ticks=5):
            x = sx(value)
            c.create_line(x, bottom, x, bottom + 4, fill="#475569")
            c.create_text(x, bottom + 15, text=f"{value:.0f}", fill="#475569", font=("TkDefaultFont", 8))
        c.create_text(17, top + plot_h / 2, text="Counts", angle=90, fill="#334155", font=("TkDefaultFont", 9))
        c.create_text(left + plot_w / 2, h - 10, text="ADU - CMN", fill="#334155", font=("TkDefaultFont", 9))

        points: list[float] = []
        edges = histogram.edges
        visible_bins = np.flatnonzero((edges[:-1] < x_max) & (edges[1:] > x_min))
        for bin_index in visible_bins:
            x0 = max(float(edges[bin_index]), x_min)
            x1 = min(float(edges[bin_index + 1]), x_max)
            y = float(counts[bin_index])
            if not points:
                points.extend([sx(x0), sy(y)])
            else:
                points.extend([sx(x0), points[-1]])
                points.extend([sx(x0), sy(y)])
            points.extend([sx(x1), sy(y)])
        if len(points) >= 4:
            c.create_line(*points, fill="#2563eb", width=1)

        legend_items = [("Data", "#2563eb")]
        if draw_fit and result is not None and np.isfinite(result.peak_adu) and np.isfinite(result.sigma_adu):
            fit_min = result.fit_min_adu if np.isfinite(result.fit_min_adu) else result.peak_adu - 5.0
            fit_max = result.fit_max_adu if np.isfinite(result.fit_max_adu) else result.peak_adu + 5.0
            x0, x1 = max(x_min, fit_min), min(x_max, fit_max)
            if x0 < x1 and np.isfinite(result.amplitude) and np.isfinite(result.baseline):
                x = np.linspace(x0, x1, 160)
                y = tp.gaussian_with_const(x, result.amplitude, result.peak_adu, result.sigma_adu, result.baseline)
                fit_points: list[float] = []
                for xx, yy in zip(x, y):
                    fit_points.extend([sx(float(xx)), sy(float(min(max(yy, 0.0), ymax)))])
                if len(fit_points) >= 4:
                    c.create_line(*fit_points, fill="#dc2626", width=2)
                    legend_items.append(("Fit", "#dc2626"))
            if x_min <= result.peak_adu <= x_max:
                xp = sx(result.peak_adu)
                c.create_line(xp, top, xp, bottom, fill="#dc2626", dash=(4, 3))
                c.create_text(left + 8, top + 8, text=f"Mean {result.peak_adu:.1f}", anchor="nw",
                              fill="#dc2626", font=("TkDefaultFont", 9, "bold"))

        if draw_fit:
            legend_x = right - 78
            legend_y = top + 10
            c.create_rectangle(
                legend_x - 8,
                legend_y - 7,
                right - 8,
                legend_y + 17 * len(legend_items) + 2,
                fill="#ffffff",
                outline="#d6dee8",
            )
            for index, (label, color) in enumerate(legend_items):
                y = legend_y + 17 * index
                c.create_line(legend_x, y, legend_x + 18, y, fill=color, width=2)
                c.create_text(legend_x + 24, y, text=label, anchor="w", fill="#334155", font=("TkDefaultFont", 8))


class TestPulseReviewTkGUI:
    def __init__(self, cfg: tp.FitConfig):
        ctk.set_appearance_mode("light")
        ctk.set_default_color_theme("blue")
        self.cfg = cfg
        self.time_ids = tp.discover_time_ids(cfg)
        if not self.time_ids:
            raise RuntimeError("No test-pulse input files were found.")
        self.all_items = [(time_id, fec) for time_id in self.time_ids for fec in cfg.fec_ids]
        self.states = load_review_states(cfg)
        self.values_cache: dict[tuple[str, int], np.ndarray] = {}
        self.histogram_cache: dict[tuple[str, int], HistogramData] = {}
        self.cached_time_ids: OrderedDict[str, None] = OrderedDict()
        self.max_cached_times = 3
        self.initial_guesses: dict[tuple[str, int], InitialGuess] = {}
        self.time_index = 0
        self.selected_fec = cfg.fec_ids[0]

        self.root = ctk.CTk()
        self.root.title("NanoGRAMS test-pulse review")
        self.root.geometry("1500x980")
        self.root.minsize(1260, 820)
        self.root.configure(fg_color="#edf2f7")

        self.peak_vars: dict[int, tk.BooleanVar] = {}
        self.entries: dict[int, dict[str, tk.StringVar]] = {}
        self.wide_panels: dict[int, HistogramPanel] = {}
        self.zoom_panels: dict[int, HistogramPanel] = {}
        self.message_var = tk.StringVar()
        self.title_var = tk.StringVar()
        self.path_var = tk.StringVar()
        self.input_file_var = tk.StringVar()
        self.review_file_var = tk.StringVar()
        self.summary_file_var = tk.StringVar()
        self.entry_widgets: dict[int, dict[str, ctk.CTkEntry]] = {}
        self.first_button: ctk.CTkButton | None = None
        self.prev_button: ctk.CTkButton | None = None
        self.next_button: ctk.CTkButton | None = None
        self.last_button: ctk.CTkButton | None = None
        self.save_button: ctk.CTkButton | None = None

        self.build_layout()
        self.set_guess_boxes_for_current_time(force=True)
        self.draw()

    def build_layout(self) -> None:
        top = ctk.CTkFrame(self.root, fg_color="#ffffff", corner_radius=0)
        top.pack(fill="x", padx=0, pady=0)
        ctk.CTkLabel(
            top,
            textvariable=self.title_var,
            font=ctk.CTkFont(size=18, weight="bold"),
            text_color="#111827",
        ).pack(pady=(10, 0))
        ctk.CTkLabel(
            top,
            textvariable=self.path_var,
            font=ctk.CTkFont(size=11),
            text_color="#475569",
        ).pack(pady=(2, 6))
        file_frame = ctk.CTkFrame(top, fg_color="#f8fafc", corner_radius=10)
        file_frame.pack(fill="x", padx=14, pady=(0, 8))
        for col in range(3):
            file_frame.grid_columnconfigure(col, weight=1)
        self.add_file_label(file_frame, 0, "Reading ROOT", self.input_file_var, "#2563eb")
        self.add_file_label(file_frame, 1, "Writing review CSV", self.review_file_var, "#0f766e")
        self.add_file_label(file_frame, 2, "Save summary CSV", self.summary_file_var, "#7c3aed")

        body = ctk.CTkFrame(self.root, fg_color="transparent")
        body.pack(fill="both", expand=True, padx=14, pady=(8, 8))
        body.grid_columnconfigure(0, weight=1)
        body.grid_columnconfigure(1, weight=1)
        body.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(
            body,
            text="Full range",
            anchor="center",
            font=ctk.CTkFont(size=24, weight="bold"),
            text_color="#172033",
        ).grid(row=0, column=0, sticky="ew", padx=(6, 12), pady=(0, 6))
        ctk.CTkLabel(
            body,
            text="Peak +/- 30 ADU",
            anchor="center",
            font=ctk.CTkFont(size=24, weight="bold"),
            text_color="#172033",
        ).grid(row=0, column=1, sticky="ew", padx=(12, 6), pady=(0, 6))

        wide = ctk.CTkFrame(body, fg_color="#ffffff", corner_radius=10)
        zoom = ctk.CTkFrame(body, fg_color="#ffffff", corner_radius=10)
        wide.grid(row=1, column=0, sticky="nsew", padx=(0, 7))
        zoom.grid(row=1, column=1, sticky="nsew", padx=(7, 0))

        for index, fec in enumerate(self.cfg.fec_ids):
            row = index // 2
            col = index % 2
            self.wide_panels[fec] = HistogramPanel(wide, f"FEC{fec}")
            self.zoom_panels[fec] = HistogramPanel(zoom, f"FEC{fec}")
            self.wide_panels[fec].frame.grid(row=row, column=col, sticky="nsew", padx=4, pady=4)
            self.zoom_panels[fec].frame.grid(row=row, column=col, sticky="nsew", padx=4, pady=4)
            self.wide_panels[fec].bind_click(lambda fec=fec: self.select_fec(fec))
            self.zoom_panels[fec].bind_click(lambda fec=fec: self.select_fec(fec))
        for container in (wide, zoom):
            for i in range(2):
                container.columnconfigure(i, weight=1)
                container.rowconfigure(i, weight=1)

        controls = ctk.CTkFrame(self.root, fg_color="#ffffff", corner_radius=12, height=146)
        controls.pack(fill="x", padx=14, pady=(0, 8))
        controls.pack_propagate(False)
        controls.grid_columnconfigure(0, weight=0)
        controls.grid_columnconfigure(1, weight=0)
        controls.grid_columnconfigure(2, weight=1)
        controls.grid_rowconfigure(0, weight=1)

        action_frame = ctk.CTkFrame(controls, fg_color="transparent")
        action_frame.grid(row=0, column=0, sticky="nsw", padx=(12, 8), pady=10)
        for col in range(5):
            action_frame.grid_columnconfigure(col, weight=0)

        self.first_button = ctk.CTkButton(action_frame, text="First", command=self.on_first, width=72, height=42, font=ctk.CTkFont(size=18))
        self.prev_button = ctk.CTkButton(action_frame, text="Prev", command=self.on_prev, width=72, height=42, font=ctk.CTkFont(size=18))
        self.prev_button.grid(row=0, column=0, padx=(0, 6), pady=(20, 2), sticky="ew")
        self.first_button.grid(row=1, column=0, padx=(0, 6), pady=(2, 20), sticky="ew")

        switch_frame = ctk.CTkFrame(action_frame, fg_color="#f8fafc", corner_radius=10)
        switch_frame.grid(row=0, column=1, rowspan=2, padx=(0, 12), pady=0, sticky="ns")
        for fec in self.cfg.fec_ids:
            var = tk.BooleanVar(value=False)
            self.peak_vars[fec] = var
            ctk.CTkSwitch(
                switch_frame,
                text=f"FEC{fec}",
                variable=var,
                command=lambda fec=fec: self.on_toggle_peak(fec),
                button_color="#94a3b8",
                button_hover_color="#64748b",
                progress_color="#2563eb",
                font=ctk.CTkFont(size=14, weight="bold"),
                text_color="#172033",
            ).pack(anchor="w", padx=10, pady=(6 if fec == self.cfg.fec_ids[0] else 2, 6))

        ctk.CTkButton(action_frame, text="Fit", command=self.on_fit, width=72, height=42, font=ctk.CTkFont(size=18)).grid(
            row=0, column=2, padx=(0, 7), pady=(20, 2), sticky="ew"
        )
        self.save_button = ctk.CTkButton(
            action_frame,
            text="Save summary",
            command=self.on_save_summary,
            width=118,
            height=42,
            fg_color="#0f766e",
            hover_color="#115e59",
            font=ctk.CTkFont(size=18)
        )
        self.save_button.grid(row=0, column=3, padx=(0, 7), pady=(20, 2), sticky="ew")
        self.next_button = ctk.CTkButton(action_frame, text="Next", command=self.on_next, width=62, height=42, font=ctk.CTkFont(size=18))
        self.last_button = ctk.CTkButton(action_frame, text="Last", command=self.on_last, width=62, height=42, font=ctk.CTkFont(size=18))
        self.next_button.grid(row=0, column=4, padx=(0, 0), pady=(20, 2), sticky="ew")
        self.last_button.grid(row=1, column=4, padx=(0, 0), pady=(2, 20), sticky="ew")

        entry_frame = ctk.CTkFrame(controls, fg_color="#f8fafc", corner_radius=10)
        entry_frame.grid(row=0, column=1, sticky="nsew", padx=(0, 8), pady=5)
        labels = ("Mean", "Height", r"Width", "Half window")
        for col, label in enumerate(("", *labels)):
            ctk.CTkLabel(
                entry_frame,
                text=label,
                font=ctk.CTkFont(size=14, weight="bold"),
                text_color="#334155",
            ).grid(row=0, column=col, padx=2, pady=(4, 1), sticky="wens")
        for row, fec in enumerate(self.cfg.fec_ids, start=1):
            ctk.CTkLabel(
                entry_frame,
                text=f"FEC{fec}",
                font=ctk.CTkFont(size=14, weight="bold"),
                text_color="#172033",
            ).grid(row=row, column=0, padx=(8, 4), pady=1, sticky="wens")
            self.entries[fec] = {}
            self.entry_widgets[fec] = {}
            for col, name in enumerate(("mean", "height", "width", "fit_window"), start=1):
                var = tk.StringVar()
                self.entries[fec][name] = var
                entry = ctk.CTkEntry(
                    entry_frame,
                    textvariable=var,
                    width=76,
                    height=24,
                    font=ctk.CTkFont(size=11),
                    fg_color="#ffffff",
                    border_color="#cbd5e1",
                )
                entry.grid(row=row, column=col, padx=4, pady=1, sticky="w")
                self.entry_widgets[fec][name] = entry

        log_frame = ctk.CTkFrame(controls, fg_color="#0f172a", corner_radius=12)
        log_frame.grid(row=0, column=2, sticky="nsew", padx=(0, 12), pady=10)
        ctk.CTkLabel(
            log_frame,
            text="Log",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color="#93c5fd",
        ).pack(anchor="w", padx=12, pady=(8, 0))
        ctk.CTkLabel(
            log_frame,
            textvariable=self.message_var,
            justify="left",
            anchor="nw",
            wraplength=560,
            font=ctk.CTkFont(size=11),
            text_color="#e5e7eb",
        ).pack(fill="both", expand=True, padx=12, pady=(2, 10))

    def add_file_label(self, parent, column: int, title: str, variable: tk.StringVar, color: str) -> None:
        frame = ctk.CTkFrame(parent, fg_color="#ffffff", corner_radius=8)
        frame.grid(row=0, column=column, sticky="ew", padx=6, pady=6)
        ctk.CTkLabel(
            frame,
            text=title,
            anchor="w",
            font=ctk.CTkFont(size=16, weight="bold"),
            text_color=color,
        ).pack(fill="x", padx=10, pady=(2, 2))
        ctk.CTkLabel(
            frame,
            textvariable=variable,
            anchor="w",
            font=ctk.CTkFont(size=12),
            text_color="#334155",
        ).pack(fill="x", padx=10, pady=(0, 0))

    def short_path(self, path: Path | str, keep_parts: int = 4) -> str:
        path = Path(path)
        parts = path.parts
        if len(parts) <= keep_parts + 1:
            return str(path)
        return ".../" + "/".join(parts[-keep_parts:])

    def current_time_id(self) -> str:
        return self.time_ids[self.time_index]

    def current_key(self) -> tuple[str, int]:
        return self.current_time_id(), self.selected_fec

    def current_file_path(self) -> Path:
        return self.cfg.data_root / self.current_time_id() / self.cfg.file_name

    def state_for(self, time_id: str, fec: int) -> ReviewState:
        key = (time_id, fec)
        self.states.setdefault(key, ReviewState())
        return self.states[key]

    def current_state(self) -> ReviewState:
        return self.state_for(*self.current_key())

    def note_cached_time(self, time_id: str) -> None:
        self.cached_time_ids[time_id] = None
        self.cached_time_ids.move_to_end(time_id)
        while len(self.cached_time_ids) > self.max_cached_times:
            old_time_id, _ = self.cached_time_ids.popitem(last=False)
            for key in list(self.values_cache):
                if key[0] == old_time_id:
                    del self.values_cache[key]
            for key in list(self.histogram_cache):
                if key[0] == old_time_id:
                    del self.histogram_cache[key]

    def values_for(self, time_id: str, fec: int) -> np.ndarray:
        key = (time_id, fec)
        if key not in self.values_cache:
            file_path = self.cfg.data_root / time_id / self.cfg.file_name
            for loaded_fec, values in tp.read_fec_channel_values_by_fec(file_path, self.cfg).items():
                self.values_cache[(time_id, loaded_fec)] = values
        self.note_cached_time(time_id)
        return self.values_cache[key]

    def histogram_edges(self) -> np.ndarray:
        return tp.integer_centered_histogram_edges(self.cfg)

    def histogram_cache_path(self, time_id: str, fec: int) -> Path:
        safe_time = time_id.replace("/", "_")
        return self.cfg.outdir / "gui_hist_cache" / f"{safe_time}_fec{fec}_ch{self.cfg.test_pulse_channel}.npz"

    def histogram_cache_metadata(self, time_id: str, fec: int) -> dict[str, float | int]:
        file_path = self.cfg.data_root / time_id / self.cfg.file_name
        try:
            source_mtime_ns = file_path.stat().st_mtime_ns
        except FileNotFoundError:
            source_mtime_ns = -1
        return {
            "source_mtime_ns": source_mtime_ns,
            "fec": fec,
            "channel": self.cfg.test_pulse_channel,
            "histogram_min_adu": self.cfg.histogram_min_adu,
            "histogram_max_adu": self.cfg.histogram_max_adu,
            "histogram_bin_width_adu": self.cfg.histogram_bin_width_adu,
            "histogram_edge_offset_adu": -0.5 * self.cfg.histogram_bin_width_adu,
            "fit_window_adu": self.cfg.fit_window_adu,
            "sigma_min_adu": self.cfg.sigma_min_adu,
            "sigma_max_adu": self.cfg.sigma_max_adu,
        }

    def load_histogram_from_disk(self, time_id: str, fec: int) -> HistogramData | None:
        path = self.histogram_cache_path(time_id, fec)
        if not path.exists():
            return None
        expected = self.histogram_cache_metadata(time_id, fec)
        try:
            with np.load(path) as cached:
                for name, value in expected.items():
                    if name not in cached or cached[name].item() != value:
                        return None
                return HistogramData(
                    counts=cached["counts"],
                    edges=cached["edges"],
                    centers=cached["centers"],
                    peak_center=float(cached["peak_center"].item()),
                    peak_count=float(cached["peak_count"].item()),
                    peak_width=float(cached["peak_width"].item()),
                )
        except Exception:
            return None

    def save_histogram_to_disk(self, time_id: str, fec: int, histogram: HistogramData) -> None:
        path = self.histogram_cache_path(time_id, fec)
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(path, counts=histogram.counts, edges=histogram.edges, centers=histogram.centers,
                                peak_center=histogram.peak_center, peak_count=histogram.peak_count,
                                peak_width=histogram.peak_width, **self.histogram_cache_metadata(time_id, fec))
        except Exception:
            pass

    def histogram_for(self, time_id: str, fec: int) -> HistogramData:
        key = (time_id, fec)
        if key in self.histogram_cache:
            self.note_cached_time(time_id)
            return self.histogram_cache[key]
        cached = self.load_histogram_from_disk(time_id, fec)
        if cached is not None:
            self.histogram_cache[key] = cached
            self.note_cached_time(time_id)
            return cached
        values = self.values_for(time_id, fec)
        edges = self.histogram_edges()
        counts, edges = np.histogram(values, bins=edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        histogram = HistogramData(counts=counts, edges=edges, centers=centers)
        if counts.size > 0 and np.max(counts) > 0:
            peak_index = int(np.argmax(counts))
            histogram.peak_center = float(centers[peak_index])
            histogram.peak_count = float(counts[peak_index])
            near_peak = values[np.abs(values - histogram.peak_center) <= self.cfg.fit_window_adu]
            width = float(np.std(near_peak)) if near_peak.size > 1 else math.nan
            if not np.isfinite(width) or width <= 0.0:
                width = max(1.0, self.cfg.sigma_min_adu)
            histogram.peak_width = float(np.clip(width, self.cfg.sigma_min_adu, self.cfg.sigma_max_adu))
        self.histogram_cache[key] = histogram
        self.save_histogram_to_disk(time_id, fec, histogram)
        return histogram

    def histogram_guess_for(self, time_id: str, fec: int) -> InitialGuess:
        histogram = self.histogram_for(time_id, fec)
        if not np.isfinite(histogram.peak_center):
            return InitialGuess()
        return InitialGuess(histogram.peak_center, histogram.peak_count,
                            histogram.peak_width, self.cfg.fit_window_adu)

    def default_initial_guess_for(self, time_id: str, fec: int) -> InitialGuess:
        state = self.state_for(time_id, fec)
        if state.result is not None and np.isfinite(state.result.peak_adu):
            fit_window = self.cfg.fit_window_adu
            if np.isfinite(state.result.fit_min_adu) and np.isfinite(state.result.fit_max_adu):
                fit_window = 0.5 * (state.result.fit_max_adu - state.result.fit_min_adu)
            return InitialGuess(state.result.peak_adu, state.result.amplitude,
                                state.result.sigma_adu, fit_window)
        return self.histogram_guess_for(time_id, fec)

    def set_guess_boxes_for_current_time(self, force: bool = False) -> None:
        time_id = self.current_time_id()
        for fec in self.cfg.fec_ids:
            key = (time_id, fec)
            if force or key not in self.initial_guesses:
                self.initial_guesses[key] = self.default_initial_guess_for(time_id, fec)
            guess = self.initial_guesses[key]
            for name in ("mean", "height", "width", "fit_window"):
                value = getattr(guess, name)
                self.entries[fec][name].set("" if not np.isfinite(value) else f"{value:.1f}")
        self.update_entry_states()

    def initial_guess_from_entries(self, fec: int) -> InitialGuess:
        values = self.entries[fec]
        return InitialGuess(
            mean=parse_optional_float(values["mean"].get()),
            height=parse_optional_float(values["height"].get()),
            width=parse_optional_float(values["width"].get()),
            fit_window=parse_optional_float(values["fit_window"].get()),
        )

    def sync_current_time_initial_guesses_from_entries(self) -> None:
        time_id = self.current_time_id()
        for fec in self.cfg.fec_ids:
            self.initial_guesses[(time_id, fec)] = self.initial_guess_from_entries(fec)

    def select_fec(self, fec: int) -> None:
        self.selected_fec = fec
        self.update_status()

    def draw(self) -> None:
        time_id = self.current_time_id()
        self.title_var.set(
            f"Test-pulse ch{self.cfg.test_pulse_channel} | {time_id} "
            f"({self.time_index + 1}/{len(self.time_ids)} files)"
        )
        self.path_var.set(self.short_path(self.current_file_path(), keep_parts=5))
        self.input_file_var.set(self.short_path(self.current_file_path(), keep_parts=4))
        self.review_file_var.set(self.short_path(self.cfg.review_csv, keep_parts=4))
        self.summary_file_var.set(self.short_path(self.cfg.summary_csv, keep_parts=4))
        self.root.title(f"Test-pulse ch{self.cfg.test_pulse_channel}: {time_id}")
        for fec in self.cfg.fec_ids:
            state = self.state_for(time_id, fec)
            self.peak_vars[fec].set(state.peak_ok is True)
            histogram = self.histogram_for(time_id, fec)
            if state.result is not None and np.isfinite(state.result.peak_adu):
                zoom_center = state.result.peak_adu
            elif np.isfinite(histogram.peak_center):
                zoom_center = histogram.peak_center
            else:
                zoom_center = 0.0
            self.wide_panels[fec].draw(histogram, (self.cfg.histogram_min_adu, self.cfg.histogram_max_adu),
                                       state.result, selected=False, draw_fit=False)
            self.zoom_panels[fec].draw(histogram, (zoom_center - 30.0, zoom_center + 30.0),
                                       state.result, selected=state.peak_ok is True, draw_fit=True)
        self.update_navigation_buttons()
        self.update_entry_states()
        self.update_status()

    def status_message(self) -> str:
        time_id = self.current_time_id()
        active = [fec for fec in self.cfg.fec_ids if self.state_for(time_id, fec).peak_ok is True]
        fitted = [fec for fec in active if self.state_for(time_id, fec).result is not None]
        accepted = sum(1 for item in self.all_items
                       if self.states.get(item, ReviewState()).peak_ok is True
                       and self.states.get(item, ReviewState()).result is not None
                       and self.states.get(item, ReviewState()).result.alive)
        decided = sum(1 for item in self.all_items if self.states.get(item, ReviewState()).peak_ok is not None)
        state = self.current_state()
        lines = [
            f"selected=FEC{self.selected_fec}; Peak ON={', '.join(f'FEC{x}' for x in active) if active else 'none'}; "
            f"fitted={len(fitted)}/{len(active)}; decided={decided}/{len(self.all_items)}; accepted={accepted}",
        ]
        if state.result is not None:
            lines.append(
                f"fit: peak={state.result.peak_adu:.3g} +/- {state.result.peak_adu_err:.2g}, "
                f"sigma={state.result.sigma_adu:.3g}, chi2/ndf={state.result.chi2_ndf:.3g}, reason={state.result.reason}"
            )
        if state.message:
            lines.append(state.message)
        return "\n".join(lines[:3])

    def update_status(self) -> None:
        self.message_var.set(self.status_message())

    def update_navigation_buttons(self) -> None:
        at_first = self.time_index == 0
        at_last = self.time_index == len(self.time_ids) - 1
        if self.first_button is not None:
            self.first_button.configure(state="disabled" if at_first else "normal")
        if self.prev_button is not None:
            self.prev_button.configure(state="disabled" if at_first else "normal")
        if self.next_button is not None:
            self.next_button.configure(state="disabled" if at_last else "normal")
        if self.last_button is not None:
            self.last_button.configure(state="disabled" if at_last else "normal")

    def update_entry_states(self) -> None:
        time_id = self.current_time_id()
        for fec in self.cfg.fec_ids:
            enabled = self.state_for(time_id, fec).peak_ok is True
            for entry in self.entry_widgets.get(fec, {}).values():
                entry.configure(
                    state="normal" if enabled else "disabled",
                    fg_color="#ffffff" if enabled else "#e5e7eb",
                    border_color="#94a3b8" if enabled else "#cbd5e1",
                    text_color="#111827" if enabled else "#94a3b8",
                )

    def persist(self) -> None:
        save_review_states(self.cfg, self.all_items, self.states)

    def on_toggle_peak(self, fec: int) -> None:
        self.selected_fec = fec
        state = self.current_state()
        state.peak_ok = bool(self.peak_vars[fec].get())
        state.use_fit = None if state.peak_ok else False
        state.message = f"FEC{fec} peak marked {'OK' if state.peak_ok else 'NG'}."
        self.persist()
        self.draw()

    def move(self, step: int) -> None:
        self.sync_current_time_initial_guesses_from_entries()
        self.time_index = max(0, min(len(self.time_ids) - 1, self.time_index + step))
        self.set_guess_boxes_for_current_time(force=False)
        self.draw()

    def jump_to(self, index: int) -> None:
        self.sync_current_time_initial_guesses_from_entries()
        self.time_index = max(0, min(len(self.time_ids) - 1, index))
        self.set_guess_boxes_for_current_time(force=False)
        self.draw()

    def on_first(self) -> None:
        self.jump_to(0)

    def on_prev(self) -> None:
        self.move(-1)

    def on_next(self) -> None:
        self.move(1)

    def on_last(self) -> None:
        self.jump_to(len(self.time_ids) - 1)

    def on_fit(self) -> None:
        time_id = self.current_time_id()
        active = [fec for fec in self.cfg.fec_ids if self.state_for(time_id, fec).peak_ok is True]
        if not active:
            self.current_state().message = "Turn on at least one Peak switch before fitting."
            self.update_status()
            return
        file_path = self.current_file_path()
        messages = []
        for fec in active:
            self.selected_fec = fec
            state = self.state_for(time_id, fec)
            guess = self.initial_guess_from_entries(fec)
            self.initial_guesses[(time_id, fec)] = guess
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
            state.use_fit = bool(state.result.alive)
            if state.result is not None and np.isfinite(state.result.peak_adu):
                self.initial_guesses[(time_id, fec)] = self.default_initial_guess_for(time_id, fec)
            messages.append(f"FEC{fec}: {state.result.reason}")
        self.current_state().message = "Fit done: " + "; ".join(messages)
        self.set_guess_boxes_for_current_time(force=False)
        self.persist()
        self.draw()
        detail_lines = [f"time_id: {time_id}"]
        for fec in active:
            result = self.state_for(time_id, fec).result
            if result is None:
                continue
            detail_lines.append(
                f"FEC{fec}: {result.reason}, peak={result.peak_adu:.3g}, "
                f"sigma={result.sigma_adu:.3g}, chi2/ndf={result.chi2_ndf:.3g}"
            )
        self.current_state().message = "Fit done:\n" + "\n".join(detail_lines)
        self.update_status()

    def on_save_summary(self) -> None:
        if self.save_button is not None:
            self.save_button.configure(fg_color="#f59e0b", hover_color="#d97706")
            self.root.update_idletasks()
        time_id = self.current_time_id()
        missing = [fec for fec in self.cfg.fec_ids
                   if self.state_for(time_id, fec).peak_ok is True
                   and self.state_for(time_id, fec).result is None]
        if missing:
            self.current_state().message = "Save blocked: missing fit for " + ", ".join(f"FEC{fec}" for fec in missing)
            if self.save_button is not None:
                self.save_button.configure(fg_color="#dc2626", hover_color="#b91c1c")
            self.update_status()
            return
        failed = [fec for fec in self.cfg.fec_ids
                  if self.state_for(time_id, fec).peak_ok is True
                  and self.state_for(time_id, fec).result is not None
                  and not self.state_for(time_id, fec).result.alive]
        if failed:
            self.current_state().message = "Save blocked: failed fit for " + ", ".join(f"FEC{fec}" for fec in failed)
            if self.save_button is not None:
                self.save_button.configure(fg_color="#dc2626", hover_color="#b91c1c")
            self.update_status()
            return
        self.persist()
        write_summary_from_review(self.cfg, self.all_items, self.states)
        self.current_state().message = f"Saved summary CSV: {self.cfg.summary_csv}"
        if self.save_button is not None:
            self.save_button.configure(fg_color="#16a34a", hover_color="#15803d")
        self.update_status()

    def show(self) -> None:
        self.root.mainloop()
