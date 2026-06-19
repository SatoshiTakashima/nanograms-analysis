"""Sky-coordinate back projection utilities for Compton cones.

This module mirrors the geometry used by ComptonSoft's BackProjectionSky:
sample directions on each Compton cone, project them either to a local sky
image whose axes are specified by latitude/longitude directions, or fill them
directly into a HEALPix map.
"""

from __future__ import annotations

import numpy as np


def _normalize_vector(vector, name):
    vector = np.asarray(vector, dtype=float)
    norm = np.linalg.norm(vector)
    if norm == 0.0:
        raise ValueError(f"{name} must be a non-zero vector")
    return vector / norm


def _calc_cone_axis_ortho(cone_axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ref = np.tile(np.array([0.0, 0.0, 1.0]), (cone_axis.shape[0], 1))
    mask = np.abs(cone_axis[:, 2]) > 0.9
    ref[mask] = np.array([1.0, 0.0, 0.0])

    u1 = np.cross(cone_axis, ref)
    u1 /= np.linalg.norm(u1, axis=1, keepdims=True)
    u2 = np.cross(cone_axis, u1)
    u2 /= np.linalg.norm(u2, axis=1, keepdims=True)
    return u1, u2


def _prepare_backprojection_geometry(dataDict):
    pos_hit1 = np.stack([dataDict[f"hit1_pos{xyz}"] for xyz in ["x", "y", "z"]], axis=1)
    pos_hit2 = np.stack([dataDict[f"hit2_pos{xyz}"] for xyz in ["x", "y", "z"]], axis=1)
    cone_vertex = pos_hit1

    dr12 = pos_hit1 - pos_hit2
    cone_axis = dr12 / np.linalg.norm(dr12, axis=1, keepdims=True)
    u1, u2 = _calc_cone_axis_ortho(cone_axis)
    return cone_vertex, cone_axis, u1, u2


def _as_radians(value, degrees: bool) -> np.ndarray:
    value = np.asarray(value, dtype=float)
    return np.deg2rad(value) if degrees else value


def _to_output_angle(value: np.ndarray, degrees: bool) -> np.ndarray:
    return np.rad2deg(value) if degrees else value


def unit_vector_from_latlon(lat, lon, degrees: bool = True) -> np.ndarray:
    """Return unit vector(s) for astronomical latitude and longitude."""
    lat = _as_radians(lat, degrees)
    lon = _as_radians(lon, degrees)
    clat = np.cos(lat)
    return np.stack([clat * np.cos(lon), clat * np.sin(lon), np.sin(lat)], axis=-1)


def unit_vector_from_theta_phi(theta, phi, degrees: bool = True) -> np.ndarray:
    """Return unit vector(s) for ComptonSoft theta/phi angles.

    theta is colatitude measured from +z, phi is azimuth around +z.
    """
    theta = _as_radians(theta, degrees)
    phi = _as_radians(phi, degrees)
    return np.stack([
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta),
    ], axis=-1)


def sky_axes_from_latlon(
    image_center_lat,
    image_center_lon,
    image_yaxis_lat,
    image_yaxis_lon,
    degrees: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build local sky x/y/z axes from center and proposed y-axis directions.

    The returned axes follow BackProjectionSky.cc:
    z is the image center direction, y is the proposed y direction projected
    onto the plane perpendicular to z, and x = y cross z.
    """
    zaxis = _normalize_vector(
        unit_vector_from_latlon(image_center_lat, image_center_lon, degrees),
        "image_center",
    )
    yaxis_proposed = unit_vector_from_latlon(image_yaxis_lat, image_yaxis_lon, degrees)
    yaxis = yaxis_proposed - np.dot(yaxis_proposed, zaxis) * zaxis
    yaxis = _normalize_vector(yaxis, "image_yaxis")
    xaxis = _normalize_vector(np.cross(yaxis, zaxis), "image_xaxis")
    return xaxis, yaxis, zaxis


def sky_axes_from_theta_phi(
    image_center_theta,
    image_center_phi,
    image_yaxis_theta,
    image_yaxis_phi,
    degrees: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build local sky axes using the same theta/phi parameters as C++."""
    center = unit_vector_from_theta_phi(image_center_theta, image_center_phi, degrees)
    yaxis = unit_vector_from_theta_phi(image_yaxis_theta, image_yaxis_phi, degrees)
    center_lat = np.arcsin(center[2])
    center_lon = np.arctan2(center[1], center[0])
    yaxis_lat = np.arcsin(yaxis[2])
    yaxis_lon = np.arctan2(yaxis[1], yaxis[0])
    return sky_axes_from_latlon(center_lat, center_lon, yaxis_lat, yaxis_lon, degrees=False)


def local_sky_xy_from_directions(
    directions: np.ndarray,
    xaxis: np.ndarray,
    yaxis: np.ndarray,
    zaxis: np.ndarray,
    degrees: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Project global unit direction vectors to BackProjectionSky local x/y."""
    ux = np.sum(directions * xaxis, axis=-1)
    uy = np.sum(directions * yaxis, axis=-1)
    uz = np.sum(directions * zaxis, axis=-1)
    local_x = np.arctan2(ux, uz)
    local_y = np.arcsin(np.clip(uy, -1.0, 1.0))
    return _to_output_angle(local_x, degrees), _to_output_angle(local_y, degrees)


def _event_weights(dataDict, event_weights, num_events: int) -> np.ndarray:
    if event_weights is None:
        for key in ("reconstruction_fraction", "fraction", "weight", "weights"):
            if key in dataDict:
                event_weights = dataDict[key]
                break
    if event_weights is None:
        return np.ones(num_events, dtype=float)
    event_weights = np.asarray(event_weights, dtype=float)
    if event_weights.shape[0] != num_events:
        raise ValueError("event_weights must have the same length as the events")
    return event_weights


def iterate_sky_backprojection_chunks(
    dataDict,
    spread_arm,
    num_points: int = 1000,
    chunk_size: int = 5000,
    rng=None,
    event_weights=None,
):
    """Yield sampled Compton-cone directions in chunks.

    Each yielded direction has unit length. The per-sample weight is
    event_weight / num_points, matching BackProjectionSky.cc's normalization.
    """
    if rng is None:
        rng = np.random.default_rng()

    cone_vertex, cone_axis, u1, u2 = _prepare_backprojection_geometry(dataDict)
    del cone_vertex
    num_events = cone_axis.shape[0]
    weights = _event_weights(dataDict, event_weights, num_events)

    for start in range(0, num_events, chunk_size):
        end = min(start + chunk_size, num_events)

        axis = cone_axis[start:end]
        uu1 = u1[start:end]
        uu2 = u2[start:end]
        chunk_events = end - start

        phi = rng.uniform(0.0, 2.0 * np.pi, size=(chunk_events, num_points))
        theta = (
            dataDict["theta_k"][start:end, None]
            + rng.normal(loc=0.0, scale=spread_arm, size=(chunk_events, num_points))
        )

        cphi = np.cos(phi)[:, :, None]
        sphi = np.sin(phi)[:, :, None]
        ctheta = np.cos(theta)[:, :, None]
        stheta = np.sin(theta)[:, :, None]

        perp_dir = cphi * uu1[:, None, :] + sphi * uu2[:, None, :]
        directions = ctheta * axis[:, None, :] + stheta * perp_dir
        directions /= np.linalg.norm(directions, axis=2, keepdims=True)

        sample_weights = np.repeat(weights[start:end] / float(num_points), num_points)
        yield {
            "start": start,
            "end": end,
            "directions": directions.astype(np.float32, copy=False),
            "weights": sample_weights.astype(np.float32, copy=False),
        }


def calc_sky_backprojection_histogram(
    dataDict,
    image_center_lat,
    image_center_lon,
    image_yaxis_lat,
    image_yaxis_lon,
    spread_arm,
    num_points: int = 1000,
    xRange=(-10.0, 10.0),
    yRange=(-10.0, 10.0),
    numBinsXY: int = 256,
    numBinsX: int | None = None,
    numBinsY: int | None = None,
    degrees: bool = True,
    chunk_size: int = 5000,
    rng=None,
    event_weights=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fill a local-sky x/y histogram from sampled Compton-cone directions."""
    if numBinsX is None:
        numBinsX = numBinsXY
    if numBinsY is None:
        numBinsY = numBinsXY

    x_edges = np.linspace(*xRange, int(numBinsX) + 1)
    y_edges = np.linspace(*yRange, int(numBinsY) + 1)
    hist = np.zeros((int(numBinsX), int(numBinsY)), dtype=float)
    axes = sky_axes_from_latlon(
        image_center_lat,
        image_center_lon,
        image_yaxis_lat,
        image_yaxis_lon,
        degrees=degrees,
    )

    for chunk in iterate_sky_backprojection_chunks(
        dataDict,
        spread_arm,
        num_points=num_points,
        chunk_size=chunk_size,
        rng=rng,
        event_weights=event_weights,
    ):
        x, y = local_sky_xy_from_directions(chunk["directions"], *axes, degrees=degrees)
        weights = chunk["weights"]
        chunk_hist, _, _ = np.histogram2d(
            x.ravel(),
            y.ravel(),
            bins=[x_edges, y_edges],
            weights=weights,
        )
        hist += chunk_hist

    return hist, x_edges, y_edges


def plot_sky_backprojection(
    dataDict,
    image_center_lat,
    image_center_lon,
    image_yaxis_lat,
    image_yaxis_lon,
    spread_arm,
    num_points: int = 1000,
    xRange=(-10.0, 10.0),
    yRange=(-10.0, 10.0),
    numBinsXY: int = 256,
    numBinsX: int | None = None,
    numBinsY: int | None = None,
    degrees: bool = True,
    chunk_size: int = 5000,
    rng=None,
    event_weights=None,
    title: str | None = None,
    figName: str | None = None,
    cmap="viridis",
):
    """Plot a BackProjectionSky-like local sky image."""
    import matplotlib.pyplot as plt

    hist, x_edges, y_edges = calc_sky_backprojection_histogram(
        dataDict,
        image_center_lat,
        image_center_lon,
        image_yaxis_lat,
        image_yaxis_lon,
        spread_arm,
        num_points=num_points,
        xRange=xRange,
        yRange=yRange,
        numBinsXY=numBinsXY,
        numBinsX=numBinsX,
        numBinsY=numBinsY,
        degrees=degrees,
        chunk_size=chunk_size,
        rng=rng,
        event_weights=event_weights,
    )

    unit_name = "deg" if degrees else "rad"
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    mesh = ax.pcolormesh(x_edges, y_edges, hist.T, cmap=cmap)
    fig.colorbar(mesh, ax=ax, label="Weighted counts")
    ax.set_xlim(*xRange)
    ax.set_ylim(*yRange)
    ax.set_xlabel(f"Sky X ({unit_name})")
    ax.set_ylabel(f"Sky Y ({unit_name})")
    ax.set_aspect("equal")
    ax.set_title(title or "Sky Back Projection")
    if figName is not None:
        fig.savefig(figName)
    return fig, ax, hist, x_edges, y_edges


def calc_healpix_backprojection(
    dataDict,
    spread_arm,
    nside: int = 64,
    num_points: int = 1000,
    nest: bool = False,
    chunk_size: int = 5000,
    rng=None,
    event_weights=None,
) -> np.ndarray:
    """Fill a HEALPix array with sampled Compton-cone directions."""
    try:
        import healpy as hp
    except ImportError as exc:
        raise ImportError("calc_healpix_backprojection requires healpy") from exc

    healpix_map = np.zeros(hp.nside2npix(nside), dtype=float)
    for chunk in iterate_sky_backprojection_chunks(
        dataDict,
        spread_arm,
        num_points=num_points,
        chunk_size=chunk_size,
        rng=rng,
        event_weights=event_weights,
    ):
        directions = chunk["directions"].reshape(-1, 3)
        pix = hp.vec2pix(
            nside,
            directions[:, 0],
            directions[:, 1],
            directions[:, 2],
            nest=nest,
        )
        np.add.at(healpix_map, pix, chunk["weights"])
    return healpix_map


def plot_healpix_backprojection_mollweide(
    dataDict,
    spread_arm,
    nside: int = 64,
    num_points: int = 1000,
    nest: bool = False,
    chunk_size: int = 5000,
    rng=None,
    event_weights=None,
    title: str = "HEALPix Back Projection",
    unit: str = "Weighted counts",
    coord=None,
    rot=None,
    flip: str = "geo",
    cmap="viridis",
    figName: str | None = None,
    return_map: bool = True,
    graticule_interval_deg: float | None = None,
    longitude_tick_interval_deg: float | None = 60.0,
    latitude_tick_interval_deg: float | None = 30.0,
    show_tick_labels: bool = True,
    longitude_axis_label: str = "y",
    latitude_axis_label: str = "z",
    tick_label_color: str = "#111827",
    tick_label_fontsize: int = 10,
    axis_label_fontsize: int = 11,
    graticule_color: str = "white",
    graticule_alpha: float = 0.55,
    graticule_linewidth: float = 0.6,
    colorbar: bool = True,
    colorbar_pad: float = 0.03,
    colorbar_fraction: float = 0.046,
    colorbar_shrink: float = 0.78,
):
    """Fill a HEALPix map and display it with healpy's Mollweide projection."""
    try:
        import healpy as hp
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("plot_healpix_backprojection_mollweide requires healpy and matplotlib") from exc

    healpix_map = calc_healpix_backprojection(
        dataDict,
        spread_arm,
        nside=nside,
        num_points=num_points,
        nest=nest,
        chunk_size=chunk_size,
        rng=rng,
        event_weights=event_weights,
    )
    hp.mollview(healpix_map, nest=nest, title=title, unit=unit,
                coord=coord, rot=rot, flip=flip, cmap=cmap, cbar=False)
    fig = plt.gcf()
    ax = plt.gca()

    if graticule_interval_deg is not None and graticule_interval_deg > 0.0:
        longitude_tick_interval_deg = graticule_interval_deg
        latitude_tick_interval_deg = graticule_interval_deg

    if (
        longitude_tick_interval_deg is not None
        and longitude_tick_interval_deg > 0.0
        and latitude_tick_interval_deg is not None
        and latitude_tick_interval_deg > 0.0
    ):
        hp.graticule(
            dpar=latitude_tick_interval_deg,
            dmer=longitude_tick_interval_deg,
            coord=coord,
            color=graticule_color,
            alpha=graticule_alpha,
            linewidth=graticule_linewidth,
        )

    if show_tick_labels:
        def degree_label(value: float) -> str:
            value = float(value)
            if abs(value) < 1.0e-10:
                value = 0.0
            return f"{value:g}\N{DEGREE SIGN}"

        if longitude_tick_interval_deg is not None and longitude_tick_interval_deg > 0.0:
            longitude_values = np.arange(
                -180.0 + longitude_tick_interval_deg,
                180.0,
                longitude_tick_interval_deg,
            )
            for lon in longitude_values:
                label_value = lon % 360.0
                hp.projtext(
                    lon,
                    0.0,
                    degree_label(label_value),
                    lonlat=True,
                    ha="center",
                    va="bottom",
                    color=tick_label_color,
                    fontsize=tick_label_fontsize,
                )

        if latitude_tick_interval_deg is not None and latitude_tick_interval_deg > 0.0:
            latitude_values = np.arange(
                -90.0 + latitude_tick_interval_deg,
                90.0,
                latitude_tick_interval_deg,
            )
            for lat in latitude_values:
                hp.projtext(
                    -178.0,
                    lat,
                    degree_label(lat),
                    lonlat=True,
                    ha="right",
                    va="center",
                    color=tick_label_color,
                    fontsize=tick_label_fontsize,
                )

    if longitude_axis_label:
        ax.text(
            0.5,
            -0.075,
            longitude_axis_label,
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=axis_label_fontsize,
            color=tick_label_color,
        )
    if latitude_axis_label:
        fig.text(
            0.03,
            0.5,
            latitude_axis_label,
            ha="center",
            va="center",
            rotation=90,
            fontsize=axis_label_fontsize,
            color=tick_label_color,
        )

    if colorbar and ax.get_images():
        cbar = fig.colorbar(
            ax.get_images()[0],
            ax=ax,
            orientation="vertical",
            pad=colorbar_pad,
            fraction=colorbar_fraction,
            shrink=colorbar_shrink,
        )
        cbar.set_label(unit)

    if figName is not None:
        fig.savefig(figName)
    if return_map:
        return fig, healpix_map
    return fig
