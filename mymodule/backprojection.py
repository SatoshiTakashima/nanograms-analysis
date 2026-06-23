from __future__ import annotations

from pathlib import Path

import numpy as np
import uproot
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams["font.size"] = 15
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13


DEFAULT_CETREE_BRANCHES = (
    "hit1_energy",
    "hit1_posx",
    "hit1_posy",
    "hit1_posz",
    "hit2_energy",
    "hit2_posx",
    "hit2_posy",
    "hit2_posz",
    "costheta",
)


def read_cetree_backprojection_data(
    file_path,
    tree_name: str = "cetree",
    cut: str | None = "num_hits==2",
    branches: tuple[str, ...] | list[str] | None = None,
    library: str = "np",
) -> dict[str, np.ndarray]:
    """Read a Compton-event tree and prepare arrays used by backprojection.

    The returned dict has ``theta_k`` even when the input tree only stores
    ``costheta``.
    """
    if branches is None:
        branches = DEFAULT_CETREE_BRANCHES

    file_path = Path(file_path)
    with uproot.open(f"{file_path}:{tree_name}") as tree:
        arrays = tree.arrays(list(branches), cut=cut, library=library)

    return _prepare_cetree_backprojection_arrays(arrays)


def _prepare_cetree_backprojection_arrays(arrays) -> dict[str, np.ndarray]:
    if "theta_k" not in arrays and "costheta" in arrays:
        arrays["theta_k"] = np.arccos(np.clip(arrays["costheta"], -1.0, 1.0))
    return dict(arrays)


def iterate_cetree_backprojection_data(
    file_path,
    tree_name: str = "cetree",
    cut: str | None = "num_hits==2",
    branches: tuple[str, ...] | list[str] | None = None,
    step_size="100 MB",
    library: str = "np",
):
    """Yield cetree arrays in chunks prepared for backprojection."""
    if branches is None:
        branches = DEFAULT_CETREE_BRANCHES

    file_path = Path(file_path)
    with uproot.open(f"{file_path}:{tree_name}") as tree:
        for arrays in tree.iterate(
            expressions=list(branches),
            cut=cut,
            step_size=step_size,
            library=library,
        ):
            yield _prepare_cetree_backprojection_arrays(arrays)


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


def build_plane_bin_array(x_range, y_range, num_bins_x, num_bins_y):
    return [
        np.linspace(*x_range, int(num_bins_x) + 1),
        np.linspace(*y_range, int(num_bins_y) + 1),
    ]


def _prepare_projection_axes(plane_normal, plane_yaxis):
    plane_normal = _normalize_vector(plane_normal, "plane_normal")
    plane_yaxis = np.asarray(plane_yaxis, dtype=float)
    plane_yaxis = plane_yaxis - np.dot(plane_yaxis, plane_normal) * plane_normal
    plane_yaxis = _normalize_vector(plane_yaxis, "plane_yaxis")
    plane_xaxis = _normalize_vector(np.cross(plane_yaxis, plane_normal), "plane_xaxis")
    return plane_normal, plane_xaxis, plane_yaxis


def section_cone_and_plane_vectorized(vertex, cone_dir, plane_normal, plane_point):
    numerator = np.sum((plane_point - vertex) * plane_normal, axis=1, keepdims=True)
    denominator = np.sum(cone_dir * plane_normal[None, None, :], axis=2)

    eps = 1.0e-12
    valid = np.abs(denominator) > eps
    t = np.where(valid, numerator / denominator, np.nan)
    cone_section = vertex[:, None, :] + t[:, :, None] * cone_dir
    positive = (t > 0.0) & valid
    return cone_section, positive


def iterate_plane_backprojection_chunks(
    dataDict,
    plane_point,
    plane_normal,
    plane_yaxis,
    spread_arm,
    num_points: int = 1000,
    chunk_size: int = 5000,
    rng=None,
):
    """Yield Compton-cone intersections projected onto a specified plane."""
    if rng is None:
        rng = np.random.default_rng()

    plane_point = np.asarray(plane_point, dtype=float)
    plane_normal, plane_xaxis, plane_yaxis = _prepare_projection_axes(plane_normal, plane_yaxis)
    cone_vertex, cone_axis, u1, u2 = _prepare_backprojection_geometry(dataDict)
    num_events = cone_axis.shape[0]
    phi = rng.uniform(0.0, 2.0 * np.pi, num_points)
    cphi = np.cos(phi)[None, :, None]
    sphi = np.sin(phi)[None, :, None]

    for start in range(0, num_events, chunk_size):
        end = min(start + chunk_size, num_events)
        chunk_events = end - start

        theta = (
            dataDict["theta_k"][start:end, None]
            + rng.normal(loc=0.0, scale=spread_arm, size=(chunk_events, num_points))
        )
        ctheta = np.cos(theta)[:, :, None]
        stheta = np.sin(theta)[:, :, None]
        perp_dir = cphi * u1[start:end, None, :] + sphi * u2[start:end, None, :]
        cone_dir = ctheta * cone_axis[start:end, None, :] + stheta * perp_dir

        cone_section, positive = section_cone_and_plane_vectorized(
            cone_vertex[start:end],
            cone_dir,
            plane_normal,
            plane_point,
        )
        plane_section = cone_section - plane_point[None, None, :]

        yield {
            "start": start,
            "end": end,
            "plane_x": np.sum(plane_section * plane_xaxis[None, None, :], axis=2).astype(np.float32, copy=False),
            "plane_y": np.sum(plane_section * plane_yaxis[None, None, :], axis=2).astype(np.float32, copy=False),
            "positive_flag": positive,
        }


def calc_plane_backprojection(
    dataDict,
    plane_point,
    plane_normal,
    plane_yaxis,
    spread_arm,
    num_points: int = 1000,
    chunk_size: int = 5000,
    rng=None,
):
    num_events = len(dataDict["theta_k"])
    plane_x = np.empty((num_events, num_points), dtype=np.float32)
    plane_y = np.empty((num_events, num_points), dtype=np.float32)
    positive = np.empty((num_events, num_points), dtype=bool)

    for chunk in iterate_plane_backprojection_chunks(
        dataDict,
        plane_point,
        plane_normal,
        plane_yaxis,
        spread_arm,
        num_points=num_points,
        chunk_size=chunk_size,
        rng=rng,
    ):
        start = chunk["start"]
        end = chunk["end"]
        plane_x[start:end] = chunk["plane_x"]
        plane_y[start:end] = chunk["plane_y"]
        positive[start:end] = chunk["positive_flag"]

    return plane_x, plane_y, positive


def calc_plane_backprojection_histogram(
    dataDict,
    plane_point,
    plane_normal,
    plane_yaxis,
    spread_arm,
    num_points: int = 1000,
    xRange=(-50.0, 50.0),
    yRange=(-50.0, 50.0),
    numBinsX: int = 200,
    numBinsY: int = 200,
    chunk_size: int = 5000,
    rng=None,
):
    x_edges, y_edges = build_plane_bin_array(xRange, yRange, numBinsX, numBinsY)
    hist = np.zeros((int(numBinsX), int(numBinsY)), dtype=float)

    for chunk in iterate_plane_backprojection_chunks(
        dataDict,
        plane_point,
        plane_normal,
        plane_yaxis,
        spread_arm,
        num_points=num_points,
        chunk_size=chunk_size,
        rng=rng,
    ):
        positive = chunk["positive_flag"]
        chunk_hist, _, _ = np.histogram2d(
            chunk["plane_x"][positive],
            chunk["plane_y"][positive],
            bins=[x_edges, y_edges],
        )
        hist += chunk_hist

    return hist, x_edges, y_edges


def plot_plane_backprojection(
    dataDict,
    plane_point,
    plane_normal,
    plane_yaxis,
    spread_arm,
    num_points: int = 10000,
    xRange=(-50.0, 50.0),
    yRange=(-50.0, 50.0),
    numBinsX: int = 200,
    numBinsY: int = 200,
    chunk_size: int = 5000,
    rng=None,
    title: str | None = None,
    figName: str | None = None,
    xlabel: str = "Plane X (cm)",
    ylabel: str = "Plane Y (cm)",
    cmap="viridis",
    return_hist: bool = True,
):
    """Plot sampled Compton cones projected onto an arbitrary plane."""
    import matplotlib.pyplot as plt

    hist, x_edges, y_edges = calc_plane_backprojection_histogram(
        dataDict,
        plane_point,
        plane_normal,
        plane_yaxis,
        spread_arm,
        num_points=num_points,
        xRange=xRange,
        yRange=yRange,
        numBinsX=numBinsX,
        numBinsY=numBinsY,
        chunk_size=chunk_size,
        rng=rng,
    )

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    mesh = ax.pcolormesh(x_edges, y_edges, hist.T, cmap=cmap)
    fig.colorbar(mesh, ax=ax, label="Counts")
    ax.set_xlim(*xRange)
    ax.set_ylim(*yRange)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_aspect("equal")
    ax.set_title(title or "Plane Back Projection")
    if figName is not None:
        fig.savefig(figName)
    if return_hist:
        return fig, ax, hist, x_edges, y_edges
    return fig, ax


def calcBackProjection_chunked(
    dataDict,
    plane_normal,
    plane_point,
    plane_yaxis,
    spread_arm,
    num_points,
    chunk_size: int = 5000,
    rng=None,
):
    """Backward-compatible wrapper for the old utilfunc argument order."""
    return calc_plane_backprojection(
        dataDict,
        plane_point,
        plane_normal,
        plane_yaxis,
        spread_arm,
        num_points,
        chunk_size=chunk_size,
        rng=rng,
    )


def iterate_backprojection_chunks(
    dataDict,
    plane_normal,
    plane_point,
    plane_yaxis,
    spread_arm,
    num_points=1000,
    chunk_size=5000,
    rng=None,
):
    """Backward-compatible wrapper for the old utilfunc argument order."""
    return iterate_plane_backprojection_chunks(
        dataDict,
        plane_point,
        plane_normal,
        plane_yaxis,
        spread_arm,
        num_points=num_points,
        chunk_size=chunk_size,
        rng=rng,
    )


calcConeAxisOrtho = _calc_cone_axis_ortho
sectionConeAndPlane_vectorized = section_cone_and_plane_vectorized


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


def _apply_axis_tick_config(ax, axis: str, tick_config):
    if tick_config is None:
        return

    if isinstance(tick_config, dict):
        ticks = tick_config.get("ticks")
        labels = tick_config.get("labels")
        minor = tick_config.get("minor", False)
        tick_params = tick_config.get("tick_params")
        text_kwargs = {
            key: value
            for key, value in tick_config.items()
            if key not in {"ticks", "labels", "minor", "tick_params"}
        }
    else:
        ticks = tick_config
        labels = None
        minor = False
        tick_params = None
        text_kwargs = {}

    if ticks is not None:
        if axis == "x":
            ax.set_xticks(ticks, labels=labels, minor=minor, **text_kwargs)
        elif axis == "y":
            ax.set_yticks(ticks, labels=labels, minor=minor, **text_kwargs)
        else:
            raise ValueError("axis must be 'x' or 'y'")

    if tick_params is not None:
        ax.tick_params(axis=axis, **tick_params)


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
        uu1 = u1[start:end] #(chunk, 3)
        uu2 = u2[start:end] #(chunk, 3)
        chunk_events = end - start

        phi = rng.uniform(0.0, 2.0 * np.pi, size=(chunk_events, num_points))
        theta = (
            dataDict["theta_k"][start:end, None]
            + rng.normal(loc=0.0, scale=spread_arm, size=(chunk_events, num_points))
        )

        cphi   = np.cos(phi)[:, :, None] #(chunk_events, num_points, 1)
        sphi   = np.sin(phi)[:, :, None] #(chunk_events, num_points, 1)
        ctheta = np.cos(theta)[:, :, None] #(chunk_events, num_points, 1)
        stheta = np.sin(theta)[:, :, None] #(chunk_events, num_points, 1)

        perp_dir    = cphi * uu1[:, None, :] + sphi * uu2[:, None, :] #(chunk_events, num_points, 3)
        directions  = ctheta * axis[:, None, :] + stheta * perp_dir  #(chunk_events, num_points, 3)
        directions /= np.linalg.norm(directions, axis=2, keepdims=True) #(chunk_events, num_points, 1)

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
    numBinsX: int = 256,
    numBinsY: int = 256,
    degrees: bool = True,
    chunk_size: int = 5000,
    rng=None,
    event_weights=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fill a local-sky x/y histogram from sampled Compton-cone directions."""
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
    numBinsX: int = 256,
    numBinsY: int = 256,
    degrees: bool = True,
    chunk_size: int = 5000,
    rng=None,
    event_weights=None,
    title: str | None = None,
    figName: str | None = None,
    cmap="viridis",
    xlabel="X",
    ylabel="Y",
    xticks=None,
    yticks=None,
    colorbar: bool = True,
    colorbar_label: str = "Weighted counts",
    colorbar_label_position: str = "top",
    colorbar_size: str = "4%",
    colorbar_pad: float = 0.15,
    colorbar_kwargs: dict | None = None,
):
    """Plot a BackProjectionSky-like local sky image."""
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

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
        numBinsX=numBinsX,
        numBinsY=numBinsY,
        degrees=degrees,
        chunk_size=chunk_size,
        rng=rng,
        event_weights=event_weights,
    )

    fig, ax = plt.subplots(1, 1, figsize=(7, 3.5))
    mesh = ax.pcolormesh(x_edges, y_edges, hist.T, cmap=cmap)
    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=colorbar_size, pad=colorbar_pad)
        cbar = fig.colorbar(mesh, cax=cax, **dict(colorbar_kwargs or {}))
        if colorbar_label_position == "top":
            cbar.ax.set_title(colorbar_label, pad=8)
        else:
            cbar.set_label(colorbar_label)
    ax.set_xlim(*xRange)
    ax.set_ylim(*yRange)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    _apply_axis_tick_config(ax, "x", xticks)
    _apply_axis_tick_config(ax, "y", yticks)
    ax.set_aspect("equal")
    if title != None:
        ax.set_title(title)
    plt.subplots_adjust(left=0.09, right=0.93, bottom=0.15, top=0.91)
    if figName is not None:
        fig.savefig(figName)
    return fig, ax, hist, x_edges, y_edges


calc_lonlat_backprojection_histogram = calc_sky_backprojection_histogram
plot_lonlat_backprojection = plot_sky_backprojection


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


def calc_healpix_backprojection_from_cetree(
    file_path,
    spread_arm,
    tree_name: str = "cetree",
    cut: str | None = "num_hits==2",
    branches: tuple[str, ...] | list[str] | None = None,
    cetree_step_size="100 MB",
    nside: int = 64,
    num_points: int = 1000,
    nest: bool = False,
    chunk_size: int = 5000,
    rng=None,
    event_weights=None,
) -> np.ndarray:
    """Fill a HEALPix array while reading a cetree in chunks."""
    try:
        import healpy as hp
    except ImportError as exc:
        raise ImportError("calc_healpix_backprojection_from_cetree requires healpy") from exc

    healpix_map = np.zeros(hp.nside2npix(nside), dtype=float)
    for dataDict in iterate_cetree_backprojection_data(
        file_path,
        tree_name=tree_name,
        cut=cut,
        branches=branches,
        step_size=cetree_step_size,
    ):
        healpix_map += calc_healpix_backprojection(
            dataDict,
            spread_arm,
            nside=nside,
            num_points=num_points,
            nest=nest,
            chunk_size=chunk_size,
            rng=rng,
            event_weights=event_weights,
        )
    return healpix_map


def plot_healpix_map_mollweide(
    healpix_map,
    nest: bool = False,
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
    tick_label_color: str = "white",
    tick_label_outline_color: str | None = "black",
    tick_label_outline_width: float = 2.0,
    tick_label_fontsize: int = 10,
    axis_label_color: str = "#111827",
    axis_label_fontsize: int = 11,
    graticule_color: str = "white",
    graticule_alpha: float = 0.55,
    graticule_linewidth: float = 0.6,
    colorbar: bool = True,
    colorbar_pad: float = 0.03,
    colorbar_fraction: float = 0.046,
    colorbar_shrink: float = 0.78,
):
    """Display a precomputed HEALPix map with healpy's Mollweide projection."""
    try:
        import healpy as hp
        import matplotlib.pyplot as plt
        import matplotlib.patheffects as path_effects
    except ImportError as exc:
        raise ImportError("plot_healpix_map_mollweide requires healpy and matplotlib") from exc

    hp.mollview(
        healpix_map,
        nest=nest,
        title=title,
        unit=unit,
        coord=coord,
        rot=rot,
        flip=flip,
        cmap=cmap,
        cbar=False,
    )
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
        tick_path_effects = None
        if tick_label_outline_color is not None and tick_label_outline_width > 0.0:
            tick_path_effects = [
                path_effects.withStroke(
                    linewidth=tick_label_outline_width,
                    foreground=tick_label_outline_color,
                )
            ]

        def style_tick_text(text):
            if text is not None and tick_path_effects is not None:
                text.set_path_effects(tick_path_effects)
            return text

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
                style_tick_text(
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
                )

        if latitude_tick_interval_deg is not None and latitude_tick_interval_deg > 0.0:
            latitude_values = np.arange(
                -90.0 + latitude_tick_interval_deg,
                90.0,
                latitude_tick_interval_deg,
            )
            for lat in latitude_values:
                style_tick_text(
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
            color=axis_label_color,
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
            color=axis_label_color,
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
    tick_label_color: str = "white",
    tick_label_outline_color: str | None = "black",
    tick_label_outline_width: float = 2.0,
    tick_label_fontsize: int = 10,
    axis_label_color: str = "#111827",
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
        import matplotlib.patheffects as path_effects
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
        tick_path_effects = None
        if tick_label_outline_color is not None and tick_label_outline_width > 0.0:
            tick_path_effects = [
                path_effects.withStroke(
                    linewidth=tick_label_outline_width,
                    foreground=tick_label_outline_color,
                )
            ]

        def style_tick_text(text):
            if text is not None and tick_path_effects is not None:
                text.set_path_effects(tick_path_effects)
            return text

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
                style_tick_text(
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
                )

        if latitude_tick_interval_deg is not None and latitude_tick_interval_deg > 0.0:
            latitude_values = np.arange(
                -90.0 + latitude_tick_interval_deg,
                90.0,
                latitude_tick_interval_deg,
            )
            for lat in latitude_values:
                style_tick_text(
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
            color=axis_label_color,
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
            color=axis_label_color,
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


def calc_plane_backprojection_histogram_from_cetree(
    file_path,
    plane_point,
    plane_normal,
    plane_yaxis,
    spread_arm,
    tree_name: str = "cetree",
    cut: str | None = "num_hits==2",
    branches: tuple[str, ...] | list[str] | None = None,
    cetree_step_size="100 MB",
    **kwargs,
):
    kwargs = dict(kwargs)
    xRange = kwargs.pop("xRange", (-50.0, 50.0))
    yRange = kwargs.pop("yRange", (-50.0, 50.0))
    numBinsX = kwargs.pop("numBinsX", 200)
    numBinsY = kwargs.pop("numBinsY", 200)
    x_edges, y_edges = build_plane_bin_array(xRange, yRange, numBinsX, numBinsY)
    hist = np.zeros((int(numBinsX), int(numBinsY)), dtype=float)

    for dataDict in iterate_cetree_backprojection_data(
        file_path,
        tree_name=tree_name,
        cut=cut,
        branches=branches,
        step_size=cetree_step_size,
    ):
        chunk_hist, _, _ = calc_plane_backprojection_histogram(
            dataDict,
            plane_point,
            plane_normal,
            plane_yaxis,
            spread_arm,
            xRange=xRange,
            yRange=yRange,
            numBinsX=numBinsX,
            numBinsY=numBinsY,
            **kwargs,
        )
        hist += chunk_hist

    return hist, x_edges, y_edges


def plot_plane_backprojection_from_cetree(
    file_path,
    plane_point,
    plane_normal,
    plane_yaxis,
    spread_arm,
    tree_name: str = "cetree",
    cut: str | None = "num_hits==2",
    branches: tuple[str, ...] | list[str] | None = None,
    cetree_step_size="100 MB",
    num_points: int = 10000,
    xRange=(-50.0, 50.0),
    yRange=(-50.0, 50.0),
    numBinsX: int = 200,
    numBinsY: int = 200,
    chunk_size: int = 5000,
    rng=None,
    title: str | None = None,
    figName: str | None = None,
    xlabel: str = "Plane X (cm)",
    ylabel: str = "Plane Y (cm)",
    cmap="viridis",
    return_hist: bool = True,
):
    """Plot plane backprojection while reading a cetree in chunks."""
    import matplotlib.pyplot as plt

    hist, x_edges, y_edges = calc_plane_backprojection_histogram_from_cetree(
        file_path,
        plane_point,
        plane_normal,
        plane_yaxis,
        spread_arm,
        tree_name=tree_name,
        cut=cut,
        branches=branches,
        cetree_step_size=cetree_step_size,
        num_points=num_points,
        xRange=xRange,
        yRange=yRange,
        numBinsX=numBinsX,
        numBinsY=numBinsY,
        chunk_size=chunk_size,
        rng=rng,
    )

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    mesh = ax.pcolormesh(x_edges, y_edges, hist.T, cmap=cmap)
    fig.colorbar(mesh, ax=ax, label="Counts")
    ax.set_xlim(*xRange)
    ax.set_ylim(*yRange)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_aspect("equal")
    ax.set_title(title or "Plane Back Projection")
    if figName is not None:
        fig.savefig(figName)
    if return_hist:
        return fig, ax, hist, x_edges, y_edges
    return fig, ax


def calc_lonlat_backprojection_histogram_from_cetree(
    file_path,
    image_center_lat,
    image_center_lon,
    image_yaxis_lat,
    image_yaxis_lon,
    spread_arm,
    tree_name: str = "cetree",
    cut: str | None = "num_hits==2",
    branches: tuple[str, ...] | list[str] | None = None,
    cetree_step_size="100 MB",
    **kwargs,
):
    kwargs = dict(kwargs)
    xRange = kwargs.pop("xRange", (-10.0, 10.0))
    yRange = kwargs.pop("yRange", (-10.0, 10.0))
    numBinsX = kwargs.pop("numBinsX", 256)
    numBinsY = kwargs.pop("numBinsY", 256)
    x_edges = np.linspace(*xRange, int(numBinsX) + 1)
    y_edges = np.linspace(*yRange, int(numBinsY) + 1)
    hist = np.zeros((int(numBinsX), int(numBinsY)), dtype=float)

    for dataDict in iterate_cetree_backprojection_data(
        file_path,
        tree_name=tree_name,
        cut=cut,
        branches=branches,
        step_size=cetree_step_size,
    ):
        chunk_hist, _, _ = calc_lonlat_backprojection_histogram(
            dataDict,
            image_center_lat,
            image_center_lon,
            image_yaxis_lat,
            image_yaxis_lon,
            spread_arm,
            xRange=xRange,
            yRange=yRange,
            numBinsX=numBinsX,
            numBinsY=numBinsY,
            **kwargs,
        )
        hist += chunk_hist

    return hist, x_edges, y_edges


def plot_lonlat_backprojection_from_cetree(
    file_path,
    image_center_lat,
    image_center_lon,
    image_yaxis_lat,
    image_yaxis_lon,
    spread_arm,
    tree_name: str = "cetree",
    cut: str | None = "num_hits==2",
    branches: tuple[str, ...] | list[str] | None = None,
    cetree_step_size="100 MB",
    **kwargs,
):
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    kwargs = dict(kwargs)
    title = kwargs.pop("title", None)
    figName = kwargs.pop("figName", None)
    cmap = kwargs.pop("cmap", "viridis")
    xlabel = kwargs.pop("xlabel", "X")
    ylabel = kwargs.pop("ylabel", "Y")
    xticks = kwargs.pop("xticks", None)
    yticks = kwargs.pop("yticks", None)
    colorbar = kwargs.pop("colorbar", True)
    colorbar_label = kwargs.pop("colorbar_label", "Weighted counts")
    colorbar_label_position = kwargs.pop("colorbar_label_position", "top")
    colorbar_size = kwargs.pop("colorbar_size", "4%")
    colorbar_pad = kwargs.pop("colorbar_pad", 0.15)
    colorbar_kwargs = kwargs.pop("colorbar_kwargs", None)

    hist, x_edges, y_edges = calc_lonlat_backprojection_histogram_from_cetree(
        file_path,
        image_center_lat,
        image_center_lon,
        image_yaxis_lat,
        image_yaxis_lon,
        spread_arm,
        tree_name=tree_name,
        cut=cut,
        branches=branches,
        cetree_step_size=cetree_step_size,
        **kwargs,
    )

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    mesh = ax.pcolormesh(x_edges, y_edges, hist.T, cmap=cmap)
    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=colorbar_size, pad=colorbar_pad)
        cbar = fig.colorbar(mesh, cax=cax, **dict(colorbar_kwargs or {}))
        if colorbar_label_position == "top":
            cbar.ax.set_title(colorbar_label, pad=8)
        else:
            cbar.set_label(colorbar_label)
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    _apply_axis_tick_config(ax, "x", xticks)
    _apply_axis_tick_config(ax, "y", yticks)
    ax.set_aspect("equal")
    ax.set_title(title or "Sky Back Projection")
    if figName is not None:
        fig.savefig(figName)
    return fig, ax, hist, x_edges, y_edges


def plot_healpix_backprojection_mollweide_from_cetree(
    file_path,
    spread_arm,
    tree_name: str = "cetree",
    cut: str | None = "num_hits==2",
    branches: tuple[str, ...] | list[str] | None = None,
    cetree_step_size="100 MB",
    **kwargs,
):
    kwargs = dict(kwargs)
    nside = kwargs.pop("nside", 64)
    num_points = kwargs.pop("num_points", 1000)
    nest = kwargs.pop("nest", False)
    chunk_size = kwargs.pop("chunk_size", 5000)
    rng = kwargs.pop("rng", None)
    event_weights = kwargs.pop("event_weights", None)

    healpix_map = calc_healpix_backprojection_from_cetree(
        file_path,
        spread_arm,
        tree_name=tree_name,
        cut=cut,
        branches=branches,
        cetree_step_size=cetree_step_size,
        nside=nside,
        num_points=num_points,
        nest=nest,
        chunk_size=chunk_size,
        rng=rng,
        event_weights=event_weights,
    )
    return plot_healpix_map_mollweide(healpix_map, nest=nest, **kwargs)
