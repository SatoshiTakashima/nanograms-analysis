import numpy as np
import uproot
import matplotlib.pyplot as plt
from os.path import join

def build_plane_bin_array(x_range, y_range, num_bins_x, num_bins_y):
    return [
        np.linspace(*x_range, int(num_bins_x) + 1),
        np.linspace(*y_range, int(num_bins_y) + 1),
    ]

def calcConeAxisOrtho(cone_axis):
    """
    cone_axis: (N, 3), 各行が単位ベクトル
    returns
      u1: (N, 3)  cone_axis に垂直な単位ベクトル
      u2: (N, 3)  cone_axis, u1 の両方に垂直な単位ベクトル
    """
    N = cone_axis.shape[0]

    ref = np.tile(np.array([0.0, 0.0, 1.0]), (N, 1))
    mask = np.abs(cone_axis[:, 2]) > 0.9
    ref[mask] = np.array([1.0, 0.0, 0.0])

    u1 = np.cross(cone_axis, ref)
    u1 /= np.linalg.norm(u1, axis=1, keepdims=True)

    u2 = np.cross(cone_axis, u1)
    u2 /= np.linalg.norm(u2, axis=1, keepdims=True)

    return u1, u2

def sectionConeAndPlane_vectorized(vertex, cone_dir, plane_normal, plane_point):
    """
    vertex   : (N, 3)
    cone_dir : (N, P, 3)

    returns
      cone_section        : (N, P, 3)
      positive_flag_array : (N, P)
    """
    numerator = np.sum((plane_point - vertex) * plane_normal, axis=1, keepdims=True)   # (N,1)
    denominator = np.sum(cone_dir * plane_normal[None, None, :], axis=2)                # (N,P)

    eps = 1e-12
    valid = np.abs(denominator) > eps

    t = np.where(valid, numerator / denominator, np.nan)   # (N,P)

    cone_section = vertex[:, None, :] + t[:, :, None] * cone_dir
    #compton coneと平面が重なるかどうか
    positive_flag_array = (t > 0.0) & valid

    return cone_section, positive_flag_array

def _prepare_backprojection_geometry(dataDict):
    pos_hit1 = np.stack([dataDict[f"hit1_pos{xyz}"] for xyz in ["x", "y", "z"]], axis=1)
    pos_hit2 = np.stack([dataDict[f"hit2_pos{xyz}"] for xyz in ["x", "y", "z"]], axis=1)
    cone_vertex = pos_hit1

    dr12 = pos_hit1 - pos_hit2
    cone_axis = dr12 / np.linalg.norm(dr12, axis=1, keepdims=True)
    u1, u2 = calcConeAxisOrtho(cone_axis)

    return cone_vertex, cone_axis, u1, u2

def iterate_backprojection_chunks(
    dataDict,
    plane_normal,
    plane_point,
    spread_arm,
    num_points=1000,
    chunk_size=5000,
    rng=None
):
    if rng is None:
        rng = np.random.default_rng()

    cone_vertex, cone_axis, u1, u2 = _prepare_backprojection_geometry(dataDict)
    numEvents = cone_axis.shape[0]
    phi_array = rng.uniform(0.0, 2.0 * np.pi, num_points)

    cphi = np.cos(phi_array)[None, :, None]
    sphi = np.sin(phi_array)[None, :, None]

    for start in range(0, numEvents, chunk_size):
        end = min(start + chunk_size, numEvents)

        vtx = cone_vertex[start:end]
        axis = cone_axis[start:end]
        uu1 = u1[start:end]
        uu2 = u2[start:end]

        #shape: (end-start, num_points)
        theta = (
            dataDict["theta_k"][start:end, None]
            + rng.normal(loc=0.0, scale=spread_arm, size=(end - start, num_points))
        )

        ctheta = np.cos(theta)[:, :, None]
        stheta = np.sin(theta)[:, :, None]

        perp_dir = cphi * uu1[:, None, :] + sphi * uu2[:, None, :]
        cone_dir = ctheta * axis[:, None, :] + stheta * perp_dir

        cone_section, positive = sectionConeAndPlane_vectorized(vtx, cone_dir, plane_normal, plane_point)

        yield {
            "start": start,
            "end": end,
            "plane_x": cone_section[:, :, 0].astype(np.float32, copy=False),
            "plane_y": cone_section[:, :, 1].astype(np.float32, copy=False),
            "positive_flag": positive,
        }

def calcBackProjection_chunked(dataDict, plane_normal, plane_point,
                               spread_arm, num_points=1000, chunk_size=5000, rng=None):
    numEvents = len(dataDict["hit1_energy"])
    plane_x_array = np.empty((numEvents, num_points), dtype=np.float32)
    plane_y_array = np.empty((numEvents, num_points), dtype=np.float32)
    positive_flag_array = np.empty((numEvents, num_points), dtype=bool)

    for chunk in iterate_backprojection_chunks(
        dataDict,
        plane_normal,
        plane_point,
        spread_arm,
        num_points=num_points,
        chunk_size=chunk_size,
        rng=rng,
    ):
        start = chunk["start"]
        end = chunk["end"]
        plane_x_array[start:end] = chunk["plane_x"]
        plane_y_array[start:end] = chunk["plane_y"]
        positive_flag_array[start:end] = chunk["positive_flag"]

    return plane_x_array, plane_y_array, positive_flag_array


def plotBackProjection(
    dataDict,
    plane_point,
    plane_normal=np.array([0,0,1]),
    xRange=(-50,50),
    yRange=(-50,50),
    numBinsXY=200,
    numBinsX=None,
    numBinsY=None,
    title=None,
    figName=None,
    spread_arm=np.deg2rad(3),
    num_points=1000,
):
    plane_x_array, plane_y_array, positive_flag_array = calcBackProjection_chunked(
        dataDict,
        plane_normal=plane_normal,
        plane_point=plane_point,
        spread_arm=spread_arm,
        num_points=num_points,
    )

    plane_x_array_masked = plane_x_array[positive_flag_array]
    plane_y_array_masked = plane_y_array[positive_flag_array]

    if numBinsX is None:
        numBinsX = numBinsXY
    if numBinsY is None:
        numBinsY = numBinsXY

    planeBinArrayXY = build_plane_bin_array(xRange, yRange, numBinsX, numBinsY)

    H = np.histogram2d(
        plane_x_array_masked.flatten(),
        plane_y_array_masked.flatten(),
        bins=planeBinArrayXY,
    )

    threshold = 0.5 / num_points

    H_masked = np.ma.masked_less(H[0].T, threshold)

    cmap = plt.cm.viridis.copy()
    cmap.set_bad("white")

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))

    pcm = ax.pcolormesh(
        H[1],
        H[2],
        H_masked,
        cmap=cmap,
        vmin=0,
    )

    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label("Counts")

    ax.set_xlim(*xRange)
    ax.set_ylim(*yRange)
    ax.set_xlabel("X (cm)")
    ax.set_ylabel("Y (cm)")
    ax.set_aspect("equal")

    if title is None:
        ax.set_title(f"z={plane_point[2]:.1f} cm")
    else:
        ax.set_title(title)

    if figName is not None:
        fig.savefig(figName)

    return fig, ax
