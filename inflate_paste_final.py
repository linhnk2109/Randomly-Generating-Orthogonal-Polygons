"""
Implementation of Inflate-Paste algorithm
from paper: "Quadratic-Time Linear-Space Algorithms for Generating
Orthogonal Polygons with a Given Number of Vertices"
by Ana Paula Tomás & António Leslie Bajuelos (EWCG 2004)

Key idea:
  At each step, pick a convex vertex vi of the current polygon P.
  Find a cell C in the Free Staircase Neighbourhood (FSN) of vi.
  Apply Inflate(P, p, q)  — shifts coords to open up space for C.
  Apply Paste            — glues a rectangle bump to P̃ at ṽi,
                           introducing +2 vertices per iteration.

Paste never fails when the selected cell C belongs to FSN(vi); this implementation explicitly constructs the grid FSN by marking exterior cells in the extended grid and extracting the maximal monotone free staircase, while keeping a final safety validity check.

This version also includes result validation, runtime/memory statistics,
CSV export, and polygon drawing utilities.
"""

import random
import time
import tracemalloc
import bisect
import csv
import math
import os
from typing import List, Tuple, Set, Optional, Dict, Any

import numpy as np

Point = Tuple[int, int]


# ──────────────────────────────────────────────────────────────────
# Minimal SortedDict replacement using bisect (no external libs)
# ──────────────────────────────────────────────────────────────────

class _SortedDict:
    """Lightweight sorted-key dict backed by bisect for O(log n) ops."""
    def __init__(self):
        self._keys: List = []
        self._vals: dict = {}

    def __setitem__(self, k, v):
        if k not in self._vals:
            bisect.insort(self._keys, k)
        self._vals[k] = v

    def __delitem__(self, k):
        if k in self._vals:
            idx = bisect.bisect_left(self._keys, k)
            if idx < len(self._keys) and self._keys[idx] == k:
                self._keys.pop(idx)
            del self._vals[k]

    def __contains__(self, k):
        return k in self._vals

    def __getitem__(self, k):
        return self._vals[k]

    def irange(self, lo, hi, inclusive=(True, True)):
        """Yield keys in [lo, hi] (both inclusive by default)."""
        if inclusive[0]:
            start = bisect.bisect_left(self._keys, lo)
        else:
            start = bisect.bisect_right(self._keys, lo)
        if inclusive[1]:
            end = bisect.bisect_right(self._keys, hi)
        else:
            end = bisect.bisect_left(self._keys, hi)
        return iter(self._keys[start:end])


# ──────────────────────────────────────────────────────────────────
# Shared helpers (mirrored from inflate_cut_new.py)
# ──────────────────────────────────────────────────────────────────

def remove_duplicates(polygon: List[Point]) -> List[Point]:
    if not polygon:
        return []
    result = [polygon[0]]
    for p in polygon[1:]:
        if p != result[-1]:
            result.append(p)
    if len(result) > 1 and result[-1] == result[0]:
        result.pop()
    return result


def find_convex_vertices(polygon: List[Point]) -> List[int]:
    """Return convex vertex indices, using the polygon orientation."""
    pts = np.asarray(polygon, dtype=np.int64)
    if len(pts) < 4:
        return []

    prev_pts = np.roll(pts, 1, axis=0)
    next_pts = np.roll(pts, -1, axis=0)
    v1 = pts - prev_pts
    v2 = next_pts - pts
    cross = v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]
    area2 = np.sum(pts[:, 0] * next_pts[:, 1] - next_pts[:, 0] * pts[:, 1])
    if area2 == 0:
        return []
    return np.flatnonzero(cross * area2 > 0).astype(int).tolist()


def is_simple_polygon(polygon: List[Point]) -> bool:
    """Sweep-line simplicity check – O(n log n)."""
    n = len(polygon)
    sp = sorted(polygon)
    for i in range(1, len(sp)):
        if sp[i] == sp[i - 1]:
            return False

    h_segs, v_segs = [], []
    adj = set()

    for i in range(n):
        p1, p2 = polygon[i], polygon[(i + 1) % n]
        if p1[0] > p2[0] or (p1[0] == p2[0] and p1[1] > p2[1]):
            p1, p2 = p2, p1
        if p1[0] == p2[0]:
            v_segs.append((p1[0], min(p1[1], p2[1]), max(p1[1], p2[1]), i))
        else:
            h_segs.append((min(p1[0], p2[0]), max(p1[0], p2[0]), p1[1], i))
        adj |= {(i, (i+1)%n), ((i+1)%n, i),
                (i, (i-1+n)%n), ((i-1+n)%n, i)}

    v_segs.sort()
    events = []
    for x1, x2, y, i in h_segs:
        events.append((x1, 0, y, i))
        events.append((x2, 1, y, i))
    for x, y1, y2, i in v_segs:
        events.append((x, 2, y1, y2, i))
    events.sort(key=lambda e: (e[0], e[1]))

    active = _SortedDict()
    for ev in events:
        if ev[1] == 0:
            active[ev[2]] = ev[3]
        elif ev[1] == 1:
            if ev[2] in active:
                del active[ev[2]]
        else:
            _, _, y1, y2, i = ev
            try:
                for y in active.irange(y1, y2):
                    j = active[y]
                    if (i, j) not in adj and (j, i) not in adj:
                        return False
            except Exception:
                pass
    return True


def is_grid_ogon(polygon: List[Point]) -> bool:
    """Check the grid n-ogon general-position condition used in the paper.

    A valid grid n-ogon has only axis-parallel nonzero edges, n/2 horizontal
    edges and n/2 vertical edges, and no two horizontal/vertical edges lie on
    the same grid line.
    """
    polygon = remove_duplicates(polygon)
    n = len(polygon)
    if n < 4 or n % 2 != 0:
        return False

    h_lines: Dict[int, int] = {}
    v_lines: Dict[int, int] = {}
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        if y1 == y2 and x1 != x2:
            h_lines[y1] = h_lines.get(y1, 0) + 1
        elif x1 == x2 and y1 != y2:
            v_lines[x1] = v_lines.get(x1, 0) + 1
        else:
            return False

    return (
        len(h_lines) == n // 2 and
        len(v_lines) == n // 2 and
        all(c == 1 for c in h_lines.values()) and
        all(c == 1 for c in v_lines.values())
    )


def polygon_area(polygon: List[Point]) -> float:
    polygon = remove_duplicates(polygon)
    n = len(polygon)
    if n < 3:
        return 0.0
    area2 = 0
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        area2 += x1 * y2 - x2 * y1
    return abs(area2) / 2.0


def polygon_perimeter(polygon: List[Point]) -> int:
    polygon = remove_duplicates(polygon)
    n = len(polygon)
    if n < 2:
        return 0
    perimeter = 0
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        perimeter += abs(x2 - x1) + abs(y2 - y1)
    return perimeter


def count_reflex_vertices(polygon: List[Point]) -> int:
    polygon = remove_duplicates(polygon)
    pts = np.asarray(polygon, dtype=np.int64)
    n = len(pts)
    if n < 4:
        return 0
    prev_pts = np.roll(pts, 1, axis=0)
    next_pts = np.roll(pts, -1, axis=0)
    v1 = pts - prev_pts
    v2 = next_pts - pts
    cross = v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]
    area2 = np.sum(pts[:, 0] * next_pts[:, 1] - next_pts[:, 0] * pts[:, 1])
    if area2 == 0:
        return 0
    return int(np.count_nonzero(cross * area2 < 0))


def convex_hull(points: List[Point]) -> List[Point]:
    """Monotone-chain convex hull."""
    pts = sorted(set(points))
    if len(pts) <= 1:
        return pts

    def cross(o: Point, a: Point, b: Point) -> int:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: List[Point] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper: List[Point] = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    return lower[:-1] + upper[:-1]


def draw_polygon(
    polygon: List[Point],
    save_path: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    label_vertices: bool = True,
    invert_y: bool = True,
) -> Optional[str]:
    """Draw and optionally save the generated orthogonal polygon."""
    polygon = remove_duplicates(polygon)
    if not polygon:
        raise ValueError("Cannot draw an empty polygon")

    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon

    closed = polygon + [polygon[0]]
    xs = [p[0] for p in closed]
    ys = [p[1] for p in closed]

    fig, ax = plt.subplots(figsize=(8, 8))
    patch = MplPolygon(polygon, closed=True, fill=True, alpha=0.25, edgecolor="black", linewidth=2)
    ax.add_patch(patch)
    ax.plot(xs, ys, marker="o", linewidth=1.5)

    if label_vertices:
        for idx, (x, y) in enumerate(polygon):
            ax.text(x, y, f" {idx}", fontsize=8, ha="left", va="bottom")

    min_x = min(x for x, _ in polygon)
    max_x = max(x for x, _ in polygon)
    min_y = min(y for _, y in polygon)
    max_y = max(y for _, y in polygon)
    pad = 1
    ax.set_xlim(min_x - pad, max_x + pad)
    ax.set_ylim(min_y - pad, max_y + pad)
    if invert_y:
        ax.invert_yaxis()

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title or f"Inflate-Paste polygon, n={len(polygon)}")
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)

    if save_path:
        directory = os.path.dirname(save_path)
        if directory:
            os.makedirs(directory, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return save_path


def inflate(polygon: List[Point], p: int, q: int) -> List[Point]:
    """Inflate transformation – O(n)."""
    pts = np.asarray(polygon, dtype=np.int64)
    inflated = pts.copy()
    inflated[:, 0] += inflated[:, 0] > p
    inflated[:, 1] += inflated[:, 1] > q
    return [tuple(map(int, point)) for point in inflated]


def replace_vertex(polygon: List[Point], idx: int,
                   chain: List[Point]) -> List[Point]:
    """Replace vertex at idx with chain – O(n)."""
    if not chain:
        return polygon
    result = list(polygon[:idx]) + list(chain) + list(polygon[idx + 1:])
    return remove_duplicates(result)


def check_chain_validity(chain: List[Point], point_set: Set[Point]) -> bool:
    if any(p in point_set for p in chain):
        return False
    return len(set(chain)) == len(chain)



def _edge_arrays_np(polygon: List[Point]):
    """Return vectorized horizontal and vertical edge arrays."""
    pts = np.asarray(polygon, dtype=np.int64)
    nxt = np.roll(pts, -1, axis=0)

    horizontal = pts[:, 1] == nxt[:, 1]
    vertical = pts[:, 0] == nxt[:, 0]

    hy = pts[horizontal, 1]
    hxlo = np.minimum(pts[horizontal, 0], nxt[horizontal, 0])
    hxhi = np.maximum(pts[horizontal, 0], nxt[horizontal, 0])

    vx = pts[vertical, 0]
    vylo = np.minimum(pts[vertical, 1], nxt[vertical, 1])
    vyhi = np.maximum(pts[vertical, 1], nxt[vertical, 1])

    return hy, hxlo, hxhi, vx, vylo, vyhi


def _rect_exterior_ok_np(
    xi: int,
    yi: int,
    p: int,
    q: int,
    hy: np.ndarray,
    hxlo: np.ndarray,
    hxhi: np.ndarray,
    vx: np.ndarray,
    vylo: np.ndarray,
    vyhi: np.ndarray,
) -> bool:
    """Check whether the open rectangle defined by vi and c=(p+1,q+1)
    crosses the current polygon boundary.

    This is the practical FSN test: a cell C is accepted only when the rectangle
    that will be pasted from vi to the center of the inflated cell lies in the
    exterior of P, except for allowed boundary contact.
    """
    cx, cy = p + 1, q + 1
    xlo, xhi = min(xi, cx), max(xi, cx)
    ylo, yhi = min(yi, cy), max(yi, cy)

    if xlo == xhi or ylo == yhi:
        return False

    if len(hy):
        h_bad = (hy > ylo) & (hy < yhi) & (hxlo < xhi) & (hxhi > xlo)
        if h_bad.any():
            return False

    if len(vx):
        v_bad = (vx > xlo) & (vx < xhi) & (vylo < yhi) & (vyhi > ylo)
        if v_bad.any():
            return False

    return True


def _paste_rectangle_is_free(
    inflated: List[Point],
    vi_idx: int,
    center: Point,
) -> bool:
    """Final safety check for the Paste rectangle in inflated coordinates.

    If C was selected from FSN(vi), Paste should not fail. This check is kept
    to protect the optimized implementation from numerical/indexing mistakes.
    """
    pts = np.asarray(inflated, dtype=np.int64)
    vi = pts[vi_idx]
    x1, x2 = sorted((int(vi[0]), int(center[0])))
    y1, y2 = sorted((int(vi[1]), int(center[1])))
    if x1 == x2 or y1 == y2:
        return False

    strictly_inside = (
        (x1 < pts[:, 0]) & (pts[:, 0] < x2)
        & (y1 < pts[:, 1]) & (pts[:, 1] < y2)
    )
    if np.any(strictly_inside):
        return False

    starts = pts
    ends = np.roll(pts, -1, axis=0)
    edge_indices = np.arange(len(pts))
    incident = (edge_indices == vi_idx) | (edge_indices == (vi_idx - 1) % len(pts))

    vertical = starts[:, 0] == ends[:, 0]
    vx = starts[:, 0]
    vy_low = np.minimum(starts[:, 1], ends[:, 1])
    vy_high = np.maximum(starts[:, 1], ends[:, 1])
    vertical_cross = (
        vertical & ~incident
        & (x1 < vx) & (vx < x2)
        & (np.maximum(vy_low, y1) < np.minimum(vy_high, y2))
    )
    if np.any(vertical_cross):
        return False

    horizontal = starts[:, 1] == ends[:, 1]
    hy = starts[:, 1]
    hx_low = np.minimum(starts[:, 0], ends[:, 0])
    hx_high = np.maximum(starts[:, 0], ends[:, 0])
    horizontal_cross = (
        horizontal & ~incident
        & (y1 < hy) & (hy < y2)
        & (np.maximum(hx_low, x1) < np.minimum(hx_high, x2))
    )
    return not bool(np.any(horizontal_cross))




def _outside_cell_mask_np(polygon: List[Point], max_x: int, max_y: int) -> np.ndarray:
    """Return a boolean mask for cells outside the polygon interior.

    Cell (p, q) is the unit grid cell with northwest corner (p, q), where
    p = 0..max_x and q = 0..max_y in the extended grid.  The test is done at
    cell centers.  Since all polygon edges lie on grid lines, no cell center
    lies on the boundary.
    """
    pts = np.asarray(polygon, dtype=np.float64)
    x = np.arange(max_x + 1, dtype=np.float64) + 0.5
    y = np.arange(max_y + 1, dtype=np.float64) + 0.5
    xx, yy = np.meshgrid(x, y, indexing="xy")

    inside = np.zeros_like(xx, dtype=bool)
    x1 = pts[:, 0]
    y1 = pts[:, 1]
    x2 = np.roll(x1, -1)
    y2 = np.roll(y1, -1)

    for a_x, a_y, b_x, b_y in zip(x1, y1, x2, y2):
        if a_y == b_y:
            continue
        y_low = min(a_y, b_y)
        y_high = max(a_y, b_y)
        crosses_y = (yy >= y_low) & (yy < y_high)
        x_intersect = a_x + (yy - a_y) * (b_x - a_x) / (b_y - a_y)
        inside ^= crosses_y & (xx < x_intersect)

    return ~inside


def _staircase_closure_ok(cell_mask: np.ndarray, anchor_p: int, anchor_q: int, p: int, q: int) -> bool:
    """Check that the grid rectangle from the anchor cell to (p, q) is free.

    In a grid orthogonal polygon, the free staircase neighbourhood of a convex
    vertex is exactly the maximal monotone staircase of free cells on the
    exterior side of its horizontal edge.  A candidate cell belongs to that
    staircase only if every cell in the monotone rectangle between the anchor
    and the candidate is exterior to P.
    """
    p0, p1 = sorted((anchor_p, p))
    q0, q1 = sorted((anchor_q, q))
    if p0 < 0 or q0 < 0 or q1 >= cell_mask.shape[0] or p1 >= cell_mask.shape[1]:
        return False
    return bool(np.all(cell_mask[q0:q1 + 1, p0:p1 + 1]))


def _fsn_candidate_cells(
    polygon: List[Point],
    vi_idx: int,
    h_type: str,
    horizontal_other: Point,
    vertical_neighbor: Point,
    max_candidates: Optional[int] = None,
) -> List[Tuple[int, int]]:
    """Construct the Free Staircase Neighbourhood FSN(vi) from the paper.

    The paper defines FSN(vi) as the largest staircase polygon in the extended
    grid that has vi as a vertex, does not intersect the interior of P, and has
    a base edge containing the horizontal edge eH(vi).  This implementation
    constructs that grid staircase explicitly:

      1. Embed P in the extended grid with one free line on each side.
      2. Mark all unit cells whose centers are outside the interior of P.
      3. Keep only the exterior side of eH(vi), with p-values such that the
         base of the staircase contains eH(vi).
      4. Accept a cell exactly when the monotone rectangle from the base-anchor
         cell to that cell is free; the accepted cells form the maximal free
         staircase anchored at vi.

    Numpy is used for the cell-interior mask and for the free-rectangle tests.
    """
    pts = np.asarray(polygon, dtype=np.int64)
    xi, yi = map(int, pts[vi_idx])
    xo, _ = horizontal_other
    _, yv = vertical_neighbor
    xo = int(xo)
    yv = int(yv)

    max_x = int(np.max(pts[:, 0]))
    max_y = int(np.max(pts[:, 1]))
    p_min, p_max = 0, max_x
    q_min, q_max = 0, max_y

    # Select the side of vi on which the horizontal edge eH(vi) extends.
    if xo > xi:
        p_values = np.arange(max(xo - 1, p_min), p_max + 1, dtype=np.int64)
        anchor_p = max(xo - 1, p_min)
    else:
        p_values = np.arange(p_min, min(xo, p_max) + 1, dtype=np.int64)
        anchor_p = min(xo, p_max)

    # Select the exterior side of the horizontal edge eH(vi).  The vertical
    # neighbour lies on the interior side, so FSN is on the opposite side.
    if yv < yi:
        q_values = np.arange(max(yi, q_min), q_max + 1, dtype=np.int64)
        anchor_q = max(yi, q_min)
    else:
        q_values = np.arange(q_min, min(yi - 1, q_max) + 1, dtype=np.int64)
        anchor_q = min(yi - 1, q_max)

    if p_values.size == 0 or q_values.size == 0:
        return []

    exterior = _outside_cell_mask_np(polygon, max_x, max_y)
    hy, hxlo, hxhi, vx, vylo, vyhi = _edge_arrays_np(polygon)
    occupied_vertices = set(polygon)

    candidates: List[Tuple[int, int]] = []
    for q in q_values:
        for p in p_values:
            pp = int(p)
            qq = int(q)
            c_tilde = (pp + 1, qq + 1)

            if c_tilde in occupied_vertices:
                continue
            if not exterior[qq, pp]:
                continue
            if not _staircase_closure_ok(exterior, anchor_p, anchor_q, pp, qq):
                continue
            if not _rect_exterior_ok_np(xi, yi, pp, qq, hy, hxlo, hxhi, vx, vylo, vyhi):
                continue

            candidates.append((pp, qq))

    random.shuffle(candidates)
    if max_candidates is not None and len(candidates) > max_candidates:
        candidates = candidates[:max_candidates]
    return candidates


# ──────────────────────────────────────────────────────────────────
# Inflate-Paste core
# ──────────────────────────────────────────────────────────────────

def _paste_step(polygon: List[Point]) -> Optional[List[Point]]:
    """
    One Inflate-Paste iteration.

    Algorithm (paper §2, Inflate-Paste):
      1. Pick a convex vertex vi of P.
      2. Determine its horizontal edge eH(vi) and Free Staircase Neighbourhood.
      3. Select a cell C = (p, q)  in FSN(vi).
      4. Apply Inflate(P, p, q)  → P̃, ṽi, c̃ = (p+1, q+1).
      5. Apply Paste:
           If eH(vi) = vi·vi+1 : replace ṽi with (x̃i, q+1), c̃, (p+1, ỹi)
           If eH(vi) = vi-1·vi  : replace ṽi with (p+1, ỹi), c̃, (x̃i, q+1)
    """
    n = len(polygon)
    convex = find_convex_vertices(polygon)
    if not convex:
        return None

    random.shuffle(convex)

    for vi_idx in convex:
        vi = polygon[vi_idx]
        xi, yi = vi
        prev_v = polygon[(vi_idx - 1) % n]
        next_v = polygon[(vi_idx + 1) % n]

        # ── Identify horizontal edge ────────────────────────────────
        if prev_v[1] == yi:          # eH = vi-1 * vi  (prev is horizontal)
            h_type   = 'prev'
            h_other  = prev_v
            vert_nbr = next_v
        elif next_v[1] == yi:        # eH = vi * vi+1  (next is horizontal)
            h_type   = 'next'
            h_other  = next_v
            vert_nbr = prev_v
        else:
            continue                 # shouldn't happen for a convex vertex

        candidates = _fsn_candidate_cells(polygon, vi_idx, h_type, h_other, vert_nbr, max_candidates=None)

        for p, q in candidates:
            c_tilde = (p + 1, q + 1)

            inflated      = inflate(polygon, p, q)
            xi_t, yi_t    = inflated[vi_idx]
            inflated_set  = set(inflated)

            if not _paste_rectangle_is_free(inflated, vi_idx, c_tilde):
                continue

            # Build chain according to which edge is horizontal
            if h_type == 'next':
                # eH = vi * vi+1
                chain = [(xi_t, q + 1), c_tilde, (p + 1, yi_t)]
            else:
                # eH = vi-1 * vi
                chain = [(p + 1, yi_t), c_tilde, (xi_t, q + 1)]

            if not check_chain_validity(chain, inflated_set):
                continue

            candidate = replace_vertex(inflated, vi_idx, chain)

            if len(candidate) != len(inflated) + 2:
                continue

            candidate = remove_duplicates(candidate)
            if len(candidate) != len(inflated) + 2:
                continue

            if not is_simple_polygon(candidate):
                continue

            return candidate

    return None   # no valid step found this round


# ──────────────────────────────────────────────────────────────────
# Generator class
# ──────────────────────────────────────────────────────────────────

class InflatePasteGenerator:
    """
    Generates a random orthogonal polygon with exactly n vertices
    using the Inflate-Paste algorithm (Tomás & Bajuelos, EWCG 2004).
    """

    def __init__(self):
        self.polygon:          List[Point] = []
        self.execution_time:   float = 0.0   # ms
        self.peak_memory_kb:   float = 0.0

    def generate(self, n: int) -> List[List[int]]:
        if n % 2 != 0 or n < 4:
            raise ValueError("n must be even and ≥ 4")

        tracemalloc.start()
        t0 = time.time()

        # Start from the unit square (4-ogon)
        polygon = [(1, 1), (1, 2), (2, 2), (2, 1)]
        consecutive_failures = 0
        max_failures = max(200, 4 * n)

        while len(polygon) < n and consecutive_failures < max_failures:
            result = _paste_step(polygon)
            if result is not None:
                polygon = result
                consecutive_failures = 0
            else:
                consecutive_failures += 1

        if len(polygon) != n:
            tracemalloc.stop()
            raise RuntimeError(
                f"Inflate-Paste failed to generate exactly {n} vertices; got {len(polygon)}. "
                "Try another random seed or increase the FSN candidate search."
            )

        t1 = time.time()
        _, peak_mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        self.execution_time = (t1 - t0) * 1000      # ms
        self.peak_memory_kb = peak_mem / 1024

        self.polygon = remove_duplicates(polygon)
        return [[int(x), int(y)] for x, y in self.polygon]


    def draw(
        self,
        save_path: Optional[str] = "inflate_paste_polygon.png",
        show: bool = False,
        label_vertices: bool = True,
    ) -> Optional[str]:
        """Draw the polygon generated most recently by `generate`."""
        if not self.polygon:
            raise RuntimeError("No polygon has been generated yet. Call generate(n) first.")
        return draw_polygon(
            self.polygon,
            save_path=save_path,
            show=show,
            title=f"Inflate-Paste polygon, n={len(self.polygon)}",
            label_vertices=label_vertices,
            invert_y=True,
        )

    def calculate_metrics(self) -> Optional[Dict[str, Any]]:
        """Calculate geometric, validity, runtime, and memory statistics."""
        if not self.polygon:
            return None

        pts = remove_duplicates(self.polygon)
        n = len(pts)
        area = polygon_area(pts)
        perimeter = polygon_perimeter(pts)
        reflex = count_reflex_vertices(pts)
        expected_reflex = (n - 4) // 2 if n >= 4 else 0

        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        width = max(xs) - min(xs) if xs else 0
        height = max(ys) - min(ys) if ys else 0
        aspect_ratio = width / height if height != 0 else float("inf")

        hull = convex_hull(pts)
        hull_area = polygon_area(hull) if hull else 0.0
        convexity = area / hull_area if hull_area > 0 else 0.0
        per_vertex_bytes = (self.peak_memory_kb * 1024.0) / n if n else 0.0

        return {
            "n": n,
            "time_ms": self.execution_time,
            "peak_memory_kb": self.peak_memory_kb,
            "per_vertex_bytes": per_vertex_bytes,
            "area": area,
            "perimeter": perimeter,
            "width": width,
            "height": height,
            "aspect_ratio": aspect_ratio,
            "convexity": convexity,
            "reflex": reflex,
            "expected_reflex": expected_reflex,
            "simple": is_simple_polygon(pts),
            "grid_ogon": is_grid_ogon(pts),
            "valid_inflate_paste_output": (
                n >= 4
                and n % 2 == 0
                and is_simple_polygon(pts)
                and is_grid_ogon(pts)
                and reflex == expected_reflex
            ),
        }


# ──────────────────────────────────────────────────────────────────
# Quick self-test
# ──────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    NUM_RUNS = 1
    TEST_SIZES = [10]
    SUMMARY_CSV = "inflate_paste_runtime_memory.csv"
    SELECTED_STATS_CSV = "inflate_paste_selected_statistics.csv"
    SAMPLE_PLOT_N = 300
    SAMPLE_PLOT_FILE = "inflate_paste_sample_polygon.png"

    if SAMPLE_PLOT_N % 2 != 0 or SAMPLE_PLOT_N < 4:
        raise ValueError("SAMPLE_PLOT_N must be even and at least 4")

    print("Runtime and memory statistics:")
    print("╔════════╦══════════════════╦══════════════════╗")
    print("║   n    ║ time ± std (ms)  ║ memory ± std (KB)║")
    print("╠════════╬══════════════════╬══════════════════╣")

    summary_rows = []
    selected_rows = []
    ok_all = True

    for n in TEST_SIZES:
        times = []
        mems = []
        last_metrics = None
        successes = 0

        for run_id in range(NUM_RUNS):
            gen = InflatePasteGenerator()
            try:
                poly = gen.generate(n)
                pts = [tuple(p) for p in poly]
                metrics = gen.calculate_metrics()
                valid = (
                    len(pts) == n
                    and metrics is not None
                    and metrics["simple"]
                    and metrics["grid_ogon"]
                    and metrics["reflex"] == metrics["expected_reflex"]
                )
                ok_all = ok_all and valid
                successes += int(valid)
                times.append(gen.execution_time)
                mems.append(gen.peak_memory_kb)
                last_metrics = metrics
            except Exception as exc:
                ok_all = False
                print(f"║ {n:6d} ║ ERROR: {type(exc).__name__}: {exc}")

        if times:
            mean_t = sum(times) / len(times)
            std_t = math.sqrt(sum((x - mean_t) ** 2 for x in times) / len(times))
            mean_m = sum(mems) / len(mems)
            std_m = math.sqrt(sum((x - mean_m) ** 2 for x in mems) / len(mems))
            print(f"║ {n:6d} ║ {mean_t:8.2f} ± {std_t:<6.2f} ║ {mean_m:8.2f} ± {std_m:<6.2f} ║")

            summary_rows.append({
                "n": n,
                "runs": NUM_RUNS,
                "successes": successes,
                "mean_time_ms": mean_t,
                "std_time_ms": std_t,
                "mean_memory_kb": mean_m,
                "std_memory_kb": std_m,
            })

            if last_metrics is not None:
                selected_rows.append(last_metrics)

    print("╚════════╩══════════════════╩══════════════════╝")
    print()
    print("All tests:", "PASSED ✓" if ok_all else "FAILED ✗")

    if summary_rows:
        with open(SUMMARY_CSV, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
            writer.writeheader()
            writer.writerows(summary_rows)
        print(f"Saved runtime/memory summary to {SUMMARY_CSV}")

    if selected_rows:
        with open(SELECTED_STATS_CSV, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(selected_rows[0].keys()))
            writer.writeheader()
            writer.writerows(selected_rows)
        print(f"Saved selected polygon statistics to {SELECTED_STATS_CSV}")

    sample_g = InflatePasteGenerator()
    sample_g.generate(SAMPLE_PLOT_N)
    sample_g.draw(SAMPLE_PLOT_FILE, show=False, label_vertices=True)
    print(f"Saved sample polygon plot to {SAMPLE_PLOT_FILE}")
