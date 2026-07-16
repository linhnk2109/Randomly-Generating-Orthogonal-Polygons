"""
Implementation of the Inflate-Cut algorithm from:

    Ana Paula Tomas and Antonio Leslie Bajuelos,
    "Generating Random Orthogonal Polygons", LNAI 3040, 2004.

The implementation follows Fig. 5 and Section 3 of the paper:
start with the unit square, repeatedly choose a unit grid cell C in the
interior of the current polygon, apply Inflate using the north-west corner
(p, q) of C, then apply Cut by removing a rectangle from the inflated polygon.
Each successful iteration adds exactly two vertices.
"""

import random
import ctypes
import subprocess
import time
import tracemalloc
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

Point = Tuple[int, int]
FloatPoint = Tuple[float, float]
_KERNEL = None


def _load_c_kernel():
    global _KERNEL
    if _KERNEL is not None:
        return _KERNEL

    base_dir = Path(__file__).resolve().parent
    source = base_dir / "inflate_cut_kernel.c"
    library = base_dir / "inflate_cut_kernel.so"

    if not source.exists():
        _KERNEL = False
        return None

    try:
        if (not library.exists()) or source.stat().st_mtime > library.stat().st_mtime:
            subprocess.run(
                [
                    "gcc",
                    "-O3",
                    "-std=c99",
                    "-shared",
                    "-fPIC",
                    str(source),
                    "-o",
                    str(library),
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

        kernel = ctypes.CDLL(str(library))
        kernel.ic_find_cut_candidate.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.int64, ndim=2, flags="C_CONTIGUOUS"),
            ctypes.c_int,
            np.ctypeslib.ndpointer(dtype=np.int64, ndim=2, flags="C_CONTIGUOUS"),
            ctypes.c_int,
            np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(dtype=np.int64, ndim=2, flags="C_CONTIGUOUS"),
        ]
        kernel.ic_find_cut_candidate.restype = ctypes.c_int
        kernel.ic_collect_interior_cells.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.int64, ndim=2, flags="C_CONTIGUOUS"),
            ctypes.c_int,
            np.ctypeslib.ndpointer(dtype=np.int64, ndim=2, flags="C_CONTIGUOUS"),
            ctypes.c_int,
        ]
        kernel.ic_collect_interior_cells.restype = ctypes.c_int
        _KERNEL = kernel
        return kernel
    except Exception:
        _KERNEL = False
        return None


def remove_duplicates(polygon: Sequence[Point]) -> List[Point]:
    """Remove consecutive duplicate vertices, including a repeated closing one."""
    if not polygon:
        return []

    result = [tuple(polygon[0])]
    for point in polygon[1:]:
        point = tuple(point)
        if point != result[-1]:
            result.append(point)

    if len(result) > 1 and result[0] == result[-1]:
        result.pop()
    return result


def signed_area(polygon: Sequence[Point]) -> float:
    area2 = 0
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        area2 += x1 * y2 - x2 * y1
    return area2 / 2.0


def point_in_polygon(point: FloatPoint, polygon: Sequence[Point]) -> bool:
    """Return True when point is strictly inside polygon."""
    x, y = point
    inside = False
    n = len(polygon)

    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]

        if ((y1 > y) != (y2 > y)):
            x_cross = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            if x_cross > x:
                inside = not inside

    return inside


def _orientation(a: Point, b: Point, c: Point) -> int:
    value = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return (value > 0) - (value < 0)


def _on_segment(a: Point, b: Point, p: Point) -> bool:
    return (
        min(a[0], b[0]) <= p[0] <= max(a[0], b[0])
        and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])
        and _orientation(a, b, p) == 0
    )


def _segments_intersect(a: Point, b: Point, c: Point, d: Point) -> bool:
    o1 = _orientation(a, b, c)
    o2 = _orientation(a, b, d)
    o3 = _orientation(c, d, a)
    o4 = _orientation(c, d, b)

    if o1 != o2 and o3 != o4:
        return True
    if o1 == 0 and _on_segment(a, b, c):
        return True
    if o2 == 0 and _on_segment(a, b, d):
        return True
    if o3 == 0 and _on_segment(c, d, a):
        return True
    if o4 == 0 and _on_segment(c, d, b):
        return True
    return False


def is_simple_polygon(polygon: Sequence[Point]) -> bool:
    """Validate that the vertex cycle is a simple orthogonal polygon."""
    polygon = remove_duplicates(polygon)
    n = len(polygon)
    if n < 4 or n % 2 != 0 or len(set(polygon)) != n:
        return False

    for i in range(n):
        a, b = polygon[i], polygon[(i + 1) % n]
        if a == b or (a[0] != b[0] and a[1] != b[1]):
            return False

    for i in range(n):
        a, b = polygon[i], polygon[(i + 1) % n]
        for j in range(i + 1, n):
            if j == i or j == (i + 1) % n or i == (j + 1) % n:
                continue
            c, d = polygon[j], polygon[(j + 1) % n]
            if _segments_intersect(a, b, c, d):
                return False

    return abs(signed_area(polygon)) > 0


def inflate(polygon: Sequence[Point], p: int, q: int) -> List[Point]:
    """
    Apply the paper's Inflate transformation using cell C with north-west
    corner (p, q).
    """
    return [
        (x if x <= p else x + 1, y if y <= q else y + 1)
        for x, y in polygon
    ]


def _as_np_polygon(polygon: Sequence[Point]) -> np.ndarray:
    return np.asarray(polygon, dtype=np.int64).reshape((-1, 2))


def _remove_consecutive_duplicates_np(polygon: np.ndarray) -> np.ndarray:
    if len(polygon) == 0:
        return polygon.reshape((0, 2))

    keep = np.ones(len(polygon), dtype=bool)
    keep[1:] = np.any(polygon[1:] != polygon[:-1], axis=1)
    result = polygon[keep]
    if len(result) > 1 and np.array_equal(result[0], result[-1]):
        result = result[:-1]
    return result


def _inflate_np(polygon: np.ndarray, p: int, q: int) -> np.ndarray:
    """Vectorized Inflate transformation from the paper."""
    inflated = polygon.copy()
    inflated[:, 0] += inflated[:, 0] > p
    inflated[:, 1] += inflated[:, 1] > q
    return inflated


def replace_vertex(polygon: Sequence[Point], index: int, chain: Sequence[Point]) -> List[Point]:
    return remove_duplicates(list(polygon[:index]) + list(chain) + list(polygon[index + 1:]))


def _interior_cells(polygon: Sequence[Point]) -> List[Tuple[int, int]]:
    """
    Enumerate unit grid cells inside P. A cell is represented by the coordinates
    (p, q) of its north-west corner, exactly as in the paper.
    """
    xs = [p[0] for p in polygon]
    ys = [p[1] for p in polygon]
    cells: List[Tuple[int, int]] = []

    for q in range(min(ys), max(ys)):
        for p in range(min(xs), max(xs)):
            if point_in_polygon((p + 0.5, q + 0.5), polygon):
                cells.append((p, q))

    return cells


def _interior_cells_np(polygon: np.ndarray) -> np.ndarray:
    """
    Scanline enumeration of unit cells inside P.

    For each row q, test the horizontal line y=q+0.5 against vertical edges.
    Consecutive pairs of intersections delimit interior cell spans. This is the
    same cell set as Fig. 5 requires, but avoids point-in-polygon work for every
    cell independently.
    """
    if len(polygon) < 4:
        return np.empty((0, 2), dtype=np.int64)

    starts = polygon
    ends = np.roll(polygon, -1, axis=0)
    vertical = starts[:, 0] == ends[:, 0]
    vx = starts[vertical, 0]
    vy1 = np.minimum(starts[vertical, 1], ends[vertical, 1])
    vy2 = np.maximum(starts[vertical, 1], ends[vertical, 1])

    cells = []
    for q in range(int(polygon[:, 1].min()), int(polygon[:, 1].max())):
        y_mid = q + 0.5
        crosses = vx[(vy1 < y_mid) & (y_mid < vy2)]
        if len(crosses) < 2:
            continue
        crosses = np.sort(crosses)
        for left, right in zip(crosses[0::2], crosses[1::2]):
            if right > left:
                ps = np.arange(int(left), int(right), dtype=np.int64)
                qs = np.full(len(ps), q, dtype=np.int64)
                cells.append(np.column_stack((ps, qs)))

    if not cells:
        return np.empty((0, 2), dtype=np.int64)
    return np.vstack(cells)


def _ray_intersections(
    polygon: Sequence[Point],
    center: Point,
) -> List[Tuple[int, Point, Point, Point, int]]:
    """
    Compute the first boundary hit for the four rays from center.

    Returns tuples:
        (direction_index, intersection s, edge_start, edge_end, edge_start_index)
    """
    cx, cy = center
    hits: List[Tuple[int, int, Point, Point, Point, int]] = []
    n = len(polygon)

    for i in range(n):
        a = polygon[i]
        b = polygon[(i + 1) % n]

        if a[0] == b[0]:
            x = a[0]
            y_low, y_high = sorted((a[1], b[1]))
            if y_low < cy < y_high:
                if x > cx:
                    hits.append((0, x - cx, (x, cy), a, b, i))
                elif x < cx:
                    hits.append((1, cx - x, (x, cy), a, b, i))
        elif a[1] == b[1]:
            y = a[1]
            x_low, x_high = sorted((a[0], b[0]))
            if x_low < cx < x_high:
                if y > cy:
                    hits.append((2, y - cy, (cx, y), a, b, i))
                elif y < cy:
                    hits.append((3, cy - y, (cx, y), a, b, i))

    first_by_direction = {}
    for direction, distance, s, a, b, i in hits:
        old = first_by_direction.get(direction)
        if old is None or distance < old[0]:
            first_by_direction[direction] = (distance, s, a, b, i)

    intersections = [
        (direction, data[1], data[2], data[3], data[4])
        for direction, data in first_by_direction.items()
    ]
    random.shuffle(intersections)
    return intersections


def _ray_intersections_np(
    polygon: np.ndarray,
    center: Tuple[int, int],
) -> List[Tuple[int, Point, Point, Point, int]]:
    cx, cy = center
    starts = polygon
    ends = np.roll(polygon, -1, axis=0)
    result: List[Tuple[int, Point, Point, Point, int]] = []

    vertical = starts[:, 0] == ends[:, 0]
    if np.any(vertical):
        idx = np.flatnonzero(vertical)
        x = starts[idx, 0]
        y_low = np.minimum(starts[idx, 1], ends[idx, 1])
        y_high = np.maximum(starts[idx, 1], ends[idx, 1])
        spans_center = (y_low < cy) & (cy < y_high)

        right = np.flatnonzero(spans_center & (x > cx))
        if len(right):
            local = right[np.argmin(x[right] - cx)]
            edge_i = int(idx[local])
            sx = int(x[local])
            result.append((0, (sx, cy), tuple(starts[edge_i]), tuple(ends[edge_i]), edge_i))

        left = np.flatnonzero(spans_center & (x < cx))
        if len(left):
            local = left[np.argmin(cx - x[left])]
            edge_i = int(idx[local])
            sx = int(x[local])
            result.append((1, (sx, cy), tuple(starts[edge_i]), tuple(ends[edge_i]), edge_i))

    horizontal = starts[:, 1] == ends[:, 1]
    if np.any(horizontal):
        idx = np.flatnonzero(horizontal)
        y = starts[idx, 1]
        x_low = np.minimum(starts[idx, 0], ends[idx, 0])
        x_high = np.maximum(starts[idx, 0], ends[idx, 0])
        spans_center = (x_low < cx) & (cx < x_high)

        up = np.flatnonzero(spans_center & (y > cy))
        if len(up):
            local = up[np.argmin(y[up] - cy)]
            edge_i = int(idx[local])
            sy = int(y[local])
            result.append((2, (cx, sy), tuple(starts[edge_i]), tuple(ends[edge_i]), edge_i))

        down = np.flatnonzero(spans_center & (y < cy))
        if len(down):
            local = down[np.argmin(cy - y[down])]
            edge_i = int(idx[local])
            sy = int(y[local])
            result.append((3, (cx, sy), tuple(starts[edge_i]), tuple(ends[edge_i]), edge_i))

    random.shuffle(result)
    return result


def _rectangle_contains_only_endpoint(
    polygon: Sequence[Point],
    center: Point,
    endpoint: Point,
) -> bool:
    x1, x2 = sorted((center[0], endpoint[0]))
    y1, y2 = sorted((center[1], endpoint[1]))

    for vertex in polygon:
        if vertex == endpoint:
            continue
        if x1 <= vertex[0] <= x2 and y1 <= vertex[1] <= y2:
            return False
    return True


def _rectangle_contains_only_endpoint_np(
    polygon: np.ndarray,
    center: Point,
    endpoint: Point,
) -> bool:
    x1, x2 = sorted((center[0], endpoint[0]))
    y1, y2 = sorted((center[1], endpoint[1]))
    inside = (
        (x1 <= polygon[:, 0]) & (polygon[:, 0] <= x2)
        & (y1 <= polygon[:, 1]) & (polygon[:, 1] <= y2)
    )
    is_endpoint = (polygon[:, 0] == endpoint[0]) & (polygon[:, 1] == endpoint[1])
    return not bool(np.any(inside & ~is_endpoint))


def _replace_vertex_np(polygon: np.ndarray, index: int, chain: Sequence[Point]) -> np.ndarray:
    chain_array = np.asarray(chain, dtype=np.int64).reshape((-1, 2))
    result = np.concatenate((polygon[:index], chain_array, polygon[index + 1:]), axis=0)
    return _remove_consecutive_duplicates_np(result)


def _cut_candidates(
    polygon: Sequence[Point],
    center: Point,
) -> List[List[Point]]:
    """
    Generate all valid results of applying Cut to the inflated polygon.
    The randomization mirrors the paper: intersection points and rectangle
    endpoints are tried in random order.
    """
    candidates: List[List[Point]] = []

    for _, s, edge_start, edge_end, edge_index in _ray_intersections(polygon, center):
        endpoint_options = [(edge_start, edge_index), (edge_end, (edge_index + 1) % len(polygon))]
        random.shuffle(endpoint_options)

        for endpoint, endpoint_index in endpoint_options:
            if not _rectangle_contains_only_endpoint(polygon, center, endpoint):
                continue

            s_prime = (
                center[0] + (endpoint[0] - s[0]),
                center[1] + (endpoint[1] - s[1]),
            )

            if endpoint_index == edge_index:
                chain = [s_prime, center, s]
            else:
                chain = [s, center, s_prime]

            if len(set(chain)) != 3:
                continue

            candidate = replace_vertex(polygon, endpoint_index, chain)
            if len(candidate) == len(polygon) + 2 and is_simple_polygon(candidate):
                candidates.append(candidate)

    random.shuffle(candidates)
    return candidates


def inflate_cut_step(polygon: Sequence[Point]) -> Optional[List[Point]]:
    """
    One successful Fig. 5 iteration, or None if all current interior cells fail.
    """
    cells = _interior_cells(polygon)
    random.shuffle(cells)

    for p, q in cells:
        inflated = inflate(polygon, p, q)
        center = (p + 1, q + 1)
        cut_results = _cut_candidates(inflated, center)
        if cut_results:
            return cut_results[0]

    return None


def _inflate_cut_step_np(polygon: np.ndarray) -> Optional[np.ndarray]:
    """
    NumPy-backed version of one Fig. 5 iteration.

    It keeps the same random choices as the paper: random interior cell, random
    ray intersection, then random endpoint of that boundary edge when possible.
    """
    cells = _interior_cells_np(polygon)
    if len(cells) == 0:
        return None

    for cell_index in np.random.permutation(len(cells)):
        p, q = map(int, cells[cell_index])
        inflated = _inflate_np(polygon, p, q)
        center = (p + 1, q + 1)

        for _, s, edge_start, edge_end, edge_index in _ray_intersections_np(inflated, center):
            endpoint_options = [
                (edge_start, edge_index),
                (edge_end, (edge_index + 1) % len(inflated)),
            ]
            random.shuffle(endpoint_options)

            for endpoint, endpoint_index in endpoint_options:
                if not _rectangle_contains_only_endpoint_np(inflated, center, endpoint):
                    continue

                s_prime = (
                    center[0] + (endpoint[0] - s[0]),
                    center[1] + (endpoint[1] - s[1]),
                )
                chain = [s_prime, center, s] if endpoint_index == edge_index else [s, center, s_prime]
                if len(set(chain)) != 3:
                    continue

                candidate = _replace_vertex_np(inflated, endpoint_index, chain)
                if len(candidate) == len(inflated) + 2:
                    return candidate

    return None


def _inflate_cut_step_c(polygon: np.ndarray) -> Optional[np.ndarray]:
    """
    C-backed version of one Fig. 5 iteration.

    Python still performs the paper-level random choice of interior cells.
    The C kernel only evaluates Inflate and Cut candidates for that shuffled
    cell order.
    """
    kernel = _load_c_kernel()
    if not kernel:
        return _inflate_cut_step_np(polygon)

    polygon_c = np.ascontiguousarray(polygon, dtype=np.int64)
    xs = polygon_c[:, 0]
    ys = polygon_c[:, 1]
    max_cells = max(int((xs.max() - xs.min()) * (ys.max() - ys.min())), 1)
    cells_buf = np.empty((max_cells, 2), dtype=np.int64)
    cell_count = kernel.ic_collect_interior_cells(polygon_c, len(polygon_c), cells_buf, max_cells)
    if cell_count <= 0:
        return None

    cells = np.ascontiguousarray(cells_buf[:cell_count][np.random.permutation(cell_count)], dtype=np.int64)
    dir_order = np.ascontiguousarray(np.random.permutation(4), dtype=np.int32)
    endpoint_order = np.ascontiguousarray(np.random.permutation(2), dtype=np.int32)
    out = np.empty((len(polygon) + 2, 2), dtype=np.int64)

    ok = kernel.ic_find_cut_candidate(
        polygon_c,
        len(polygon_c),
        cells,
        len(cells),
        dir_order,
        endpoint_order,
        out,
    )
    if ok:
        return _remove_consecutive_duplicates_np(out)

    return None


class InflateCutGenerator:
    """Generate random grid n-ogons using the Inflate-Cut algorithm."""

    def __init__(self):
        self.polygon: List[Point] = []
        self.execution_time = 0.0
        self.peak_memory_kb = 0.0

    def generate(self, n: int) -> List[List[int]]:
        if n % 2 != 0 or n < 4:
            raise ValueError("n must be even and >= 4")

        if tracemalloc.is_tracing():
            tracemalloc.stop()
        tracemalloc.start()
        start = time.perf_counter()

        # Unit square, in the same grid coordinates used by the paper.
        polygon = np.asarray([(1, 1), (1, 2), (2, 2), (2, 1)], dtype=np.int64)
        remaining_reflex_vertices = n // 2 - 2

        while remaining_reflex_vertices > 0:
            result = _inflate_cut_step_c(polygon)
            if result is None:
                raise RuntimeError("Inflate-Cut failed to find a cuttable interior cell")
            polygon = result
            remaining_reflex_vertices -= 1

        end = time.perf_counter()
        _, peak_mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        self.execution_time = (end - start) * 1000.0
        self.peak_memory_kb = peak_mem / 1024.0
        self.polygon = [tuple(map(int, point)) for point in _remove_consecutive_duplicates_np(polygon)]
        return [[int(x), int(y)] for x, y in self.polygon]

    def visualize(self, save_path: Optional[str] = None):
        import matplotlib.pyplot as plt

        if not self.polygon:
            print("No polygon to visualize!")
            return

        xs = [p[0] for p in self.polygon] + [self.polygon[0][0]]
        ys = [p[1] for p in self.polygon] + [self.polygon[0][1]]

        plt.figure(figsize=(10, 10))
        plt.fill(xs, ys, alpha=0.4, color="lightcoral", edgecolor="darkred", linewidth=2)
        plt.plot(xs[:-1], ys[:-1], "ro", markersize=6)

        for i, (x, y) in enumerate(self.polygon):
            plt.annotate(str(i), (x, y), xytext=(5, 5), textcoords="offset points", fontsize=8)

        plt.axis("equal")
        plt.grid(True, alpha=0.3, linestyle="--")
        plt.title(f"Inflate-Cut Algorithm: {len(self.polygon)} Vertices")

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            print(f"Saved to {save_path}")
        plt.show()


if __name__ == "__main__":
    for n_value in [50000, 100000]:
        generator = InflateCutGenerator()
        polygon_value = generator.generate(n_value)
        print(
            f"n={n_value:4d} achieved={len(polygon_value):4d} "
            f"simple={is_simple_polygon([tuple(p) for p in polygon_value])} "
            f"time={generator.execution_time:.2f} ms "
            f"mem={generator.peak_memory_kb:.2f} KB"
        )
