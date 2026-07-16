#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>

typedef struct {
    int dir;
    int sx;
    int sy;
    int edge_index;
} Hit;

static int cmp_int64(const void *a, const void *b) {
    int64_t av = *(const int64_t *)a;
    int64_t bv = *(const int64_t *)b;
    return (av > bv) - (av < bv);
}

int ic_collect_interior_cells(
    const int64_t *polygon,
    int n,
    int64_t *out_cells,
    int max_cells
) {
    if (n < 4 || max_cells <= 0) {
        return 0;
    }

    int64_t min_y = polygon[1];
    int64_t max_y = polygon[1];
    for (int i = 1; i < n; ++i) {
        int64_t y = polygon[2 * i + 1];
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }

    int64_t *crosses = (int64_t *)malloc((size_t)n * sizeof(int64_t));
    if (!crosses) {
        return 0;
    }

    int count = 0;
    for (int64_t q = min_y; q < max_y; ++q) {
        int cross_count = 0;
        for (int i = 0; i < n; ++i) {
            int j = (i + 1) % n;
            int64_t ax = polygon[2 * i];
            int64_t ay = polygon[2 * i + 1];
            int64_t bx = polygon[2 * j];
            int64_t by = polygon[2 * j + 1];

            if (ax != bx) {
                continue;
            }

            int64_t y_low = ay < by ? ay : by;
            int64_t y_high = ay < by ? by : ay;
            if (y_low <= q && q < y_high) {
                crosses[cross_count++] = ax;
            }
        }

        if (cross_count < 2) {
            continue;
        }

        qsort(crosses, (size_t)cross_count, sizeof(int64_t), cmp_int64);
        for (int k = 0; k + 1 < cross_count; k += 2) {
            int64_t left = crosses[k];
            int64_t right = crosses[k + 1];
            for (int64_t p = left; p < right; ++p) {
                if (count >= max_cells) {
                    free(crosses);
                    return count;
                }
                out_cells[2 * count] = p;
                out_cells[2 * count + 1] = q;
                count++;
            }
        }
    }

    free(crosses);
    return count;
}

static int rect_contains_only_endpoint(
    const int64_t *poly,
    int n,
    int cx,
    int cy,
    int ex,
    int ey
) {
    int x1 = cx < ex ? cx : ex;
    int x2 = cx < ex ? ex : cx;
    int y1 = cy < ey ? cy : ey;
    int y2 = cy < ey ? ey : cy;

    for (int i = 0; i < n; ++i) {
        int vx = (int)poly[2 * i];
        int vy = (int)poly[2 * i + 1];
        if (vx == ex && vy == ey) {
            continue;
        }
        if (x1 <= vx && vx <= x2 && y1 <= vy && vy <= y2) {
            return 0;
        }
    }
    return 1;
}

static int collect_ray_hits(
    const int64_t *poly,
    int n,
    int cx,
    int cy,
    Hit hits[4]
) {
    for (int h = 0; h < 4; ++h) {
        hits[h].dir = -1;
        hits[h].sx = 0;
        hits[h].sy = 0;
        hits[h].edge_index = -1;
    }

    int best_dist[4] = {0, 0, 0, 0};
    int found[4] = {0, 0, 0, 0};

    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        int ax = (int)poly[2 * i];
        int ay = (int)poly[2 * i + 1];
        int bx = (int)poly[2 * j];
        int by = (int)poly[2 * j + 1];

        if (ax == bx) {
            int y_low = ay < by ? ay : by;
            int y_high = ay < by ? by : ay;
            if (y_low < cy && cy < y_high) {
                int dir = -1;
                int dist = 0;
                if (ax > cx) {
                    dir = 0;
                    dist = ax - cx;
                } else if (ax < cx) {
                    dir = 1;
                    dist = cx - ax;
                }
                if (dir >= 0 && (!found[dir] || dist < best_dist[dir])) {
                    found[dir] = 1;
                    best_dist[dir] = dist;
                    hits[dir].dir = dir;
                    hits[dir].sx = ax;
                    hits[dir].sy = cy;
                    hits[dir].edge_index = i;
                }
            }
        } else if (ay == by) {
            int x_low = ax < bx ? ax : bx;
            int x_high = ax < bx ? bx : ax;
            if (x_low < cx && cx < x_high) {
                int dir = -1;
                int dist = 0;
                if (ay > cy) {
                    dir = 2;
                    dist = ay - cy;
                } else if (ay < cy) {
                    dir = 3;
                    dist = cy - ay;
                }
                if (dir >= 0 && (!found[dir] || dist < best_dist[dir])) {
                    found[dir] = 1;
                    best_dist[dir] = dist;
                    hits[dir].dir = dir;
                    hits[dir].sx = cx;
                    hits[dir].sy = ay;
                    hits[dir].edge_index = i;
                }
            }
        }
    }

    return found[0] + found[1] + found[2] + found[3];
}

static int build_candidate(
    const int64_t *inflated,
    int n,
    int endpoint_index,
    int sx,
    int sy,
    int cx,
    int cy,
    int ex,
    int ey,
    int endpoint_is_edge_start,
    int64_t *out_poly
) {
    int spx = cx + (ex - sx);
    int spy = cy + (ey - sy);

    if ((spx == cx && spy == cy) || (sx == cx && sy == cy) || (spx == sx && spy == sy)) {
        return 0;
    }

    if (!((spx == cx && cy == sy) || (spy == cy && cx == sx))) {
        return 0;
    }

    int out_i = 0;
    for (int i = 0; i < n; ++i) {
        if (i != endpoint_index) {
            out_poly[2 * out_i] = inflated[2 * i];
            out_poly[2 * out_i + 1] = inflated[2 * i + 1];
            out_i++;
            continue;
        }

        if (endpoint_is_edge_start) {
            out_poly[2 * out_i] = spx;
            out_poly[2 * out_i + 1] = spy;
            out_i++;
            out_poly[2 * out_i] = cx;
            out_poly[2 * out_i + 1] = cy;
            out_i++;
            out_poly[2 * out_i] = sx;
            out_poly[2 * out_i + 1] = sy;
            out_i++;
        } else {
            out_poly[2 * out_i] = sx;
            out_poly[2 * out_i + 1] = sy;
            out_i++;
            out_poly[2 * out_i] = cx;
            out_poly[2 * out_i + 1] = cy;
            out_i++;
            out_poly[2 * out_i] = spx;
            out_poly[2 * out_i + 1] = spy;
            out_i++;
        }
    }

    return out_i == n + 2;
}

int ic_find_cut_candidate(
    const int64_t *polygon,
    int n,
    const int64_t *cells,
    int cell_count,
    const int *dir_order,
    const int *endpoint_order,
    int64_t *out_poly
) {
    if (n < 4 || cell_count <= 0) {
        return 0;
    }

    int64_t *inflated = (int64_t *)malloc((size_t)n * 2 * sizeof(int64_t));
    if (!inflated) {
        return 0;
    }

    int result = 0;
    for (int c = 0; c < cell_count; ++c) {
        int p = (int)cells[2 * c];
        int q = (int)cells[2 * c + 1];
        int cx = p + 1;
        int cy = q + 1;

        for (int i = 0; i < n; ++i) {
            int64_t x = polygon[2 * i];
            int64_t y = polygon[2 * i + 1];
            inflated[2 * i] = x + (x > p);
            inflated[2 * i + 1] = y + (y > q);
        }

        Hit hits[4];
        int found_count = collect_ray_hits(inflated, n, cx, cy, hits);
        if (found_count == 0) {
            continue;
        }

        for (int di = 0; di < 4; ++di) {
            int dir = dir_order[di];
            if (dir < 0 || dir > 3 || hits[dir].dir != dir || hits[dir].edge_index < 0) {
                continue;
            }

            int edge_index = hits[dir].edge_index;
            int next_index = (edge_index + 1) % n;
            int sx = hits[dir].sx;
            int sy = hits[dir].sy;

            for (int eo = 0; eo < 2; ++eo) {
                int endpoint_is_edge_start = endpoint_order[eo] == 0;
                int endpoint_index = endpoint_is_edge_start ? edge_index : next_index;
                int ex = (int)inflated[2 * endpoint_index];
                int ey = (int)inflated[2 * endpoint_index + 1];

                if ((ex == sx && ey == sy) || (ex == cx && ey == cy)) {
                    continue;
                }

                if (!rect_contains_only_endpoint(inflated, n, cx, cy, ex, ey)) {
                    continue;
                }

                if (build_candidate(
                        inflated, n, endpoint_index, sx, sy, cx, cy, ex, ey,
                        endpoint_is_edge_start, out_poly)) {
                    result = 1;
                    break;
                }
            }
            if (result) {
                break;
            }
        }
        if (result) {
            break;
        }
    }

    free(inflated);
    return result;
}
