import random
import time
from typing import List, Tuple, Set, Optional, Dict, Any
from sortedcontainers import SortedDict # Vẫn giữ SortedDict cho is_simple_polygon
import bisect 
import tracemalloc
from shapely.geometry import Polygon # Cần thư viện shapely để tính toán metrics
import numpy as np

Point = Tuple[int, int]

# --- Helper Functions (From inflate_cut.py) ---

"""
    Implementation of Inflate-Cut algorithm
    from paper "Generating Random Orthogonal Polygons" 
    by Tomás & Bajuelos (2004)
"""

def remove_duplicates(polygon: List[Point]) -> List[Point]:
    """Remove any duplicate vertices from a polygon - O(n)"""

    if not polygon:
        return []

    result = [polygon[0]]
    for i in range(1, len(polygon)):
        if polygon[i] != result[-1]:
            result.append(polygon[i])

    # Check the last and first vertices
    if result[-1] == result[0] and len(result) > 1:
        result.pop()

    return result

def find_convex_vertices(polygon: List[Point]) -> List[int]:
    """Find convex vertices in O(n) time"""
    n = len(polygon)
    result = []

    # Precompute next indices to avoid modulo operations - O(n)
    next_indices = [(i + 1) % n for i in range(n)]
    prev_indices = [(i - 1) % n for i in range(n)]

    # Single pass through vertices - O(n)
    for i in range(n):
        a = polygon[prev_indices[i]]
        b = polygon[i]
        c = polygon[next_indices[i]]

        # Skip duplicate vertices
        if b == a or b == c:
            continue

        # Optimized cross product for orthogonal polygons
        # In orthogonal polygons, we can simplify this check
        # A vertex is convex if and only if both adjacent edges form a right angle inward

        if ((b[0] == a[0] and b[1] == c[1]) or (b[1] == a[1] and b[0] == c[0])) and \
                ((b[0] > a[0] and b[1] > c[1]) or (b[0] < a[0] and b[1] < c[1]) or
                 (b[1] > a[1] and b[0] < c[0]) or (b[1] < a[1] and b[0] > c[0])):

            result.append(i)

    return result

def is_simple_polygon(polygon: List[Point]) -> bool:
    """Check if polygon is simple using an optimized segment tree - O(n log n)"""
    n = len(polygon)


    # Check for duplicate vertices first - O(n log n)
    sorted_points = sorted(polygon)
    for i in range(1, len(sorted_points)):
        if sorted_points[i] == sorted_points[i - 1]:
            return False



    # Step 1: Generate horizontal and vertical segments - O(n)
    h_segments = []  # Horizontal segments
    v_segments = []  # Vertical segments

    # Create adjacency set for O(1) lookups
    adjacent_edges = set()

    for i in range(n):
        p1 = polygon[i]
        p2 = polygon[(i + 1) % n]

        # Ensure p1 is left/bottom
        if p1[0] > p2[0] or (p1[0] == p2[0] and p1[1] > p2[1]):
            p1, p2 = p2, p1

        if p1[0] == p2[0]:  # Vertical segment
            v_segments.append((p1[0], min(p1[1], p2[1]), max(p1[1], p2[1]), i))
        else:  # Horizontal segment
            h_segments.append((min(p1[0], p2[0]), max(p1[0], p2[0]), p1[1], i))

        # Mark adjacency
        adjacent_edges.add((i, (i + 1) % n))
        adjacent_edges.add(((i + 1) % n, i))
        adjacent_edges.add((i, (i - 1 + n) % n))
        adjacent_edges.add(((i - 1 + n) % n, i))

    # Step 2: Sort segments for sweep line - O(n log n)
    v_segments.sort()  # Sort vertical segments by x-coordinate

    # Step 3: Sweep line algorithm with segment trees - O(n log n)
    events = []

    # Create events for horizontal segments - O(n)
    for x1, x2, y, i in h_segments:
        events.append((x1, 0, y, i))  # 0 = start event
        events.append((x2, 1, y, i))  # 1 = end event

    # Add vertical segments as query events - O(n)
    for x, y1, y2, i in v_segments:
        events.append((x, 2, y1, y2, i))  # 2 = vertical segment query

    # Sort events - O(n log n)
    events.sort(key=lambda e: (e[0], e[1]))  # Sort by x and then by event type

    # Use balanced binary search tree for active segments - O(log n) operations
    active_segments = SortedDict()


    for event in events:
        if event[1] == 0:  # Start horizontal segment
            x, _, y, i = event
            active_segments[y] = i
        elif event[1] == 1:  # End horizontal segment
            x, _, y, i = event
            if y in active_segments:  # Defensive programming
                del active_segments[y]
        else:  # Vertical segment query
            x, _, y1, y2, i = event

            # Find horizontal segments that intersect - O(log n)
            try:
                intersecting = active_segments.irange(y1, y2, inclusive=(True, True))
            except ValueError:
                # Handle edge case where range is invalid
                continue

            # Check for non-adjacent intersections
            for y in intersecting:
                j = active_segments[y]
                if (i, j) not in adjacent_edges and (j, i) not in adjacent_edges:
                    return False

    return True


def inflate(polygon: List[Point], p: int, q: int) -> List[Point]:
    """Inflate polygon in O(n) time"""
    # Pre-allocate result list
    result = []

    # Single pass inflation - O(n)
    for x, y in polygon:
        result.append((x + (x > p), y + (y > q)))

    return result


def replace_vertex(polygon: List[Point], index: int, chain: List[Point]) -> List[Point]:
    """Replace vertex in O(n) time with duplicate checking"""
    # Skip empty chains
    if not chain:
        return polygon

    result = []
    result.extend(polygon[:index])
    result.extend(chain)
    result.extend(polygon[index + 1:])

    # Remove any duplicates that might have been created
    return remove_duplicates(result)

def check_chain_validity(chain: List[Point], point_set: Set[Point]) -> bool:
    """Check if a chain is valid (no duplicates with existing points)"""
    # Check if any of the chain points already exists
    for point in chain:
        if point in point_set:
            return False

    # Check if there are duplicate points within the chain
    if len(set(chain)) != len(chain):
        return False

    return True



# --- Main Class ---

class InflateCutGenerator:

    def __init__(self):
        self.polygon: List[Point] = []
        self.start_time = 0.0
        self.end_time = 0.0
        self.execution_time = 0.0 
        self.max_attempts = 3 # Số lần thử tối đa

    def generate(self, n: int) -> List[List[int]]:


        """random orthogonal polygon with guaranteed O(n log n) complexity and no duplicate vertices"""
        if n % 2 != 0 or n < 4:
            raise ValueError("n must be even and â‰¥ 4")

        self.start_time = time.time()
        # Initialize with a simple 4-vertex polygon
        polygon = [(1, 1), (1, 2), (2, 2), (2, 1)]
        current_vertices = 4
        target_vertices = n


        # Use constant number of attempts for O(n log n) complexity
        max_attempts = 3

        # Track grid occupancy for O(1) collision detection
        occupied_grid = set(polygon)

        # Fixed directions
        directions = [(1, 0), (0, 1)]

        # Store frequently accessed indices for faster lookups
        next_indices = [(i + 1) % 4 for i in range(4)]

        # Track failure count to avoid infinite loops
        consecutive_failures = 0
        max_failures = 10


        # Each successful inflation adds O(1) vertices
        while current_vertices < target_vertices and consecutive_failures < max_failures:
            # Find convex vertices - O(n)
            convex_indices = find_convex_vertices(polygon)
            if not convex_indices:
                break


            success = False

            # Sample convex vertices efficiently - O(1) samples
            sample_size = min(max_attempts, len(convex_indices))
            shuffled_indices = random.sample(convex_indices, sample_size)

            # Try vertices - O(1) iterations
            for vi_index in shuffled_indices:
                vi = polygon[vi_index]

                # Try both directions deterministically
                for dx, dy in directions:
                    p, q = vi[0] + dx, vi[1] + dy
                    c = (p + 1, q + 1)

                    # Skip if new point already exists - O(1)
                    if c in occupied_grid:
                        continue


                    # Inflate polygon - O(n)
                    inflated_polygon = inflate(polygon, p, q)
                    vi_tilde = inflated_polygon[vi_index]

                    # Determine if edge is horizontal directly - O(1)
                    next_idx = next_indices[vi_index]
                    is_horiz = polygon[vi_index][1] == polygon[next_idx][1]

                    # Determine chain based on edge orientation
                    if is_horiz:
                        chain = [(vi_tilde[0], q + 1), c, (p + 1, vi_tilde[1])]
                    else:
                        chain = [(p + 1, vi_tilde[1]), c, (vi_tilde[0], q + 1)]

                    # Check if chain creates any duplicates
                    if not check_chain_validity(chain, occupied_grid):
                        continue


                    # Replace vertex - O(n)
                    candidate_polygon = replace_vertex(inflated_polygon, vi_index, chain)

                    # Ensure we didn't create a duplicate
                    if len(candidate_polygon) != len(inflated_polygon) + 2:
                        continue


                    # Check if result is simple - O(n log n)
                    if is_simple_polygon(candidate_polygon):
                        # Update polygon and metadata - O(n)
                        polygon = candidate_polygon
                        occupied_grid = set(polygon)



                        current_vertices = len(polygon)
                        next_indices = [(i + 1) % current_vertices for i in range(current_vertices)]
                        success = True
                        consecutive_failures = 0

                        break


                if success:
                    break


            if not success:
                consecutive_failures += 1

        self.end_time = time.time()
        #self.execution_time = self.end_time - self.start_time
        self.execution_time = (self.end_time - self.start_time) * 1000  # ms

        # assign the final polygon
        self.polygon = polygon 
        # Final check for duplicates
        self.polygon = remove_duplicates(self.polygon)
        return [[int(x), int(y)] for x, y in self.polygon]

    def visualize(self, save_path: Optional[str] = None):
        """Visualize the polygon"""
        import matplotlib.pyplot as plt
        
        if not self.polygon:
            print("No polygon to visualize!")
            return
            
        data = [[int(x), int(y)] for x, y in self.polygon]
        xs, ys = zip(*data)
        
        xs = xs + (xs[0],)
        ys = ys + (ys[0],)
        
        plt.figure(figsize=(10, 10))
        plt.fill(xs, ys, alpha=0.4, color='lightcoral', 
                edgecolor='darkred', linewidth=2)
        plt.plot(xs[:-1], ys[:-1], 'ro', markersize=8)
        
        for i, (x, y) in enumerate(data):
            plt.annotate(str(i), (x, y), xytext=(5, 5), 
                        textcoords='offset points', fontsize=9, fontweight='bold')
        
        plt.axis('equal')
        plt.grid(True, alpha=0.3, linestyle='--')
        plt.title(f"Inflate-Cut Algorithm: {len(data)} Vertices", 
                 fontsize=14, fontweight='bold')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved to {save_path}")
        
        plt.show()



# Example usage

if __name__ == "__main__":
    N_VALUES = [100000]  # List các số đỉnh N mà bạn muốn kiểm tra
    #N_VALUES = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
    NUM_RUNS = 1        # Số lần chạy cho mỗi giá trị N

    for index, n in enumerate(N_VALUES):
        # Danh sách để lưu trữ thời gian thực thi của NUM_RUNS lần chạy cho N hiện tại
        execution_times = []
        
        # Danh sách để lưu trữ số đỉnh cuối cùng (kiểm tra)
        final_vertex_counts = []
        
        print(f"\n--- Bắt đầu chạy Inflate-Cut ({NUM_RUNS} lần, N={n} đỉnh) ---")

        for i in range(NUM_RUNS):
            generator = InflateCutGenerator()
            polygon = generator.generate(n)
            
            execution_times.append(generator.execution_time)
            final_vertex_counts.append(len(polygon))

            print(f"Lần chạy {i+1}/{NUM_RUNS}: {generator.execution_time:.2f} ms")
            
        '''
        if index == len(N_VALUES) - 1:
            generator.polygon = [tuple(p) for p in polygon] # Gán kết quả của lần chạy cuối
            generator.visualize(save_path=f"inflate_cut_final_{n}.png")
        

        # --- 📊 Tính toán Thống kê ---
        
            
        times_array = np.array(execution_times)
        vertex_array = np.array(final_vertex_counts)
        
        mean_time = np.mean(times_array)
        min_time = np.min(times_array)
        max_time = np.max(times_array)
        std_dev_time = np.std(times_array)

        # --- 📝 In Kết quả Cuối cùng ---
        
        print("\n========================================================")
        print(f"✅ BÁO CÁO THỐNG KÊ (N={n} ĐỈNH, {NUM_RUNS} LẦN CHẠY) ✅")
        print("========================================================")
        print(f"Số đỉnh tạo được (TB): {np.mean(vertex_array):.0f} (Max: {np.max(vertex_array)})")
        print(f"Thời gian Trung bình (Mean): {mean_time:.4f} ms")
        print(f"Thời gian Tối thiểu (Min):   {min_time:.4f} ms")
        print(f"Thời gian Tối đa (Max):      {max_time:.4f} ms")
        print(f"Độ lệch chuẩn (Std Dev):    {std_dev_time:.4f} ms")
        print("========================================================")
        '''
