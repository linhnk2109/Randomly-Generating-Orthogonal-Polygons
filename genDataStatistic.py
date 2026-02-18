import random
import matplotlib.pyplot as plt
import time
import numpy as np
from shapely.geometry import Polygon
import tracemalloc
import sys


class OrthogonalPolygonGenerator:
    def __init__(self):
        self.d = {}  # Dictionary để đánh dấu tọa độ đã sử dụng
        self.data = [[0, 0], [1, 0], [1, 1], [0, 1]]  # Khởi tạo hình ban đầu
        self.a = [[1, 0], [-1, 0], [0, 1], [0, -1]]  # Các hướng di chuyển
        
        # Khởi tạo các metrics
        self.start_time = 0
        self.end_time = 0
        self.iteration_count = 0
        self.successful_additions = 0
        
        # Khởi tạo tọa độ ban đầu
        self.d[(0, 0)] = 1
        self.d[(1, 0)] = 1
        self.d[(1, 1)] = 1
        self.d[(0, 1)] = 1

    def solve1(self, x, y, x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i):
        f = [[x0, y0], [x2, y2], [x6, y6], [x5, y5], [x3, y3]]
        self.d[(x1, y1)] = 0; self.d[(x2, y2)] = 1; self.d[(x3, y3)] = 1
        self.d[(x4, y4)] = 0; self.d[(x5, y5)] = 1; self.d[(x6, y6)] = 1
        self.d[(x2, y2)] = 1; self.d[(x, y)] = 0
        self.data = self.data[:i] + f + self.data[i + 1:]
        self.successful_additions += 1

    def solve2(self, x, y, x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i):
        f = [[x3, y3], [x5, y5], [x6, y6], [x2, y2], [x0, y0]]
        self.d[(x1, y1)] = 0; self.d[(x2, y2)] = 1; self.d[(x3, y3)] = 1
        self.d[(x4, y4)] = 0; self.d[(x5, y5)] = 1; self.d[(x6, y6)] = 1
        self.d[(x2, y2)] = 1; self.d[(x, y)] = 0
        self.data = self.data[:i] + f + self.data[i + 1:]
        self.successful_additions += 1
        
    def solve3(self, x, y, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i):
        f = [[x6, y6], [x5, y5], [x2, y2], [x4, y4], [x3, y3]]
        self.d[(x1, y1)] = 0; self.d[(x2, y2)] = 1; self.d[(x3, y3)] = 1
        self.d[(x4, y4)] = 1; self.d[(x5, y5)] = 1; self.d[(x6, y6)] = 1
        self.d[(x, y)] = 0
        self.data = self.data[:i] + f + self.data[i + 1:]
        self.successful_additions += 1

    def solve4(self, x, y, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i):
        f = [[x3, y3], [x4, y4], [x2, y2], [x5, y5], [x6, y6]]
        self.d[(x1, y1)] = 0; self.d[(x2, y2)] = 1; self.d[(x3, y3)] = 1
        self.d[(x4, y4)] = 1; self.d[(x5, y5)] = 1; self.d[(x6, y6)] = 1
        self.d[(x, y)] = 0
        self.data = self.data[:i] + f + self.data[i + 1:]
        self.successful_additions += 1

    def generate(self, n):
        """Sinh đa giác với n đỉnh"""

        # Bắt đầu tracking memory
        tracemalloc.start()
                
        self.start_time = time.time()
        self.iteration_count = 0
        self.successful_additions = 0
        
        if n > 4:
            if n % 4 != 0:
                self.data = [[0, 0], [1, 0], [1, 1], [2, 1], [2, 2], [0, 2]]
                self.d[(2, 1)] = 1; self.d[(2, 2)] = 1; self.d[(1, 2)] = 0
                self.d[(0, 2)] = 1; self.d[(0, 1)] = 0
                
            while True:
                if len(self.data) == n:
                    break
                    
                self.iteration_count += 1
                m = random.randint(0, len(self.data) - 1)
                
                for i in range(m, len(self.data)):
                    if len(self.data) == n:
                        break
                        
                    x = self.data[i][0]
                    y = self.data[i][1]
                    
                    if ((x, y + 1) in self.d) and ((x, y - 1) in self.d) and \
                       ((x + 1, y) in self.d) and ((x - 1, y) in self.d):
                        continue
                        
                    j = random.randint(0, 3)
                    x0 = x + self.a[j][0]
                    y0 = y + self.a[j][1]
                    
                    if ((x0, y0 + 1) in self.d) and ((x0, y0 - 1) in self.d) and \
                       ((x0 + 1, y0) in self.d) and ((x0 - 1, y0) in self.d):
                        continue
                        
                    if (x0, y0) in self.d:
                        if len(self.data) == n:
                            break
                        self._process_vertex(x, y, x0, y0, j, i)
                        
        self.end_time = time.time()

        # Lấy peak memory
        current_mem, peak_mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        # Memory từ tracemalloc 
        self.peak_memory_kb = peak_mem / 1024  # KB
        

        return self.data

    def _process_vertex(self, x, y, x0, y0, j, i):
        '''
        if j == 0:
            x1 = x; y1 = y + 1; x2 = x0; y2 = y0 + 1
            x22 = x0 + 1; y22 = y0
            if ((x22, y22) in self.d) and (self.d[(x0, y0)] == 0):
                if ((x1, y1) not in self.d) and ((x2, y2) not in self.d):
                    x3 = x1; y3 = y1 + 1; x4 = x2; y4 = y2 + 1
                    x5 = x3 + 2; y5 = y3; x6 = x2 + 1; y6 = y2
                    if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and \
                       ((x5, y5) not in self.d) and ((x6, y6) not in self.d):
                        self.solve1(x, y, x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i)
        '''
        if j == 0: #Điểm (x0, y0)đang xét ở bên phải điểm (x, y) và xét cạnh nằm ngang của đa giác
            x1 = x; y1 = y + 1; x2 = x0; y2 = y0 + 1 #(x1, y1), (x2, y2) lần lượt là các điểm nằm bên trên (x, y) và (x0, x0)
            x22 = x0 + 1; y22 = y0 #(x22, y22) là các điểm bên phải (x0, x0)
            if ((x22, y22) in self.d) and (self.d[(x0, y0)] == 0): #Nếu điểm (x22, y22) chưa được sử dụng và (x0, y0) không phải là đỉnh của đa giác
                if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                    x3 = x1; y3 = y1 + 1; x4 = x2; y4 = y2 + 1 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên trên (x1, y1) và (x2, y2)
                    x5 = x3 + 2; y5 = y3; x6 = x2 + 1; y6 = y2 # (x5, y5), (x6, y6) là các điểm nằm bên phải (x3, y3) và (x2, y2)
                    if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                        self.solve1(x, y, x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
                else:
                    x1 = x; y1 = y - 1; x2 = x0; y2 = y0 - 1 #(x1, y1), (x2, y2) lần lượt là các điểm nằm bên self.dưới (x, y) và (x0, x0)
                    if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #Nếu các điểm chưa được sử dụng
                        x3 = x1; y3 = y1 - 1; x4 = x2; y4 = y2 - 1 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên self.dưới (x1, y1) và (x2, y2)
                        x5 = x3 + 2; y5 = y3; x6 = x2 + 1; y6 = y2 #(x5, y5), (x6, y6) là các điểm nằm bên phải (x3, y3) và (x2, y2)
                        if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                            self.solve2(x, y, x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
            else:
                if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                    x3 = x1; y3 = y1 + 1; x4 = x2; y4 = y2 + 1  # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên trên (x1, y1) và (x2, y2)
                    x5 = x1 + 2; y5 = y1; x6 = x0 + 1; y6 = y0 #(x5, y5), (x6, y6) là các điểm nằm bên phải (x1, y1) và (x0, y0)
                    if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                        self.solve3(x, y, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
                else:
                    x1 = x; y1 = y - 1; x2 = x0; y2 = y0 - 1 #(x1, y1), (x2, y2) lần lượt là các điểm nằm bên self.dưới (x, y) và (x0, x0)
                    if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                        x3 = x1; y3 = y1 - 1; x4 = x2; y4 = y2 - 1 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên self.dưới (x1, y1) và (x2, y2)
                        x5 = x1 + 2; y5 = y1; x6 = x0 + 1; y6 = y0  #(x5, y5), (x6, y6) là các điểm nằm bên phải (x1, y1) và (x0, y0)
                        if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                            self.solve4(x, y, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
        if j == 1: #Điểm (x0, y0) đang xét ở bên trái điểm (x, y) và xét cạnh nằm ngang của đa giác
            x1 = x; y1 = y + 1; x2 = x0; y2 = y0 + 1 #(x1, y1), (x2, y2) lần lượt là các điểm nằm bên trên (x, y) và (x0, x0)
            x22 = x0 - 1; y22 = y0 #(x22, y22) là điểm bên trái (x0, x0)
            if ((x22, y22) in self.d) and (self.d[(x0, y0)] == 0): #Nếu điểm (x22, y22) chưa được sử dụng và (x0, y0) không phải là đỉnh của đa giác
                if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                    x3 = x1; y3 = y1 + 1; x4 = x2; y4 = y2 + 1 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên trên (x1, y1) và (x2, y2)
                    x5 = x3 - 2; y5 = y3; x6 = x2 - 1; y6 = y2 #(x5, y5), (x6, y6) là các điểm nằm bên trái (x3, y3) và (x2, y2)
                    if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d):
                        self.solve2(x, y, x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
                else:
                    x1 = x; y1 = y - 1; x2 = x0; y2 = y0 - 1 #(x1, y1), (x2, y2) lần lượt là các điểm nằm bên self.dưới (x, y) và (x0, x0)
                    if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                        x3 = x1; y3 = y1 - 1; x4 = x2; y4 = y2 - 1 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên self.dưới (x1, y1) và (x2, y2)
                        x5 = x3 - 2; y5 = y3; x6 = x2 - 1; y6 = y2 #(x5, y5), (x6, y6) là các điểm nằm bên trái (x3, y3) và (x2, y2)
                        if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d):
                            self.solve1(x, y, x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
            else:
                if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                    x3 = x1; y3 = y1 + 1; x4 = x2; y4 = y2 + 1 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bêntrên (x1, y1) và (x2, y2)
                    x5 = x1 - 2; y5 = y1; x6 = x0 - 1; y6 = y0 #(x5, y5), (x6, y6) lần lượt là các điểm nằm bên trái (x1, y1) và (x0, y0)
                    if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d):
                        self.solve4(x, y, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
                else:
                    x1 = x; y1 = y - 1; x2 = x0; y2 = y0 - 1 #(x1, y1), (x2, y2) lần lượt là các điểm nằm bên self.dưới (x, y) và (x0, x0)
                    if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                        x3 = x1; y3 = y1 - 1; x4 = x2; y4 = y2 - 1 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên self.dưới (x1, y1) và (x2, y2)
                        x5 = x1 - 2; y5 = y1; x6 = x0 - 1; y6 = y0 #(x5, y5), (x6, y6) lần lượt là các điểm nằm bên trái (x1, y1) và (x0, y0)
                        if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                            self.solve3(x, y, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
                                
        if j == 2: #Điểm (x0, y0) đang xét ở trên bên trên điểm (x, y) và xét cạnh nằm self.dọc của đa giác
            x1 = x + 1; y1 = y; x2 = x0 + 1; y2 = y0 #(x1, y1), (x2, y2) lần lượt là các điểm nằm bên phải (x, y) và (x0, x0)
            x22 = x0; y22 = y0 + 1 #(x22, y22) là điểm bên trên (x0, x0)
            if ((x22, y22) in self.d) and (self.d[(x0, y0)] == 0): #Nếu điểm (x22, y22) chưa được sử dụng và (x0, y0) không phải là đỉnh của đa giác
                if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                    x3 = x1 + 1; y3 = y1; x4 = x2 + 1; y4 = y2 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên phải (x1, y1) và (x2, y2)
                    x5 = x3; y5 = y3 + 2; x6 = x2; y6 = y2 + 1 #(x5, y5), (x6, y6) lần lượt là các điểm nằm bên trên (x3, y3) và (x2, y2)
                    if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                        self.solve2(x, y, x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
                else:
                    x1 = x - 1; y1 = y; x2 = x0 - 1; y2 = y0 #(x1, y1), (x2, y2) lần lượt là các điểm nằm bên trái (x, y) và (x0, x0)
                    if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                        x3 = x1 - 1; y3 = y1; x4 = x2 - 1; y4 = y2 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên trái (x1, y1) và (x2, y2)
                        x5 = x3; y5 = y3 + 2; x6 = x2; y6 = y2 + 1 #(x5, y5), (x6, y6) lần lượt là các điểm nằm bên trên (x3, y3) và (x2, y2)
                        if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                            self.solve1(x, y, x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
            else:
                if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                    x3 = x1 + 1; y3 = y1; x4 = x2 + 1; y4 = y2 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên phải (x1, y1) và (x2, y2)
                    x5 = x1; y5 = y1 + 2; x6 = x0 ; y6 = y0 + 1 #(x5, y5), (x6, y6) lần lượt là các điểm nằm bên trên (x1, y1) và (x0, y0)
                    if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                        self.solve4(x, y, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
                else:
                    x1 = x - 1; y1 = y; x2 = x0 - 1; y2 = y0 #(x1, y1), (x2, y2) lần lượt là các điểm nằm bên trái (x, y) và (x0, x0)
                    if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                        x3 = x1 - 1; y3 = y1; x4 = x2 - 1; y4 = y2 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên trái (x1, y1) và (x2, y2)
                        x5 = x1; y5 = y1 + 2; x6 = x0; y6 = y0 + 1 #(x5, y5), (x6, y6) lần lượt là các điểm nằm bên trên (x1, y1) và (x0, y0)
                        if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                            self.solve3(x, y, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
        if j == 3: #Điểm (x0, y0) đang xét ở self.dưới bên self.dưới điểm (x, y) và xét cạnh nằm self.dọc của đa giác
            x1 = x + 1; y1 = y; x2 = x0 + 1; y2 = y0 #(x1, y1), (x2, y2) lần lượt là các điểm nằm bên phải (x, y) và (x0, x0)
            x22 = x0; y22 = y0 - 1 #(x22, y22) là điểm bên self.dưới và (x0, x0)
            if ((x22, y22) in self.d) and (self.d[(x0, y0)] == 0): #Nếu điểm (x22, y22) chưa được sử dụng và (x0, y0) không phải là đỉnh của đa giác
                if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                    x3 = x1 + 1; y3 = y1; x4 = x2 + 1; y4 = y2 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên phải (x1, y1) và (x2, y2)
                    x5 = x3; y5 = y3 - 2; x6 = x2; y6 = y2 - 1 #(x5, y5), (x6, y6) lần lượt là các điểm nằm bên self.dưới (x3, y3) và (x2, y2)
                    if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                        self.solve1(x, y, x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
                else:
                    x1 = x - 1; y1 = y; x2 = x0 - 1; y2 = y0 #(x1, y1), (x2, y2) lần lượt là các điểm nằm bên trái (x, y) và (x0, x0)
                    if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                        x3 = x1 - 1; y3 = y1; x4 = x2 - 1; y4 = y2 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên trái (x1, y1) và (x2, y2)
                        x5 = x3; y5 = y3 - 2; x6 = x2; y6 = y2 - 1 #(x5, y5), (x6, y6) lần lượt là các điểm nằm bên self.dưới (x3, y3) và (x2, y2)
                        if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                            self.solve2(x, y, x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
            else:
                if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                    x3 = x1 + 1; y3 = y1; x4 = x2 + 1; y4 = y2 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên phải (x1, y1) và (x2, y2)
                    x5 = x1; y5 = y1 - 2; x6 = x0 ; y6 = y0 - 1 #(x5, y5), (x6, y6) lần lượt là các điểm nằm bên self.dưới (x1, y1) và (x0, y0)
                    if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                        self.solve3(x, y, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác
                else:
                    x1 = x - 1; y1 = y; x2 = x0 - 1; y2 = y0 #(x1, y1), (x2, y2) lần lượt là các điểm nằm bên trái (x, y) và (x0, x0)
                    if ((x1, y1) not in self.d) and ((x2, y2) not in self.d): #(x1, y1), (x2, y2) chưa được sử dụng
                        x3 = x1 - 1; y3 = y1; x4 = x2 - 1; y4 = y2 # (x3, y3), (x4, y4) lần lượt là các điểm nằm bên trái (x1, y1) và (x2, y2)
                        x5 = x1; y5 = y1 - 2; x6 = x0; y6 = y0 - 1#(x5, y5), (x6, y6) lần lượt là các điểm nằm bên self.dưới (x1, y1) và (x0, y0)
                        if ((x3, y3) not in self.d) and ((x4, y4) not in self.d) and ((x5, y5) not in self.d) and ((x6, y6) not in self.d): #Nếu các điểm này chưa được sử dụng
                            self.solve4(x, y, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, i) #Tiến hành thêm bớt các đỉnh của đa giác



    def calculate_metrics(self):
        """Tính toán các thông số của đa giác"""
        if len(self.data) < 3:
            return None
            
        # Tạo Shapely polygon để tính các thông số hình học
        poly = Polygon(self.data)
        
        # 1. Thời gian thực thi
        execution_time = (self.end_time - self.start_time) * 1000  # ms
        # execution_time = (self.end_time - self.start_time)  # second
        
        # 2. Số đỉnh
        num_vertices = len(self.data)
        
        # 3. Diện tích
        area = poly.area
        
        # 4. Chu vi
        perimeter = poly.length
        
        # 5. aspect ratio
        xs = [p[0] for p in self.data]
        ys = [p[1] for p in self.data]
        bbox_width = max(xs) - min(xs)
        bbox_height = max(ys) - min(ys)
        #aspect_ratio = max(bbox_width, bbox_height) / min(bbox_width, bbox_height) if min(bbox_width, bbox_height) > 0 else 0
        
        aspect_ratio = poly.length / poly.area if poly.area > 0 else 0
        
        # 6. Convexity measure (Area / Convex Hull Area)
        convex_hull = poly.convex_hull
        convexity = area / convex_hull.area if convex_hull.area > 0 else 0
        
        # 7. Đếm reflex vertices (góc 270 độ)
        reflex_count = self._count_reflex_vertices()
        
        # 8. Đếm cạnh ngang và dọc
        horizontal_edges, vertical_edges = self._count_edge_orientations()
        
        # 9. Độ phức tạp (số lần lặp / số đỉnh)
        complexity_ratio = self.iteration_count / num_vertices if num_vertices > 0 else 0
        
        peak_memory_kb = self.peak_memory_kb
        per_vertex_memory_bytes = (peak_memory_kb * 1024) / num_vertices if num_vertices > 0 else 0
        
        return {
            "execution_time_ms": round(execution_time, 3),
            "num_vertices": num_vertices,
            "area": round(area, 2),
            "perimeter": round(perimeter, 2),
            "bbox_width": bbox_width,
            "bbox_height": bbox_height,
            "aspect_ratio": round(aspect_ratio, 2),
            "convexity_measure": round(convexity, 3),
            "reflex_vertices": reflex_count,
            "horizontal_edges": horizontal_edges,
            "vertical_edges": vertical_edges,
            "iteration_count": self.iteration_count,
            "successful_additions": self.successful_additions,
            "complexity_ratio": round(complexity_ratio, 2),
            "peak_memory_kb": round(peak_memory_kb, 2),
            "per_vertex_memory_bytes": round(per_vertex_memory_bytes, 2)
        }
        


    
    def _count_reflex_vertices(self):
        """Đếm số đỉnh lõm (reflex vertices - góc 270°)"""
        n = len(self.data)
        reflex_count = 0
        
        for i in range(n):
            prev = self.data[(i - 1) % n]
            curr = self.data[i]
            next_v = self.data[(i + 1) % n]
            
            # Tính cross product để xác định góc
            v1 = (curr[0] - prev[0], curr[1] - prev[1])
            v2 = (next_v[0] - curr[0], next_v[1] - curr[1])
            cross = v1[0] * v2[1] - v1[1] * v2[0]
            
            # Góc reflex khi cross product < 0 (góc lớn hơn 180°)
            if cross < 0:
                reflex_count += 1
                
        return reflex_count
    
    def _count_edge_orientations(self):
        """Đếm số cạnh ngang và dọc"""
        horizontal = 0
        vertical = 0
        
        n = len(self.data)
        for i in range(n):
            curr = self.data[i]
            next_v = self.data[(i + 1) % n]
            
            if curr[1] == next_v[1]:  # Cùng y -> cạnh ngang
                horizontal += 1
            elif curr[0] == next_v[0]:  # Cùng x -> cạnh dọc
                vertical += 1
                
        return horizontal, vertical
    
    def visualize(self, save_path=None):
        """Vẽ đa giác"""
        xs, ys = zip(*self.data)
        
        plt.figure(figsize=(10, 10))
        plt.fill(xs, ys, alpha=0.4, color='lightblue', edgecolor='blue', linewidth=2)
        plt.plot(xs, ys, 'ro', markersize=6)
        
        # Đánh số các đỉnh
        for i, (x, y) in enumerate(self.data):
            plt.annotate(str(i), (x, y), xytext=(5, 5), textcoords='offset points', fontsize=8)
        
        plt.axis('equal')
        plt.grid(True, alpha=0.3)
        plt.title(f"Orthogonal Polygon with {len(self.data)} vertices", fontsize=14, fontweight='bold')
        plt.xlabel("X coordinate", fontsize=12)
        plt.ylabel("Y coordinate", fontsize=12)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()


def run_experiment(vertex_counts, trials_per_size=10):
    all_results = []
    
    print("=" * 80)
    print("STARTING EXPERIMENTAL EVALUATION")
    print("=" * 80)
    
    for n in vertex_counts:
        print(f"\n{'='*60}")
        print(f"Testing with n = {n} vertices ({trials_per_size} trials)")
        print(f"{'='*60}")
        
        trial_results = []
        
        for trial in range(trials_per_size):
            generator = OrthogonalPolygonGenerator()
            generator.generate(n)
            metrics = generator.calculate_metrics()
            
            if metrics:
                trial_results.append(metrics)
                print(f"  Trial {trial+1}/{trials_per_size}: {metrics['execution_time_ms']:.2f} ms, "
                      f"Area: {metrics['area']:.1f}, Memory: {metrics['peak_memory_kb']:.2f} KB")
        
        if trial_results:
            stats = calculate_statistics(trial_results)
            stats['n'] = n
            all_results.append(stats)
            
            print(f"\n  Summary for n={n}:")
            print(f"    Mean time: {stats['mean_time']:.2f} ± {stats['std_time']:.2f} ms")
            print(f"    Mean area: {stats['mean_area']:.2f} ± {stats['std_area']:.2f}")
            print(f"    Mean memory: {stats['mean_peak_memory_kb']:.2f} KB")
    
    return all_results


def calculate_statistics(results):
    times = [r['execution_time_ms'] for r in results]
    areas = [r['area'] for r in results]
    perimeters = [r['perimeter'] for r in results]
    aspect_ratios = [r['aspect_ratio'] for r in results]
    convexities = [r['convexity_measure'] for r in results]
    reflex_counts = [r['reflex_vertices'] for r in results]
    peak_memories = [r['peak_memory_kb'] for r in results]
    per_vertex_memories = [r['per_vertex_memory_bytes'] for r in results]
    
    return {
        'mean_time': np.mean(times),
        'median_time': np.median(times),
        'std_time': np.std(times),
        'min_time': np.min(times),
        'max_time': np.max(times),
        'mean_area': np.mean(areas),
        'std_area': np.std(areas),
        'min_area': np.min(areas),
        'max_area': np.max(areas),
        'mean_perimeter': np.mean(perimeters),
        'mean_aspect_ratio': np.mean(aspect_ratios),
        'mean_convexity': np.mean(convexities),
        'mean_reflex': np.mean(reflex_counts),
        'mean_peak_memory_kb': np.mean(peak_memories),
        'mean_per_vertex_bytes': np.mean(per_vertex_memories)
    }


def print_summary_table(results, output_file=None):
    lines = []
    
    lines.append("\n" + "=" * 100)
    lines.append("PERFORMANCE STATISTICS TABLE")
    lines.append("=" * 100)
    lines.append(f"{'Vertices':<10} {'Mean (ms)':<12} {'Median (ms)':<13} {'Std Dev':<10} {'Min (ms)':<10} {'Max (ms)':<10}")
    lines.append("-" * 100)
    
    for r in results:
        lines.append(f"{r['n']:<10} {r['mean_time']:<12.2f} {r['median_time']:<13.2f} "
                    f"{r['std_time']:<10.2f} {r['min_time']:<10.2f} {r['max_time']:<10.2f}")
    
    lines.append("\n" + "=" * 100)
    lines.append("MEMORY USAGE STATISTICS TABLE")
    lines.append("=" * 100)
    lines.append(f"{'Vertices':<10} {'Peak Mem (KB)':<20} {'Per-Vertex (bytes)':<25}")
    lines.append("-" * 100)
    
    for r in results:
        lines.append(f"{r['n']:<10} {r['mean_peak_memory_kb']:<20.2f} {r['mean_per_vertex_bytes']:<25.2f}")
    
    lines.append("\n" + "=" * 100)
    lines.append("GEOMETRIC PROPERTIES TABLE")
    lines.append("=" * 100)
    lines.append(f"{'Vertices':<10} {'Mean Area':<12} {'Std Area':<12} {'Mean Aspect':<13} {'Mean Convex':<12} {'Mean Reflex':<12}")
    lines.append("-" * 100)
    
    for r in results:
        lines.append(f"{r['n']:<10} {r['mean_area']:<12.2f} {r['std_area']:<12.2f} "
                    f"{r['mean_aspect_ratio']:<13.2f} {r['mean_convexity']:<12.3f} {r['mean_reflex']:<12.1f}")
    
    for line in lines:
        print(line)
    
    if output_file:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        print(f"\n✓ Summary tables saved to '{output_file}'")


if __name__ == "__main__":
    VERTEX_COUNTS = [10, 20, 50, 100, 200, 500, 1000, 100000]
    #VERTEX_COUNTS = [1000000]
    TRIALS_PER_SIZE = 50
    
    results = run_experiment(VERTEX_COUNTS, TRIALS_PER_SIZE)
    
    print_summary_table(results, output_file='experiment_summary.txt')
    
    with open('experiment_results_detailed.txt', 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("DETAILED EXPERIMENTAL RESULTS\n")
        f.write("=" * 80 + "\n\n")
        
        for r in results:
            f.write(f"\n{'='*60}\n")
            f.write(f"RESULTS FOR n = {r['n']} VERTICES\n")
            f.write(f"{'='*60}\n\n")
            
            f.write("Performance Metrics:\n")
            f.write(f"  Mean execution time:     {r['mean_time']:.3f} ms\n")
            f.write(f"  Median execution time:   {r['median_time']:.3f} ms\n")
            f.write(f"  Std deviation:           {r['std_time']:.3f} ms\n")
            f.write(f"  Min time:                {r['min_time']:.3f} ms\n")
            f.write(f"  Max time:                {r['max_time']:.3f} ms\n\n")
            
            f.write("Memory Usage:\n")
            f.write(f"  Peak memory:             {r['mean_peak_memory_kb']:.2f} KB\n")
            f.write(f"  Per-vertex memory:       {r['mean_per_vertex_bytes']:.2f} bytes\n\n")
            
            f.write("Geometric Properties:\n")
            f.write(f"  Mean area:               {r['mean_area']:.2f} sq. units\n")
            f.write(f"  Std area:                {r['std_area']:.2f} sq. units\n")
            f.write(f"  Min area:                {r['min_area']:.2f} sq. units\n")
            f.write(f"  Max area:                {r['max_area']:.2f} sq. units\n")
            f.write(f"  Mean perimeter:          {r['mean_perimeter']:.2f} units\n")
            f.write(f"  Mean aspect ratio:       {r['mean_aspect_ratio']:.2f}\n")
            f.write(f"  Mean convexity:          {r['mean_convexity']:.3f}\n")
            f.write(f"  Mean reflex vertices:    {r['mean_reflex']:.1f}\n")
            f.write("\n")
    
    print("✓ Detailed results saved to 'experiment_results_detailed.txt'")
    
    
    print("\n" + "=" * 80)
    print("GENERATING SAMPLE POLYGON")
    print("=" * 80)
    
    generator = OrthogonalPolygonGenerator()
    generator.generate(1000)
    metrics = generator.calculate_metrics()
    
    print("\nSample Polygon Metrics:")
    for key, value in metrics.items():
        print(f"  {key}: {value}")
    
    generator.visualize(save_path='sample_polygon.png')
    print("\n✓ Sample polygon visualization saved to 'sample_polygon.png'")
    
