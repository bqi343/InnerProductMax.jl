#%%
import time

import numpy as np
import scipy
from scipy.spatial import ConvexHull

print(scipy.__version__)


def gen_on_unit_sphere(N):
    points = np.random.randn(int(N), 2)
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    return points

def benchmark(N):
    for trial in range(1):
        points = gen_on_unit_sphere(N)
        start = time.time()
        hull = ConvexHull(points)
        print(vars(hull).keys())
        end = time.time()
        assert len(hull.vertices) == N
        print(hull.vertices)
        print(hull.simplices)
        print(hull.equations)
        # print(type(hull.simplices), hull.simplices.dtype)
        # print(end - start) # ~ 0.1s

benchmark(10)

# %%
