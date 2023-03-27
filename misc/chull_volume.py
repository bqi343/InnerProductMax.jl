import numpy as np
from scipy.spatial import ConvexHull

points = np.array([[-1.0, 1.0, 1.0],
    [-1.0, -1.0, 1.0],
    [1.0, 0.0, -1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, -1.0]])

print(ConvexHull(points).volume) # 2.6666666666666665
print(ConvexHull(points).simplices + 1) 