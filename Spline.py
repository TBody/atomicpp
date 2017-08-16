from interpolate import RectBivariateSpline
import numpy as np

x = np.array([1,2,3,4,5])
y = np.array([1,2.5,3,4,5])

z = np.array([[1,2,3,4,5], [1,2,3,4,5], [2,2,3,4,5], [1,2,3,5,5],[1,2,3,4,5]])

interp = RectBivariateSpline(x, y, z)

print(interp(1,2))
