import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

c_w_arr_1 = np.arange(.1, .015, -.001)
c_h_arr_1 = np.zeros(len(c_w_arr_1))
for i in np.arange(len(c_w_arr_1)):
    c_w = c_w_arr_1[i]
    func = lambda c_h: (c_w + c_h)**0.2 * (c_w + 2*c_h) / (c_w * c_h) - 100
    c_h_final = fsolve(func, 0.01)
    c_h_arr_1[i] = c_h_final


c_h_arr_2 = np.arange(.015, .1, .001)
c_w_arr_2 = np.zeros(len(c_h_arr_2))
for i in np.arange(len(c_h_arr_2)):
    c_h = c_h_arr_2[i]
    func = lambda c_w: (c_w + c_h)**0.2 * (c_w + 2*c_h) / (c_w * c_h) - 100
    c_w_final = fsolve(func, .016)
    c_w_arr_2[i] = c_w_final


c_h_arr = np.concatenate([c_h_arr_1, c_h_arr_2])
c_w_arr = np.concatenate([c_w_arr_1, c_w_arr_2])
dP_arr = np.zeros(len(c_w_arr))
for i in np.arange(len(c_w_arr)):
    dP_arr[i] = (c_w_arr[i] + c_h_arr[i]) / ((c_w_arr[i]*c_h_arr[i])**3)



plt.figure()
plt.plot(c_w_arr, c_h_arr)
plt.xlabel('Channel Width')
plt.ylabel('Channel Height')
plt.title('Heat Flux Constraint')


color_map = np.arange(len(c_w_arr))
plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(c_w_arr, c_h_arr, dP_arr/1000000000, c = dP_arr)
ax.set_xlabel('Channel Width')
ax.set_ylabel('Channel Height')
ax.set_zlabel('dP (no units)')
plt.title('dP Minimization')
plt.show()


