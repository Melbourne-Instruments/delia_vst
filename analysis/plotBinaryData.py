from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import data
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack

filestring = "./param_smooth_test2.bin"
array = np.fromfile(filestring, dtype="<f")
array_1_d = np.gradient(array)
fig, ax = plt.subplots()

ax.plot(array)
filestring = "./param_smooth_test.bin"
array = np.fromfile(filestring, dtype="<f")
array_2_d = np.gradient(array)
ax.plot(array)
max1  = np.max(array_1_d)
max2  = np.max(array_2_d)
max1 = max(max1,max2)
ax.plot(array_1_d/max1)
ax.plot(array_2_d/max1)
ax.legend(["1", "2", "1d", "2d"], prop={'size': 20})
plt.show()
#filestring = "./float_filter_dump_in.bin"