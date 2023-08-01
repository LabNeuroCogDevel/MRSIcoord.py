#!/usr/bin/env python3
import cv2  # resize: si integral to match rorig
import matplotlib.pyplot as plt  # for color scale
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg  # numpy plots
from matplotlib.figure import Figure
from PIL import Image, ImageTk

from siarray import Scout, SIArray, adjust
si_array_file = "/Volumes/Hera/Projects/7TBrainMech/subjs/11868_20211025/slice_PFC/MRSI_roi/raw/siarray.1.1"
pixdim = (216, 216, 99)

siarray = SIArray(si_array_file)
intgrl = siarray.integrateSI(0)
cm = plt.get_cmap("viridis")
res = cv2.resize(intgrl, pixdim[0:2], interpolation=cv2.INTER_NEAREST)
colored = cm(res / res.max()) * 255
im = Image.fromarray(colored[:, :, 0:3].astype(np.uint8))

# things to hold onto
#imgtk = ImageTk.PhotoImage(image=im)

bone = plt.get_cmap("bone").copy()
bone.set_bad(color="black")
plt.imshow(np.ma.masked_where(abs(intgrl) < 1,
                              adjust(intgrl, 0, 100, 5)),
           cmap=bone, vmax=30)
