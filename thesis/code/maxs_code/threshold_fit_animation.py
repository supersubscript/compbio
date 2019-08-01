'''
Created on 19 Jul 2017

@author: maxbrambach
'''
import numpy as np
import mayavi.mlab as mlab
import  moviepy.editor as mpy
import automated_phenotyping as ap
from tissueviewer.tvtiff import tiffread
from handy_functions import *


duration = 6 # duration of the animation in seconds (it will loop)

# MAKE A FIGURE WITH MAYAVI

fig_myv = mlab.figure(size=(512,512), bgcolor=(1,1,1))
# mlab.view(distance='auto')#azimuth=200.)


# INITIALISE DATA
data_large,_ = tiffread('/home/maxbrambach/workspace/automated_phenotyping/test_images/nt1.tif')
A = ap.AutoPhenotype(data = data_large)
A.reduce(2,spline=False)
A.contour_fit_threshold(.7, 3)

# view3d(A.contour, True)

# src = mlab.pipeline.scalar_field(A.data)
 
# ANIMATE THE FIGURE WITH MOVIEPY, WRITE AN ANIMATED GIF
 
def make_frame(t):
    mlab.clf() # clear the figure (to reset the colors)
    mlab.contour3d(A.contour, contours=[.5],transparent = False)
#     mlab.pipeline.image_plane_widget(src, plane_orientation='x_axes', colormap='gray',slice_index=0)
    mlab.view(azimuth=t/duration*360.,distance=1000)
    return mlab.screenshot(antialiased=True)
     
 
animation = mpy.VideoClip(make_frame, duration=duration)
animation.write_gif("test.gif", fps=30)
