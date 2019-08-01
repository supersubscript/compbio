'''
Created on 19 Jul 2017

@author: maxbrambach
'''
import numpy as np
import mayavi.mlab as mlab
import  moviepy.editor as mpy
import automated_phenotyping as ap
from tissueviewer.tvtiff import tiffread


duration = 6 # duration of the animation in seconds (it will loop)

# MAKE A FIGURE WITH MAYAVI




# INITIALISE DATA
data_large,_ = tiffread('/home/maxbrambach/workspace/automated_phenotyping/test_images/nt1.tif')
A = ap.AutoPhenotype(data = data_large)
A.reduce(2)
A.set_contour_to_box()
A.load('/home/maxbrambach/workspace/automated_phenotyping/anim_ACWE')
A.contour_fit(50, 1., 3)
# A.save('/home/maxbrambach/workspace/automated_phenotyping/anim_ACWE')
# src = mlab.pipeline.scalar_field(A.data)

# ANIMATE THE FIGURE WITH MOVIEPY, WRITE AN ANIMATED GIF

fig_myv = mlab.figure(size=(512,512), bgcolor=(1,1,1))
# mlab.view(azimuth=200.)

def make_frame(t):
    mlab.clf() # clear the figure (to reset the colors)
    mlab.contour3d(A.contour, contours=[.5],transparent = False)
#     mlab.pipeline.image_plane_widget(src, plane_orientation='x_axes', colormap='gray',slice_index=0)
    mlab.view(azimuth=t/duration*360.,distance=1000)
    return mlab.screenshot(antialiased=True)

animation = mpy.VideoClip(make_frame, duration=duration)
animation.write_gif("test.gif", fps=30)
