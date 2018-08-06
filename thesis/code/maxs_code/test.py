'''
Created on 6 Jul 2017

@author: maxbrambach
'''
import automated_phenotyping as ap
from tissuelab.tltiff import tiffread
from handy_functions import *
from scipy.ndimage import zoom

# tic()
# data_large,_ = tiffread('/home/maxbrambach/workspace/automated_phenotyping/test_images/nt1.tif')
# A = ap.AutoPhenotype(data = data_large)
# A.set_contour_to_box()
# # A.contour_fit_threshold(.8, 5)
# # A.contour_fit_two_stage(100, 1., 1, 30, 0.2, 3, 0.1)
# # # data_mid = zoom(data_large,0.5)
# # data_small = zoom(data_large,zoom=0.1)
# # A.data = data_small
# # A.set_contour_to_box()
# # A.contour_fit(100, 1., 1)
# # 
# # # view3d(A.contour, True)
# # 
# # 
# # # data_mid_shape = np.shape(data_mid)
# # # contour_mid = zoom(A.contour,zoom=5.05)
# # # contour_mid = contour_mid[0:data_mid_shape[0],0:data_mid_shape[1],0:data_mid_shape[2]]
# # # 
# # # A  = ap.AutoPhenotype(data=data_mid,contour=contour_mid)
# # # A.contour_fit(20, .2, 3)
# # # 
# # # view3d(A.contour, True)
# # 
# # data_large_shape = np.shape(data_large)
# # contour_large = zoom(A.contour,zoom=10.05)
# # contour_large = contour_large[0:data_large_shape[0],0:data_large_shape[1],0:data_large_shape[2]]
# # 
# # 
# # A  = ap.AutoPhenotype(data=data_large,contour=contour_large)
# # A.contour_fit(20, .2, 3)
# # 
# # # print test
# toc()
# 
# view3d(A.contour, True)


data_large,_ = tiffread('/home/maxbrambach/workspace/automated_phenotyping/test_images/nt1.tif')
A = ap.AutoPhenotype(data = data_large)
A.reduce(2,spline=False)
A.contour_fit_threshold(.7, 3)

A.mesh_conversion()
A.show_mesh()
A.clean_mesh()
A.show_mesh()