import numpy as np
from tissueviewer.tvtiff import tiffread, tiffsave
import automated_phenotyping as ap
from mayavi import mlab
from handy_functions import *
import os


# tic('t_tot')
# directory_name = '/home/maxbrambach/workspace/automated_phenotyping/henriks_files/plant1/processed_tiffs/acylYFP'
# # print directory_name+'/'+os.listdir(directory_name)[0]
# #os.listdir(directory_name)[0][:-4]
# filenames = os.listdir(directory_name)
# for i in range(len(filenames)):
#     tic()
#     print filenames[i]
#     print '------------'
#     A = ap.AutoPhenotype()
#     A.data,_ = tiffread(directory_name+'/'+filenames[i])
#     toc()
#     print 'contour fit'
#     A.contour_fit_threshold(.8, 3)
#     toc()
#     print 'mesh conversion'
#     A.mesh_conversion()
#     A.smooth_mesh(300, .5)
#     A.clean_mesh()
#     A.curvature_slice(0,'mean')
#     A.feature_extraction(20)
#     toc()
#     print 'sphere fitting'
#     A.sphere_fit()
#     A.sphere_evaluation()
#     toc()
#     print 'paraboloid fit'
#     A.paraboloid_fit_mersitem()
#     toc()
#     print 'save all'
#     A.save('/home/maxbrambach/workspace/automated_phenotyping/henriks_files/plant1/done/'+filenames[i][:-4])
#     toc()
#     print 'done'
# print 'total time'   
# toc('t_tot')

directory_name = '/home/maxbrambach/plants_henrik/acylYFP'
filenames = os.listdir(directory_name)
files = []
for i in range(len(filenames)):
    files.append(directory_name+'/'+filenames[i])


def henriks_process(file_loc):
    A = ap.AutoPhenotype()
    A.data,_ = tiffread(file_loc)
    A.contour_fit_threshold(.8, 3)
    A.mesh_conversion()
    A.smooth_mesh(300, .5)
    A.clean_mesh()
    A.curvature_slice(0,'mean')
    A.feature_extraction(20)
    A.sphere_fit()
    A.sphere_evaluation()
    A.paraboloid_fit_mersitem()
    A.save(file_loc[:-4])

from multiprocessing import Pool


p = Pool(1)
p.map(henriks_process, files)

directory_name = '/home/maxbrambach/workspace/automated_phenotyping/images_from_yi'
filenames = os.listdir(directory_name)
files = []
for i in range(len(filenames)):
    files.append(directory_name+'/'+filenames[i])
    
p = Pool(1)
p.map(henriks_process, files)