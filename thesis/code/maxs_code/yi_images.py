'''
Created on 17 Jul 2017

@author: maxbrambach
'''
import automated_phenotyping as ap
from tissuelab.tltiff import tiffread


A = ap.AutoPhenotype()
A.data, _ = tiffread('/home/maxbrambach/workspace/automated_phenotyping/images_from_yi/Col0 (0-5um)-1.lsm')

# A.reduce(2)

A.contour_fit_threshold(.8, 3)

print 'contour done'

A.mesh_conversion()
# A.clear('contour')

# A.show_mesh()

A.smooth_mesh(300, 0.3)
 
print 'mesh done'
 
A.curvature_slice(threshold=-0.0,curvature_type='mean')
A.show_mesh()
 
A.feature_extraction(100)
  
  
print 'features done'
  
# A.clear('mesh')
A.sphere_fit()
  
print 'fit done done'
  
# A.clear('features')
A.sphere_evaluation()
print A.results
print A.get_div_angle('sphere_R')
A.show_spheres(meristem_first=False)