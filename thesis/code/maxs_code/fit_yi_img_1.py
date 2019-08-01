'''
Created on 17 Jul 2017

@author: maxbrambach
'''
import automated_phenotyping as ap
from tissuelab.tltiff import tiffread
from handy_functions import *

A = ap.AutoPhenotype()
A.data, _ = tiffread('/home/maxbrambach/workspace/automated_phenotyping/images_from_yi/Col0 (0-5um)-1.lsm')
A.set_contour_to_box()
A.contour_fit(250, .8, 3)
A.save('/home/maxbrambach/workspace/automated_phenotyping/yi_contour/Col0')

A = ap.AutoPhenotype()
A.data, _ = tiffread('/home/maxbrambach/workspace/automated_phenotyping/images_from_yi/cu4 (0-5um)-4.lsm')
A.set_contour_to_box()
A.contour_fit(250, .8, 3)
A.save('/home/maxbrambach/workspace/automated_phenotyping/yi_contour/cu4')

A = ap.AutoPhenotype()
A.data, _ = tiffread('/home/maxbrambach/workspace/automated_phenotyping/images_from_yi/det3 (0-5um)-2.lsm')
A.set_contour_to_box()
A.contour_fit(250, .8, 3)
A.save('/home/maxbrambach/workspace/automated_phenotyping/yi_contour/det3')
