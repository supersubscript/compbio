from tissueviewer.tvtiff import tiffread
import pickle
import numpy as np
import time
import scipy.optimize as opt
import vtk

def tic(name='time1'):
    globals()[name] = time.time()
    
def toc(name='time1',print_it=True):
    total_time = time.time() - globals()[name]
    if print_it == True:
        if total_time < 0.001:
            print '--- ', round(total_time*1000.,2), 'ms', ' ---'
        elif total_time >= 0.001 and total_time < 60:
            print '--- ', round(total_time,2), 's', ' ---'
        elif total_time >= 60 and total_time/3600. < 1:
            print '--- ', round(total_time/60.,2), 'min', ' ---'
        else:
            print '--- ', round(total_time/3600.,2), 'h', ' ---'
    else:
        return total_time
    
def readImages(imageFileName):
    image, tags = tiffread(imageFileName)
    return image, tags

def circle_levelset(shape, center, sqradius):
    """Build a binary function with a circle as the 0.5-levelset."""
    grid = np.mgrid[list(map(slice, shape))].T - center
    phi = sqradius - np.sqrt(np.sum((grid.T)**2, 0))
    u = np.float_(phi > 0)
    return u

def fit_sphere(data,init=[0,0,0,10]):
    def fitfunc(p, coords):
        x0, y0, z0, _ = p
        x, y, z = coords.T
        return ((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
    errfunc = lambda p, x: fitfunc(p, x) - p[3]**2.
    p1, _ = opt.leastsq(errfunc, init, args=(np.array(np.nonzero(data)).T,))
    p1[3] = abs(p1[3])
    return p1

def view3d(data,contour=False):
    from mayavi import mlab
    data = np.array(data)
    if data.dtype == 'bool':
        data = np.array(data,dtype='int')
    mlab.gcf()
    mlab.clf()
    if contour == False:
        mlab.points3d(np.nonzero(data)[0],np.nonzero(data)[1],np.nonzero(data)[2],scale_factor=.5)
    else:
        mlab.contour3d(data, contours=[0.5])
    mlab.show()

def save_var(variables,path,confirm=False):
    with open(path,'w') as f:
        pickle.dump(variables, f)
    if confirm != False: 
        print 'all saved'
        
def load_var(path):
    with open(path) as f:
        return pickle.load(f)
    
def shake(array):
    msk = np.array(array)
    msk[1::,:,:] = msk[:-1:,:,:] + msk[1::,:,:]
    msk[:-1:,:,:] = msk[:-1:,:,:] + msk[1::,:,:]
    msk[:,1::,:] = msk[:,:-1:,:] + msk[:,1::,:]
    msk[:,:-1:,:] = msk[:,:-1:,:] + msk[:,1::,:]
    msk[:,:,1::] = msk[:,:,1:] + msk[:,:,:-1:]
    msk[:,:,:-1:] = msk[:,:,1:] + msk[:,:,:-1:]
    return np.array(msk,dtype='bool')


def sort_a_along_b(b,a):
    return np.array(sorted(zip(a,b)))[:,1]
# 
def view_polydata(poly):
    if np.shape(poly) == ():
        numel = 1
        poly = [poly]
    else:
        numel = np.shape(poly)[0]
    colors = [(1.,1.,1.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.),(1.,1.,0.),
              (1.,0.,1.),(0.,1.,1.),(1.,0.5,0.5),(0.5,1.,0.5),(0.5,0.5,1.),
              (1.,1.,0.5),(1.,0.5,1.),(0.5,1.,1.)] # this is not very nice 
    Mappers = []
    Actors = []
    render = vtk.vtkRenderer()
    for i in range(numel):
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(poly[i])
        mapper.ScalarVisibilityOff()
        mapper.Update()
        Mappers.append(mapper)
        actor = vtk.vtkActor()
        actor.SetMapper(Mappers[i])
        actor.GetProperty().SetColor(colors[i])
        Actors.append(actor)
        render.AddActor(Actors[i])
    renderwindow = vtk.vtkRenderWindow()
    renderwindow.AddRenderer(render)
    renderwindow.SetSize(600,600)
    interactrender = vtk.vtkRenderWindowInteractor()
    interactrender.SetRenderWindow(renderwindow)
    interactrender.Initialize()
    axes = vtk.vtkAxesActor()
    widget = vtk.vtkOrientationMarkerWidget()
    widget.SetOutlineColor( 0.9300, 0.5700, 0.1300 )
    widget.SetOrientationMarker( axes )
    widget.SetInteractor( interactrender )
    widget.SetViewport( 0.0, 0.0, 0.2, 0.2 )
    widget.SetEnabled( 1 )
    widget.InteractiveOn()
    render.ResetCamera()
    renderwindow.Render()
    interactrender.Start()

    
# def array_from_vtk_polydata(poly,size=[]):
#     if np.shape(size) == np.shape([]):
#         size = np.array(poly.GetPoints().GetBounds(),dtype='int')[1::2]
#     indices = np.array(vtk_to_numpy(poly.GetPoints().GetData()),dtype='int')
#     out = np.zeros(size)
#     out[indices[:,0]-1,indices[:,1]-1,indices[:,2]-1] = 1
#     return np.array(out)
 
 

# def vtk_polydata_from_array(array):
#     out = vtk.vtkPolyData()
#     longarray = (numpy_to_vtk(np.nonzero(array)))
#     out = vtk.vtkPointData()
#     out.SetInputData(longarray)
#     NumPy_data_shape = array.shape
#     VTK_data = numpy_to_vtk(num_array=array.ravel(), deep=True, array_type=vtk.VTK_POINT_DATA)

        
def spherefit_results(spheres):
    """
    gives several results from an array of spheres, such as distance between the first sphere (mersitem) and the other spheres (organs).
    Input:
        np.array[[x_center_meristem, y_center_meristem, z_center_meristem, radius_mersitem],
                [x_center_organ1, y_center_organ1, z_center_organ1, radius_organ1]
                ...]
    Output:
        np.array[[voulme_meristem, 0,0,0,0,0,0]
                [volume_organ1, location_organ1_realtive_to_meristem_x, y, z, r, theta, phi, projected_theta]
                ...]
    note: for spherical coordinates:  xyz -> yzx 
    angels in rad, distances in voxel
    """
    
    num_obj = np.shape(spheres)[0]
    out = np.zeros((num_obj,8))
    def sphere_voulume(radius):
        return  4./3.*np.pi*radius**3.
    
    out[:,0] = sphere_voulume(spheres[:,-1]) # voulumes
    out[1:,1] = spheres[1:,0]-spheres[0,0] #x relative to meristem
    out[1:,2] = spheres[1:,1]-spheres[0,1] #y
    out[1:,3] = spheres[1:,2]-spheres[0,2] #z
    out[1:,4] = np.sqrt(out[1:,1]**2. + out[1:,2]**2. + out[1:,3]**2.) #r
    out[1:,5] = np.arccos(out[1:,1]/out[1:,4]) #theta
    out[1:,6] = np.arctan(out[1:,3]/out[1:,2]) #phi
    out[1:,7] = np.arctan2(out[1:,2],out[1:,3])
    for i in range(1,num_obj):
        if out[i,7] < 0:
            out[i,7] = out[i,7] + 2.*np.pi
    
    return out        
        
        
        
        
        
        
        
        