'''
Created on 6 Jul 2017

@author: Max Brambach
'''
from scipy.ndimage.morphology import binary_dilation, binary_erosion
from scipy.ndimage.interpolation import zoom
from itertools import cycle
import numpy as np
import scipy.optimize as opt
from tissueviewer.mesh import tvMeshImageVTK
import vtk
import pandas as pd
from vtk.util.numpy_support import vtk_to_numpy
import colorsys
import sys, os
from tissueviewer.tvtiff import tiffread, tiffsave

"""
===================
AutoPhenotype class
===================
"""

"""
Definition & Methds
===================
"""
class AutoPhenotype(object):
    """A python class for the automated extraction of SAM features in three 
    dimensions from a stack of confocal microscope images. 
    
    AutoPhenotype is a python class for the automated extraction of features of 
    a shoot apical meristem (SAM) from a stack of confocal microscope images.
    Its main features are:
    
    *Extraction of the SAMs surface using active contour without edges (ACWE).
    *Separation of the SAMs features into primordiae and the meristem.
    *Evaluation of the features in terms of location, size and divergence angle.
    
    [MORE INFO NEEDED]
    
    
    Attributes
    ----------
    data : 3D intensity array
        Image data e.g. from .tif file.
        
    contour : 3D boolean array
        Contour generated from data.
        
    mesh : vtkPolyData
        Mesh generated from contour. Can contain curvature values.
        
    features : list of vtkPolyData
        The extracted features from mesh.
    
    results : panda.DataFrame
        The results of the evaluation operations (e.g. sphere fit). Also
        contains the number of points contained in each feature.
        The number of rows is the number of primordea, whereas the  number of
        culumns depends on the number of evaluation operations performed. New
        results are appended.  
    """
    def __init__(self, data=None, tags=None, contour=None, mesh=None, features=None):
        """Builder"""
        self.data = data
        self.tags = tags
        self.contour = contour
        self.mesh = mesh
        self.features = features
        self.results = pd.DataFrame([],columns=['points_in_feature'])
        
    def contour_fit_threshold(self, threshold=.8, smooth_iterate=3):
        """Generate contour fit via a threshold. 
        
        Uses a threshold based approach to generate the initial model of 
        self.data. Threshold default is half of the images mean intensity; 
        can be adjusted.
        Subsequently, the contour is smoothed
        
        Parameters
        ----------
        threshold : float
            Sets the threshold for a voxel to be considered in- or outside the
            contour.
                *Inside: intensity(voxel) > mean_intensity/2*threshold
                *Outside: else
        
        iterate_smooth : int
            Number of times, the smoothing operation is performed.
            
        Return
        ------
        no return :
            Overwrites contour.
        """
        self.contour = np.array(self.data < 
                                threshold*0.5*np.mean(self.data.flatten()),
                                dtype = 'int')
        self.smooth_contour(smooth_iterate)
     
    def step(self,weighting):
        """Perform a single step of the morphological Chan-Vese evolution.
        
        Implementation of the first two parts of Equation 32 in Marquez-Neila
        et al. 2014. (DOI: 10.1109/TPAMI.2013.106). The third part (smoothing) 
        is implemented in the smooth_contour method.
        
        Parameters
        ----------
        weighting : float
            Ratio of the weightings of the inside to outside intensity. 
            (lambda_1/lambda_2 in Equation 32, where lambda_2 is chosen to be 1
        
        Return
        ------
        no return : 
            Overwrite self.contour.
        """
        if type(self.contour) == type(None):
            raise ValueError("the levelset function is not set (use .contour=)")
        inside = self.contour>0
        outside = self.contour<=0
        c0 = self.data[outside].sum() / float(outside.sum())
        c1 = self.data[inside].sum() / float(inside.sum())
        dres = np.array(np.gradient(self.contour))
        abs_dres = np.abs(dres).sum(0)
        aux = abs_dres * (weighting*(self.data - c1)**2 - 
                          1./float(weighting)*(self.data - c0)**2)
        self.contour[aux < 0] = 1
        self.contour[aux > 0] = 0
    
    def smooth_contour(self,iterate=1):
        """Smooth contour with ISoSI, SIoIS operators.
        
        Uses the SIoIS and ISoSI operator of Marquez-Neila et al. 2014 
        (Section 3.5) to smooth the contour. At least one iteration after step 
        is needed to get a good contour.
        """
        curvop = fcycle([SIoIS, ISoSI])
        res = self.contour
        for _i in range(iterate):
            res = curvop(res)
        self.contour = np.array(res,dtype='int')
        
    def reduce(self,factor=2,spline=False):
        """Reduce the size of contour and data by the specified factor using 
        slicing. 
        
        The method deletes the unused planes (keeps only every nth - n=factor).
        Can also use spline filters for up and downsampling.
        
        Parameter
        ---------
        factor : int if spline = False, float else
            Specify the factor of reduction.
            
        spline : boolean
            *if False (default): use numpy slicing to reduce number of planes.
            *if True: use spline interpolation for reducing / enlarging data and
              contour.
            
        Return
        ------
        no return :
            overwrite self.contour and self.data.
        """
        if spline == False:
            if type(self.data) != type(None):
                self.data = self.data[::factor,::factor,::factor]
            if type(self.contour) != type(None):
                self.contour = self.contour[::factor,::factor,::factor]
        if spline == True:
            if type(self.data) != type(None):
                self.data = zoom(self.data, zoom=1./float(factor))
            if type(self.contour) != type(None):
                self.contour = zoom(self.contour, zoom=1./float(factor))
                
    def contour_fit(self, iterations, weighting, iterate_smooth):
        """Run several full iterations of the morphological Chan-Vese method.
        
        Implementation of Equation 32 in Marquez-Neila et al. 2014. 
        (DOI: 10.1109/TPAMI.2013.106).
        
        Parameter
        ---------
        iterations : int
            Number of times, the Chan-Vese Method is iterated.
            
        weighting : float
            Ratio of the weightings of the outside and inside contour.
            
        iterate_smooth : int
            Number of times, the smoothing operation is performed during one
            iteration.
            
        Return
        ------
        no return : 
            overwrite contour.
        
        """
        for i in range(iterations):
            self.step(weighting)
            self.smooth_contour(iterate_smooth)
            print("Contour fit: iteration %s/%s..." % (i + 1, iterations))
            
    def contour_fit_two_stage(self, iterations1, weighting1, iterate_smooth1,
                              iterations2, weighting2, iterate_smooth2, 
                              zoom_factor):
        """Run a two staged contour fit. First fit is on down sampled image.
        
        The generated contour from the first fit is then up sampled and is used
        as initial contour for the regular fit.  

        Parameter
        ---------
        iterations1,2 : int
            Number of times, the Chan-Vese Method is iterated.
            1 ->
            
        weighting1,2 : float
            Ratio of the weightings of the outside and inside contour.
            
        iterate_smooth1,2 : int
            Number of times, the smoothing operation is performed during one
            iteration.
            
        Return
        ------
        no return : 
            overwrite contour.
        
        """
        data_large = self.data
        data_small = zoom(data_large,zoom=float(zoom_factor))
        self.data = data_small
        self.set_contour_to_box()
        self.contour_fit(iterations1, weighting1, iterate_smooth1)
        data_large_shape = np.shape(data_large)
        contour_large = zoom(self.contour, zoom = 1./float(zoom_factor) + .05)
        contour_large = contour_large[0:data_large_shape[0],
                                      0:data_large_shape[1],
                                      0:data_large_shape[2]]
        self.data = data_large
        self.contour = contour_large
        self.contour_fit(iterations2, weighting2, iterate_smooth2)
        
        

    def set_contour_to_box(self):
        """Set the contour attribute to a box. 
        
        The box consists of ones on the surfaces except the edges and the first
        plane on the first axis (axis = 0). The rest is zeros.
        
        Return
        ------
        no return : 
            Overwrites contour.
        """
        contour = getplanes(np.shape(self.data))
        contour = setedges(contour, 0)
        contour[0,:,:] = 0
        self.contour = contour
        
    def mesh_conversion(self):
        """Convert contour to mesh using marching cubes. 

        Uses the marching cubes algorithm from vtk and the tvMeshImageVTK 
        function is from tissuviewer. 
        
        Return
        ------
        no return : 
            Overwrites mesh.
        """
        fitval = 122
        fit = self.contour
        fit[fit==0] = fitval
        blockPrint()
        mesh = tvMeshImageVTK(fit, removeBoundaryCells = False, reduction = 0, 
                              smooth_steps = 0)
        enablePrint()
        self.mesh = mesh[fitval]

    def clean_mesh(self):
        """Extract largest connected set in self.mesh.
        
        Can be used to reduce residues of the contour fit inside the meristem.
        Works only if residues are not connected (share at least one point with)
        the meristems surface.
        
        Return
        ------
        no return : 
            Overwrites mesh.
        """
        connect = vtk.vtkConnectivityFilter()
        connect.SetInput(self.mesh)
        connect.SetExtractionModeToLargestRegion()
        connect.Update()
        geofilter = vtk.vtkGeometryFilter()
        geofilter.SetInput(connect.GetOutput())
        geofilter.Update()
        self.mesh = geofilter.GetOutput()

    def smooth_mesh(self, iterations=500, relaxation_factor=.5):
        """Smooth mesh.
        
        Uses vtk methods to clean and smooth the mesh. See documentation of 
        vtkSmoothPolyData and vtkCleanPolyData for more info.
        
        Parameters
        ----------
        iterations : int
            Number of iterations of the smooth algorithm.
            
        relaxation_factor: float
            Relaxation factor of the smooth algorithm.
        
        Return
        ------
        no return : 
            Overwrites mesh.
        """
        cleanPolyData = vtk.vtkCleanPolyData()
        cleanPolyData.SetInput(self.mesh)
        smoothFilter = vtk.vtkSmoothPolyDataFilter()
        smoothFilter.SetInputConnection(cleanPolyData.GetOutputPort())
        smoothFilter.SetNumberOfIterations(iterations)
        smoothFilter.SetRelaxationFactor(relaxation_factor)
        smoothFilter.Update()
        self.mesh = smoothFilter.GetOutput()
        
    def curvature_slice(self, threshold=0., curvature_type='mean'):
        """Slice the mesh along negative curvature.
        
        Computes the curvature of the mesh and then uses a threshold filter to
        remove the parts with negative curvature.
        
        Return
        ------
        no return : 
            Overwrites mesh.
        """
        curvature = vtk.vtkCurvatures()
        if curvature_type == 'max':
            curvature.SetCurvatureTypeToMaximum()
        elif curvature_type == 'mean':
            curvature.SetCurvatureTypeToMean()
        elif curvature_type == 'gauss':
            curvature.SetCurvatureTypeToGaussian()
        curvature.SetInput(self.mesh)
        curvature.Update()
        borders = vtk.vtkThreshold()
        if curvature_type == 'max' or 'mean':
            borders.ThresholdByLower(threshold)
        if curvature_type == 'gauss':
            borders.ThresholdByUpper(threshold)
        borders.SetInputConnection(curvature.GetOutputPort())
        geofilter = vtk.vtkGeometryFilter()
        geofilter.SetInputConnection(borders.GetOutputPort())
        geofilter.Update()
        self.mesh = geofilter.GetOutput()

    def feature_extraction(self,res=50):
        """Extract the SAM features from the sliced mesh.
        
        Uses a vtk connectivity filter for the feature extraction.
        The points of the mesh are numbered. The method steps through the points
        (stepsize = total#ofPoints/res) and selects the connected surface in 
        which the point lies. The selected surface is saved if:
            *Its #ofPoints is larger than stepsize
            *It has not already been selected 
                (#ofPoints != #ofPoints(previous_iterations) )
        
        Parameters
        ----------
        res : int
            'resolution' of the algorithm. Determines the step size:
            stepsize = total#ofPoints/res
        
        Return
        ------
        no return : 
            Overwrites features.
        """
        connect = vtk.vtkPolyDataConnectivityFilter()
        connect.SetExtractionModeToPointSeededRegions()
        connect.SetInput(self.mesh)
        connect.Update()
        points = connect.GetInput().GetNumberOfPoints()
        thresh = int(points/float(res))
        objects_vtk = []
        lastobj = []
        for i in range(2,res):
            connect.InitializeSeedList()
            connect.AddSeed(int(points/res*float(i))-1)
            connect.Update()
            temp = []
            temp = connect.GetOutput()
            temp.Update()
            if temp.GetNumberOfPoints() > thresh:
                if not(temp.GetNumberOfPoints() in lastobj):
                    objects_vtk.append(vtk.vtkPolyData())
                    objects_vtk[-1].DeepCopy(temp)
                    lastobj.append(temp.GetNumberOfPoints())
        self.features = objects_vtk
        npoints = pd.DataFrame(lastobj, columns=['points_in_feature'])
        self.results = self.results.append(npoints)
                
    def sphere_fit(self):
        """Fit spheres onto the features.
        
        Iterates over the vtkPolyData objects in features and performs a least
        square fit of a sphere to each feature. 
        
        Return
        ------
        self.results : pandas.DataFrame
            Four additional rows in self.results:
                *sphere_x_abs: absolute x-coordinate of the center of the 
                    fitted sphere
                *sphere_y_abs: absolute y-coordinate of the center of the 
                    fitted sphere
                *sphere_z_abs: absolute z-coordinate of the center of the 
                    fitted sphere
            Note: absolute means in units of the coordinate system used in 
            self.features.
        """
        out = []
        for i in range(np.shape(self.features)[0]):
            out.append(fit_sphere(array_from_vtk_polydata(self.features[i])))
        fitval = pd.DataFrame(np.array(out),columns=['sphere_x_abs',
                                                     'sphere_y_abs',
                                                     'sphere_z_abs',
                                                     'sphere_radius',
                                                     'sphere_res_var'])
        self.results = pd.concat([self.results,fitval],axis=1)
    
    def sort_results(self,column='index',ascending=False,reset_index=False):
        """Sort results by specified column.
        
        The sorting direction can be adjusted.
        
        Parameters
        ----------
        column : str
            Name of the column by which the array should be sorted.
            If column = 'index', the array will be sorted by its index
        
        ascending : bool
            Specifies the sorting direction:
                *False: high->low
                *True: low->high
        
        reset_index : bool
            Reset the row index after sorting the array. If True, the order of
            self.features is also changed accordingly.
        
        Return
        ------
        no return : 
            Overwrites results
        """
        if column == 'index':
            self.results.sort_index(ascending=True, inplace=True)
        else:
            self.results.sort_values(column, ascending=ascending, inplace=True)
        if reset_index == True:
            self.features = list(np.array(self.features)
                                 [self.results.index.values.tolist()])
            self.results.reset_index(inplace=True,drop=True)
        
    def sphere_evaluation(self):
        """Add several results to self.results.
        
        Added results are:
            *'sphere_x_rel': x-location of the primordium relative to the 
                meristem.
            *'sphere_y_rel': y-location of the primordium relative to the 
                meristem.
            *'sphere_z_rel': z-location of the primordium relative to the 
                meristem.
            *'sphere_volume': Volume of the sphere
            *'sphere_R': Distance of the primordium to the meristem.
            *'sphere_angle_raw': Angle between the primordia in the y,z plane. 
                Zero is chosen to be z = 0 and y > 0. 
            *--results missing--
        Note: List of results is newly sorted in a way that the row with index 
        0 is the meristem. The row label is changed accordingly.
        Return
        ------
        no return : 
            Results are updated (see description).
        """
        if 'sphere_radius' not in self.results.columns:
            raise ValueError('Perform sphere fit first. Use self.sphere_fit()')
        self.sort_results('sphere_radius', reset_index=True)
        num_obj = self.results.shape[0]
        out = np.zeros((num_obj,6))
        out[:,0] = sphere_voulume(self.results['sphere_radius']) # volumes
        out[1:,1] = self.results['sphere_x_abs'][1:]-self.results[
            'sphere_x_abs'][0] #x relative to meristem
        out[1:,2] = self.results['sphere_y_abs'][1:]-self.results[
            'sphere_y_abs'][0] #y relative to meristem
        out[1:,3] = self.results['sphere_z_abs'][1:]-self.results[
            'sphere_z_abs'][0] #z relative to meristem
        out[1:,4] = np.sqrt(out[1:,1]**2. + out[1:,2]**2. + out[1:,3]**2.) #R
        out[1:,5] = np.arctan2(out[1:,2],out[1:,3]) #theta'
        for i in range(1,num_obj):
            if out[i,5] < 0:
                out[i,5] = out[i,5] + 2.*np.pi
        out[1:,5] = out[1:,5]/2./np.pi*360.
        out_pd = pd.DataFrame(np.array(out),columns=['sphere_volume',
                                                     'sphere_x_rel',
                                                     'sphere_y_rel',
                                                     'sphere_z_resl',
                                                     'sphere_R',
                                                     'sphere_angle_raw'])
        self.results = pd.concat([self.results,out_pd],axis=1)

    def paraboloid_fit_mersitem(self):
        """Fit the mersitem with a paraboloid.
        
        Uses the first entry of self.features as meristem.
        See fit_paraboloid() function for more info.
        
        Return
        ------
        no return :
            New results are added to self.results
            *para_p1 ... para_p5: Parameters of the paraboloid fit.
            *para_alpha,beta,gamma: Rotation of the paraboloid relative to
              image.
            *para_apex_x,y,z: Location of the paraboloids apex in coordinates
              of the image.
        """
        out = np.zeros([self.results.shape[0],11])
        indices = np.array(np.nonzero(
            array_from_vtk_polydata(self.features[0]))).T
        popt = fit_paraboloid(indices)
        apex = get_paraboloid_apex(popt)
        out[0,:] = np.array(list(popt)+list(apex))
        fitval = pd.DataFrame(np.array(out),columns=['para_p1',
                                                     'para_p2',
                                                     'para_p3',
                                                     'para_p4',
                                                     'para_p5',
                                                     'para_alpha',
                                                     'para_beta',
                                                     'para_gamma',
                                                     'para_apex_x',
                                                     'para_apex_y',
                                                     'para_apex_z'])
        self.results = pd.concat([self.results,fitval],axis=1)
     
        
    def save(self, where):
        """Saves the all data stored in an AutoPhenotype Object. 
        
        Recognises whether or not data is available and saves only available. 
        Note: it is advised to create a new folder for every save, since 
        multiple files are created which always have the same name.
        The following files are saved (if avaliable):
            *data.tif : processed input data (e.g. reduced) as .tif stack
            *contour.tif : contour of data as .tif stack
            *mesh.vtp : mesh as vtk .vtp data (readable with 
                vtk.vtkXMLPolyDataReader() )
            *featuresX.vtp : features in same format as mesh.vtp. Each feature
                is a separate file and X is a running number starting from 0.
            *results.csv : the results from e.g. the spherical fit as .csv 
                data generated with pandas DataFrame.to_csv().
        
        Parameters
        ----------
        where : string, path
            Path to the folder in which the data should be saved. If folder 
            does not exist, it is created.
        """
        logfile = pd.DataFrame(np.zeros((1,6)),columns=['data',
                                                    'contour',
                                                    'mesh',
                                                    'features',
                                                    'results',
                                                    'tags'])
        if not os.path.exists(where):  # checks if specified directory exists
            os.makedirs(where)         # creates one if not
        if type(self.data) != type(None):
            logfile['data'] = 1
            tiffsave(np.array(self.data,'uint16'),where+'/data.tif')
        if type(self.contour) != type(None):
            logfile['contour'] = 1
            tiffsave(np.array(self.contour,'uint16'),where+'/contour.tif')
        if type(self.mesh) != type(None):
            logfile['mesh'] = 1
            meshwriter = vtk.vtkXMLPolyDataWriter()
            meshwriter.SetInput(self.mesh)
            meshwriter.SetFileName(where+'/mesh.vtp')
            meshwriter.Write()
        if type(self.features) != type(None):
            os.makedirs(where+'/features')
            logfile['features'] = 1
            for i in range(np.shape(self.features)[0]):
                featurewriter = vtk.vtkXMLPolyDataWriter()    
                featurewriter.SetInput(self.features[i])
                featurewriter.SetFileName(where+'/features/feature%s.vtp' 
                                          % str(i))
                featurewriter.Write()
        if self.results.shape[0] != 0:
            logfile['results'] = 1  
            self.results.to_csv(where+'/results.csv')
        logfile.to_csv(where+'/logfile.csv')
            
    def load(self,where,what=None):
        """Load data to existing AutoPhenotype Object.
    
        What is loaded can be specified.
        
        Parameters
        ----------
        where : string, path
            Path to the folder from which the data should be loaded from.
        
        what : either None, logfile or list of strings
            *None: Logfile from specified folder is used.
            *logfile: No logfile is imported, instead the specified logfile is
                used. E.g. generated by create_logfile() function
            *list of strings: List of keywords specifying which data should be 
                loaded. See logfile_from_str() for more information.
        Return
        ------
        no return :
            Attributes are overwritten.
        """
        if type(what) == type(None):
            logfile = pd.read_csv(where+'/logfile.csv')
        elif type(what) == type(create_logfile([1,1,1,1,1,1])):
            logfile = what
        elif type(what) == type(['1','2']):
            logfile = logfile_from_str(what)
        else:
            print ' - parameter what is unknown - '
            print 'If logfile.csv exists in directory (where) use what=None'
            print 'Else: either create logfile with create_logfile() function'
            print 'or specify files to load with strings:'
            print 'all, data, contour, mesh, features, results, tags'
            raise ValueError('parameter what is unknown')
        if logfile['data'][0] != 0:
            self.data,_ = tiffread(where+'/data.tif')
        if logfile['contour'][0] != 0:
            self.contour,_ = tiffread(where+'/contour.tif')
        if logfile['mesh'][0] != 0:
            meshreader = vtk.vtkXMLPolyDataReader()
            meshreader.SetFileName(where+'/mesh.vtp')
            meshreader.Update()
            self.mesh = meshreader.GetOutput()
        if logfile['features'][0] != 0:
            self.features = []
            number_of_features = len(next(os.walk(where+'/features'))[2])
            for i in range(number_of_features):
                featurereader = vtk.vtkXMLPolyDataReader()
                featurereader.SetFileName(where+'/features/feature%s.vtp' 
                                          % str(i))
                featurereader.Update()
                self.features.append(vtk.vtkPolyData())
                self.features[-1].DeepCopy(featurereader.GetOutput())
        if logfile['results'][0] != 0:
            self.results = pd.read_csv(where+'/results.csv')   

    def reset_results(self,keep=['points_in_feature']):
        """Reset the result attribute.
    
        Results to be kept can be specified. Keeps points_in_feature by default.
        If all results should be deleted use keep = []
        
        Parameters
        ----------
        keep : list of strings
            List specifying the results to be kept. 
            E.g. ['points_in_feature', 'sphere_radius']
            Note: Has to be list even if it only has one entry.
            Use [] for no results to be kept. 
        
        Return
        ------
        no return :
            Overwrites self.results.
        """
        if type(keep) == type('string'):
            keep = [keep]
        self.results = self.results[keep]
        
    def clear(self, what):
        """Clear specified attributes.
    
        Attributes to be specified can be:
            *'all'
            *'data'
            *'tags'
            *'contour'
            *'mesh'
            *'features'
            *'results'
        Specified attributes are reset into initial condition.
                
        Parameters
        ----------
        what : list of strings
            List of attributes to be cleared. See description.
            
        Return
        ------
        no return :
            Clears attributes.
        """
        if type(what) == type('string'):
            what = [what]
        if any(t == 'all' for t in what):
            self.data = None
            self.tags = None
            self.contour = None
            self.mesh = None
            self.features = None
            self.results = pd.DataFrame([],columns=['points_in_feature'])
        if any(t == 'data' for t in what):
            self.data = None
        if any(t == 'tags' for t in what):
            self.tags = None
        if any(t == 'contour' for t in what):
            self.contour = None
        if any(t == 'mesh' for t in what):
            self.mesh = None
        if any(t == 'features' for t in what):
            self.features = None
        if any(t == 'results' for t in what):
            self.results = pd.DataFrame([],columns=['points_in_feature'])
        
    def get_div_angle(self,sort_by='sphere_radius', sort_results=False):
        """Return divergence angles for angles sorted by sort_by.
    
        Can only be used after self.sphere_evaluation() has been used.
        Note: all angles must be in degree. Output is also in degree.
        
        Parameters
        ----------
        sort_by : str, index in self.results
            specify the order in which the organs are. 
            E.g. sort by 'sphere_radius'
            Default is 'sphere_radius'.
            
        sort_results : bool
            Specify, if self.results should be sorted accordingly.
            Default is no sorting.
                
        Return
        ------
        clockwise : float list
            Difference between angles of the primordiae, 
            clockwise orientation.
            
        counterclockwise : float list
            Difference between angles of the primordiae, 
            counterclockwise orientation.
        """
        results = self.results[['sphere_angle_raw',sort_by]]
        results = results.drop('meristem')
        results.sort_values(sort_by, ascending=False, inplace=True)
        if sort_results == True:
            self.sort_results(sort_by)
        return angle_difference(results['sphere_angle_raw'])

    def show_spheres(self,meristem_first=False):
        """"3D visualisation of the sphere fit.
        
        Uses vtk to show the fitted spheres. 
        Color coding:
            *White: first entry in self.results
            *R->G->B->P: following results
        Note: This makes the script pause at the position of the call. Closing
        the render window lets the script continue. 
        
        Return
        ------
        no return : 
            Opens a render window.
        """
        if meristem_first == True:
            firstcolor = (1.,1.,1.)
            lastcolor = ()
        if meristem_first == False:
            lastcolor = (1.,1.,1.)
            firstcolor = ()
        spheres = self.results[['sphere_x_abs','sphere_y_abs','sphere_z_abs',
                                'sphere_radius']].as_matrix()
        sphereResolution = 50
        spheresSources = []
        for i in range(np.shape(spheres)[0]):
            spherevtk = vtk.vtkSphereSource()
            spherevtk.SetCenter(spheres[i,0], spheres[i,1], spheres[i,2])
            spherevtk.SetRadius(spheres[i,3])
            spherevtk.SetThetaResolution(sphereResolution)
            spherevtk.SetPhiResolution(sphereResolution)
            spherevtk.Update()
            spheresSources.append(spherevtk.GetOutput())
        view_polydata(spheresSources,firstcolor,lastcolor)
        
    def show_spheres_and_features(self):
        """"//3D visualisation of the sphere fit.
        
        Uses vtk to show the fitted spheres. 
        Color coding:
            *White: first entry in self.results
            *R->G->B->P: following results
        Note: This makes the script pause at the position of the call. Closing
        the render window lets the script continue. 
        
        Return
        ------
        no return : 
            Opens a render window.
        """
        features = self.features
        numel = len(features)
        spheres = self.results[['sphere_x_abs','sphere_y_abs','sphere_z_abs',
                                'sphere_radius']].as_matrix()
        sphereResolution = 50
        spheresSources = []
        for i in range(np.shape(spheres)[0]):
            spherevtk = vtk.vtkSphereSource()
            spherevtk.SetCenter(spheres[i,0], spheres[i,1], spheres[i,2])
            spherevtk.SetRadius(spheres[i,3])
            spherevtk.SetThetaResolution(sphereResolution)
            spherevtk.SetPhiResolution(sphereResolution)
            spherevtk.Update()
            spheresSources.append(spherevtk.GetOutput())
        sphereMappers = []
        featureMappers = []
        Actors = []
        render = vtk.vtkRenderer()
        s_colors = rgb_list(numel, (1.,1.,1.))
        f_colors = rgb_list(numel, (.7,.7,.7),v=.7)
        for i in range(numel):
            s_mapper = vtk.vtkPolyDataMapper()
            s_mapper.SetInput(spheresSources[i])
            s_mapper.ScalarVisibilityOff()
            s_mapper.Update()
            sphereMappers.append(s_mapper)
            f_mapper = vtk.vtkPolyDataMapper()
            f_mapper.SetInput(features[i])
            f_mapper.ScalarVisibilityOff()
            f_mapper.Update()
            featureMappers.append(f_mapper)
            s_actor = vtk.vtkActor()
            s_actor.SetMapper(sphereMappers[i])
            s_actor.GetProperty().SetColor(s_colors[i])
            Actors.append(s_actor)
            f_actor = vtk.vtkActor()
            f_actor.SetMapper(featureMappers[i])
            f_actor.GetProperty().SetColor(f_colors[i])
            Actors.append(f_actor)
            render.AddActor(Actors[-1])
            render.AddActor(Actors[-2])
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
    
    def show_mesh(self):
        """"Visualise the mesh.
        
        Uses vtk to visualise the mesh.
        This can be done before or after the segmentation (using 
        self.curvatures_slice() ).
        Note: This makes the script pause at the position of the call. Closing
        the render window lets the script continue. 

        
        Return
        ------
        no return : 
            Opens a render window.
        """
        view_polydata(self.mesh,(1.,1.,1.),())
        
    def show_features(self):
        """"___.
        
        Uses vtk to visualise the mesh.
        This can be done before or after the segmentation (using 
        self.curvatures_slice() ).
        Note: This makes the script pause at the position of the call. Closing
        the render window lets the script continue. 

        
        Return
        ------
        no return : 
            Opens a render window.
        """
        view_polydata(self.features,(1.,1.,1.),())
###############################################################################
"""
Used functions
==============
"""
def setedges(array,value=0):
    """Set the edges of a 3D array to a specified value.
    
    Parameter
    ---------
    array : numpy array with shape() = (x,y,z)
        Array which edges are to be set to value.
        
    value : int, float
        Desired value for the edges of the array.
        
    Return:
    -------
    array : numpy array with shape() = (x,y,z)
        Input array with edges set to value.
    
    """
    array[[0,0,-1,-1],[0,-1,0,-1],:] = value
    array[:,[0,0,-1,-1],[0,-1,0,-1]] = value
    array[[0,-1,0,-1],:,[0,0,-1,-1]] = value
    return array

def getedges(shape):
    """Create a 3D numpy array with zeros on the edges and ones else.
    
    Used to suppress the fitting of edges by the smooth() function 
    (ISoSI operator mask).
    
    Parameter
    ---------
    shape : tuple with three components
        Shape of the returned array
        
    Return:
    -------
    array : numpy array with shape() = shape
        Three dimensional numpy array with zeros on the edges and ones else.    
    """
    array = np.ones(shape)
    setedges(array)
    return array

def setplanes(array,value=0):
    """Set the surface of a 3D array to a specified value.
    
    Parameter
    ---------
    array : numpy array with shape() = (x,y,z)
        Array which surface are to be set to value.
        
    value : int, float
        Desired value for the surface of the array.
        
    Return:
    -------
    array : numpy array with shape() = (x,y,z)
        Input array with surface set to value.
    
    """
    array[[0,-1],:,:] = value
    array[:,[0,-1],:] = value
    array[:,:,[0,-1]] = value
    
def getplanes(shape,invert=True):
    """Create a 3D numpy array with zeros on the surface and ones else or the
    other way around.
    
    Used as initial contour for the AutoPhenotype.step() method. 
    
    Parameter
    ---------
    shape : tuple with three components
        Shape of the returned array
    
    invert : bool
        *True: Ones on surface, zeros else
        *False: Zeros on surface, ones else 
        
    Return:
    -------
    array : numpy array with shape() = shape
        Three dimensional numpy array with zeros on the surface and ones else 
        or the other way around.    
    """
    if invert == True:
        array = np.zeros(shape)
        setplanes(array,value=1)
    if invert == False:
        array = np.ones(shape)
        setplanes(array,value=0)
    return array

class fcycle(object):
    """Call functions from the iterable each time it is called."""
    def __init__(self, iterable):
        self.funcs = cycle(iterable)
    def __call__(self, *args, **kwargs):
        f = next(self.funcs)
        return f(*args, **kwargs)   

def SI(u,iterate=1):
    """SI operator.
    
    Marquez-Neila et al. 2014. (DOI: 10.1109/TPAMI.2013.106)
    """
    P = [np.zeros((3,3,3)) for i in range(9)]
    P[0][:,:,1] = 1
    P[1][:,1,:] = 1
    P[2][1,:,:] = 1
    P[3][:,[0,1,2],[0,1,2]] = 1
    P[4][:,[0,1,2],[2,1,0]] = 1
    P[5][[0,1,2],:,[0,1,2]] = 1
    P[6][[0,1,2],:,[2,1,0]] = 1
    P[7][[0,1,2],[0,1,2],:] = 1
    P[8][[0,1,2],[2,1,0],:] = 1
    _aux = np.zeros((0))
    if u.shape != _aux.shape[1:]:
        _aux = np.zeros((len(P),) + u.shape)    
    for i in range(len(P)):
        _aux[i] = binary_erosion(u, P[i],iterations=iterate,
                                 mask=getedges(np.shape(u)))    
    return _aux.max(0)

def IS(u,iterate=1):
    """IS operator.
    
    Marquez-Neila et al. 2014. (DOI: 10.1109/TPAMI.2013.106)
    """
    P = [np.zeros((3,3,3)) for i in range(9)]
    P[0][:,:,1] = 1
    P[1][:,1,:] = 1
    P[2][1,:,:] = 1
    P[3][:,[0,1,2],[0,1,2]] = 1
    P[4][:,[0,1,2],[2,1,0]] = 1
    P[5][[0,1,2],:,[0,1,2]] = 1
    P[6][[0,1,2],:,[2,1,0]] = 1
    P[7][[0,1,2],[0,1,2],:] = 1
    P[8][[0,1,2],[2,1,0],:] = 1
    _aux = np.zeros((0))
    if u.shape != _aux.shape[1:]:
        _aux = np.zeros((len(P),) + u.shape)
    
    for i in range(len(P)):
        _aux[i] = binary_dilation(u, P[i],iterations=iterate,
                                  mask=getedges(np.shape(u)))
    return _aux.min(0)

def SIoIS(u):
    """SIoIS operator.
    
    Marquez-Neila et al. 2014. (DOI: 10.1109/TPAMI.2013.106)
    """
    return SI(IS(u))

def ISoSI(u):
    """ISoSI operator.
    
    Marquez-Neila et al. 2014. (DOI: 10.1109/TPAMI.2013.106)
    """
    return IS(SI(u))

def sort_a_along_b(b,a):
    """Return list 'a' sorted following the sorting of list 'b'. 
    
    List 'b' is sorted from low to high values. Elements in 'a' follow the
    sorting in list 'b'.
    Example: 
        a = [p,c,e,s,e,i,l]; b = [5,2,7,6,1,4,3]
        -> sort_a_along_b(a,b) = [e,c,l,i,p,s,e]; 
        # and b would be [1,2,3,4,5,6,7]
    
    Parameters
    ----------
    a : list (same length as b)
        List to be sorted.
        
    b : list (same length as a)
        Reference list for sorting.
    
    Return
    ------
    a' : list
        List 'a' sorted with respect to 'b'.
    """
    return np.array(sorted(zip(a,b)))[:,1]

def fit_sphere(data,init=[0,0,0,10]):
    """Fit a sphere to specified data.
    
    Uses a least square fit for optimisation. 
    Return coordinates of the sphere center and its radius as well as the 
    residual variance of the fit.
    
    Parameters
    ----------
    data : 3D numpy array
        Data to be fit with a sphere
        
    init : list
        List of initial parameters: 
        [x0, y0, z0, r0]
    
    Return
    ------
    parameter : list
        list of fitted parameters and residual variance:
        [x,y,z,r,res_var]
    """
    def fitfunc(p, coords):
        x0, y0, z0, _ = p
        x, y, z = coords.T
        return ((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
    errfunc = lambda p, x: fitfunc(p, x) - p[3]**2.
    index = np.array(np.nonzero(data)).T
    p1,_ = opt.leastsq(errfunc, init, args=(index,))
    p1[3] = abs(p1[3])
    p1 = list(p1)
    p1.append(np.var(np.sqrt(np.square(index - p1[:3]).sum(1))-p1[3])) # res_var
    return p1

def fit_paraboloid(data, init=[1,1,1,1,1,0,0,0]):
    """Fit a paraboloid to arbitrarily oriented 3D data.
    
    The paraboloid data can by oriented along an arbitrary axis. Not neccesary
    x,y,z. The function rotates the data points and returns the rotation angles
    along the x,y,z axis.
    Returns the parameters for a praboloid along the z-axis. The angles can be
    used to correct the paraboloid for rotation. 
    Paraboloid equation : p1*x**2.+p2*y**2.+p3*x+p4*y+p5 = z 
        
    Parameters
    ----------
    data : array or list of cartesian coordinates
        Cartesian coordinates of the data points.
    
    init : list of 8 scalars
        Initial parameter set.
        [p1, p2, p3, p4, p5, alpha, beta, gamma]
        *alpha: rotation around x axis
        *beta: rotation around y axis
        *gamma: rotation around z axis
    
    Return
    ------
    popt : list of 8 scalars
        Optimised parameters.
        [p1, p2, p3, p4, p5, alpha, beta, gamma]
        *alpha: rotation around x axis
        *beta: rotation around y axis
        *gamma: rotation around z axis
    """
    def errfunc(p,coord):
        p1, p2, p3, p4, p5, alpha, beta, gamma = p
        coord = rot_coord(coord, [alpha,beta,gamma])
        x, y, z = np.array(coord).T
        return p1*x**2.+p2*y**2.+p3*x+p4*y+p5 - z
#         return p1*x**2.+p2*y**2.+p3*x*y+p4*x+p5*y+p6 - z
    popt,_ = opt.leastsq(errfunc, init, args=(data,))
    return popt

def get_paraboloid_apex(p):
    """Return the apex of a paraboloid.
    
    Use the return of fit_paraboloid() to compute the apex of the paraboloid.
    The return is in the coordinate system of the data, meaning that the 
    coordinates have been corected for the rotation angles.
    
    Parameters
    ----------
    p : list of 8 scalars
        Optimised parameters.
        [p1, p2, p3, p4, p5, alpha, beta, gamma]
        *alpha: rotation around x axis
        *beta: rotation around y axis
        *gamma: rotation around z axis
    
    Return
    ------
    coord : list
        List of the apex' [x,y,z] coordinates.
    """
    p1, p2, p3, p4, p5, alpha, beta, gamma = p
    x = -p3/(2.*p1)
    y = -p4/(2.*p2)
    z = p1*x**2.+p2*y**2.+p3*x+p4*y+p5
    return rot_coord(np.array([[x,y,z],]), [alpha, beta, gamma], True)[0]

def rot_coord(coord, angles, invert=False):
    """Rotate given coordinates by specified angles.
    
    Use rotation matrices to rotate a list of coordinates around the x,y,z axis
    by specified angles alpha,beta,gamma.
    
    Parameters
    ----------
    coord : list of coordinates
        List of catesian coordinates
        
    angles : list of three scalars
        List specifying the rotation angles. See description.
    
    invert : boolean
        True inverts the used rotation matrix. This can be used for undoing a
        rotation.
    
    Return
    ------
    rotated_coord : list of cordinates
        List of rotated cartesian coordinates
    """
    alpha, beta, gamma = angles
    xyz = np.zeros(np.shape(coord))
    Rx = np.array([[1,0,0],
                   [0,np.cos(alpha),-np.sin(alpha)],
                   [0,np.sin(alpha),np.cos(alpha)]])
    Ry = np.array([[np.cos(beta),0,np.sin(beta)],
                   [0,1,0],
                   [-np.sin(beta),0,np.cos(beta)]])
    Rz = np.array([[np.cos(gamma),-np.sin(gamma),0],
                   [np.sin(gamma),np.cos(gamma),0],
                   [0,0,1]])
    if invert == True:
        R = np.linalg.inv(Rx.dot(Ry.dot(Rz)))
    elif invert == False:
        R = Rx.dot(Ry.dot(Rz))
    for i in range(np.shape(coord)[0]):
        xyz [i,:] = R.dot(np.array(coord[i,:]))
    return xyz

def cart2sphere(xyz):
    """Convert cartesian coordinates into spherical coordinates.
    
    Convert a list of cartesian coordinates x,y,z to spherical coordinates
    r,theta,phi. theta is defined as 0 along z-axis.
    
    Parameters
    ----------
    xyz : list
        List of cartesian coordinates
    
    Return
    ------
    rtp : list
        List of spherical coordinates
    """
    rtp = np.zeros(xyz.shape)
    xy = xyz[:, 0] ** 2 + xyz[:, 1] ** 2
    rtp[:, 0] = np.sqrt(xy + xyz[:, 2] ** 2)
    rtp[:, 1] = np.arctan2(np.sqrt(xy), xyz[:, 2])  # for elevation angle defined from Z-axis down
    rtp[:, 2] = np.arctan2(xyz[:, 1], xyz[:, 0])
    return rtp

def sphere2cart(rtp):
    """Convert spherical coordinates into cartesian coordinates.
    
    Convert a list of spherical coordinates r,theta,phi to cartesian coordinates
    x,y,z. theta is defined as 0 along z-axis.
    
    Parameters
    ----------
    rtp : list
        List of spherical coordinates
    
    Return
    ------
    xyz : list
        List of cartesian coordinates
    """
    xyz = np.zeros(rtp.shape)
    xyz[:,0] = rtp[:,0]*np.sin(rtp[:,1])*np.cos(rtp[:,2])
    xyz[:,1] = rtp[:,0]*np.sin(rtp[:,1])*np.sin(rtp[:,2])
    xyz[:,2] = rtp[:,0]*np.cos(rtp[:,1])
    return xyz


def array_from_vtk_polydata(poly,size=[]):
    """Create a boolean numpy array from vtkPolyData.
    
    The size of the created numpy array can be specified. It has to be equal or 
    larger than the bounds of the polydata. Scalar values in the polyData are
    ignored. All nonzero values in poly are transformed into ones in the return.
    Note: works only in 3D.
    
    Parameters
    ----------
    poly : vtkPolyData (3D)
        polyData to be transformed.
        
    size : tuple / list of int
        Size of the created numpy array. By default set to the boundd of the 
        polyData. e.g. (3,4,5)
    
    Return
    ------
    array : numpy array
        Boolean numpy array corresponding to the nonzero points in poly.
    """
    if np.shape(size) == np.shape([]):
        size = np.array(poly.GetPoints().GetBounds(),dtype='int')[1::2]
    indices = np.array(vtk_to_numpy(poly.GetPoints().GetData()),dtype='int')
    out = np.zeros(size)
    out[indices[:,0]-1,indices[:,1]-1,indices[:,2]-1] = 1
    return np.array(out)

def view_polydata(poly,firstcolor=(),lastcolor=()):
    """Display vtkPolyData. Can show superposition of many vtkPolyData.
    
    If input is a list of vtkPolyData, displays all of them in one viewer.
    
    Parameters
    ----------
    poly : vtkPolyData (3D) / list of vtkPolyData (3D)
        polyData to be displayed. List or single polydata.
        
    Return
    ------
    no return : 
        Opens render window.

    """
    if np.shape(poly) == ():
        numel = 1
        poly = [poly]
    else:
        numel = np.shape(poly)[0]
    if np.shape(firstcolor) != np.shape(()) and np.shape(lastcolor) != np.shape(()):
        colors = rgb_list(numel,firstcolor=firstcolor,lastcolor=lastcolor)
    elif np.shape(firstcolor) != np.shape(()):
        colors = rgb_list(numel,firstcolor=firstcolor)
    elif np.shape(lastcolor) != np.shape(()):
        colors = rgb_list(numel,lastcolor=lastcolor)
    else:
        colors = rgb_list(numel)
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
   
def rgb_list(N,firstcolor=(),lastcolor=(),s=1.,v=1.):
    """Generate a list of N distinct RGB tuples.
    
    The first and last entry of the list can be specified. The list will still
    have N entries. The range of each tuple entry is between 0. and 1. The list
    goes from red over green to blue and purple.
    
    Parameters
    ----------
    N : int
        Number of returned RGB-tuples.
    
    firstcolor : list, three components
        First entry of the returned RGB-list.
        
    lastcolor : list, three components
        Last entry of the returned RGB-list.
        
    s : float between 0 and 1
        Saturation value of the list
        
    v : float between 0 and 1
        Value value of the list
    
    Return
    ------
    RGB-list : list of tuples with three components
        List with N RGB-tuples
    """
    if len(firstcolor) == 3 and len(lastcolor) == 3:
        HSV_tuples = [(float(x)/float(N), s, v) for x in range(N-2)]
        RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
        RGB_tuples.insert(0, firstcolor)
        RGB_tuples.append(lastcolor)
    elif len(firstcolor) == 3:
        HSV_tuples = [(float(x)/float(N), s, v) for x in range(N-1)]
        RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
        RGB_tuples.insert(0, firstcolor)
    elif len(lastcolor) == 3:
        HSV_tuples = [(float(x)/float(N), s, v) for x in range(N-1)]
        RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
        RGB_tuples.append(lastcolor)
    else:
        HSV_tuples = [(float(x)/float(N), s, v) for x in range(N)]
        RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    return RGB_tuples

def sphere_voulume(radius):
    """Return the volume of a sphere.
    
    Is able to process numpy arrays.
    
    Parameters
    ----------
    R : ndarray of floats
        Radi of spheres.
        
    Return
    ------
    Volumes : ndarray of floats
        Volumes of spheres.
    """
    return  4./3.*np.pi*radius**3.

def angle_difference(array):
    """Return the differences between consecutive angles in an array.
    
    Computes both clockwise and counterclockwise angle differences.
    Angles need to be in degree.
    
    Parameters
    ----------
    array : list with n entries
        List of angles.
        
    Return
    ------
    angle differences : list with two entries each with n-1 entries
        *return[0]: clockwise angle differences
        *return[1]: counterclockwise angle differences
    """
    clockwise = np.ediff1d(array) % 360.
    counterclockwise = abs(360.-clockwise)
    return clockwise, counterclockwise

def create_logfile(logs,path = None):
    """Create a logfile for saving and loading AutoPhenotype data from boolean 
    list.
    
    Return logfile and optionally save it as .csv file.
    
    Parameters
    ----------
    logs : list of 1,0
        Has to have six entries. 0: deactivated, 1: activated
            *Entry 0: data
            *Entry 1: contour
            *Entry 2: mesh
            *Entry 3: features
            *Entry 4: results
            *Entry 5: tags
        
    path = string
        If specified, the logfile is created at the specified location as .csv
        file. Path has to include file name.
    
    Return
    ------
    logfile : logfile
        Logfile which can be used in saving, loading AutoPhenotype data. Format
        is pandas.DataFrame()
    """
    logfile = pd.DataFrame([logs],columns=['data',
                                           'contour',
                                           'mesh',
                                           'features',
                                           'results',
                                           'tags'])
    if type(path) == type(None):
        return logfile
    if type(path) == type('string'):
        logfile.to_csv(path)
        return logfile

def logfile_from_str(what, path = None):
    """Create a logfile for saving and loading AutoPhenotype data from keywords
    
    Uses keywords to generate a logfile. Keywords can be:
        *'all': everything below
        *'data': input data as .tif
        *'contour': contour fit as .tif
        *'mesh': mesh as vtk data .vtp
        *'features': features as vtk data .vtp
        *'results': results as .csv
    
    Parameters
    ----------
    what : list of strings
        Use specified strings from description.
        Note: has to be list, even if only one keyword is specified.
    
    path = string
        If specified, the logfile is created at the specified location as .csv
        file. Path has to include file name.
    
    Return
    ------
    logfile : logfile
        Logfile which can be used in saving, loading AutoPhenotype data. Format
        is pandas.DataFrame()
    """
    logs = [0,0,0,0,0,0]
    if any(t == 'all' for t in what):
        logs = [1,1,1,1,1,1]
    if any(t == 'data' for t in what):
        logs[0] = 1
    if any(t == 'contour' for t in what):
        logs[1] = 1
    if any(t == 'mesh' for t in what):
        logs[2] = 1
    if any(t == 'features' for t in what):
        logs[3] = 1
    if any(t == 'results' for t in what):
        logs[4] = 1
    if any(t == 'tags' for t in what):
        logs[5] = 1
    return create_logfile(logs, path)

def blockPrint():
    """Suppress all print output after call."""
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    """Enable all print output after call, if it was blocked before."""
    sys.stdout = sys.__stdout__

#Docstring template
"""[short description]

[Detailed description]

Parameters
----------
parameter_1 : dtype
    description

Return
------
return_1 : dtype
    description
"""

