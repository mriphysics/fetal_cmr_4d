#### path/variable admin
fcmrNum = 202
foldExt = '_hrh_fullRecon'
# foldExt = ''
pvStateExt = ''
velDir = '\\vel_vol_4d'
# velDir = '\\vel_vol_trans_4d'

# pvFold = '\\paraview' + '\\'
# pvFold = '\\paraview_polyCorr' + '\\'
pvFold = '\\paraview_polyCorr_aorta_LV_RV_LOT_ROT_LA_RA_IVC_SVC_PA_DA' + '\\'
# pvFold = '\\paraview_polyCorr_aorta_LV_RV_LOT_ROT_LA_RA_IVC_SVC' + '\\'

fcmrDir = 'E:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\4D_Flow_Paper\\fcmr'
path = fcmrDir + str(fcmrNum) + foldExt + velDir + pvFold

pvStateFilename = 'fcmr' + str(fcmrNum) + '_' + pvStateExt + 'paraview.pvsm'

filenameCine = 'cine_vol_masked_t'
filenameVelVol = 'vel_vol_masked_VxVy-Vz_t'
numFrames = 25



#### functions

def makeArrayFilenames( path, filename, numFiles ):
    
    arrayFilenames = []
    for f in range(numFiles):
        currFilename = path + filename + str(f) + '.vtk'
        arrayFilenames.append(currFilename)
    
    return arrayFilenames



#____________________ Paraview Script ____________________#

# import
from paraview.simple import * #### import the simple module from the paraview

#### Load data

# disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# load in cine_vol .vtk files
arrayFilenamesCine = makeArrayFilenames(path, filenameCine, numFrames)
cine_vol = LegacyVTKReader(FileNames=arrayFilenamesCine, registrationName="cine_vol_masked")

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# load in vel_vol .vtk files
arrayFilenamesVelVol = makeArrayFilenames(path, filenameVelVol, numFrames)
vel_vol = LegacyVTKReader(FileNames=arrayFilenamesVelVol, registrationName="vel_vol_masked")

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [955, 778]



#### Press Apply

# show data in view
cine_vol_Display = Show(cine_vol, renderView1)
# trace defaults for the display properties.
cine_vol_Display.Representation = 'Outline'
cine_vol_Display.AmbientColor = [0.0, 0.0, 0.0]
cine_vol_Display.ColorArrayName = ['POINTS', '']
cine_vol_Display.OSPRayScaleArray = 'magnitude_intensity'
cine_vol_Display.OSPRayScaleFunction = 'PiecewiseFunction'
cine_vol_Display.SelectOrientationVectors = 'None'
cine_vol_Display.ScaleFactor = 7.2
cine_vol_Display.SelectScaleArray = 'magnitude_intensity'
cine_vol_Display.GlyphType = 'Arrow'
cine_vol_Display.GlyphTableIndexArray = 'magnitude_intensity'
cine_vol_Display.DataAxesGrid = 'GridAxesRepresentation'
cine_vol_Display.PolarAxes = 'PolarAxesRepresentation'
cine_vol_Display.ScalarOpacityUnitDistance = 1.7383632950303343
#cine_vol_Display.InputVectors = [None, '']

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
cine_vol_Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
cine_vol_Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
cine_vol_Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
cine_vol_Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
cine_vol_Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
cine_vol_Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
cine_vol_Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
cine_vol_Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
cine_vol_Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
cine_vol_Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
cine_vol_Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# show data in view
vel_vol_Display = Show(vel_vol, renderView1)
# trace defaults for the display properties.
vel_vol_Display.Representation = 'Outline'
vel_vol_Display.AmbientColor = [0.0, 0.0, 0.0]
vel_vol_Display.ColorArrayName = ['POINTS', '']
vel_vol_Display.OSPRayScaleArray = 'velocity_magnitude'
vel_vol_Display.OSPRayScaleFunction = 'PiecewiseFunction'
vel_vol_Display.SelectOrientationVectors = 'vector_field'
vel_vol_Display.ScaleFactor = 7.2
vel_vol_Display.SelectScaleArray = 'velocity_magnitude'
vel_vol_Display.GlyphType = 'Arrow'
vel_vol_Display.GlyphTableIndexArray = 'velocity_magnitude'
vel_vol_Display.DataAxesGrid = 'GridAxesRepresentation'
vel_vol_Display.PolarAxes = 'PolarAxesRepresentation'
vel_vol_Display.ScalarOpacityUnitDistance = 1.7383632950303343
#vel_vol_Display.InputVectors = ['POINTS', 'vector_field']

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
vel_vol_Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
vel_vol_Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
vel_vol_Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
vel_vol_Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
vel_vol_Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
vel_vol_Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
vel_vol_Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
vel_vol_Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
vel_vol_Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
vel_vol_Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
vel_vol_Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()



#### Update colormap/transparency for cine_vol

# set active source
SetActiveSource(cine_vol)

# set scalar coloring
ColorBy(cine_vol_Display, ('POINTS', 'magnitude_intensity'))

# rescale color and/or opacity maps used to include current data range
cine_vol_Display.RescaleTransferFunctionToDataRange(True, True)

# change representation type
cine_vol_Display.SetRepresentationType('Volume')

# Properties modified on cine_vol_Display
cine_vol_Display.SelectMapper = 'Resample To Image'

# get color transfer function/color map for 'magnitude_intensity'
magnitude_intensityLUT = GetColorTransferFunction('magnitude_intensity')
magnitude_intensityLUT.RGBPoints = [-7.28000020980835, 0.0, 0.0, 0.0, 87.65499806404114, 0.0, 0.0, 0.0, 182.58999633789062, 1.0, 1.0, 1.0]
magnitude_intensityLUT.NanColor = [1.0, 0.0, 0.0]
magnitude_intensityLUT.ScalarRangeInitialized = 1.0

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
magnitude_intensityLUT.ApplyPreset('Grayscale', True)

# get opacity transfer function/opacity map for 'magnitude_intensity'
magnitude_intensityPWF = GetOpacityTransferFunction('magnitude_intensity')
magnitude_intensityPWF.Points = [-7.28000020980835, 0.0, 0.5, 0.0, 87.65499806404114, 0.0, 0.5, 0.0, 182.58999633789062, 1.0, 0.5, 0.0]
magnitude_intensityPWF.AllowDuplicateScalars = 1
magnitude_intensityPWF.ScalarRangeInitialized = 1

# update the view to ensure updated data information
renderView1.Update()




#### Make velocity vectors for vel_vol

# set active source
SetActiveSource(vel_vol)

# create a new 'Glyph'
glyph1 = Glyph(Input=vel_vol,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'velocity_magnitude']
glyph1.Vectors = ['POINTS', 'vector_field']
glyph1.ScaleMode = 'scalar'
glyph1.ScaleFactor = 7.2
glyph1.GlyphMode = 'All Points'
glyph1.MaximumNumberOfSamplePoints = 10000000
glyph1.GlyphTransform = 'Transform2'

# hide data in view
Hide(vel_vol, renderView1)

# Properties modified on glyph1
glyph1.ScaleFactor = 0.1

# get color transfer function/color map for 'velocity_magnitude'
velocity_magnitudeLUT = GetColorTransferFunction('velocity_magnitude')
velocity_magnitudeLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 11.645000457763672, 0.865003, 0.865003, 0.865003, 23.290000915527344, 0.705882, 0.0156863, 0.14902]
velocity_magnitudeLUT.ScalarRangeInitialized = 1.0

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.AmbientColor = [0.0, 0.0, 0.0]
glyph1Display.ColorArrayName = ['POINTS', 'velocity_magnitude']
glyph1Display.LookupTable = velocity_magnitudeLUT
glyph1Display.OSPRayScaleArray = 'velocity_magnitude'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'GlyphVector'
glyph1Display.ScaleFactor = 7.2
glyph1Display.SelectScaleArray = 'velocity_magnitude'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'velocity_magnitude'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'
glyph1Display.GaussianRadius = 3.6
glyph1Display.SetScaleArray = ['POINTS', 'velocity_magnitude']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'velocity_magnitude']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
#glyph1Display.InputVectors = ['POINTS', 'GlyphVector']

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
glyph1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
glyph1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on renderView1
renderView1.OrientationAxesVisibility = 0

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# Rescale transfer function
velocity_magnitudeLUT.RescaleTransferFunction(0.0, 30.0)

# get opacity transfer function/opacity map for 'velocity_magnitude'
velocity_magnitudePWF = GetOpacityTransferFunction('velocity_magnitude')
velocity_magnitudePWF.Points = [0.0, 0.0, 0.5, 0.0, 23.290000915527344, 1.0, 0.5, 0.0]
velocity_magnitudePWF.ScalarRangeInitialized = 1

# Rescale transfer function
velocity_magnitudePWF.RescaleTransferFunction(0.0, 30.0)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
velocity_magnitudeLUT.ApplyPreset('Rainbow Desaturated', True)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [955, 778]

# reset camera centre of rotation
ResetCamera()
renderView1.CenterOfRotation = GetActiveCamera().GetFocalPoint()
Render()

# current camera placement for renderView1
renderView1.CameraPosition = [50.0, 100.0, 30.0]
renderView1.CameraFocalPoint = [0.0, -110.0, 30.0]
renderView1.CameraViewUp = [0.0, 0.0, 0.0]
renderView1.CameraParallelScale = 55

# save screenshot - for debugging...
# SaveScreenshot(path + '/test.png', renderView1, ImageResolution=[955, 778])

# save state
SaveState(fcmrDir + str(fcmrNum) + foldExt + velDir + pvFold + pvStateFilename)