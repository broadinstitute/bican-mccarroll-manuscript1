selectAnnotations()

runPlugin('qupath.imagej.detect.cells.WatershedCellDetection',
    '{"detectionImage": "DAPI",' +
    ' "requestedPixelSizeMicrons": 0.5,' +
    ' "backgroundRadiusMicrons": 8.0,' +
    ' "backgroundByReconstruction": true,' +
    ' "medianRadiusMicrons": 0.0,' +
    ' "sigmaMicrons": 1.5,' +
    ' "minAreaMicrons": 15.0,' +
    ' "maxAreaMicrons": 400.0,' +
    ' "threshold": 100.0,' +
    ' "watershedPostProcess": true,' +
    ' "cellExpansionMicrons": 5.0,' +
    ' "includeNuclei": true,' +
    ' "smoothBoundaries": true,' +
    ' "makeMeasurements": true}')

def n = getDetectionObjects().size()
def um = getCurrentServer().getPixelCalibration().getPixelWidthMicrons()
double mm2 = getAnnotationObjects().sum { it.getROI().getArea() } * um * um / 1e6
println String.format("%-50s  %7d cells   %6.1f mm²   %6.0f cells/mm²",
    getProjectEntry().getImageName(), n, mm2, n / mm2)