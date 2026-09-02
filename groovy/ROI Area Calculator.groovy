def um = getCurrentServer().getPixelCalibration().getPixelWidthMicrons()
double area = 0
int n = 0
for (ann in getAnnotationObjects()) {
    if (ann.getROI() == null || !ann.getROI().isArea()) continue
    area += ann.getROI().getArea() * um * um / 1e6
    n++
}
println String.format("%-58s  %d ann  %8.2f mm²", getProjectEntry().getImageName(), n, area)