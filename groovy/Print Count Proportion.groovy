def POS_CLASS = ''      // must match your class name exactly

def dets = getDetectionObjects()
int total = dets.size()
int pos = dets.count { it.getPathClass()?.toString() == POS_CLASS }
int unclassified = dets.count { it.getPathClass() == null }

double um = getCurrentServer().getPixelCalibration().getPixelWidthMicrons()
double areaPx = 0
for (a in getAnnotationObjects()) {
    if (a.getROI() != null && a.getROI().isArea()) areaPx += a.getROI().getArea()
}
double mm2 = areaPx * um * um / 1e6

println String.format("%s\t%.2f\t%d\t%d\t%d\t%.3f",
    getProjectEntry().getImageName(), mm2, total, pos, unclassified,
    total > 0 ? 100.0 * pos / total : 0)