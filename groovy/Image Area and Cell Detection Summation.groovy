// ============================================================
//  Project totals — DAPI detections and analysed area
//  Plain Run (walks the project itself, not Run for project)
//
//  Ensure accuracy of image file names if including or excluding. 
// ============================================================

def INCLUDE = []                        // whitelist, e.g. ['Donor1_', 'Donor2_']
def EXCLUDE = []                        // ignored when INCLUDE is non-empty

def keep = { String name ->
    if (!INCLUDE.isEmpty()) return INCLUDE.any { name.contains(it) }
    return !EXCLUDE.any { name.contains(it) }
}

def rows = []
def skipped = []
double totalArea = 0
long totalCells = 0

for (entry in getProject().getImageList()) {
    def name = entry.getImageName()
    if (!keep(name)) { skipped << name; continue }

    def imageData = entry.readImageData()
    def hier = imageData.getHierarchy()
    double um = imageData.getServer().getPixelCalibration().getPixelWidthMicrons()

    int n = hier.getDetectionObjects().size()
    double areaPx = 0
    for (a in hier.getAnnotationObjects()) {
        if (a.getROI() != null && a.getROI().isArea()) areaPx += a.getROI().getArea()
    }
    double mm2 = areaPx * um * um / 1e6

    rows << [name, mm2, n, mm2 > 0 ? n / mm2 : 0]
    totalArea += mm2
    totalCells += n
    imageData.getServer().close()
}

// ---------- formatted table ----------
println String.format("\n%-52s %10s %12s %10s", "image", "area_mm2", "nuclei", "nuc/mm2")
println "-".multiply(88)
rows.each {
    println String.format("%-52s %10.2f %12s %10.0f",
        it[0], it[1], String.format("%,d", it[2]), it[3])
}
println "-".multiply(88)
println String.format("%-52s %10.2f %12s %10.0f",
    "TOTAL (${rows.size()} sections)", totalArea, String.format("%,d", totalCells),
    totalArea > 0 ? totalCells / totalArea : 0)

// ---------- excluded ----------
println ""
if (skipped) {
    println "EXCLUDED (${skipped.size()}):"
    skipped.each { println "   ${it}" }
} else {
    println "EXCLUDED: none"
}
println "Included ${rows.size()} of ${getProject().getImageList().size()} images in project"

// ---------- raw values being summed ----------
println "\n--- individual values summed (tab-separated, paste into Excel) ---"
println "image\tarea_mm2\tnuclei\tnuc_per_mm2"
rows.each { println String.format("%s\t%.4f\t%d\t%.1f", it[0], it[1], it[2], it[3]) }

println "\nareas summed (${rows.size()} values):"
println "   " + rows.collect { String.format("%.4f", it[1]) }.join(" + ")
println "   = ${String.format('%.4f', totalArea)} mm²"

println "\nnuclei summed (${rows.size()} values):"
println "   " + rows.collect { it[2] as String }.join(" + ")
println "   = ${totalCells}"

// ---------- summary ----------
def dens = rows.collect { it[3] }.toSorted()
double med = dens.isEmpty() ? 0 :
    (dens.size() % 2 == 1 ? dens[dens.size().intdiv(2)]
                          : (dens[dens.size().intdiv(2) - 1] + dens[dens.size().intdiv(2)]) / 2.0)

println "\nFOR METHODS:  ${String.format('%,d', totalCells)} nuclei across " +
        "${String.format('%.1f', totalArea)} mm² (${rows.size()} sections)"
println "median nuclear density: ${String.format('%.0f', med)} nuclei/mm²"