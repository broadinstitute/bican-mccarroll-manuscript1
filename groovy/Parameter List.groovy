def ml = getDetectionObjects()[0].getMeasurementList()
def names = []
try { names = ml.getMeasurementNames() } catch (Exception e) { names = ml.getNames() }
println "\n${names.size()} measurements:"
names.each { println "   ${it}" }