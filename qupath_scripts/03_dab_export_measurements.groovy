import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathDetectionObject

// Get the list of all images in the current project
def project = getProject()
def imagesToExport = project.getImageList()

// Separate each measurement value in the output file with a tab ("\t")
def separator = "\t"

// Choose the columns that will be included in the export
// Note: if 'columnsToInclude' is empty, all columns will be included
def columnsToInclude = new String[]{"Image", "Object ID", "ROI", "Centroid X µm", 
"Centroid Y µm", "Area µm^2", "Length µm", "Circularity", "Solidity", "Max diameter µm",
"Min diameter µm", "Hematoxylin: Mean", "Hematoxylin: Median", "Hematoxylin: Min", 
"Hematoxylin: Max", "Hematoxylin: Std.Dev.", "DAB: Mean", "DAB: Median", "DAB: Min", 
"DAB: Max", "DAB: Std.Dev."}

// Choose the type of objects that the export will process
// Other possibilities include:
//    1. PathAnnotationObject
//    2. PathDetectionObject
//    3. PathRootObject
// Note: import statements should then be modified accordingly
def exportType = PathDetectionObject.class

// Create the output file path
def projectName = project.toString().replaceFirst(/^Project:\s*/, "").replaceFirst(/-project$/, "") // removes 'Project: ' prefix and '-project' suffix
def outputPath = buildFilePath(PROJECT_BASE_DIR, "${projectName}_measurements.txt")
def outputFile = new File(outputPath)

// Create the measurementExporter and start the export
def exporter  = new MeasurementExporter()
                  .imageList(imagesToExport)            // Images from which measurements will be exported
                  .separator(separator)                 // Character that separates values
                  .includeOnlyColumns(columnsToInclude) // Columns are case-sensitive
                  .exportType(exportType)               // Type of objects to export
                  //.filter(obj -> obj.getPathClass() == getPathClass("Tumor"))    // Keep only objects with class 'Tumor'
                  .exportMeasurements(outputFile)        // Start the export process

print "File saved: ${outputPath}"
print "Done!"