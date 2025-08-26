
// Import required QuPath and Java classes
import qupath.lib.images.servers.ImageServers
import javax.imageio.ImageIO
import qupath.lib.regions.RegionRequest
import qupath.lib.gui.QuPathGUI
import qupath.lib.projects.ProjectImageEntry
import javafx.application.Platform


// Path to the main TMA image and the name to use in the project
def imagePath = "/home/aboffelli/PACC/TMA/Data/Images/SweBcg_CD8/CD8_20A.npdi"
def imageName = "main_tma"


// The script needs a project to run, so check that first
def project = getProject()
if (project == null) {
    println "A project must be opened before running this script"
    return
}


// Check if the main TMA image is already in the project
def existingEntry = project.getImageList().find { entry -> entry.getImageName() == imageName }

if (!existingEntry) {
    // If not, add the image to the project and set its name
    def imageServer = ImageServers.buildServer(imagePath)
    def imageEntry = project.addImage(imageServer.getBuilder())
    imageEntry.setImageName(imageName)
    // Load the image entry as the current image data
    getQuPath().openImageEntry(imageEntry)
} else {
    // If already present, just open it
    getQuPath().openImageEntry(existingEntry)
}

// Sync and refresh project to ensure changes are visible
project.syncChanges()
getQuPath().refreshProject()


// Define export resolution (1.0 = full size)
double downsample = 1.0

// Check if TMA dearraying has been performed; if not, run the plugin and exit
if (!isTMADearrayed()) {
    print "Running TMA Dearrayer"
    runPlugin('qupath.imagej.detect.dearray.TMADearrayerPluginIJ', '{"coreDiameterMM":1.3,"labelsHorizontal":"A-M","labelsVertical":"1-10","labelOrder":"Row first","densityThreshold":5,"boundsScale":105}')
    println "Please check the TMA grid, fix if needed and rerun the script."
    return
} else {
    print "TMA Dearrayer already performed."
}


// Create output directory for core images inside the project
def dirPath = buildFilePath(PROJECT_BASE_DIR)
def folderName = dirPath.split("/")[-1]  // Linux/Mac
// def folderName = dirPath.split("\\\\")[-1]  // Windows (uncomment if needed)

if (folderName) {
    // Build the full output path and create the output directory for core images
    dirPath = buildFilePath(dirPath, folderName)
    def dirOutput = dirPath + "_cores"
    mkdirs(dirOutput)

    // Export each TMA core as a separate TIFF image
    def server = getCurrentImageData().getServer()
    def path = server.getPath()
    def cores = getTMACoreList()
    
    cores.each { core ->
        def coreName = core.getName()
        def outputFile = new File(dirOutput, core.getName() + '.tif')
        // Read the region for the current core and save as TIFF
        img = server.readRegion(RegionRequest.createInstance(path, downsample, core.getROI()))
        ImageIO.write(img, 'TIF', outputFile)
        println "Image saved for core: ${coreName}"
    }
    print("Images saved on: " + dirOutput)

    // Add exported core images to the project, replacing any with the same name
    def folder = new File(dirOutput)
    def imageFiles = folder.listFiles().findAll { it.name.endsWith('.tif') }  // Adjust the extension as needed
    
    // Flag to track if any existing images were found
    def existingImageFound = false      
    
    imageFiles.each { file ->
        // Get the image path for the current file
        def currentImagePath = file.getAbsolutePath()
        // Build the image server for each image file
        def currentImageServer = ImageServers.buildServer(currentImagePath)
        // Set a unique name based on the file name (without extension)
        def currentImageName = file.getName().replaceAll(/\.\w+/, "")  // Removes extension
        // Check if the image with this name is already in the project
        def existingImage = project.getImageList().find { entry -> entry.getImageName() == currentImageName }
        if (existingImage) {
            println "Removing old image '${existingImage}' from project."
            project.removeImage(existingImage, false)  // false avoids moving to trash
            project.syncChanges()
            getQuPath().refreshProject()
            existingImageFound = true // Set the flag to true
        }
        // Add image to the project
        def currentImageEntry = project.addImage(currentImageServer.getBuilder())      
        currentImageEntry.setImageName(currentImageName)
        println "Added image: ${currentImageName}"
    }
    // Sync and refresh project after adding images
    project.syncChanges()
    getQuPath().refreshProject()
    
    // Print message only if existing images were found
    if (existingImageFound) {
        println "All images added to the project successfully. Since one or more images changed, you may need to restart QuPath."
    } else {
        println "All images added to the project successfully."
    }
}



// End of automatic upload of images.
