import javax.imageio.ImageIO
import qupath.lib.regions.RegionRequest


// Define resolution - 1.0 means full size
double downsample = 1.0

if (!isTMADearrayed()) {
	runPlugin('qupath.imagej.detect.dearray.TMADearrayerPluginIJ', '{"coreDiameterMM":1.3,"labelsHorizontal":"A-M","labelsVertical":"1-10","labelOrder":"Row first","densityThreshold":5,"boundsScale":105}')
	}

// Create output directory inside the project
def dirPath = buildFilePath(PROJECT_BASE_DIR)
def folderName = dirPath.split("/")[-1]  // Linux/Mac
// def folderName = dirPath.split("\\\\")[-1]  \\ Windows 


if (folderName) {
    dirPath = buildFilePath(dirPath, folderName)
    def dirOutput = dirPath + "_cores"
    mkdirs(dirOutput)

    // Write the cores
    def server = getCurrentImageData().getServer()
    def path = server.getPath()
    getTMACoreList().parallelStream().forEach({core ->
        img = server.readRegion(RegionRequest.createInstance(path, downsample, core.getROI()))
        ImageIO.write(img, 'TIF', new File(dirOutput, core.getName() + '.tif'))
})
    print("Images saved on: " + dirOutput)
}

print('Done!')
