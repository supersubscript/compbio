dir=getDirectory("Choose a Directory"); 
print(dir); 
splitDir=dir; 
print(splitDir);
fileList = getFileList(dir); 

for (i=0; i<fileList.length; i++) { 
	if (endsWith(fileList[i], "_trim_stretch-clv3.tif")){ 
		print(i + ": " + dir + fileList[i]); 
		open(dir + fileList[i]); 
		nameIntensityImage =getTitle(); 
		baseNameEndIntImage=indexOf(nameIntensityImage, "_trim_stretch-clv3.tif"); 
		baseNameIntImage=substring(nameIntensityImage, 0, baseNameEndIntImage); 
		selectWindow(nameIntensityImage); 
                run("Brightness/Contrast...");
                setMinAndMax(0, 163);
                call("ij.ImagePlus.setDefault16bitRange", 8); 
	} 
}

for (i=0; i<fileList.length; i++) { 

        if (endsWith(fileList[i], "_trim_stretch-clv3_seg.tif")){
                print(i + ": " + dir+fileList[i]); 
		open(dir+fileList[i]); 
		nameWallsImage = getTitle(); 
              //  hminInd = indexOf(nameWallsImage, "_hmin"); 
               // cleanInd = indexOf(nameWallsImage, "clean_"); 
                //parametersSeg =substring(nameWallsImage, hminInd, cleanInd); 
                run("8-bit");  
                run("Merge Channels...", "c4=" + nameWallsImage + " c7=" + nameIntensityImage + " create keep ");
                selectWindow("Composite");
               // saveAs("Tiff", splitDir + baseNameIntImage + parametersSeg + "under_segmentation.tif");
                saveAs("Tiff", splitDir + baseNameIntImage + "under_segmentation.tif");
                close();

                selectWindow(nameWallsImage);
                run("Invert LUT");
                setMinAndMax(-833, 1343);
                call("ij.ImagePlus.setDefault16bitRange", 8);
                run("Merge Channels...", "c4=" + nameWallsImage + " c7=" +nameIntensityImage + " create keep");
                selectWindow("Composite"); 
               // saveAs("Tiff", splitDir + baseNameIntImage + parametersSeg + "over_segmentation.tif");
                saveAs("Tiff", splitDir + baseNameIntImage + "over_segmentation.tif");
                close(); 
                selectWindow(nameWallsImage);
                close(); 
        } 
} 

run("Close All"); 
