dir=getDirectory("Choose a Directory"); 
print(dir); 
dirList = getFileList(dir);

for (j=0; j<dirList.length; j++) { 
   if (endsWith(dirList[j], "hrs/")){  
       fileList = getFileList(dir+ dirList[j]); 

       for (i=0; i<fileList.length; i++) { 
	    if (endsWith(fileList[i], "_trim-acylYFP.tif")){ 
		    print(i + ": " + dir+dirList[j] + fileList[i]); 
		    open(dir+dirList[j] + fileList[i]); 
		    nameIntensityImage =getTitle(); 
		    baseNameEndIntImage=indexOf(nameIntensityImage, "_trim-acylYFP.tif"); 
		    baseNameIntImage=substring(nameIntensityImage, 0, baseNameEndIntImage); 
		    selectWindow(nameIntensityImage); 
                    run("Brightness/Contrast...");
                    setMinAndMax(0, 163);
                    call("ij.ImagePlus.setDefault16bitRange", 8); 
	    } 
       }

       for (i=0; i<fileList.length; i++) { 
            if (endsWith(fileList[i], "_walls.tif")){
                    print(i + ": " + dir + dirList[j] + fileList[i]); 
		    open(dir + dirList[j] + fileList[i]); 
		    nameWallsImage = getTitle(); 
                    hminInd = indexOf(nameWallsImage, "_hmin"); 
                    cleanInd = indexOf(nameWallsImage, "clean_"); 
                    parametersSeg =substring(nameWallsImage, hminInd, cleanInd); 

                    run("Merge Channels...", "c4=" + nameWallsImage + " c7=" + nameIntensityImage + " keep");
                    selectWindow("RGB");
                    saveAs("Tiff", dir + dirList[j] + baseNameIntImage + parametersSeg + "under_segmentation.tif");
                    close();
                    selectWindow(nameWallsImage);
                    run("Invert LUT");
                    setMinAndMax(-833, 1343);
                    call("ij.ImagePlus.setDefault16bitRange", 8);
                    run("Merge Channels...", "c4=" + nameWallsImage + " c7=" + nameIntensityImage + " keep");
                    selectWindow("RGB");
                    saveAs("Tiff", dir + dirList[j] + baseNameIntImage + parametersSeg + "over_segmentation.tif");
                    close();
                    selectWindow(nameWallsImage);
                    close();
            } 
       }   
       run("Close All"); 
    }
}
