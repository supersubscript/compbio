// CLV3 image processing + Costanza analysis
// 28 November 2014
// For use with pCLV3-dsRED-pWUS-GFP-OutNPA time series

dir1=getDirectory("Choose Destination Directory ");

run("Split Channels");
selectImage(1);
close();
selectImage(2);
close();
selectImage(1);
getVoxelSize(xp,yp,zp,unit);
title = getTitle();
run("Fire");
fileName1="raw-"+title
saveAs("Tiff", dir1+fileName1);


run("Costanza Plugin");
close("Costanza Plugin");
selectWindow("Costanza - Basins of attractions (BOA)");
close();
selectWindow("Costanza - Cell centers");
close();
selectWindow("Costanza - Basins of attractions (BOA)-intensity");
setVoxelSize(xp,yp,zp,unit);
fileName2="BOAi-"+title
saveAs("Tiff", dir1+fileName2);
close()
selectWindow("Results");
run(Save)
close();

