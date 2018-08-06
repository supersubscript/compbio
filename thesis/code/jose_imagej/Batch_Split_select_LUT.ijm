dir1 = getDirectory("Choose Source Directory ");
dir2 = getDirectory("Choose Destination Directory ");

list = getFileList(dir1);

setBatchMode(true);

for (i=0; i<list.length; i++) {
	
	showProgress(i+1, list.length);
	open(dir1+list[i]);

	
	run("Split Channels");
	selectImage(1);
	close();
	selectImage(2);
	close();
	selectImage(1);
	getVoxelSize(xp,yp,zp,unit);
	title = getTitle();
	run("Fire");
	fileName1="raw-"+title;
	saveAs("Tiff", dir2+fileName1);
	run("Costanza Plugin");
	close("Costanza Plugin");
	selectWindow("Costanza - Basins of attractions (BOA)");
	close();
	selectWindow("Costanza - Cell centers");
	close();
	selectWindow("Costanza - Basins of attractions (BOA)-intensity");
	setVoxelSize(xp,yp,zp,unit);
	fileName2="BOAi-"+title;
	saveAs("Tiff", dir2+fileName2);
	close();
	selectWindow("Results");

}