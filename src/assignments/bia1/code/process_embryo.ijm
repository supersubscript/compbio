// ImageJ macro for command-line usage. Remember to change the paths if 
// necessary. Usage: 
// 
// path-to/ImageJ-linux64 --allow-multiple --headless -macro path-to/macro.ijm \
//    'Threshold_method useOutliers Gaussian_sigma \
//    Maximum_filter_radius XY_radius Z_radius useGaussian' 

path="/home/hpa22/project_out/"
filelist = getArgument();
file = split(filelist,' ');
argssuff = file[0]+"_"+file[1]+"_"+file[2]+"_"+file[3]+"_"+file[4]+"_"+file[5] +"_"+ file[6];

method=file[0];
open("project_data/Embryo.ome.tif");
run("Duplicate...", "duplicate range=1-100");
run("Properties...", "channels=1 slices=100 frames=1 unit=pixel pixel_width=0.5 pixel_height=0.5 voxel_depth=1.98 global");
run("Duplicate...", "duplicate range=1-100");
Stack.getDimensions(width, height, channels, slices, frames);
Stack.setSlice(50);

// Applying filters
print("Applying filters...");
run("Enhance Contrast", "saturated=0.35");
run("Apply LUT", "stack");
if(file[6] == 1){
  run("Gaussian Blur...", "sigma="+file[2]+" stack");
}
run("Maximum...", "radius="+file[3]+" stack");
//run("Despeckle", "stack");
if(file[1] == 1){
  run("Remove Outliers...", "radius=10 threshold=0 which=Bright stack");
}
rename("Embryo.ome-2.tif");
run("Duplicate...", "duplicate range=1-100");

// Set threshold and apply mask
print("Setting threshold and applying mask...");
setAutoThreshold(""+method+" dark");
run("Convert to Mask", "method="+method+" background=Dark");
saveAs("Embryo.ome-3.tif",path+argssuff+".tif");
rename(path+argssuff+".tif");
run("Divide...", "value=255 stack");
imageCalculator("Multiply create stack", "Embryo.ome-2.tif", path+argssuff+".tif");
rename("Result of Embryo.ome-2.tif");
selectWindow("Result of Embryo.ome-2.tif");
//run("8-bit");
run("Invert LUT");
saveAs("Result of Embryo.ome-2.tif", path+argssuff+"_multiplied.tif");

// Find maxima
print("Doing 3D Maxima finding");
xyrad=file[4];
zrad=file[5];
nnoise=0;

side_chunks = 1;
for(i=0; i<side_chunks; i++){
  for(j=0; j<side_chunks; j++){
  print("(" + i + ", " + j + ")" + "...");
  selectWindow("Result of Embryo.ome-2.tif");
  makeRectangle(i*width/side_chunks, j*height/side_chunks, width/side_chunks, height/side_chunks);
  run("3D Maxima Finder", "radiusxy="+xyrad+"  radiusz="+zrad+" noise="+nnoise);
  saveAs("Results", path+argssuff+"_res"+i+j+".dat"); 
}

exit;
