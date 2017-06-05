// CLV3 image processing + Costanza analysis
// 28 November 2014
// For use with pCLV3-dsRED-pWUS-GFP-OutNPA time series

run("Split Channels");
selectImage(1);
close();
selectImage(2);
close();
selectImage(1);
run("Fire");
run("Save");

run("Costanza Plugin");
close("Costanza Plugin");
selectWindow("Costanza - Basins of attractions (BOA)");
close();
selectWindow("Costanza - Cell centers");
close();
selectWindow("Costanza - Basins of attractions (BOA)-intensity");
run("Properties...", "unit=pixel pixel_width=0.1757176 pixel_height=0.1757176 voxel_depth=1.0000000");
run("Save");
selectWindow("Results");
run("Save");
close();

