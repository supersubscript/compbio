//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");
setMinAndMax(32935, 33265);
run("Apply LUT", "stack");
run("Remove Outliers...", "radius=30 threshold=0 which=Dark stack");
run("Remove Outliers...", "radius=30 threshold=0 which=Dark stack");
run("Remove Outliers...", "radius=1 threshold=0 which=Bright stack");
rename("N2");
saveAs("N2", "/home/henrik/n2_spotremoved.tif");
//run("Fire");
//run("TrackMate");

