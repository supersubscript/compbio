
function action(input, output, filename) {
  open(input + filename);
	run("Flip Z");
	run("Slice Keeper", "first=1 last=25 increment=1");
	run("Flip Z");
	strlen = lengthOf(filename);
	newstr = substring(filename, 0, strlen - 4);
  saveAs("Tiff", output + newstr + "_top25.tif");
  close();
  close();

  open(input + filename);
	run("Flip Z");
	run("Slice Keeper", "first=1 last=50 increment=1");
	run("Flip Z");
	strlen = lengthOf(filename);
	newstr = substring(filename, 0, strlen - 4);
  saveAs("Tiff", output + newstr + "_top50.tif");
  close();
  close();

  open(input + filename);
	run("Flip Z");
	run("Slice Keeper", "first=1 last=75 increment=1");
	run("Flip Z");
	strlen = lengthOf(filename);
	newstr = substring(filename, 0, strlen - 4);
  saveAs("Tiff", output + newstr + "_top75.tif");
  close();
  close();

  open(input + filename);
	run("Flip Z");
	run("Slice Keeper", "first=1 last=100 increment=1");
	run("Flip Z");
	strlen = lengthOf(filename);
	newstr = substring(filename, 0, strlen - 4);
  saveAs("Tiff", output + newstr + "_top100.tif");
  close();
  close();
  
  open(input + filename);
	run("Flip Z");
	run("Slice Keeper", "first=1 last=125 increment=1");
	run("Flip Z");
	strlen = lengthOf(filename);
	newstr = substring(filename, 0, strlen - 4);
  saveAs("Tiff", output + newstr + "_top125.tif");
  close();
  close();

  open(input + filename);
	run("Flip Z");
	run("Slice Keeper", "first=1 last=150 increment=1");
	run("Flip Z");
	strlen = lengthOf(filename);
	newstr = substring(filename, 0, strlen - 4);
  saveAs("Tiff", output + newstr + "_top150.tif");
  close();
  close();
}

input = "/home/henrik/to_slice/";
output = "/home/henrik/sliced/";

setBatchMode(true); 
list = getFileList(input);
for (i = 0; i < list.length; i++)
        action(input, output, list[i]);
setBatchMode(false);

