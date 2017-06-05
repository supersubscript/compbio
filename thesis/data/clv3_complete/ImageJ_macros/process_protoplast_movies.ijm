// This is a basic Macro to open the timecourses in a certain folder obtained from the different Nikon and spinning disk microscopes. 
// This macro reads the different nd files in a folder to extract the information, so it automatically goes across the different generated positions, time points and z-stacks, if any. 

// The generetaed output is stored in a subfolder called "processed".

// The macro puts together the different time points in a multitiff file and shows stores it in Fire colour code. 
// For z-stacks, the macro collapses the z-stacks both summing the z and/or averaging its fluorescense.

// Last modification Monday 11th May 2015.


setBatchMode(true);

// Macro parameters 
quickloop=0; 		// This is a parameter to avoid certain scrip functions
channel_set=0; 		// This is to select the set of channels to process. Note different microscopes will have different channels. This can be improved it the future, and to be read from the nd files.
colum_stages_row=10; // Row number in the nd file where the number of positions is found (not sure why is 10 and not 9?)

// Macro input strings
path_in="G:\\Research\\Data\\Protoplasts\\";
path_out="M:\\teamJL_2\\Pau\\";
path_in="Y:\\"
path_out=path_in; // I want the path for the output files to be the same path for the input files. The output will be generated in the "processed" subfolder. 


// Name of the folder where data is stored. The different channel sets have to be adapted by the macro user.
datafoldername="2015-03-13_movie_protoplasts_clock_genefr";
datafoldername="2015-03-17_movie_oldprotoplasts_clock"; 
datafoldername="2015-03-11_movie_WT_GI_CCA1";
datafoldername="2015-03-09_movie_GI_Pau";
datafoldername="2015-03-13_movie_protoplasts_clock_genefr";


datafoldername="2015-02-02_Proto35S_Pau_spin";
datafoldername="2015-02-18_Protoplasts40x_tl_Pau";
datafoldername="2015-02-02_Proto35S_Pau";
datafoldername="2015-05-09_SpinDsk_Proto35SPGMM_Celldir_10x";

//datafoldername="2015-03-18_movie_oldprotoplasts"; 


fullpath_in=path_in+datafoldername+"\\";

prefullpath_out=path_out+datafoldername+"\\processed";
fullpath_out=prefullpath_out+"\\";

if (File.exists(prefullpath_out)==0)    
  {File.makeDirectory(prefullpath_out); 
  print ("folder created:", prefullpath_out );}

list=getFileList(fullpath_in);


// This loop reads and stores the nd files in the chosen data folder.
for (i=0; i<list.length; i++){
if (endsWith(list[i], ".nd"))
{write(list[i]);}
}

// This loop goes across each nd file and performs the data preprocessing.
for (i=0; i<list.length; i++) {

if (endsWith(list[i], ".nd")){

fullfiletxt=list[i];
pathfile=fullpath_in+fullfiletxt;
filestring=File.openAsString(pathfile); 
rows=split(filestring, "\n");
write(rows[colum_stages_row]);
columns=split(rows[colum_stages_row],",");
numpositions=columns[1];
print(numpositions); 
index2=indexOf(fullfiletxt,".");
filetxt=substring(fullfiletxt,0,index2);

action_preprocess(fullpath_in,fullpath_out,filetxt,numpositions,channel_set,quickloop);
};
}

setBatchMode(false);



function action_preprocess(fullpath_in,fullpath_out,filetxt,numpositions,channel_set,quickloop){

if (channel_set == 0)
    {
    channels=newArray("w1GFP Confocal - Benoit","w2DSRED Confocal");
    //channels=newArray("w2DSRED Confocal","w1GFP Confocal - Benoit");
    }
    else if (channel_set == 1) 
    channels=newArray("w1YFP Confocal","w2FM464 Confocal");
    else if (channel_set == 2) 
    channels=newArray("w1GFP Confocal FRAP","w2DSRED Confocal FRAP");

numchannels=channels.length;

if(quickloop==1)
	numchannels=1;
	
for (i = 0; i < numpositions; i++){
position=toString(i+1);

for (j = 0; j < numchannels; j++){
channel=channels[j];

filerootname=filetxt+"_"+channel+"_s"+position;
namefile=filerootname+"_t1.TIF";
relevantstring=filerootname+"_";

open(fullpath_in+namefile);
getDimensions(dummy, dummy, dummy, sliceCount, nFrames);
slices=sliceCount;
write("Z Slices "+toString(sliceCount));
close();

write("Working with "+filerootname);
run("Image Sequence...", "open=["+fullpath_in+namefile+"] file=["+relevantstring+"] sort");
write("Imported "+relevantstring);
// we need to know only how many frames there are
getDimensions(dummy, dummy, dummy, sliceCount, nFrames);

timepoints=sliceCount/slices;
write("timepts "+timepoints);
numslices=toString(slices);
numframes=toString(timepoints);
write("Z Slices "+toString(slices)+"Total frames "+toString(sliceCount)+" Timepoints "+toString(timepoints));
//run("Enhance Contrast", "saturated=0.35");


if (slices > 1){
run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices="+numslices+" frames="+numframes+" display=Color");
}

rawid = getImageID();
hyperstackid = getImageID();

if (slices > 1){
run("Z Project...", "projection=[Sum Slices] all");}

run("Fire");
if (slices > 1){
	saveAs("Tiff", fullpath_out+relevantstring+"SUMFIRE.tif");}
else 
{saveAs("Tiff", fullpath_out+relevantstring+"FIRE.tif");}

close();


if (slices > 1){
if (quickloop == 0){
	//run("Z Project...", "projection=[Max Intensity] all");
	selectImage(hyperstackid);
	run("Z Project...", "projection=[Max Intensity] all");
	run("Fire");
	saveAs("Tiff", fullpath_out+relevantstring+"MAXINTFIRE.tif");
close();}}

}
}
}