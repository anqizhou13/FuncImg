macro "Break up a lif into individual TIFF" {
// open the file manager to select a lif file to break it into TIFFs
// in this case, only the metadata specific to a series will be written

path = File.openDialog("Select a File");

run("Bio-Formats Macro Extensions");
Ext.setId(path);
Ext.getCurrentFile(file);
Ext.getSeriesCount(seriesCount);

for (s=1; s<=seriesCount; s++) {
// Bio-Formats Importer uses an argument that can be built by concatenate a set of strings
run("Bio-Formats Importer", "open=&path autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_"+s);
run("8-bit");


imageName = getTitle();
newFileName = getTitle() + "_t";
newFileName = replace(newFileName, " ", "");
newFileName = replace(newFileName, ".lif", "");
dirNewFolder = "//Volumes//TOSHIBA//IMAGING//GCaMP_data//binary//R22E09_38H09LexA_LexAopGCaMP7s//processed_" + newFileName + "//" ; // TO CHANGE
File.makeDirectory(dirNewFolder); 

run("StackReg", "transformation=[Rigid Body]"); // adjust for xy shift
run("Fire"); 
setTool(3);                                          //Freehand
selectWindow(imageName);
beep();
waitForUser("Select region of interest");
getSelectionBounds(x, y, width, height);
run("Plot Z-axis Profile");
xpoints = newArray ();
ypoints = newArray ();
Plot.getValues (xpoints, ypoints);
for (i = 0; i < xpoints.length; i++) {
    setResult ("x", i, xpoints[i]);
    setResult ("y", i, ypoints[i]);
  }
  updateResults ();
saveAs("Results", dirNewFolder + "Values.csv");
selectWindow(imageName);
run("Draw", "stack");



setTool(0);                                          //Rectangle tool 
selectWindow(imageName);
beep();
waitForUser("Select background");
getSelectionBounds(x, y, width, height);
run("Plot Z-axis Profile");
xpoints = newArray ();
ypoints = newArray ();
Plot.getValues (xpoints, ypoints);
for (i = 0; i < xpoints.length; i++) {
    setResult ("x", i, xpoints[i]);
    setResult ("y", i, ypoints[i]);
  }
  updateResults ();
saveAs("Results", dirNewFolder + "ValuesBackground.csv");

selectWindow(imageName);
run("Draw", "stack");
newFileDir = dirNewFolder + newFileName + ".avi" ;
run("AVI... ", "compression=PNG frame=1 save=["+newFileDir+"]]");
run("Close All");
    }
}
