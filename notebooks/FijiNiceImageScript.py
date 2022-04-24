from ij import IJ, WindowManager
from ij.io import DirectoryChooser
from ij.gui import WaitForUserDialog, GenericDialog
from ij.measure import Calibration
import os


dc=DirectoryChooser("Choose a folder")
platefolder=dc.getDirectory()
if platefolder is None:
    print("User Canceled")

else:
  gd=GenericDialog(platefolder)
  gd.addRadioButtonGroup("Binning",["1","2"],1,2,"2")
  gd.addRadioButtonGroup("Microscope",["Phenix","Yokogawa"],1,2,"Phenix")
  gd.addRadioButtonGroup("Make montage in addition to single crop?",["Yes","No"],1,2,"Yes")
  gd.addNumericField("Number of fields in this folder",3,0)
  gd.addNumericField("Min final brightness",0,0)
  gd.addNumericField("Max final brightness",40000,0)
  gd.addStringField("Output folder", "/your/path/to/2022_Cimini_NatureProtocols/figures/image_crops", 100)
  gd.showDialog()  
  if gd.wasCanceled():
    print("User canceled dialog!")
  else:
    bin=int(gd.getNextRadioButton())
    scope=gd.getNextRadioButton().lower()
    do_montage = gd.getNextRadioButton()
    n_fields = int(gd.getNextNumber())
    minval=int(gd.getNextNumber())
    maxval=int(gd.getNextNumber())
    output_path=gd.getNextString()
    
    if scope == 'phenix':
        plate = platefolder.split('/')[-2]
        batch = platefolder.split('/')[-4]
    else:
        plate = platefolder.split('/')[-2]
        batch = platefolder.split('/')[-3]
    
    if scope == 'phenix':
        IJ.run("Image Sequence...", "open="+platefolder+"/Images sort")
    else:
        IJ.run("Image Sequence...", "open="+platefolder+" sort")
    im = IJ.getImage()
    calibration=Calibration()
    calibration.setUnit("micron")
    im.setCalibration(calibration)
    channels = im.getImageStackSize() / n_fields
    IJ.run("Stack to Hyperstack...", "order=xytcz channels="+str(channels)+" slices="+str(n_fields)+" frames=1 display=Composite")
    if bin == 2:
        pixel_size=".59797"
    else:
        if scope == 'phenix':
            pixel_size=".29898"
        else:
            pixel_size=".325"
    IJ.run("Properties...", "channels="+str(channels)+" slices="+str(n_fields)+" frames=1 pixel_width="+pixel_size+" pixel_height="+pixel_size+" voxel_depth=1")
    IJ.run("Make Montage...", "scale=1")
    im=IJ.getImage()
    IJ.run("Duplicate...", "duplicate channels=1-5 title=Composite")
    im=IJ.getImage()
    LUT_dict = {1:"Magenta",2:"Red", 3:"Yellow",4:"Green",5:"Cyan"}
    for eachkey in LUT_dict.keys():
        im.setSlice(eachkey)
        IJ.run(LUT_dict[eachkey])
    IJ.setMinAndMax(minval, maxval)
    if bin == 2:
        boxsize=200
    else:
        boxsize=400
    IJ.makeRectangle(0, 0, boxsize, boxsize)
    myWait = WaitForUserDialog ("", "Drag the box to a single representative field, and adjust contrast if needed (but try not to). The script will crop it, apply a scale bar, and save it.")
    myWait.show()
    IJ.run("Crop")
    im = IJ.getImage()
    if do_montage=='No':
        if bin ==2:
            IJ.run("Scale Bar...", "width=25 height=4 font=14 color=White background=None location=[Lower Right] bold")
        else:
            IJ.run("Scale Bar...", "width=25 height=8 font=28 color=White background=None location=[Lower Right] bold")
        IJ.saveAs("Tiff",os.path.join(output_path,batch+'-'+plate+'.tiff'))
    else:
        #start of Magic Montage code - https://wsr.imagej.net/macros/toolsets/Magic%20Montage.txt
        b=im.bitDepth
        IJ.newImage("tempmont", "RGB", boxsize, boxsize,6)
        for i in range(1,6):
            IJ.setPasteMode("copy")
            WindowManager.setTempCurrentImage(WindowManager.getImage("Composite"))	
            IJ.run("Duplicate...", "duplicate channels="+str(i)+" title=temp"+str(i))
            WindowManager.setTempCurrentImage(WindowManager.getImage("Composite-"+str(i)))
            IJ.run("RGB Color")
            IJ.run("Copy")
            WindowManager.setTempCurrentImage(WindowManager.getImage("tempmont"))
            im2=IJ.getImage()
            im2.setSlice(i)
            IJ.run("Paste")
        WindowManager.setTempCurrentImage(WindowManager.getImage("Composite"))
        IJ.run("RGB Color")
        IJ.run("Copy")
        WindowManager.setTempCurrentImage(WindowManager.getImage("tempmont"))
        im2=IJ.getImage()
        im2.setSlice(6)
        IJ.run("Paste")
        im2.setCalibration(calibration)
        IJ.run("Properties...", "channels=1 slices=6 frames=1 pixel_width="+pixel_size+" pixel_height="+pixel_size+" voxel_depth=1")
        IJ.run("Make Montage...", "columns=3 rows=2 scale=1 border=1")
        if bin ==2:
            IJ.run("Scale Bar...", "width=25 height=4 font=14 color=White background=None location=[Lower Right] bold")
        else:
            IJ.run("Scale Bar...", "width=25 height=8 font=28 color=White background=None location=[Lower Right] bold")
        IJ.saveAs("Tiff",os.path.join(output_path,batch+'-'+plate+'_montage.tiff'))
        WindowManager.setTempCurrentImage(WindowManager.getImage("Composite (RGB)"))
        if bin ==2:
            IJ.run("Scale Bar...", "width=25 height=4 font=14 color=White background=None location=[Lower Right] bold")
        else:
            IJ.run("Scale Bar...", "width=25 height=8 font=28 color=White background=None location=[Lower Right] bold")
        IJ.saveAs("Tiff",os.path.join(output_path,batch+'-'+plate+'.tiff'))
    IJ.run("Close All")