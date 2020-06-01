/***===============================================================================
 
 SpermQ-MF Version v0.3.0
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation (http://www.gnu.org/licenses/gpl.txt )

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 See the GNU General Public License for more details.
 
 @author Jan N Hansen
 @copyright (C) 2013 - 2020: Jan N Hansen and Jan F Jikeli
   
 For any questions please feel free to contact me (jan.hansen@uni-bonn.de).

==============================================================================**/
package spermQ_mf;

import java.awt.event.*;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;

import javax.swing.UIManager;
import javax.swing.filechooser.FileSystemView;

import ij.*;
import ij.gui.*;
import ij.io.*;
import ij.measure.*;
import ij.plugin.*;
import ij.plugin.frame.RoiManager;
import ij.text.TextPanel;
import spermQ_mf.jnhsupport.*;

public class multi_focal implements PlugIn, Measurements{
	//Name
		public static final String PLUGINNAME = "SpermQ-MF";
		public static final String PLUGINVERSION = "v0.3.0";
		public static final double [] PLANEPOSITIONS20X = {31.00735294, 22.23676471, 9.266176471, 17.48970588};
		public static final double [] PLANEPOSITIONS32X = {0.0, 3.5, 5, 8.5};
	//Default Settings loader
		String [] selectionsDSL = {"Mouse 20x", "Mouse 20x (editable)", 
				"Mouse 32x", "Mouse 32x (editable)",
				"Human 20x", "Human 20x (editable)",
				"Human 32x", "Human 32x (editable)",
				"Fluorescent Tethered Human 20x", "Fluorescent Tethered Human 20x (editable)",
				"Fluorescent Tethered Human 32x", "Fluorescent Tethered Human 32x (editable)"}; 
		String selectedDSL = selectionsDSL [5];
		
	//variables		
		double [] slicePosition = PLANEPOSITIONS20X;
		
		boolean calibrationMode = false;
		boolean calibrationTestMode = false;
		boolean orientation3D = false;
		
		double xyCal = 0.55;	//20x (1.6x compared to xyCal 32x)
		
		public static final String [] TRACEDETERMINATION = {"sharpest plane", "maximum-intensity-projection", "average-intensity-projection", "sum of planes", "sharpest plane in time projection"};
		String selectedTraceDetermination = TRACEDETERMINATION [4];
		static final String [] thresholdMethods = {"Default", "IJ_IsoData", "Huang", "Intermodes", "IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"};
		String selectedThresholdMethod = thresholdMethods [15];	//Triangle
		
		int upscaleFold = 3;
		double gaussSigma = 2.0;
		int maxVectorLength = 20;
		double normalLength = 24.0;	//48 for width calibration (20x)
		boolean smoothNormal = true,
				useNormalMaxForZGauss = true;
		boolean getWidthFitZ = true;
		boolean getG4PZ = true;
		double LUTStepSize = 0.1;
		
		double acceptedZSmoothDistance = 9.6;
		int plusMinusDistanceForSmooth = 5*upscaleFold;
		double minRefDist = 0.0, maxRefDist = 6.4;
		static final String [] zSmoothingMethods = {"none", "median", "mean"}; 
		String zSmoothingMethod = zSmoothingMethods [1];
		
		boolean tethered = false,
				unifyStartPoints = false,
				saveVNRois = false, 
				addCOM = false;
		
		boolean preventHeadFromCorr = true;		
		int	preventPoints = 10;
		int hrPlusMinusRange = 10;
		double curvRefDist = 10.0;
		
		double sampleRate = 500.0; 
		int groupedTimesteps = 200;
		double neglectedInitialArclength = 20.0;
		
	//dialog	
		boolean done = false;
		ProgressDialog progress;
	
	//calibration algorithm dialog
		int minSharpest = 300, maxSharpest = 400;
		double zStepSize = 0.1;
		
	@Override
	public void run(String arg) {
		/**&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		Load Default Settings
		&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
		
		GenericDialog gdDSL = new GenericDialog(PLUGINNAME + " - Default Settings Loader");		
		gdDSL.setInsets(0,0,0);		gdDSL.addMessage(PLUGINNAME + ", version " + PLUGINVERSION + " (\u00a9 2013-" + constants.dateY.format(new Date()) + ", JN Hansen \u0026 JF Jikeli)", constants.Head1);
		gdDSL.setInsets(10,0,0);	gdDSL.addMessage("Default Settings Loader ", constants.Head2);
		gdDSL.setInsets(10,0,0);	gdDSL.addMessage("Experimental setup ", constants.BoldTxt);
		gdDSL.setInsets(0,0,0);		gdDSL.addNumericField("Sample rate [Hz]", sampleRate, 2);
		gdDSL.setInsets(0, 0, 0);	gdDSL.addCheckbox("Sperm head is tethered", tethered);
		gdDSL.setInsets(0, 0, 0);	gdDSL.addCheckbox("Orient points in 3D", orientation3D);
		
		gdDSL.setInsets(10,0,0);	gdDSL.addMessage("Experimental setup ", constants.BoldTxt);
		gdDSL.setInsets(0,0,0);		gdDSL.addChoice("Start with default settings for ", selectionsDSL, selectedDSL);
		gdDSL.setInsets(0,0,0);		gdDSL.addMessage("If you select <editable> settings, you can adjust parameters in the next step.", constants.PlTxt);
		
		gdDSL.setInsets(10,0,0);	gdDSL.addMessage("Special running modes ", constants.BoldTxt);
		gdDSL.setInsets(0,0,0);		gdDSL.addCheckbox("Width calibration mode", calibrationMode);
		gdDSL.setInsets(0,0,0);		gdDSL.addCheckbox("Alternatively: test calibration mode", calibrationTestMode);
		
		gdDSL.setInsets(20,0,0);	gdDSL.addMessage("IMPORTANT NOTE: If you analyze stacks containing lots of frames (e.g. > 4000), your RAM capacity", constants.PlTxt);
		gdDSL.setInsets(0,0,0);		gdDSL.addMessage("may be exceeded and the plugin stops working. Thus, it is recommended to split the stack into", constants.PlTxt);
		gdDSL.setInsets(0,0,0);		gdDSL.addMessage("seperate sub-stacks and analyze those sub-stacks separately.", constants.PlTxt);
		
		gdDSL.showDialog();
		
		sampleRate = gdDSL.getNextNumber();
		tethered = gdDSL.getNextBoolean();
		orientation3D = gdDSL.getNextBoolean();
		
	 	selectedDSL = gdDSL.getNextChoice();

	 	calibrationMode = gdDSL.getNextBoolean();
	 	calibrationTestMode = gdDSL.getNextBoolean();
	 	
		if (gdDSL.wasCanceled())return;
		
//		IJ.log("got");
		groupedTimesteps = (int)sampleRate;
		double speciesLength = multi_focal_tools.SPECIESLENGTH_HUMAN;
		if(selectedDSL.equals(selectionsDSL[0]) || selectedDSL.equals(selectionsDSL[1])){
			//mouse 20x
			slicePosition = PLANEPOSITIONS20X;
			xyCal = 0.55;
			unifyStartPoints = false;
			selectedThresholdMethod = thresholdMethods [5];	//Li (mouse)
			upscaleFold = 3;
			addCOM = false;
			maxVectorLength = 20;
			normalLength = 5.0;
			smoothNormal = true;
			saveVNRois = false;
			preventHeadFromCorr = true;
			preventPoints = 10;
			plusMinusDistanceForSmooth = 5*upscaleFold;
			minRefDist = 0.0;
			maxRefDist = 15;
			acceptedZSmoothDistance = 9.6;
			zSmoothingMethod = zSmoothingMethods [1];
			curvRefDist = 10.0;
			neglectedInitialArclength = 20.0;
			hrPlusMinusRange = 10;
			gaussSigma = 2.0;
			speciesLength = multi_focal_tools.SPECIESLENGTH_MOUSE;
			selectedTraceDetermination = TRACEDETERMINATION [4];
		}
		else if(selectedDSL.equals(selectionsDSL[2]) || selectedDSL.equals(selectionsDSL[3])){
			//mouse 32x
			slicePosition = PLANEPOSITIONS32X; 
			xyCal = 0.34375;
			unifyStartPoints = false;
			selectedThresholdMethod = thresholdMethods [5];	//Li (mouse)
			upscaleFold = 3;
			addCOM = false;
			maxVectorLength = 30;
			normalLength = 5.0;
			smoothNormal = true;
			saveVNRois = false;
			preventHeadFromCorr = true;
			preventPoints = 15;
			plusMinusDistanceForSmooth = 7*upscaleFold;
			minRefDist = 0.0;
			maxRefDist = 15;
			acceptedZSmoothDistance = 9.6;
			zSmoothingMethod = zSmoothingMethods [1];
			curvRefDist = 10.0;
			neglectedInitialArclength = 20.0;
			hrPlusMinusRange = 15;
			gaussSigma = 2.0;
			speciesLength = multi_focal_tools.SPECIESLENGTH_MOUSE;
			selectedTraceDetermination = TRACEDETERMINATION [4];
		}
		else if(selectedDSL.equals(selectionsDSL[4]) || selectedDSL.equals(selectionsDSL[5])){
			//human 20x
			slicePosition = PLANEPOSITIONS20X;
			xyCal = 0.55;	//20x
			unifyStartPoints = false;
			selectedThresholdMethod = thresholdMethods [15];	//Triangle (human)
			upscaleFold = 3;
			addCOM = true;
			maxVectorLength = 20;
			normalLength = 24.0;	// 36 for width calibration of human sperm (20x)
			smoothNormal = true;
			saveVNRois = false;
			preventHeadFromCorr = false;
			preventPoints = 15;			
			plusMinusDistanceForSmooth = 5*upscaleFold;
			minRefDist = 4.0;
			maxRefDist = 10.0;
			acceptedZSmoothDistance = 9.6;
			zSmoothingMethod = zSmoothingMethods [1];
			curvRefDist = 10.0;
			neglectedInitialArclength = 0.0;
			hrPlusMinusRange = 10;
			gaussSigma = 2.0;
			speciesLength = multi_focal_tools.SPECIESLENGTH_HUMAN;
			selectedTraceDetermination = TRACEDETERMINATION [4];
		}
		else if(selectedDSL.equals(selectionsDSL[6]) || selectedDSL.equals(selectionsDSL[7])){
			//human 32x
			slicePosition = PLANEPOSITIONS32X;
			xyCal = 34375;	//32x
			unifyStartPoints = false;
			selectedThresholdMethod = thresholdMethods [15];	//Triangle (human)
			upscaleFold = 3;
			addCOM = true;
			maxVectorLength = 30;
			normalLength = 24.0;	// 36 for width calibration of human sperm (20x)
			smoothNormal = true;
			saveVNRois = false;
			preventHeadFromCorr = false;
			preventPoints = 15;		
			plusMinusDistanceForSmooth = 5*upscaleFold;
			minRefDist = 0.0;
			maxRefDist = 6.4;
			acceptedZSmoothDistance = 9.6;
			zSmoothingMethod = zSmoothingMethods [1];
			curvRefDist = 10.0;
			neglectedInitialArclength = 0.0;
			hrPlusMinusRange = 15;
			gaussSigma = 2.0;
			speciesLength = multi_focal_tools.SPECIESLENGTH_HUMAN;
			selectedTraceDetermination = TRACEDETERMINATION [4];
		}
		else if(selectedDSL.equals(selectionsDSL[8]) || selectedDSL.equals(selectionsDSL[9])){
			//human 20x calcium
			slicePosition = PLANEPOSITIONS20X;
			xyCal = 0.55;	//20x
			unifyStartPoints = true;
			selectedThresholdMethod = thresholdMethods [15];	//Triangle (human)
			upscaleFold = 3;
			addCOM = true;
			maxVectorLength = 20;
			normalLength = 24.0;	// 36 for width calibration of human sperm (20x)
			smoothNormal = true;
			saveVNRois = false;
			preventHeadFromCorr = false;
			preventPoints = 15;			
			plusMinusDistanceForSmooth = 5*upscaleFold;
			minRefDist = 0.0;
			maxRefDist = 6.4;
			acceptedZSmoothDistance = 9.6;
			zSmoothingMethod = zSmoothingMethods [1];
			curvRefDist = 10.0;
			neglectedInitialArclength = 0.0;
			hrPlusMinusRange = 15;
			gaussSigma = 4.0;
			speciesLength = multi_focal_tools.SPECIESLENGTH_HUMAN;
			selectedTraceDetermination = TRACEDETERMINATION [4];
		}
		else if(selectedDSL.equals(selectionsDSL[10]) || selectedDSL.equals(selectionsDSL[11])){
			//human 32x calcium
			slicePosition = PLANEPOSITIONS32X;
			xyCal = 34375;	//32x
			unifyStartPoints = true;
			selectedThresholdMethod = thresholdMethods [15];	//Triangle (human)
			upscaleFold = 3;
			addCOM = true;
			maxVectorLength = 20;
			normalLength = 24.0;	// 36 for width calibration of human sperm (20x)
			smoothNormal = true;
			saveVNRois = false;
			preventHeadFromCorr = false;
			preventPoints = 15;
			plusMinusDistanceForSmooth = 5*upscaleFold;
			minRefDist = 0.0;
			maxRefDist = 6.4;
			acceptedZSmoothDistance = 9.6;
			zSmoothingMethod = zSmoothingMethods [1];
			curvRefDist = 10.0;
			neglectedInitialArclength = 0.0;
			hrPlusMinusRange = 15;
			gaussSigma = 4.0;
			speciesLength = multi_focal_tools.SPECIESLENGTH_MOUSE;
			selectedTraceDetermination = TRACEDETERMINATION [4];
		}
		
		
		/**&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		GenericDialog
		&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
		if(selectedDSL.equals(selectionsDSL[1]) || selectedDSL.equals(selectionsDSL[3]) || selectedDSL.equals(selectionsDSL[5])
				|| selectedDSL.equals(selectionsDSL[7]) || selectedDSL.equals(selectionsDSL[9]) || selectedDSL.equals(selectionsDSL[11])){
			GenericDialog gd = new GenericDialog(PLUGINNAME + " - detailed settings");		
//			setInsets(top, left, bottom)
			gd.setInsets(0,0,0);	gd.addMessage(PLUGINNAME + ", version " + PLUGINVERSION + " (\u00a9 2013-" + constants.dateY.format(new Date()) + ", JN Hansen \u0026 JF Jikeli)", constants.Head1);
						
			gd.setInsets(10,0,0);	gd.addNumericField("xy calibration [um]", xyCal, 5);
			gd.setInsets(0,0,0);	gd.addNumericField("Slice Positions [um]: 1 | 2", slicePosition [0], 5);
			gd.setInsets(-23,55,0);	gd.addNumericField("", slicePosition [1], 5);
			gd.setInsets(0,0,0);	gd.addNumericField("3 | 4", slicePosition [2], 5);
			gd.setInsets(-23,55,0);	gd.addNumericField("", slicePosition [3], 5);
			
			gd.setInsets(10,0,0);	gd.addMessage("Trace generation ", constants.BoldTxt);
			gd.setInsets(0,0,0);	gd.addChoice("Trace determination in ", TRACEDETERMINATION, selectedTraceDetermination);
			gd.setInsets(0,0,0);	gd.addChoice("Thresholding Method", thresholdMethods, selectedThresholdMethod);
			gd.setInsets(0,0,0);	gd.addNumericField("Gauss sigma (defines size of detected objects)", gaussSigma, 2);
			gd.setInsets(0,0,0);	gd.addNumericField("Upscaling of points (fold)", upscaleFold, 0);
			gd.setInsets(0,0,0);	gd.addCheckbox("Add head center-of-mass as first point", addCOM);
			gd.setInsets(0,0,0);	gd.addCheckbox("Unify start points (for tethered sperm only!)", unifyStartPoints);
			
			gd.setInsets(10,0,0);	gd.addMessage("Trace improvements ", constants.BoldTxt);
			gd.setInsets(0,0,0);	gd.addNumericField("Maximum vector length (points)", maxVectorLength, 0);
			gd.setInsets(0,0,0);	gd.addNumericField("Normal radius for gauss fit [um]", normalLength, 2);
			gd.setInsets(0,0,0);	gd.addCheckbox("Smooth normal for XY gauss fit and Z gauss fit", smoothNormal);
//			gd.setInsets(0,0,0);	gd.addCheckbox("Use normal maximum for z gauss fit", useNormalMaxForZGauss);
			gd.setInsets(0,0,0);	gd.addCheckbox("Exclude head from correction/deletion/Zkymographs (initial (" + preventPoints + " * upscaling factor) points)", preventHeadFromCorr);
			gd.setInsets(0,0,0);	gd.addCheckbox("Save Roi-sets of vectors and normals", saveVNRois);
			
			gd.setInsets(10,0,0);	gd.addMessage("z-determination, smoothing, kymographs", constants.BoldTxt);
			gd.setInsets(0,0,0);	gd.addCheckbox("Use Gauss-4-Point fit for z determination:", getG4PZ);
			gd.setInsets(0,0,0);	gd.addCheckbox("Use width for z-determination: LUT step size [um]", getWidthFitZ);
			gd.setInsets(-23,0,0);	gd.addNumericField("", LUTStepSize, 2);
			gd.setInsets(0,0,0);	gd.addChoice("z-smoothing ", zSmoothingMethods, zSmoothingMethod);
			gd.setInsets(0,0,0);	gd.addNumericField("Accepted xy distance of points for z smoothing [um]", acceptedZSmoothDistance, 5);
			gd.setInsets(0,0,0);	gd.addNumericField("# (+/-)-consecutive points for xy- and z-smoothing", plusMinusDistanceForSmooth, 0);
			gd.setInsets(0,0,0);	gd.addNumericField("Min | Max xy-arcus position for reference vector [um]", minRefDist, 3);
			gd.setInsets(-23,55,0);	gd.addNumericField("", maxRefDist, 3);

			gd.setInsets(10,0,0);	gd.addMessage("Additional calculations", constants.BoldTxt);
			gd.setInsets(0,0,0);	gd.addNumericField("Curvature: reference point distance", curvRefDist, 4);
			gd.setInsets(0,0,0);	gd.addNumericField("FFT: Grouped consecutive time-steps", groupedTimesteps, 0);
			gd.setInsets(0,0,0);	gd.addNumericField("FFT: Do not analyze initial ... µm from head", neglectedInitialArclength, 0);
			gd.setInsets(0,0,0);	gd.addNumericField("Head rotation matrix radius", hrPlusMinusRange, 0);
					
			gd.showDialog();
				 	
		 	xyCal = gd.getNextNumber();
		 	
			slicePosition [0] = (double) gd.getNextNumber();
			slicePosition [1] = (double) gd.getNextNumber();
			slicePosition [2] = (double) gd.getNextNumber();
			slicePosition [3] = (double) gd.getNextNumber();
			
			selectedTraceDetermination = gd.getNextChoice();
		 	selectedThresholdMethod = gd.getNextChoice();
		 	gaussSigma = (double) gd.getNextNumber();
			upscaleFold = (int) gd.getNextNumber();
			addCOM = gd.getNextBoolean();
			unifyStartPoints = gd.getNextBoolean();

			maxVectorLength = (int) gd.getNextNumber();
			normalLength = (double) gd.getNextNumber();
			smoothNormal = gd.getNextBoolean();
//			useNormalMaxForZGauss = gd.getNextBoolean();
			preventHeadFromCorr = gd.getNextBoolean();
			saveVNRois = gd.getNextBoolean();
			
			getG4PZ = gd.getNextBoolean();
			getWidthFitZ = gd.getNextBoolean();
			LUTStepSize = (double) gd.getNextNumber();
			
			zSmoothingMethod = gd.getNextChoice();
			acceptedZSmoothDistance = (double) gd.getNextNumber();
			plusMinusDistanceForSmooth = (int) gd.getNextNumber();
			minRefDist = (double) gd.getNextNumber();
			maxRefDist = (double) gd.getNextNumber();
			
			curvRefDist = (double) gd.getNextNumber();
			groupedTimesteps = (int) gd.getNextNumber();
			neglectedInitialArclength = (double) gd.getNextNumber();
			
			hrPlusMinusRange = (int) hrPlusMinusRange;
			
			if (gd.wasCanceled())return;
		}
				
//		IJ.log("slices not sort: " + slicePosition [0] + " " + slicePosition [1] + " " + slicePosition [2] + " " + slicePosition [3]);
		double [] slicesSorted = slicePosition.clone();
	  	Arrays.sort(slicesSorted);
//	  	IJ.log("slices not sort: " + slicePosition [0] + " " + slicePosition [1] + " " + slicePosition [2] + " " + slicePosition [3]);
//	  	IJ.log("slices sort: " + slicesSorted [0] + " " + slicesSorted [1] + " " + slicesSorted [2] + " " + slicesSorted [3]);
	  	
	  	preventPoints *= upscaleFold;
		if(preventHeadFromCorr == false){
			preventPoints = 0;
		}
			  	
	  	int encoding = 0;
  		if(zSmoothingMethod.equals(zSmoothingMethods[0])){	//NONE
  			encoding = multi_focal_tools.PUREZ;					  			  		
	  	}else if(zSmoothingMethod.equals(zSmoothingMethods[1])){	//MEDIAN
	  		encoding = multi_focal_tools.MEDIANZ;	
	  	}else if(zSmoothingMethod.equals(zSmoothingMethods[2])){	//MEAN
	  		encoding = multi_focal_tools.MEANZ;	
	  	}
				
  		/**&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		Initiate multi task management
		&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
		try{
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		}catch(Exception e){}
		
		//image selection
		OpenFilesDialog od = new OpenFilesDialog ();
		od.setLocation(0,0);
		od.setVisible(true);
		
		od.addWindowListener(new java.awt.event.WindowAdapter() {
	        public void windowClosing(WindowEvent winEvt) {
//	        	IJ.log("Analysis canceled!");
	        	return;
	        }
	    });
	
		while(od.done==false){
			 try{
				 Thread.currentThread().sleep(50);
		     }catch(Exception e){
		     }
		}
		
		int tasks = od.filesToOpen.size();
		String [] name = new String [tasks];
		String [] dir = new String [tasks];
		boolean tasksSuccessfull [] = new boolean [tasks];
		for(int task = 0; task < tasks; task++){
			name [task] = od.filesToOpen.get(task).getName();
			dir [task] = od.filesToOpen.get(task).getParent() + System.getProperty("file.separator");
			tasksSuccessfull [task] = false;
		}	
		
		
		//start progress dialog
		progress = new ProgressDialog(name, tasks);
		progress.setVisible(true);
		progress.addWindowListener(new java.awt.event.WindowAdapter() {
	        public void windowClosing(WindowEvent winEvt) {
	        	progress.stopProcessing();
	        	if(done==false){
	        		IJ.error("Script stopped...");
	        	}       	
	        	System.gc();
	        	return;
	        }
		});	

		
		/**&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		Processing
		&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

		// open width fit look up table
	    	OpenDialog oLUT;
	    	String nameLUT = "";
	    	String dirLUT = "";
	    	if(getWidthFitZ && !calibrationMode){
	    		oLUT = new OpenDialog("Open width fit LUT file", null);
	    		progress.replaceBarText("Open a width-LUT file");
	    		nameLUT = oLUT.getFileName();
		    	dirLUT = oLUT.getDirectory();
	    	}
		
		//Initialize
			ImagePlus imp;		
			String homePath = FileSystemView.getFileSystemView().getHomeDirectory().getAbsolutePath();
			
			//get head selections
			Roi [] selections = new Roi [tasks];
			Roi [] selectionsSD = new Roi [tasks];
			IJ.setTool("polygon");
			{
				ImagePlus maxImp;
				for(int task = 0; task < tasks; task++){
					imp = IJ.openVirtual(dir [task] + name [task]);
					IJ.run(imp, "Z Project...", "projection=[Max Intensity] all");
					imp.changes = false;
					imp.close();
					imp = WindowManager.getCurrentImage();
					imp.hide();
					IJ.run(imp, "Z Project...", "projection=[Max Intensity]");
					maxImp = WindowManager.getCurrentImage();
					imp.changes = false;
					imp.close();
									
					while(true){
						progress.replaceBarText("user interaction required... [task " + (task+1) + "/" + tasks + "]");
						new WaitForUserDialog("Set a Roi containing parts of the cell in every frame [task " + (task+1) + "/" + tasks + "]").show();
						if(maxImp.getRoi()!=null) break;
					}		
					selections [task] = new PolygonRoi(maxImp.getRoi().getPolygon(), PolygonRoi.POLYGON);
					
					if(selectedTraceDetermination.equals(TRACEDETERMINATION[0]) || selectedTraceDetermination.equals(TRACEDETERMINATION[4])){
						maxImp.deleteRoi();
						while(true){
							progress.replaceBarText("user interaction required... [task " + (task+1) + "/" + tasks + "]");
							new WaitForUserDialog("Set a Roi in which the SD shall be calculated for finding the sharpest plane [task " + (task+1) + "/" + tasks + "]").show();
							if(maxImp.getRoi()!=null) break;
						}		
						selectionsSD [task] = new PolygonRoi(maxImp.getRoi().getPolygon(), PolygonRoi.POLYGON);
					}
					
					
					maxImp.changes = false;
					maxImp.close();
					System.gc();
				}
			}
			System.gc();
		
		//Initialize Variables
		double [] widthZCorrFactor = new double [10];
    	String saveName, savePath;
    	double [][] freqSummaryTheta2D = new double [tasks][2];
    	double [][] freqSummaryTheta3D = new double [tasks][2];
    	double [][] freqSummaryHrMaxPos = new double [tasks][2];
    	double [][] freqSummaryHrMaxInt = new double [tasks][2];
    	double [][] freqSummaryHrAngle = new double [tasks][2];
		double [][][] freqSummaryX = new double [tasks][2][3];
		double [][][] freqSummaryY = new double [tasks][2][3];
		double [][][] freqSummaryZ = new double [tasks][2][3];
		double [][][] freqSummaryCurv2D = new double [tasks][2][3];
		double [][][] freqSummaryCurv3D = new double [tasks][2][3];
		double [][][] freqSummaryCAngle2D = new double [tasks][2][3];
		double [][][] freqSummaryCAngle3D = new double [tasks][2][3];
		
		ArrayList <trace> traces = null;
		RoiManager rm;		
		double [] minIntensity = {Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY},
	  			maxIntensity = {0.0,0.0,0.0,0.0}, 
	  			minPosition = {Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY},
	  			maxPosition = {0.0,0.0,0.0,0.0};
		Date startDate;
		TextPanel tp;
		
		//processing
		tasking: for(int task = 0; task < tasks; task++){
			progress.updateBarText("in progress...");
			startDate = new Date();
			progress.clearLog();
			
			running: while(true){
				// get image and open
			    	saveName = name [task].substring(0,name [task].lastIndexOf(".")) + "_mfa_" + constants.dateName.format(startDate);
			    	new File(dir [task]+ saveName).mkdirs();
			    	imp = IJ.openImage(dir [task] + name [task]);	
			    	imp.getCalibration().pixelHeight = xyCal;
			    	imp.getCalibration().pixelWidth = xyCal;
			    	
			    	savePath = dir [task] + saveName + System.getProperty("file.separator") + saveName;			    	
		    	// get image and open (end)		
				
				//Check if image ist processable
				if(imp.getBitDepth()==24){
		    		progress.notifyMessage("ERROR: Multi focal analysis cannot be used for RGB images!", ProgressDialog.ERROR);
			  		break running;
		    	}
				  	
				 //get traces
				  	progress.updateBarText("get object traces...");			 
				  	if(calibrationTestMode){
				  		traces = multi_focal_tools.getObjectTraces(imp, selectedThresholdMethod, gaussSigma, 
				  				progress, selections [task], selectionsSD [task], true, selectedTraceDetermination, addCOM,
				  				minRefDist, maxRefDist);
				  	}else{
				  		traces = multi_focal_tools.getObjectTraces(imp, selectedThresholdMethod, gaussSigma, 
				  				progress, selections [task], selectionsSD [task], calibrationMode, selectedTraceDetermination, addCOM,
				  				minRefDist, maxRefDist);
				  	}			  	
				  	if(progress.isStopped())break running;
				  	
//					  	multi_focal_tools.saveTraceImage(imp, traces, multi_focal_tools.NOZ,
//				  				slicesSorted [0], slicesSorted[slicesSorted.length-1],savePath);
				  	System.gc();
				  	
				  	
				//search for reversed traces (if tethered sperm)
				  	if(tethered){
				  		multi_focal_tools.reverseReversedTraces(traces, progress);
				  		if(progress.isStopped())break running;
				  		if(unifyStartPoints){
					  		multi_focal_tools.unifyStartPoints(traces, progress);
					  		if(progress.isStopped())break running;
					  	}
				  	}					  	
				
				//upscale point trace
				  	progress.updateBarText("upscale traces...");
				  	for(int i = 0; i < traces.size(); i++){
				  		traces.get(i).upscalePoints(upscaleFold);
				  	}
				  	
//					  	multi_focal_tools.saveTraceImage(imp, traces, multi_focal_tools.NOZ,
//				  				slicesSorted [0], slicesSorted[slicesSorted.length-1],savePath + "PU");
				  	System.gc();
				
				//improve x,y and calculate xy-gaussWidth
				  	if(!calibrationMode && !calibrationTestMode){
				  		progress.updateBarText("improve traces...");
				  			multi_focal_tools.adjustPointsViaNormalVector(imp, traces, progress, saveVNRois,
				  					dir [task] + saveName + System.getProperty("file.separator"), maxVectorLength, normalLength, smoothNormal,
				  					preventHeadFromCorr, preventPoints);
					  	if(progress.isStopped())break running;
					  	
//						  	multi_focal_tools.saveTraceImage(imp, traces, multi_focal_tools.NOZ,
//					  				slicesSorted [0], slicesSorted[slicesSorted.length-1],savePath + "pfi");
					  	System.gc();
				  	}
				
				//linear interpolate xy
				  	if(!calibrationMode && !calibrationTestMode){
				  		for(int i = 0; i < traces.size(); i++){
				  			progress.updateBarText("resorting points ... (" + (i+1) + "/" + traces.size() + ")");
				  			traces.get(i).sortPointsAndRemoveOutliers(progress, multi_focal_tools.NOZ);			  			
				  		}
//					  		multi_focal_tools.saveTraceImage(imp, traces, multi_focal_tools.NOZ,
//					  				slicesSorted [0], slicesSorted[slicesSorted.length-1],savePath + "pfiS1");
				  		
				  	 	
				  	 	for(int i = 0; i < traces.size(); i++){
				  	 		progress.updateBarText("interpolate xy ... (" + (i+1) + "/" + traces.size());
				  	 		traces.get(i).interpolateXYLinear(plusMinusDistanceForSmooth, multi_focal_tools.MEAN);
				  		}
			  	 		
//				  	 		multi_focal_tools.saveTraceImage(imp, traces, multi_focal_tools.NOZ,
//					  				slicesSorted [0], slicesSorted[slicesSorted.length-1],savePath + "pfiIp");
			  	 		
					  	for(int i = 0; i < traces.size(); i++){
					  		progress.updateBarText("resorting points ... (" + (i+1) + "/" + traces.size() + ")");
				  			traces.get(i).sortPointsAndRemoveOutliers(progress ,multi_focal_tools.NOZ);
				  			
				  		}
				  		multi_focal_tools.saveTraceImage(imp, traces, multi_focal_tools.NOZ,
				  			slicesSorted [0], slicesSorted[slicesSorted.length-1], savePath + "", preventPoints);			  		
				  		System.gc();
				  	}
				 
				  
				//obtain better width fit data (no xy correction)
					if(!calibrationMode){
						progress.updateBarText("get width fit data...");
						multi_focal_tools.updateWidthFitData(imp, traces, progress, saveVNRois, dir [task] + saveName + System.getProperty("file.separator"),
								maxVectorLength, normalLength, smoothNormal, preventHeadFromCorr, preventPoints);
						if(progress.isStopped())break running;
						
						System.gc();
					}
			  				
				//calculate arc lengths
				  	progress.updateBarText("calculate arc lengths...");
				  	for(int i = 0; i < traces.size(); i++){
				  		progress.updateBarText("calculate arc lengths... " + (i+1) + "/" + traces.size());
				  		traces.get(i).calculateArcLengthPositions();
				  	}
				  	if(progress.isStopped())break running;  	
				  	
				//in calibration mode: remove non-fit points post arc length calculation
				  	if(calibrationMode){
				  		progress.updateBarText("get width fit data...");
						multi_focal_tools.updateWidthFitData(imp, traces, progress, saveVNRois, dir [task] + saveName + System.getProperty("file.separator"),
								maxVectorLength, normalLength, smoothNormal, preventHeadFromCorr, preventPoints);
						if(progress.isStopped())break running;
				  	}
				 	
				//save width results and trace
				  	multi_focal_tools.saveXYWidthGraph(traces, (1.0 / xyCal), dir [task] + saveName + System.getProperty("file.separator") + saveName);
//					  	multi_focal_tools.saveTraceImage(imp, traces, multi_focal_tools.NOZ, 0.0, 0.0,savePath);
				  	multi_focal_tools.saveGaussMaxGraphs(traces, xyCal, dir [task] + saveName + System.getProperty("file.separator") + saveName);
			  	
				  	System.gc();
				  				
				  	if(!calibrationMode){
				  		//z Gaussfit + filtering
				  			progress.updateBarText("find out z...");
						  	if(getWidthFitZ && getG4PZ){
						  		multi_focal_tools.determineZInfoIncWidth(imp, traces, slicePosition, progress, useNormalMaxForZGauss, dirLUT+nameLUT, LUTStepSize);
						  	}else if(getWidthFitZ && !getG4PZ){
						  		multi_focal_tools.determineZInfoInWidthsOnly(imp, traces, slicePosition, progress, dirLUT+nameLUT, LUTStepSize);
						  	}else if(!getWidthFitZ && getG4PZ){
						  		multi_focal_tools.determineZInfo(imp, traces, slicePosition, progress, useNormalMaxForZGauss);
						  	}					  	
						  	if(progress.isStopped())break running;
					  	
						//calibrate width fit and z fit
						  	if(getWidthFitZ && getG4PZ){
						  		widthZCorrFactor = multi_focal_tools.getWidthZCorrectionFactor(traces, progress);
							  	System.gc();
						  	}					  	
						  	
					  	//interpolate z fit
					  		progress.updateBarText("interpolate z...");
					  		multi_focal_tools.interpolateZLinear(traces, progress, acceptedZSmoothDistance, plusMinusDistanceForSmooth);
		//			  		for(int i = 0; i < traces.size(); i++){
		//			  			traces.get(i).sortPointsAndRemoveOutliers(progress, multi_focal_tools.G4PZMEDIAN);
		//			  		}			  		
					  		if(progress.isStopped())break running;
				  		
				  		//determine orientation vector and oriented points				  	
					  		for(int i = 0; i < traces.size(); i++){
								progress.updateBarText("determining arc length and orientation vector... " + (i+1) + "/" + traces.size());
								traces.get(i).calculateArcLengthPositions();
								traces.get(i).updateAllOrientationVectors(minRefDist, maxRefDist);					
								traces.get(i).calculateOrientedPoints2D();
								traces.get(i).calculateOrientedPoints3D();
							}
							if(progress.isStopped())break running;
							
						//filter traces
							multi_focal_tools.removeProblematicTraces(progress, traces);
							if(progress.isStopped())break running;
				  	}
				  	
				  //determine curvature and angles	
				  	if(!calibrationMode){
				  		for(int i = 0; i < traces.size(); i++){
							progress.updateBarText("determining arc length and orientation vector... (" + (i+1) + "/" + traces.size() + ")");
							traces.get(i).calculateDAnglesAndCurvatures(curvRefDist, preventPoints, multi_focal_tools.NOZ);
							traces.get(i).calculateDAnglesAndCurvatures(curvRefDist, preventPoints, encoding);
						}
				  		if(progress.isStopped())break running;
				  		System.gc();
				  	}				  			
				  	
				  //calculate head rotation parameters and save data
				  	Arrays.fill(minIntensity, Double.POSITIVE_INFINITY);
				  	Arrays.fill(maxIntensity, Double.NEGATIVE_INFINITY);
				  	Arrays.fill(minIntensity, Double.POSITIVE_INFINITY);
				  	Arrays.fill(maxIntensity, 0.0);
				  	
				  	if(!calibrationMode){					
				  		IJ.log("sp " + slicePosition [0] + " - " + slicePosition [1] + " - " + slicePosition [2] + " - " + slicePosition [3] + "");
				  		for(int i = 0; i < traces.size(); i++){
				  			progress.updateBarText("calculating head rotation paramters... " + (i+1) + "/" + traces.size());
				  			if(!traces.get(i).oriented)	continue;
				  			traces.get(i).determineHeadRotation2D(imp, slicePosition, hrPlusMinusRange);
				  			if(!traces.get(i).hrVset)	continue;
				  			
				  			for(int s = 0; s < 4; s++){
				  				if(traces.get(i).getHeadRotationMaximumPositions() [s] > maxPosition [s]){
					  				maxPosition [s] = traces.get(i).getHeadRotationMaximumPositions() [s];  
					  			}
				  				if(traces.get(i).getHeadRotationMaximumPositions() [s] < minPosition [s]){
					  				minPosition [s] = traces.get(i).getHeadRotationMaximumPositions() [s];  
					  			}
				  				if(traces.get(i).getHeadRotationMaximumIntensities() [s] > maxIntensity [s]){
					  				maxIntensity [s] = traces.get(i).getHeadRotationMaximumIntensities() [s];  
					  			}
				  				if(traces.get(i).getHeadRotationMaximumIntensities() [s] < minIntensity [s]){
					  				minIntensity [s] = traces.get(i).getHeadRotationMaximumIntensities() [s];  
					  			}
				  			}
				  		}
				  		if(progress.isStopped())break running;
				  		
				  		multi_focal_tools.saveHeadRotationMatrixImage(traces, slicePosition, imp.getCalibration().pixelWidth, savePath);
				  		if(progress.isStopped())break running;
				  	}			
				  	System.gc();
				  	
				  	//Frequencies: analysis
				  	if(!calibrationMode && !calibrationTestMode){
				  		progress.updateBarText("calculate X frequencies");
					  	freqSummaryX [task] = multi_focal_tools.getAndSaveFrequencies(traces, savePath, xyCal, groupedTimesteps,
					  			multi_focal_tools.KYMOX, multi_focal_tools.NOZ, sampleRate, neglectedInitialArclength, orientation3D, speciesLength);
					  	
					  	progress.updateBarText("calculate Y frequencies");
					  	freqSummaryY [task] = multi_focal_tools.getAndSaveFrequencies(traces, savePath, xyCal, groupedTimesteps,
					  			multi_focal_tools.KYMOY, multi_focal_tools.NOZ, sampleRate, neglectedInitialArclength, orientation3D, speciesLength);
					  	
					  	progress.updateBarText("calculate Z frequencies ...");
				  		freqSummaryZ [task] = multi_focal_tools.getAndSaveFrequencies(traces, savePath, xyCal, groupedTimesteps,
				  				multi_focal_tools.KYMOZ, encoding, sampleRate, neglectedInitialArclength, orientation3D, speciesLength);
				  		
				  		progress.updateBarText("calculate curvature frequencies ...");
				  		freqSummaryCurv2D [task] = multi_focal_tools.getAndSaveFrequencies(traces, savePath, xyCal, groupedTimesteps,
				  				multi_focal_tools.KYMOCURV2D, multi_focal_tools.NOZ, sampleRate, neglectedInitialArclength, orientation3D, speciesLength);
						System.gc();
				  		freqSummaryCurv3D [task] = multi_focal_tools.getAndSaveFrequencies(traces, savePath, xyCal, groupedTimesteps,
				  				multi_focal_tools.KYMOCURV3D, encoding, sampleRate, neglectedInitialArclength, orientation3D, speciesLength);
				  		System.gc();
				  		
				  		progress.updateBarText("calculate curvature angle frequencies ...");
				  		freqSummaryCAngle2D [task] = multi_focal_tools.getAndSaveFrequencies(traces, savePath, xyCal, groupedTimesteps,
				  				multi_focal_tools.KYMOCANGLE2D, multi_focal_tools.NOZ, sampleRate, neglectedInitialArclength, orientation3D, speciesLength);
						System.gc();
				  		freqSummaryCAngle3D [task] = multi_focal_tools.getAndSaveFrequencies(traces, savePath, xyCal, groupedTimesteps,
				  				multi_focal_tools.KYMOCANGLE3D, encoding, sampleRate, neglectedInitialArclength, orientation3D, speciesLength);
				  		System.gc();
				  		
				  		progress.updateBarText("calculate theta freqs ...");
				  		freqSummaryTheta2D [task] = multi_focal_tools.getAndSaveGroupedFrequenciesTraceParam(traces,
				  				multi_focal_tools.TRACE_THETA, multi_focal_tools.NOZ, groupedTimesteps, sampleRate, savePath, 
				  				minIntensity, maxIntensity, minPosition, maxPosition);
				  		freqSummaryTheta3D [task] = multi_focal_tools.getAndSaveGroupedFrequenciesTraceParam(traces,
				  				multi_focal_tools.TRACE_THETA, encoding, groupedTimesteps, sampleRate, savePath, 
				  				minIntensity, maxIntensity, minPosition, maxPosition);
						System.gc();
				  		
						progress.updateBarText("calculate head rotation freqs ...");
				  		freqSummaryHrMaxPos [task] = multi_focal_tools.getAndSaveGroupedFrequenciesTraceParam(traces,
				  				multi_focal_tools.TRACE_HRMAXPOSITION, multi_focal_tools.NOZ, groupedTimesteps, sampleRate, savePath, 
				  				minIntensity, maxIntensity, minPosition, maxPosition);
				  		freqSummaryHrMaxInt [task] = multi_focal_tools.getAndSaveGroupedFrequenciesTraceParam(traces,
				  				multi_focal_tools.TRACE_HRMAXINTENSITY, multi_focal_tools.NOZ, groupedTimesteps, sampleRate, savePath, 
				  				minIntensity, maxIntensity, minPosition, maxPosition);
				  		freqSummaryHrAngle [task] = multi_focal_tools.getAndSaveGroupedFrequenciesTraceParam(traces,
				  				multi_focal_tools.TRACE_HRANGLE, multi_focal_tools.NOZ, groupedTimesteps, sampleRate, savePath, 
				  				minIntensity, maxIntensity, minPosition, maxPosition);					  		
						System.gc();
				  		
				  	}
				  					  		
				//save into file
				  	tp = new TextPanel("Results");
				  	tp.append("Saving date:	" + constants.dateTab.format(new Date()) + "	Analysis started:	" + constants.dateTab.format(startDate));
					tp.append("Processed image:	" + name [task]);
//						tp.append("");
//						tp.append("Registration file:	" + om.getFileName());
					tp.append("");
					tp.append("Settings: ");
					tp.append("	" + "sample rate [Hz]:	" + constants.df6US.format(sampleRate));
					tp.append("	" + "Orient points in 3D:	" + orientation3D);
					tp.append("	" + "Sperm head is tethered:	" + tethered);
					tp.append("	" + "xy calibration [um]:	" + constants.df6US.format(xyCal));
					tp.append("	" + "Slice Positions [um] 1|2|3|4:	" 
							+ constants.df3US.format(slicePosition [0]) 
							+ "	" + constants.df3US.format(slicePosition [1])
							+ "	" + constants.df3US.format(slicePosition [2]) 
							+ "	" + constants.df3US.format(slicePosition [3]));
					tp.append("	" + "trace determined in " + selectedTraceDetermination);
					tp.append("	" + "selected thresholding method:	" + selectedThresholdMethod);
					tp.append("	" + "gaussian blur for trace generation - sigma:	" + constants.df3US.format(gaussSigma));
					tp.append("	" + "Upscale trace (fold):	" + upscaleFold);
					tp.append("	" + "Head center-of-mass added as first point of trace:	" + addCOM);
					tp.append("	" + "Unify start points (e.g. for tethered sperm):	" + unifyStartPoints);
					tp.append("	" + "max vector length [points]:	" + constants.df0.format(maxVectorLength));
					tp.append("	" + "normal vector radius [um]:	" + constants.df6US.format(normalLength));
					tp.append("	" + "Smooth normal for xyz fit:	" + smoothNormal);
//					tp.append("	" + "Use normal maximum for z fit:	" + useNormalMaxForZGauss);
					tp.append("	" + "Exclude head points from correction/deletion/zkymographs (initial " + preventPoints + " points):	" + preventHeadFromCorr);
					tp.append("	" + "Use Gauss-4-Point fit for z determination:	" + getG4PZ);
					tp.append("	" + "Use width fit for z determination:	" + getWidthFitZ);
					tp.append("	" + "Width fit - LUT step size [um]:	" + constants.df3US.format(LUTStepSize));
					tp.append("	" + "Z-smoothing method:	" + zSmoothingMethod);
					tp.append("	" + "Accepted xy distance for z-smoothing:	" + constants.df6US.format(acceptedZSmoothDistance));					
					tp.append("	" + "# points before or after individual point included in xy- and z-smoothing:	" + constants.df0.format(plusMinusDistanceForSmooth));
					tp.append("	" + "Minimum xy-arcus position of points for reference vector:	" + constants.df6US.format(minRefDist));
					tp.append("	" + "Maximum xy-arcus position of points for reference vector:	" + constants.df6US.format(maxRefDist));
					tp.append("	" + "Curvature and dAngle Calcultion - distance of upstream reference point:	" + constants.df6US.format(curvRefDist));
					tp.append("	" + "GroupedTimesteps for fourier transform:	" + constants.df0.format(groupedTimesteps));
					tp.append("	" + "initial arclengths neglected for forier transform (µm):	" + constants.df6US.format(neglectedInitialArclength));
					tp.append("	" + "head rotation matrix radius (points):	" + constants.df0.format(hrPlusMinusRange));
					tp.append("");
					
					if(!calibrationMode && getWidthFitZ){
						tp.append("	" + "Width fit used for z determination: " + getWidthFitZ); 
//							tp.append("	" + "	" + "Z position of image 1 in width fit LUT:	" + widthLUTOffset); 
						tp.append("	" + "	" + "Width fit LUT step size:	" + LUTStepSize);
						tp.append("	" + "	" + "LUT file name:	" + nameLUT);
//							tp.append("	" + "	" + "Width corresponding to intensity of 1:	" + LUTStepSize);
//							tp.append("	" + "	" + "Width corresponding to intensity of 254:	" + LUTStepSize);
					}else{
						tp.append("");
//							tp.append("");
						tp.append("");
						tp.append("");
//							tp.append("");
//							tp.append("");
					}
					if(!calibrationMode && getG4PZ){
						tp.append("	" + "Gauss 4-point-fit used for z determination: " + getG4PZ); 
					}else{
						tp.append("");
					}
					tp.append("");
					
					if(!calibrationMode){
						String add;
						if(getWidthFitZ && getG4PZ){
							add = "Width to Z correction factor:";
							for(int j = 0; j < 10; j++){
								add += "	" + constants.df6US.format(widthZCorrFactor [j]);
							}						 
							tp.append(add);
						}					
						tp.append("Angle Theta [°]:		2D	3D no smoothing	3D mean smoothing	3D median smoothing");
						for(int i = 0; i < traces.size(); i++){
							tp.append("	" + traces.get(i).getFrame()
									+ "	" + constants.df6US.format(traces.get(i).getThetaDegree(multi_focal_tools.NOZ))
									+ "	" + constants.df6US.format(traces.get(i).getThetaDegree(multi_focal_tools.PUREZ))
									+ "	" + constants.df6US.format(traces.get(i).getThetaDegree(multi_focal_tools.MEANZ))
									+ "	" + constants.df6US.format(traces.get(i).getThetaDegree(multi_focal_tools.MEDIANZ)));
						}
						if(!calibrationTestMode){
							add = "	" + "found average primary freq.:" 
									+ "	" + constants.df6US.format(freqSummaryTheta2D[task][0]);
							add += "	"; if(encoding == multi_focal_tools.PUREZ)	add += constants.df6US.format(freqSummaryTheta3D[task][0]);
							add += "	"; if(encoding == multi_focal_tools.MEANZ)	add += constants.df6US.format(freqSummaryTheta3D[task][0]);
							add += "	"; if(encoding == multi_focal_tools.MEDIANZ)	add += constants.df6US.format(freqSummaryTheta3D[task][0]);
							tp.append(add);
							add = "	" + "found average secondary freq.:" 
									+ "	" + constants.df6US.format(freqSummaryTheta2D[task][1]);
							add += "	"; if(encoding == multi_focal_tools.PUREZ)	add += constants.df6US.format(freqSummaryTheta3D[task][1]);
							add += "	"; if(encoding == multi_focal_tools.MEANZ)	add += constants.df6US.format(freqSummaryTheta3D[task][1]);
							add += "	"; if(encoding == multi_focal_tools.MEDIANZ)	add += constants.df6US.format(freqSummaryTheta3D[task][1]);
							tp.append(add);	
							add = "	" + "found average COM freq.:" 
									+ "	" + constants.df6US.format(freqSummaryTheta2D[task][2]);
							add += "	"; if(encoding == multi_focal_tools.PUREZ)	add += constants.df6US.format(freqSummaryTheta3D[task][2]);
							add += "	"; if(encoding == multi_focal_tools.MEANZ)	add += constants.df6US.format(freqSummaryTheta3D[task][2]);
							add += "	"; if(encoding == multi_focal_tools.MEDIANZ)	add += constants.df6US.format(freqSummaryTheta3D[task][2]);
							tp.append(add);	
						}
						
						tp.append("");
						tp.append("head position		x [µm]	y [µm]	z (no interpol) [µm]	z (mean) [µm]	z (median) [µm]");
						for(int i = 0; i < traces.size(); i++){
							tp.append("	" + traces.get(i).getFrame()
									+ "	" + constants.df6US.format(traces.get(i).getTracePoints().get(0).getX())
									+ "	" + constants.df6US.format(traces.get(i).getTracePoints().get(0).getY())
									+ "	" + constants.df6US.format(traces.get(i).getTracePoints().get(0).getZ(multi_focal_tools.PUREZ))
									+ "	" + constants.df6US.format(traces.get(i).getTracePoints().get(0).getZ(multi_focal_tools.MEANZ))
									+ "	" + constants.df6US.format(traces.get(i).getTracePoints().get(0).getZ(multi_focal_tools.MEDIANZ)));
						}
						
						tp.append("");
						tp.append("head rotation		COM-angle	COM x	COM y	Maximum-Line-Angle	Maximum position sorted by plane position				"
								+ "Maximum position (normalized)				"
								+ "Maximum intensity sorted by plane position				"
								+ "Maximum intensity (normalized)				"
								+ "position-intensity vector angle");							
						for(int i = 0; i < traces.size(); i++){								
							String appendTxt = "";
							appendTxt += "	" + traces.get(i).getFrame();
							appendTxt += "	" + constants.df6US.format(traces.get(i).getHeadRotationCOMAngleInDegree());
							appendTxt += "	" + constants.df6US.format(traces.get(i).getHeadRotationMatrixCOM() [0]);
							appendTxt += "	" + constants.df6US.format(traces.get(i).getHeadRotationMatrixCOM() [1]);
							appendTxt += "	" + constants.df6US.format(traces.get(i).getHeadRotationMaximaAngleInDegree());
							
							int index = 0;
							for(int s = 0; s < 4; s++){
								index = tools.getIndexOfClosestValue(slicePosition, slicesSorted [s]);
								appendTxt += "	" + constants.df6US.format(traces.get(i).getHeadRotationMaximumPositions() [index]);
							}												
							for(int s = 0; s < 4; s++){
								index = tools.getIndexOfClosestValue(slicePosition, slicesSorted [s]);
								appendTxt += "	" + constants.df6US.format(tools.getNormalizedValue(
										traces.get(i).getHeadRotationMaximumPositions() [index],
										minPosition [index],
										maxPosition [index]) - 0.5);
							}
							for(int s = 0; s < 4; s++){
								index = tools.getIndexOfClosestValue(slicePosition, slicesSorted [s]);
								appendTxt += "	" + constants.df6US.format(traces.get(i).getHeadRotationMaximumIntensities() [index]);
							}												
							for(int s = 0; s < 4; s++){
								index = tools.getIndexOfClosestValue(slicePosition, slicesSorted [s]);
								appendTxt += "	" + constants.df6US.format(tools.getNormalizedValue(
										traces.get(i).getHeadRotationMaximumIntensities() [index],
										minIntensity [index],
										maxIntensity [index]) - 0.5);
							}
							
							double [] angle = traces.get(i).getHeadRotationAngle(minIntensity, maxIntensity, minPosition, maxPosition);
							for(int s = 0; s < 4; s++){
								index = tools.getIndexOfClosestValue(slicePosition, slicesSorted [s]);
								appendTxt += "	" + constants.df6US.format(angle[index]);
							}
							tp.append(appendTxt);
						}	
						if(!calibrationTestMode){
							tp.append("	average primary freq.					"
									+ constants.df6US.format(freqSummaryHrMaxPos [task][0]) 
									+ "				"
									+ "				"
									+ constants.df6US.format(freqSummaryHrMaxInt [task][0]) 
									+ "				"
									+ "				"
									+ constants.df6US.format(freqSummaryHrAngle [task][0]));
							tp.append("	average secondary freq.					"
									+ constants.df6US.format(freqSummaryHrMaxPos [task][1]) 
									+ "				"
									+ "				"
									+ constants.df6US.format(freqSummaryHrMaxInt [task][1]) 
									+ "				"
									+ "				"
									+ constants.df6US.format(freqSummaryHrAngle [task][1]));
							tp.append("	average COM freq.					"
									+ constants.df6US.format(freqSummaryHrMaxPos [task][2]) 
									+ "				"
									+ "				"
									+ constants.df6US.format(freqSummaryHrMaxInt [task][2]) 
									+ "				"
									+ "				"
									+ constants.df6US.format(freqSummaryHrAngle [task][2]));		
						}
						
						
						if(!calibrationTestMode){
							tp.append("");
							tp.append("Average found frequencies for arc-length oriented points (X)");
							tp.append("	" + "	" + "1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
							tp.append("	" + "1st frequency peak	" + constants.df6US.format(freqSummaryX [task][0][0]) 
								+ "	" + constants.df6US.format(freqSummaryX [task][0][1])
								+ "	" + constants.df6US.format(freqSummaryX [task][0][2]));
							tp.append("	" + "2nd frequency peak	" + constants.df6US.format(freqSummaryX [task][1][0]) 
								+ "	" + constants.df6US.format(freqSummaryX [task][1][1])
								+ "	" + constants.df6US.format(freqSummaryX [task][1][2]));
							tp.append("	" + "COM frequency	" + constants.df6US.format(freqSummaryX [task][2][0]) 
							+ "	" + constants.df6US.format(freqSummaryX [task][2][1])
							+ "	" + constants.df6US.format(freqSummaryX [task][2][2]));
							tp.append("");
							
							tp.append("Average found frequencies for arc-length oriented points (Y)");
							tp.append("	" + "	" + "1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
							tp.append("	" + "1st frequency peak	" + constants.df6US.format(freqSummaryY [task][0][0]) 
								+ "	" + constants.df6US.format(freqSummaryY [task][0][1])
								+ "	" + constants.df6US.format(freqSummaryY [task][0][2]));
							tp.append("	" + "2nd frequency peak	" + constants.df6US.format(freqSummaryY [task][1][0]) 
								+ "	" + constants.df6US.format(freqSummaryY [task][1][1])
								+ "	" + constants.df6US.format(freqSummaryY [task][1][2]));
							tp.append("	" + "COM frequency	" + constants.df6US.format(freqSummaryY [task][2][0]) 
							+ "	" + constants.df6US.format(freqSummaryY [task][2][1])
							+ "	" + constants.df6US.format(freqSummaryY [task][2][2]));
							tp.append("");
							
							tp.append("Average found frequencies for arc-length oriented points (Z)");
							tp.append("	" + "	" + "1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
							tp.append("	" + "1st frequency peak	" + constants.df6US.format(freqSummaryZ [task][0][0]) 
								+ "	" + constants.df6US.format(freqSummaryZ [task][0][1])
								+ "	" + constants.df6US.format(freqSummaryZ [task][0][2]));
							tp.append("	" + "2nd frequency peak	" + constants.df6US.format(freqSummaryZ [task][1][0]) 
								+ "	" + constants.df6US.format(freqSummaryZ [task][1][1])
								+ "	" + constants.df6US.format(freqSummaryZ [task][1][2]));
							tp.append("	" + "COM frequency	" + constants.df6US.format(freqSummaryZ [task][2][0]) 
							+ "	" + constants.df6US.format(freqSummaryZ [task][2][1])
							+ "	" + constants.df6US.format(freqSummaryZ [task][2][2]));
							tp.append("");
							
							tp.append("Average found frequencies for arc-length oriented points (Curvature 2D)");
							tp.append("	" + "	" + "1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
							tp.append("	" + "1st frequency peak	" + constants.df6US.format(freqSummaryCurv2D [task][0][0]) 
								+ "	" + constants.df6US.format(freqSummaryCurv2D [task][0][1])
								+ "	" + constants.df6US.format(freqSummaryCurv2D [task][0][2]));
							tp.append("	" + "2nd frequency peak	" + constants.df6US.format(freqSummaryCurv2D [task][1][0]) 
								+ "	" + constants.df6US.format(freqSummaryCurv2D [task][1][1])
								+ "	" + constants.df6US.format(freqSummaryCurv2D [task][1][2]));
							tp.append("	" + "COM frequency	" + constants.df6US.format(freqSummaryCurv2D [task][2][0]) 
								+ "	" + constants.df6US.format(freqSummaryCurv2D [task][2][1])
								+ "	" + constants.df6US.format(freqSummaryCurv2D [task][2][2]));
							tp.append("");
							
							tp.append("Average found frequencies for arc-length oriented points (Curvature 3D)");
							tp.append("	" + "	" + "1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
							tp.append("	" + "1st frequency peak	" + constants.df6US.format(freqSummaryCurv3D [task][0][0]) 
								+ "	" + constants.df6US.format(freqSummaryCurv3D [task][0][1])
								+ "	" + constants.df6US.format(freqSummaryCurv3D [task][0][2]));
							tp.append("	" + "2nd frequency peak	" + constants.df6US.format(freqSummaryCurv3D [task][1][0]) 
								+ "	" + constants.df6US.format(freqSummaryCurv3D [task][1][1])
								+ "	" + constants.df6US.format(freqSummaryCurv3D [task][1][2]));
							tp.append("	" + "COM frequency	" + constants.df6US.format(freqSummaryCurv3D [task][2][0]) 
								+ "	" + constants.df6US.format(freqSummaryCurv3D [task][2][1])
								+ "	" + constants.df6US.format(freqSummaryCurv3D [task][2][2]));	
							tp.append("");
							
							tp.append("Average found frequencies for arc-length oriented points (Curvature angle 2D)");
							tp.append("	" + "	" + "1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
							tp.append("	" + "1st frequency peak	" + constants.df6US.format(freqSummaryCAngle2D [task][0][0]) 
								+ "	" + constants.df6US.format(freqSummaryCAngle2D [task][0][1])
								+ "	" + constants.df6US.format(freqSummaryCAngle2D [task][0][2]));
							tp.append("	" + "2nd frequency peak	" + constants.df6US.format(freqSummaryCAngle2D [task][1][0]) 
								+ "	" + constants.df6US.format(freqSummaryCAngle2D [task][1][1])
								+ "	" + constants.df6US.format(freqSummaryCAngle2D [task][1][2]));
							tp.append("	" + "COM frequency	" + constants.df6US.format(freqSummaryCAngle2D [task][2][0]) 
								+ "	" + constants.df6US.format(freqSummaryCAngle2D [task][2][1])
								+ "	" + constants.df6US.format(freqSummaryCAngle2D [task][2][2]));
							tp.append("");
							
							tp.append("Average found frequencies for arc-length oriented points (Curvature angle 3D)");
							tp.append("	" + "	" + "1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
							tp.append("	" + "1st frequency peak	" + constants.df6US.format(freqSummaryCAngle3D [task][0][0]) 
								+ "	" + constants.df6US.format(freqSummaryCAngle3D [task][0][1])
								+ "	" + constants.df6US.format(freqSummaryCAngle3D [task][0][2]));
							tp.append("	" + "2nd frequency peak	" + constants.df6US.format(freqSummaryCAngle3D [task][1][0]) 
								+ "	" + constants.df6US.format(freqSummaryCAngle3D [task][1][1])
								+ "	" + constants.df6US.format(freqSummaryCAngle3D [task][1][2]));	
							tp.append("	" + "COM frequency	" + constants.df6US.format(freqSummaryCAngle3D [task][2][0]) 
								+ "	" + constants.df6US.format(freqSummaryCAngle3D [task][2][1])
								+ "	" + constants.df6US.format(freqSummaryCAngle3D [task][2][2]));		
							tp.append("");
						}
						
					}else{
						tp.append("Angle Theta [°]:		2D");
						for(int i = 0; i < traces.size(); i++){
							tp.append("	" + traces.get(i).getFrame()
									+ "	" + constants.df6US.format(traces.get(i).getThetaDegree(multi_focal_tools.NOZ)));
						}
						tp.append("");
					}								
					
					multi_focal_tools.addFooter(tp);		
				  	tp.saveAs(dir [task] + saveName + System.getProperty("file.separator") + "results.txt");
				  	
				//save selection
				  	rm = RoiManager.getInstance();
					if (rm==null) rm = new RoiManager();
					rm.reset();
					selections [task].setName("trace selection");
					rm.addRoi(selections [task]);						
					if(selectedTraceDetermination.equals(TRACEDETERMINATION[0])){
						selectionsSD [task].setName("SD selection");
						rm.addRoi(selectionsSD [task]);
					}						
					rm.runCommand("Save", savePath + "selectedRoi.zip");	
				
				//save images
				  	if(!calibrationMode){	
				  		multi_focal_tools.saveTraceImage(imp, traces, encoding,
				  				slicesSorted [0], slicesSorted[slicesSorted.length-1], savePath, preventPoints);	
				  		System.gc();
				  		multi_focal_tools.saveOrientedTraceImage(imp, traces, encoding, savePath, xyCal);
				  		System.gc();
				  		
				  		multi_focal_tools.saveKymograph(traces, xyCal, savePath, encoding, multi_focal_tools.KYMOX, preventPoints, orientation3D);
				  		System.gc();
				  		multi_focal_tools.saveKymograph(traces, xyCal, savePath, encoding, multi_focal_tools.KYMOY, preventPoints, orientation3D);
				  		System.gc();
				  		multi_focal_tools.saveKymograph(traces, xyCal, savePath, encoding, multi_focal_tools.KYMOZ, preventPoints, orientation3D);
				  		System.gc();
				  		multi_focal_tools.saveKymograph(traces, xyCal, savePath, encoding, multi_focal_tools.KYMOCURV2D, preventPoints, orientation3D);
				  		System.gc();
				  		multi_focal_tools.saveKymograph(traces, xyCal, savePath, encoding, multi_focal_tools.KYMOCURV3D, preventPoints, orientation3D);
				  		System.gc();
				  		multi_focal_tools.saveKymograph(traces, xyCal, savePath, encoding, multi_focal_tools.KYMOCANGLE2D, preventPoints, orientation3D);
				  		System.gc();
				  		multi_focal_tools.saveKymograph(traces, xyCal, savePath, encoding, multi_focal_tools.KYMOCANGLE3D, preventPoints, orientation3D);
				  		System.gc();	
				  		multi_focal_tools.saveKymograph(traces, xyCal, savePath, encoding, multi_focal_tools.KYMOMAXINTENSITY, preventPoints, orientation3D);
				  		System.gc();	
				  		multi_focal_tools.saveXYZCoordinates(traces, xyCal, savePath, encoding, progress);
				  		System.gc();
				  	}			  	
				  	
				//save progress dialog log file
				  	progress.saveLog(dir [task] + saveName + System.getProperty("file.separator") +"log.txt");
				  
			  	//finish progress dialog
				  	imp.changes = false;
				  	imp.close();
				  	
				  	progress.setBar(1.0);
				  	done = true;
				  	tasksSuccessfull [task] = true;
				  	break running;		  	
			}//(end runnning)
			
			if(progress.isStopped()) break tasking;
			progress.moveTask(task);			
		}
		
		boolean success = false;
		for(int task = 0; task < tasks; task++){
			if(tasksSuccessfull [task]){
				success = true;
			}
		}		
		
		if(success && !calibrationMode && !calibrationTestMode){
			tp = new TextPanel("Results Summary");
			tp.append("Summary results for processing of the following images:");
			tp.append("	" + "image ID	directory	name");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" + (task+1) + "	" + dir [task] + "	" + name [task]);
				}else{
					tp.append("	" + (task+1) + "	" + dir [task] + "	" + name [task] + "	NOT SUCCESSFULLY PROCESSED!");
				}
			}
			tp.append("");
			tp.append("#####################################################");
			tp.append("");
			
			tp.append("Average Found 1st-Peak frequencies for arc-length oriented points (X)");
			tp.append("	image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryX [task][0][0]) 
					+ "	" + constants.df6US.format(freqSummaryX [task][0][1])
					+ "	" + constants.df6US.format(freqSummaryX [task][0][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found 2nd-Peak frequencies for arc-length oriented points (X)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryX [task][1][0]) 
					+ "	" + constants.df6US.format(freqSummaryX [task][1][1])
					+ "	" + constants.df6US.format(freqSummaryX [task][1][2]));
				}else{
					tp.append("	" +(task+1));
				}			
			}	
			tp.append("");
			
			tp.append("Average Found COM frequencies for arc-length oriented points (X)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryX [task][2][0]) 
					+ "	" + constants.df6US.format(freqSummaryX [task][2][1])
					+ "	" + constants.df6US.format(freqSummaryX [task][2][2]));
				}else{
					tp.append("	" +(task+1));
				}			
			}	
			tp.append("");
			
			tp.append("Average Found 1st-Peak frequencies for arc-length oriented points (Y)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryY [task][0][0]) 
					+ "	" + constants.df6US.format(freqSummaryY [task][0][1])
					+ "	" + constants.df6US.format(freqSummaryY [task][0][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found 2nd-Peak frequencies for arc-length oriented points (Y)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryY [task][1][0]) 
					+ "	" + constants.df6US.format(freqSummaryY [task][1][1])
					+ "	" + constants.df6US.format(freqSummaryY [task][1][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found COM frequencies for arc-length oriented points (Y)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryY [task][2][0]) 
					+ "	" + constants.df6US.format(freqSummaryY [task][2][1])
					+ "	" + constants.df6US.format(freqSummaryY [task][2][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			
			tp.append("Average Found 1st-Peak frequencies for arc-length oriented points (Z)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryZ [task][0][0]) 
					+ "	" + constants.df6US.format(freqSummaryZ [task][0][1])
					+ "	" + constants.df6US.format(freqSummaryZ [task][0][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found 2nd-Peak frequencies for arc-length oriented points (Z)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryZ [task][1][0]) 
					+ "	" + constants.df6US.format(freqSummaryZ [task][1][1])
					+ "	" + constants.df6US.format(freqSummaryZ [task][1][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found COM frequencies for arc-length oriented points (Z)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryZ [task][2][0]) 
					+ "	" + constants.df6US.format(freqSummaryZ [task][2][1])
					+ "	" + constants.df6US.format(freqSummaryZ [task][2][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			
			tp.append("Average Found 1st-Peak frequencies for arc-length oriented points (Curvature 2D)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryCurv2D [task][0][0]) 
					+ "	" + constants.df6US.format(freqSummaryCurv2D [task][0][1])
					+ "	" + constants.df6US.format(freqSummaryCurv2D [task][0][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found 2nd-Peak frequencies for arc-length oriented points (Curvature 2D)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryCurv2D [task][1][0]) 
					+ "	" + constants.df6US.format(freqSummaryCurv2D [task][1][1])
					+ "	" + constants.df6US.format(freqSummaryCurv2D [task][1][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found COM frequencies for arc-length oriented points (Curvature 2D)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryCurv2D [task][2][0]) 
					+ "	" + constants.df6US.format(freqSummaryCurv2D [task][2][1])
					+ "	" + constants.df6US.format(freqSummaryCurv2D [task][2][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
						
			tp.append("Average Found 1st-Peak frequencies for arc-length oriented points (Curvature 3D)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryCurv3D [task][0][0]) 
					+ "	" + constants.df6US.format(freqSummaryCurv3D [task][0][1])
					+ "	" + constants.df6US.format(freqSummaryCurv3D [task][0][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found 2nd-Peak frequencies for arc-length oriented points (Curvature 3D)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryCurv3D [task][1][0]) 
					+ "	" + constants.df6US.format(freqSummaryCurv3D [task][1][1])
					+ "	" + constants.df6US.format(freqSummaryCurv3D [task][1][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found COM frequencies for arc-length oriented points (Curvature 3D)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryCurv3D [task][2][0]) 
					+ "	" + constants.df6US.format(freqSummaryCurv3D [task][2][1])
					+ "	" + constants.df6US.format(freqSummaryCurv3D [task][2][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found 1st-Peak frequencies for arc-length oriented points (Curvature Angle 2D)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryCAngle2D [task][0][0]) 
					+ "	" + constants.df6US.format(freqSummaryCAngle2D [task][0][1])
					+ "	" + constants.df6US.format(freqSummaryCAngle2D [task][0][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found 2nd-Peak frequencies for arc-length oriented points (Curvature Angle 2D)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryCAngle2D [task][1][0]) 
					+ "	" + constants.df6US.format(freqSummaryCAngle2D [task][1][1])
					+ "	" + constants.df6US.format(freqSummaryCAngle2D [task][1][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found COM frequencies for arc-length oriented points (Curvature Angle 2D)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryCAngle2D [task][2][0]) 
					+ "	" + constants.df6US.format(freqSummaryCAngle2D [task][2][1])
					+ "	" + constants.df6US.format(freqSummaryCAngle2D [task][2][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
						
			tp.append("Average Found 1st-Peak frequencies for arc-length oriented points (Curvature Angle 3D)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryCAngle3D [task][0][0]) 
					+ "	" + constants.df6US.format(freqSummaryCAngle3D [task][0][1])
					+ "	" + constants.df6US.format(freqSummaryCAngle3D [task][0][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found 2nd-Peak frequencies for arc-length oriented points (Curvature Angle 3D)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryCAngle3D [task][1][0]) 
					+ "	" + constants.df6US.format(freqSummaryCAngle3D [task][1][1])
					+ "	" + constants.df6US.format(freqSummaryCAngle3D [task][1][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average Found COM frequencies for arc-length oriented points (Curvature Angle 3D)");
			tp.append("	" + "image ID	1st third of flagellum	2nd third of flagellum	3rd third of flagellum");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) + "	" + constants.df6US.format(freqSummaryCAngle3D [task][2][0]) 
					+ "	" + constants.df6US.format(freqSummaryCAngle3D [task][2][1])
					+ "	" + constants.df6US.format(freqSummaryCAngle3D [task][2][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			tp.append("");
			
			tp.append("Average found frequencies for angle theta");
			tp.append("	" + "image ID"
					+ "	" + "Primary freq. (2D)" + "	" + "Secondary freq. (2D)" 
					+ "	" + "COM freq. (2D)"
					+ "	" + "Primary freq. (3D)" + "	" + "Secondary freq. (3D)"
					+ "	" + "COM freq. (3D)");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) 
							+ "	" + constants.df6US.format(freqSummaryTheta2D [task][0])
							+ "	" + constants.df6US.format(freqSummaryTheta2D [task][1])
							+ "	" + constants.df6US.format(freqSummaryTheta2D [task][2])
							+ "	" + constants.df6US.format(freqSummaryTheta3D [task][0])
							+ "	" + constants.df6US.format(freqSummaryTheta3D [task][1])
							+ "	" + constants.df6US.format(freqSummaryTheta3D [task][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			
			tp.append("Average found head rotation frequencies");
			tp.append("	" + "image ID"
					+ "	" + "Primary freq. (max position)" + "	" + "Secondary freq. (max position)" + "	" + "COM freq. (max position)"
					+ "	" + "Primary freq. (max intensity)" + "	" + "Secondary freq. (max intensity)" + "	" + "COM freq. (max intensity)"
					+ "	" + "Primary freq. (pos.-int.-angle)" + "	" + "Secondary freq. (pos.-int.-angle)" + "	" + "COM freq. (pos.-int.-angle)");
			for(int task = 0; task < tasks; task++){
				if(tasksSuccessfull [task]){
					tp.append("	" +(task+1) 
							+ "	" + constants.df6US.format(freqSummaryHrMaxPos [task][0])
							+ "	" + constants.df6US.format(freqSummaryHrMaxPos [task][1])
							+ "	" + constants.df6US.format(freqSummaryHrMaxPos [task][2])
							+ "	" + constants.df6US.format(freqSummaryHrMaxInt [task][0])
							+ "	" + constants.df6US.format(freqSummaryHrMaxInt [task][1])
							+ "	" + constants.df6US.format(freqSummaryHrMaxInt [task][2])
							+ "	" + constants.df6US.format(freqSummaryHrAngle [task][0])
							+ "	" + constants.df6US.format(freqSummaryHrAngle [task][1])
							+ "	" + constants.df6US.format(freqSummaryHrAngle [task][2]));
				}else{
					tp.append("	" +(task+1));
				}	
			}	
			
			multi_focal_tools.addFooter(tp);
			tp.saveAs(homePath + System.getProperty("file.separator") + "SpermQ-MF_Summary_" + constants.dateName.format(new Date()) + ".txt");
			System.gc();
			done = true;
			new WaitForUserDialog("All tasks have been processed. A summary file has been saved on the desktop!").show();
		}	  	
	}
}