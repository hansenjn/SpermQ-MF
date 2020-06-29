/***===============================================================================
 
 SpermQ-MF Version v0.3.2

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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.LinkedList;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.process.*;
import ij.text.TextPanel;
import ij.measure.*;
import ij.plugin.frame.RoiManager;
import spermQ_mf.edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import spermQ_mf.jnhsupport.*;
import spermQ_mf.skeleton_analysis.*;
import spermQ_mf.skeletonize3D.Skeletonize3D_;


public class multi_focal_tools implements Measurements{
	//selections
	public static final int NOZ = 0,
			PUREZ = 1,
			MEANZ = 2,
			MEDIANZ = 3;
	
	public static final int MINDIST = 0,
			MAXDIST = 1,
			MEAN = 2,
			MEDIAN = 3;
	
	//Kymograph types
		public static final int KYMOX = 0,
				KYMOY = 1,
				KYMOZ = 2,
				KYMOMAXINTENSITY = 3,
				KYMOCURV2D = 42,
				KYMOCURV3D = 43,
				KYMOCANGLE2D = 51,
				KYMOCANGLE3D = 52,
				KYMODZ = 6,
				KYMOTANGENTANGLE = 7;
		
		//Kymograph minima and maxima
		public static final double KYMIN_X = -50.0, KYMAX_X = 130.0,
				KYMIN_Y = -130.0, KYMAX_Y = 130.0,
				KYMIN_Z = -10.0, KYMAX_Z = 100.0,
				KYMIN_MAXINTENSITY = 1.0, KYMAX_MAXINTENSITY = 65533.0,
				KYMIN_CURV = 1.0, KYMAX_CURV = 4.0,
				KYMIN_CANGLE2D = -180.0, KYMAX_CANGLE2D = 180.0,
				KYMIN_CANGLE3D = 0.0, KYMAX_CANGLE3D = 270.0,
				KYMIN_DZ = -15.0, KYMAX_DZ = 15.0,
				KYMIN_FREQ = 1.0, KYMAX_FREQ = 200.0,
				KYMIN_AMPLITUDE = 1.0, KYMAX_AMPLITUDE = 65533.0;
	
	//Trace parameters
		public static final int TRACE_THETA = 0,
			TRACE_HRMAXINTENSITY = 1,
			TRACE_HRMAXPOSITION = 2,
			TRACE_HRANGLE = 3;
			
		public static final double SPECIESLENGTH_HUMAN = 40,
			SPECIESLENGTH_MOUSE = 100;
		
	/**
	 * @deprecated
	 * */
	public static ImagePlus calibrateIntensity(ImagePlus imp, double [] calibration){
		final int slices = imp.getNSlices();		
		if(slices!=4)	return null;
		
		for(int c = 0; c < imp.getNChannels(); c++){
			for(int t = 0; t < imp.getNFrames(); t++){
				for(int x = 0; x < imp.getWidth(); x++){
					for(int y = 0; y < imp.getHeight(); y++){
						for(int z = 0; z < slices; z++){
							int s = imp.getStackIndex(c+1, z+1, t+1)-1;
							imp.getStack().setVoxel(x, y, s, imp.getStack().getVoxel(x, y, s) * calibration[z]);
						}						
					}
				}
			}
		}
		
		return imp;
	}
	
	/**
	 * @deprecated
	 * */
	public static ImagePlus registerByFile(ImagePlus imp, String filePath, ProgressDialog progress){
		if(imp.getNChannels()>1)	IJ.error("registration not implemented for > 1 channel...");
		int c = 1;
		
		ImagePlus tempImp;
		for(int t = 0; t < imp.getNFrames(); t++){
			tempImp = IJ.createImage("temp imp", imp.getWidth(), imp.getHeight(), imp.getNSlices(), imp.getBitDepth());
			
			for(int x = 0; x < imp.getWidth(); x++){
				for(int y = 0; y < imp.getHeight(); y++){
					for(int z = 0; z < imp.getNSlices(); z++){
						int s = imp.getStackIndex(c+1, z+1, t+1)-1;
						tempImp.getStack().setVoxel(x, y, z, imp.getStack().getVoxel(x, y, s));
						tempImp.getStack().setSliceLabel("slice " + (z+1), (z+1));
					}						
				}
			}
			
//			userCheck(tempImp);		
//			IJ.log(filePath);
			
			tempImp.show();		
			
			IJ.run(tempImp, "MultiStackReg", "stack_1=[temp imp] action_1=[Load Transformation File] file_1=["+filePath+"]"
//			+ " stack_2=None action_2=Ignore file_2=[] transformation=[Rigid Body]");
			+ " stack_2=None action_2=Ignore file_2=[] transformation=[Scaled Rotation]");
			
//			userCheck(tempImp);
			tempImp.hide();
			
			for(int x = 0; x < imp.getWidth(); x++){
				for(int y = 0; y < imp.getHeight(); y++){
					for(int z = 0; z < imp.getNSlices(); z++){
						int s = imp.getStackIndex(c+1, z+1, t+1)-1;
						imp.getStack().setVoxel(x, y, s, tempImp.getStack().getVoxel(x, y, z));
					}						
				}
			}
			
			tempImp.changes=false;
			tempImp.close();
			progress.addToBar(0.1*(1.0/(double)imp.getNFrames()));
			progress.updateBarText("registering t=" + (t+1) + "...");
		}
		System.gc();
		return imp;
	}
	
	public static Roi getSelection(ImagePlus imp){
		//TODO variant using Tiff virtualstack?
		ImagePlus maxImp = maxIntensityProjection(imp);
		maxImp.show();
		IJ.setTool("polygon");
		while(true){
			new WaitForUserDialog("Set a Roi containing parts of the cell in every frame!").show();
			if(maxImp.getRoi()!=null) break;
		}		
		Roi roi = new PolygonRoi(maxImp.getRoi().getPolygon(), PolygonRoi.POLYGON);
		maxImp.changes = false;
		maxImp.close();
		
		return roi;
	}
	
	public static Roi getSelectionInStack(ImagePlus imp){
		imp.show();
		IJ.setTool("polygon");
		while(true){
			new WaitForUserDialog("Set a Roi circumscribing the Roi to measure SD!").show();
			if(imp.getRoi()!=null) break;
		}		
		Roi roi = new PolygonRoi(imp.getRoi().getPolygon(), PolygonRoi.POLYGON);
		
		imp.hide();
		
		return roi;
	}
	
	public static ImagePlus maxIntensityProjection(ImagePlus imp){
		ImagePlus impMax = IJ.createImage("maximum projection", imp.getWidth(), imp.getHeight(), 1, imp.getBitDepth());
		
		for(int x = 0; x < imp.getWidth(); x++){
			for(int y = 0; y < imp.getHeight(); y++){
				impMax.getStack().setVoxel(x, y, 0, 0.0);
				for(int s = 0; s < imp.getStackSize(); s++){
					if(imp.getStack().getVoxel(x, y, s) > impMax.getStack().getVoxel(x, y, 0)){
						impMax.getStack().setVoxel(x, y, 0, imp.getStack().getVoxel(x, y, s));
					}
				}
			}
		}		
		
		impMax.setCalibration(imp.getCalibration());
		return impMax;
	}
	
	/**
	 * @returns an independent ImagePlus of the sharpest image in the ImagePlus
	 * The sharpest image is defined as the image with the lowest Standard Deviation
	 * imp: stack-image where the sharpest image shall be found
	 * */	
	public static ImagePlus getSharpestPlane(ImagePlus imp){
		double [] SD = new double [imp.getStackSize()];
		double maximum = 0.0; int maximumPos = -1;
		
		double average = 0.0;
		int counter = 0;
		for(int z = 0; z < imp.getStackSize(); z++){
			average = 0.0;
			counter = 0;
			for(int x = 0; x < imp.getWidth(); x++){
				for(int y = 0; y < imp.getHeight(); y++){
					if(imp.getStack().getVoxel(x, y, z) != 0.0){
						average += imp.getStack().getVoxel(x, y, z);
						counter++;
					}
				}
			}	
			average /= counter;
			
			for(int x = 0; x < imp.getWidth(); x++){
				for(int y = 0; y < imp.getHeight(); y++){
					if(imp.getStack().getVoxel(x, y, z) != 0.0){
						SD [z] += Math.pow(imp.getStack().getVoxel(x, y, z) - average, 2.0);
					}					
				}
			}
			SD [z] = Math.sqrt(SD [z] / (double)(counter-1));
			if(SD[z] > maximum) {
				maximum = SD [z]; 
				maximumPos = z;			
			}
		}	
		
		return getSingleImageFromStack(imp, maximumPos);
	}
	
	/**
	 * @returns an independent ImagePlus of the sharpest image in the ImagePlus
	 * The sharpest image is defined as the image with the lowest Standard Deviation
	 * imp: stack-image where the sharpest image shall be found
	 * */		
	public static ImagePlus getSharpestPlane(ImagePlus imp, Roi roi){
		double [] SD = new double [imp.getStackSize()];
		double maximum = 0.0; int maximumPos = -1;
		
		double average = 0.0;
		int counter = 0;
		for(int z = 0; z < imp.getStackSize(); z++){
			average = 0.0;
			counter = 0;
			for(int x = 0; x < imp.getWidth(); x++){
				for(int y = 0; y < imp.getHeight(); y++){
					if(imp.getStack().getVoxel(x, y, z) != 0.0 
							&& roi.contains(x,y)){
						average += imp.getStack().getVoxel(x, y, z);
						counter++;
					}
				}
			}	
			average /= counter;
			
			for(int x = 0; x < imp.getWidth(); x++){
				for(int y = 0; y < imp.getHeight(); y++){
					if(imp.getStack().getVoxel(x, y, z)!=0.0 && roi.contains(x,y)){
						SD [z] += Math.pow(imp.getStack().getVoxel(x, y, z) - average, 2.0);
					}					
				}
			}
			SD [z] = Math.sqrt(SD [z] / (double)(counter-1));
			if(SD[z] > maximum) {
				maximum = SD [z]; 
				maximumPos = z;			
			}
		}	
		
		return getSingleImageFromStack(imp, maximumPos);
	}
	
	/**
	 * @returns an independent ImagePlus of the sharpest image in the ImagePlus within the time-frame range
	 * The sharpest image is defined as the image with the lowest Standard Deviation		 
	 * */	
	public static ImagePlus getSharpestPlane(ImagePlus imp, double [] sliceDistances, double stepDistance, int rangeMin, int rangeMax){
			
		double [] SD = new double [imp.getStackSize()];
		double maximum = 0.0; int maximumPos = -1;
		
		int z, counter; 
		double average;
		int newRangeMin, newRangeMax;
		for(int s = 0; s < 4; s++){
			newRangeMin = rangeMin + (int)Math.round((sliceDistances [s] - sliceDistances [0]) / stepDistance);
			if(newRangeMin < 0)	newRangeMin = 0;
			newRangeMax = rangeMin + (int)Math.round((sliceDistances [s] - sliceDistances [0]) / stepDistance);
			
			for(int t = newRangeMin; t <= newRangeMax && t <= imp.getNFrames(); t++){
				z = imp.getStackIndex(1, s+1, t+1);
				average = 0.0;
				counter = 0;
				for(int x = 0; x < imp.getWidth(); x++){
					for(int y = 0; y < imp.getHeight(); y++){
						if(imp.getStack().getVoxel(x, y, z)!=0.0){
							average += imp.getStack().getVoxel(x, y, z);
							counter++;
						}
					}
				}
				average /= counter;
				
				for(int x = 0; x < imp.getWidth(); x++){
					for(int y = 0; y < imp.getHeight(); y++){
						if(imp.getStack().getVoxel(x, y, z)!=0.0){
							SD [z] += Math.pow(imp.getStack().getVoxel(x, y, z) - average, 2.0);
						}					
					}
				}
				SD [z] = Math.sqrt(SD [z] / (double)(counter-1));
				if(SD[z] > maximum) {
					maximum = SD [z]; 
					maximumPos = z;			
				}
			}		
		}	
		
		return getSingleImageFromStack(imp, maximumPos);
	}
	
	/**
	 * @return an independent ImagePlus of the selected stack image
	 * @param imp: image where single timepoint shall be derived from
	 * @param i: stack image 0 < i < # stack images
	 * */
	public static ImagePlus getSingleImageFromStack(ImagePlus imp, int i){
		ImagePlus newImp = IJ.createImage("Z" + i, imp.getWidth(), imp.getHeight(), 1, imp.getBitDepth());
		for(int x = 0; x < imp.getWidth(); x++){
			for(int y = 0; y < imp.getHeight(); y++){
				newImp.getStack().setVoxel(x, y, 0, imp.getStack().getVoxel(x, y, i));
			}
		}
		newImp.setCalibration(imp.getCalibration());
		return newImp;
	}
	
	public static ImagePlus getSingleTimepoint(ImagePlus imp, int t){
		/**
		 * Returns an independent ImagePlus of the selected Timepoint
		 * imp: image where single timepoint shall be derived from
		 * t: selected timestep 0 < t < #time-steps in image)
		 * */
		ImagePlus newImp = IJ.createImage("T" + t, imp.getWidth(), imp.getHeight(), imp.getNSlices(), imp.getBitDepth());
		for(int s = 0; s < imp.getNSlices(); s++){
			for(int x = 0; x < imp.getWidth(); x++){
				for(int y = 0; y < imp.getHeight(); y++){
					newImp.getStack().setVoxel(x, y, s, imp.getStack().getVoxel(x, y, imp.getStackIndex(1, s+1, t+1)-1));
				}
			}
		}	
		newImp.setCalibration(imp.getCalibration());
		return newImp;
	}
	
	/**
	 * @return ArrayList of traces, which each belong to a different time-step
	 * @param imp: ImagePlus where the object of interest is contained
	 * @param algorithm: the selected threshold algorithm for processing
	 * @param sigma: the sigma-value of the gaussian blur applied before and after binarization
	 * @param progress: the progress dialog where log/error/update messages are supposed to be displayed
	 * @param selection: ROI in which the object of interest is present in every frame
	 * @param selectionSD: ROI in which (if selected in traceDet), the sharpest plane shall be calculated
	 * @param calibrationMode: If selected the same trace is saved for every time-step of the stack
	 * @param traceDet: sets the variant of processing the 4 focal planes into one image, in which the trace shall be calculated
	 * */
	public static ArrayList<trace> getObjectTraces(ImagePlus imp, String algorithm, double sigma, ProgressDialog progress, Roi selection, Roi selectionSD,
			boolean calibrationMode, String traceDet, boolean addCOM, double minRefDist, double maxRefDist){
		//not implemented for multichannel so far... 
		if(imp.getNChannels()!=1) return null;
		
		ArrayList<trace> traceList = new ArrayList<trace>(imp.getNFrames());
		
		//calibration mode
		if(calibrationMode){
			//get sharpest image
			progress.replaceBarText("Looking for sharpest image");
			
//			ImagePlus sp = getSharpestPlane(imp, sliceDistances, sliceStep, minZPos, maxZPos);
			ImagePlus sp = getSharpestPlane(imp, selectionSD);	
			trace mipTrace = getTraceBySkeletonization(sp, 0, selection, sigma, algorithm, addCOM,  minRefDist, maxRefDist);
			sp.close();
			System.gc();
			
			ArrayList<trackPoint> points;
			for(int t = 0; t < imp.getNFrames();t++){
				if(progress.isStopped())	return null;
				progress.addToBar(0.1*(1.0/(double)imp.getNFrames()));
				progress.updateBarText("generate trace for t=" + (t+1) + "...");
				
				points = new ArrayList<trackPoint>(mipTrace.getTracePoints().size());
				//Copy list
				for(int i = 0; i < mipTrace.getTracePoints().size(); i++){
					points.add(new trackPoint(mipTrace.getTracePoints().get(i).getX(),
							mipTrace.getTracePoints().get(i).getY()));
				}
				traceList.add(new trace(points, t, minRefDist, maxRefDist));
				points.clear();
				
				if(t%100 == 0)	System.gc();
			}
			System.gc();
			return traceList;
		}
		
		// analysis mode
		
		//find out sharpest plane for all
		int maxPos = 0;
		if(traceDet.equals(multi_focal.TRACEDETERMINATION [4])){	//"sharpest plane in time projection"
			maxPos = impProcessing.getSharpestPlanePosition(impProcessing.maxIntensityProjectionOfTimesteps(imp), selectionSD);
			System.gc();
		}
		
		ImagePlus selImp;
		for(int t = 0; t < imp.getNFrames(); t++){
			if(progress.isStopped())	return null;
			progress.addToBar(0.1*(1.0/(double)imp.getNFrames()));
			progress.updateBarText("generate trace for t=" + (t+1) + "...");
			
			selImp = getSingleTimepoint(imp,t);
			if(traceDet.equals(multi_focal.TRACEDETERMINATION [0])){	//"sharpest plane"
				selImp = getSharpestPlane(selImp);
			}else if(traceDet.equals(multi_focal.TRACEDETERMINATION [1])){	//"maximum-intensity-projection"
				selImp = impProcessing.maxIntensityProjection(selImp);
			}else if(traceDet.equals(multi_focal.TRACEDETERMINATION [2])){	//"average-intensity-projection"
				selImp = impProcessing.averageIntensityProjection(selImp);
			}else if(traceDet.equals(multi_focal.TRACEDETERMINATION [3])){	//"sum of planes"
				selImp = impProcessing.stackSum(selImp);
			}else if(traceDet.equals(multi_focal.TRACEDETERMINATION [4])){	//"sharpest plane in time projection"
				selImp = impProcessing.getSingleImageFromStack(selImp, maxPos);
			}
			traceList.add(getTraceBySkeletonizationV2(selImp, t, selection, sigma, algorithm, addCOM, minRefDist, maxRefDist, progress));
			
			if(t%100 == 0)	System.gc();
		}
		System.gc();
		return traceList;
	}
	
	/**
	 * @deprecated in versions v0.3.1 and higher
	 * */
	private static trace getTraceBySkeletonization (ImagePlus impMax, int frame, Roi selection, double sigma, String thresholdAlgorithm, boolean addCOM,
			double minRefDist, double maxRefDist){	
		ImagePlus impMaxSave = impMax.duplicate();
		
		//eventually scale down before finding threshold?	TODO
		//eventually 8-bit conv before finding threshold?	TODO
					
		//gauss-filter
			impMax.getProcessor().blurGaussian(sigma);
//			IJ.log("bl1");userCheck(impMax);
		
		//threshold image
			thresholdImage(impMax,thresholdAlgorithm);
			ImagePlus impMaxSaveThresholded = impMax.duplicate();
			
		//gauss-filter
			impMax.getProcessor().blurGaussian(sigma);	
//			IJ.log("blur2");userCheck(impMax);
					
		//convert to 8-bit
			optimal8BitConversion(impMax);
					
		//skeletonize			
//			IJ.run(impMax,"Skeletonize (2D/3D)","");
			Skeletonize3D_ skelProc = new Skeletonize3D_();
			skelProc.setup("", impMax);
			skelProc.run(impMax.getProcessor());
//			IJ.log("skeletonized"); userCheck(impMax); 
//			IJ.log("ph " + impMax.getCalibration().pixelHeight);
		
		//analyze skeleton for shortest path
			AnalyzeSkeleton_ skel = new AnalyzeSkeleton_();
			skel.calculateShortestPath = true;
			skel.setup("", impMax);
			
			SkeletonResult sklRes = skel.run(AnalyzeSkeleton_.NONE, false, true, null, true, false);
			//run(int pruneIndex, boolean pruneEnds, boolean shortPath, ImagePlus origIP, boolean silent, boolean verbose)
			ArrayList<spermQ_mf.skeleton_analysis.Point>[] shortestPath = skel.getShortestPathPoints();
				
		//find longest Skeleton touching selection
			int chosenShortestPath = 0;
			if(shortestPath.length==0){
				return null;
			}else if(shortestPath.length!=1){
				double maxLength = 0.0;
				for(int i = 0; i < sklRes.getNumOfTrees(); i++){
					checking: for(int j = 0; j < shortestPath[i].size(); j++){
						if(selection.contains(shortestPath[i].get(j).x, shortestPath[i].get(j).y)){
							if(sklRes.getAverageBranchLength()[i]*sklRes.getBranches()[i]>maxLength){
								maxLength = sklRes.getAverageBranchLength()[i]*sklRes.getBranches()[i];
								chosenShortestPath = i;
							}
							break checking;
						}							
					}											
				}
			}		
			
		//count points 
			int nPoints = shortestPath[chosenShortestPath].size();
//			IJ.log("n" + nPoints);
			
		//get trace list -  new method
			ArrayList<trackPoint> list = new ArrayList<trackPoint>(nPoints);
			LinkedList<trackPoint> unsortedList = new LinkedList<trackPoint>();
							
			//find end point
//			IJ.log("searching w" + impMax.getCalibration().pixelWidth + " h" + impMax.getCalibration().pixelHeight);
			trackPoint startEnd = null; int startIndex = -1;
			searching: for(int i = 0; i < nPoints; i++){
				double counter = 0;
				trackPoint p = new trackPoint(shortestPath[chosenShortestPath].get(i).x,
						shortestPath[chosenShortestPath].get(i).y);
				for(int j = 0; j < nPoints; j++){
					if(i!=j){
						trackPoint q = new trackPoint(shortestPath[chosenShortestPath].get(j).x,
								shortestPath[chosenShortestPath].get(j).y);
//						IJ.log(i + " - " + j + " = " + get2DDistance(p,q));
						if(getDistance(p,q,NOZ) <= constants.sqrt2){
//							IJ.log(i + "-" + j);
							counter++;
						}
					}
					
				}
				if (counter==1){ 
//					IJ.log("found start" + i);
//					IJ.log("counter " + counter);
					startEnd = new trackPoint(shortestPath[chosenShortestPath].get(i).x * impMax.getCalibration().pixelWidth,
							shortestPath[chosenShortestPath].get(i).y * impMax.getCalibration().pixelHeight);
					startIndex = i;
					break searching;
				}								
			}
			
			if(startIndex == -1){
				IJ.log("Problem no start found");
			}
			
			//save unsorted list
			for(int i = 0; i < nPoints; i++){
				if(i != startIndex){
					unsortedList.add(new trackPoint(shortestPath[chosenShortestPath].get(i).x*impMax.getCalibration().pixelWidth,
							shortestPath[chosenShortestPath].get(i).y*impMax.getCalibration().pixelHeight));
				}					
			}
			
			//create sortedList (list)
			{
				list.add(startEnd);
				
//				IJ.log(unsortedList.size() + " uls");
//				IJ.log(list.size() + " ls");
				
				int index = 0;
				double distance;
				trackPoint p;
				int pIndex;
				while(!unsortedList.isEmpty()){
					distance = Double.POSITIVE_INFINITY; 
					p = null;
					pIndex = -1;
					for(int i = 0; i < unsortedList.size(); i++){
						if(getDistance(unsortedList.get(i),list.get(index),NOZ) < distance){
							p = unsortedList.get(i);
							pIndex = i;
							distance = getDistance(unsortedList.get(i),list.get(index),NOZ);
						}
					}
					if(startIndex == -1 || p.equals(null)){
						IJ.log("Problem no next point found");
					}					
					unsortedList.remove(pIndex);
					list.add(p);
					index++;
				}
			}
			list.trimToSize();				
			
			//get statistics for first point				
				OvalRoi roi0 = new OvalRoi((int)Math.round(list.get(0).getX()/impMax.getCalibration().pixelWidth) - 8,
						(int)Math.round(list.get(0).getY()/impMax.getCalibration().pixelHeight) - 8, 17, 17);
				impMaxSave.setRoi(roi0);
				ImageStatistics stats0 = impMaxSave.getStatistics();
				double intensity0 = stats0.area * stats0.mean;
				double [] COM0 = getXYCenterOfMass(impMaxSaveThresholded,roi0);
				
			//get statistics for last point
				OvalRoi roiE = new OvalRoi((int)Math.round(list.get(list.size()-1).getX()/impMax.getCalibration().pixelWidth) - 8,
						(int)Math.round(list.get(list.size()-1).getY()/impMax.getCalibration().pixelHeight) - 8, 17, 17);
				impMaxSave.setRoi(roiE);
				ImageStatistics statsE = impMaxSave.getStatistics();
				double intensityE = statsE.area*statsE.mean;
				double [] COME = getXYCenterOfMass(impMaxSaveThresholded,roiE);
				
			//close images
				impMaxSave.changes = false;
				impMaxSave.close();
				impMaxSaveThresholded.changes = false;
				impMaxSaveThresholded.close();
				impMax.changes = false;
				impMax.close();
				
			//if first point is actually last point reverse list
				ArrayList <trackPoint> newList = new ArrayList <trackPoint>(nPoints);
			if(intensity0 < intensityE){
				//add center of mass point of head (only if not equal to first skeletal point)
				if(list.get(nPoints-1).getX() == COME [0] && list.get(nPoints-1).getY() == COME [1]){
//					IJ.log("first point = COME");
				}else if(addCOM){
					newList.ensureCapacity(nPoints+1);
					newList.add(new trackPoint(COME [0], COME [1]));
				}				
								
				//invert list
				for(int i = nPoints-1; i >= 0; i--){	
					newList.add(list.get(i));	
				}					
			}else{				
				//add center of mass point of head (only if not equal to first skeletal point)
				if(list.get(0).getX() == COM0 [0] && list.get(0).getY() == COM0 [1]){
//					IJ.log("first point = COM0");
				}else if(addCOM){
					newList.ensureCapacity(nPoints+1);
					newList.add(new trackPoint(COM0 [0], COM0 [1]));
				}
								
				//invert list
				for(int i = 0; i < nPoints; i++){	
					newList.add(list.get(i));	
				}
			}	
			newList.trimToSize();
			return new trace(newList,frame, minRefDist, maxRefDist);
	}
		
	
	/*
	 * New skeleton analysis method to avoid interruptions at side branches,
	 * used from version v0.3.1 on and higher
	 * 
	 * To establish this method, code was derived from SpermQ v0.2.1, https://github.com/hansenjn/SpermQ‚
	 * 	 * */
	private static trace getTraceBySkeletonizationV2 (ImagePlus impMax, int frame, Roi selection, double sigma, String thresholdAlgorithm, boolean addCOM,
			double minRefDist, double maxRefDist, ProgressDialog progress){	
		ImagePlus impMaxSave = impMax.duplicate();
		
		//eventually scale down before finding threshold?	TODO
		//eventually 8-bit conv before finding threshold?	TODO
					
		//gauss-filter
			impMax.getProcessor().blurGaussian(sigma);
//			IJ.log("bl1");userCheck(impMax);
		
		//threshold image
			thresholdImage(impMax,thresholdAlgorithm);
			ImagePlus impMaxSaveThresholded = impMax.duplicate();
			
		//gauss-filter
			impMax.getProcessor().blurGaussian(sigma);	
//			IJ.log("blur2");userCheck(impMax);
					
		//convert to 8-bit
			optimal8BitConversion(impMax);
					
		//skeletonize
			double pixelWidth = impMax.getCalibration().pixelWidth;
			double pixelHeight = impMax.getCalibration().pixelHeight;
			
			impMax.getCalibration().pixelWidth = 1.0;
			impMax.getCalibration().pixelHeight = 1.0;
			
			
//			IJ.run(impMax,"Skeletonize (2D/3D)","");
			Skeletonize3D_ skelProc = new Skeletonize3D_();
			skelProc.setup("", impMax);
			skelProc.run(impMax.getProcessor());
//			IJ.log("skeletonized"); userCheck(impMax); 
//			IJ.log("ph " + impMax.getCalibration().pixelHeight);
		
		//analyze skeleton for shortest path
			AnalyzeSkeleton_ skel = new AnalyzeSkeleton_();
			skel.calculateShortestPath = true;
			skel.setup("", impMax);
			
			SkeletonResult sklRes = skel.run(AnalyzeSkeleton_.NONE, false, true, null, true, false);
			//run(int pruneIndex, boolean pruneEnds, boolean shortPath, ImagePlus origIP, boolean silent, boolean verbose)
			ArrayList<spermQ_mf.skeleton_analysis.Point>[] shortestPath = skel.getShortestPathPoints();
				
		//find longest Skeleton touching selection
			int chosenShortestPath = 0;
			if(shortestPath.length==0){
				return null;
			}else if(shortestPath.length!=1){
				double maxLength = 0.0;
				for(int i = 0; i < sklRes.getNumOfTrees(); i++){
					checking: for(int j = 0; j < shortestPath[i].size(); j++){
						if(selection.contains(shortestPath[i].get(j).x, shortestPath[i].get(j).y)){
							if(sklRes.getAverageBranchLength()[i]*sklRes.getBranches()[i]>maxLength){
								maxLength = sklRes.getAverageBranchLength()[i]*sklRes.getBranches()[i];
								chosenShortestPath = i;
							}
							break checking;
						}							
					}											
				}
			}			
			
		//count points 
			int nPoints = shortestPath[chosenShortestPath].size();
//			IJ.log("n" + nPoints);
			
		//get trace list -  new method
			ArrayList<trackPoint> list = new ArrayList<trackPoint>(nPoints);
			LinkedList<trackPoint> unsortedList = new LinkedList<trackPoint>();
							
			//find end point
//			trackPoint2D startEnd = null; int startIndex = -1;
			LinkedList<trackPoint> startEnds = new LinkedList <trackPoint>();
			LinkedList<Integer> startIndexes = new LinkedList <Integer>();
			
			int counter;
			trackPoint p, q;
			for(int i = 0; i < nPoints; i++){
				counter = 0;
				p = new trackPoint(shortestPath[chosenShortestPath].get(i).x,
						shortestPath[chosenShortestPath].get(i).y);
				searchingSecond: for(int j = 0; j < nPoints; j++){
					if(i!=j){
						q = new trackPoint(shortestPath[chosenShortestPath].get(j).x,
								shortestPath[chosenShortestPath].get(j).y);
//						IJ.log(i + " - " + j + " = " + get2DDistance(p,q));
						if(getDistance(p,q,NOZ) <= constants.sqrt2){
//							IJ.log(i + "-" + j);
							counter++;
							if(counter==2){
								break searchingSecond;
							}
						}
					}					
				}
				if (counter==1){ 
//					IJ.log("found start" + i);
//					IJ.log("counter " + counter);
					startEnds.add(new trackPoint(shortestPath[chosenShortestPath].get(i).x * pixelWidth,
							shortestPath[chosenShortestPath].get(i).y * pixelHeight));
					startIndexes.add(i);
//					break searching;
				}								
			}
			System.gc();
			
			if(startIndexes.isEmpty()){//Stop checking!
				progress.notifyMessage("frame " + frame + ": no start found", ProgressDialog.NOTIFICATION);
				return null;
			}
			
			//obtain best point list
			int minLeftOverPoints = Integer.MAX_VALUE, bestIndex = -2;; 
			boolean isIn;
			findingBest: for(int se = 0; se < startEnds.size(); se++){
				ArrayList<trackPoint> tempList = new ArrayList<trackPoint>(nPoints);
				//save unsorted list
				for(int i = 0; i < nPoints; i++){
					if(i != startIndexes.get(se)){
						unsortedList.add(new trackPoint(shortestPath[chosenShortestPath].get(i).x*pixelWidth,
								shortestPath[chosenShortestPath].get(i).y*pixelHeight));
					}					
				}
				
				//create sortedList (list)
				{
					tempList.add(startEnds.get(se));
					
//					IJ.log(unsortedList.size() + " uls");
//					IJ.log(list.size() + " ls");
					
					int index = 0;
					double distance;
					int pIndex;
					
					sorting: while(!unsortedList.isEmpty()){
						distance = Double.POSITIVE_INFINITY; 
						p = null;
						pIndex = -1;
						for(int i = 0; i < unsortedList.size(); i++){
							if(getDistance(unsortedList.get(i),tempList.get(index),NOZ) < distance){
								p = unsortedList.get(i);
								pIndex = i;
								distance = getDistance(unsortedList.get(i),tempList.get(index),NOZ);
							}
						}
						if(p.equals(null)){
							IJ.log("Problem no next point found");
						}					
						unsortedList.remove(pIndex);						
						if(Math.sqrt(Math.pow((int)Math.round(p.getX()/pixelWidth)-(int)Math.round(tempList.get(index).getX()/pixelWidth),2.0)+
								Math.pow((int)Math.round(p.getY()/pixelHeight)-(int)Math.round(tempList.get(index).getY()/pixelHeight),2.0)) > constants.sqrt2){
							isIn = false;
							scanJnctn: for(int jnctn = 0; jnctn < sklRes.getListOfJunctionVoxels().size(); jnctn++){
								if(sklRes.getListOfJunctionVoxels().get(jnctn).x == (int)Math.round(tempList.get(index).getX()/pixelWidth) 
										&& sklRes.getListOfJunctionVoxels().get(jnctn).y == (int)Math.round(tempList.get(index).getY()/pixelWidth)){
									isIn = true;
									break scanJnctn;
								}
							}
							if(!isIn){
//								progress.notifyMessage("frame " + frame + " - try " + (se+1) + "/" + startEnds.size() + ": " + unsortedList.size() + " points discarded. X " 
//										+ (int)Math.round(p.getX()/pixelWidth)
//										+ " _ " + (int)Math.round(tempList.get(index).getX()/pixelWidth)
//										+ " Y " + (int)Math.round(p.getY()/pixelHeight) 
//										+ " _ " + (int)Math.round(tempList.get(index).getY()/pixelHeight), ProgressDialog.LOG);
								break sorting;
							}else{
								tempList.add(p);
								index++;
							}							
						}else{
							tempList.add(p);
							index++;
						}					
					}
					tempList.trimToSize();	
					if(!unsortedList.isEmpty()){						
						if(minLeftOverPoints > unsortedList.size()){
							list = new ArrayList<trackPoint>(tempList);
							bestIndex = se;
							minLeftOverPoints = unsortedList.size();
						}else{
							tempList.clear();
							tempList = null;
							System.gc();
						}
					}else{
//						progress.notifyMessage("frame " + frame + " - try " + (se+1) + "/" + startEnds.size() + ": " + unsortedList.size() + " points discarded.", ProgressDialog.LOG);
						list = new ArrayList<trackPoint>(tempList);
						bestIndex = se;
						minLeftOverPoints = 0;
						break findingBest;
					}
				}				
			}
//			progress.notifyMessage("frame " + frame + " best index: " + (bestIndex+1) + " with " + minLeftOverPoints + " left-over points.", ProgressDialog.LOG);
			nPoints = list.size();
			
			
			//get statistics for first point				
				OvalRoi roi0 = new OvalRoi((int)Math.round(list.get(0).getX()/impMax.getCalibration().pixelWidth) - 8,
						(int)Math.round(list.get(0).getY()/impMax.getCalibration().pixelHeight) - 8, 17, 17);
				impMaxSave.setRoi(roi0);
				ImageStatistics stats0 = impMaxSave.getStatistics();
				double intensity0 = stats0.area * stats0.mean;
				double [] COM0 = getXYCenterOfMass(impMaxSaveThresholded,roi0);
				
			//get statistics for last point
				OvalRoi roiE = new OvalRoi((int)Math.round(list.get(list.size()-1).getX()/impMax.getCalibration().pixelWidth) - 8,
						(int)Math.round(list.get(list.size()-1).getY()/impMax.getCalibration().pixelHeight) - 8, 17, 17);
				impMaxSave.setRoi(roiE);
				ImageStatistics statsE = impMaxSave.getStatistics();
				double intensityE = statsE.area*statsE.mean;
				double [] COME = getXYCenterOfMass(impMaxSaveThresholded,roiE);
				
			//close images
				impMaxSave.changes = false;
				impMaxSave.close();
//				impMaxSaveThresholded.changes = false;
//				impMaxSaveThresholded.close();
				impMax.changes = false;
				impMax.close();
				
			//if first point is actually last point reverse list
				ArrayList <trackPoint> newList = new ArrayList <trackPoint>(nPoints);
			if(intensity0 < intensityE){
				//add center of mass point of head (only if not equal to first skeletal point)
				if(list.get(nPoints-1).getX() == COME [0] && list.get(nPoints-1).getY() == COME [1]){
//					IJ.log("first point = COME");
				}else if(addCOM){
					newList.ensureCapacity(nPoints+1);
					newList.add(new trackPoint(COME [0], COME [1]));
				}				
								
				//invert list
				for(int i = nPoints-1; i >= 0; i--){	
					newList.add(list.get(i));	
				}					
			}else{				
				//add center of mass point of head (only if not equal to first skeletal point)
				if(list.get(0).getX() == COM0 [0] && list.get(0).getY() == COM0 [1]){
//					IJ.log("first point = COM0");
				}else if(addCOM){
					newList.ensureCapacity(nPoints+1);
					newList.add(new trackPoint(COM0 [0], COM0 [1]));
				}
								
				//invert list
				for(int i = 0; i < nPoints; i++){	
					newList.add(list.get(i));	
				}
			}	
			newList.trimToSize();
			return new trace(newList,frame, minRefDist, maxRefDist);
	}
	
	private static void userCheck(ImagePlus impMax) {
		impMax.show();
		new WaitForUserDialog("Check").show();
		impMax.hide();
	}

	/**
	 * @returns the X- and Y-coordinate of the center-of-mass within a selection at the current hyperstack-position.
	 * @param imp: the image of which the COM shall be derived
	 * @param selection: the Roi in which calculations are performed - if null, the whole image is included into analysis
	 * 
	 * */
	public static double [] getXYCenterOfMass (ImagePlus imp, Roi selection){
		double [] COM = {0.0, 0.0};
		double intensitySum = 0.0;
		if(selection == null){
			imp.deleteRoi();
			for(int x = 0; x < imp.getWidth(); x++){
				for(int y = 0; y < imp.getHeight(); y++){
					intensitySum += imp.getStack().getVoxel(x, y, imp.getZ()-1);
					COM [0] += imp.getStack().getVoxel(x, y, imp.getZ()-1) * x * imp.getCalibration().pixelWidth;
					COM [1] += imp.getStack().getVoxel(x, y, imp.getZ()-1) * y * imp.getCalibration().pixelHeight;
				}
			}			
		}else{
			for(int x = 0; x < imp.getWidth(); x++){
				for(int y = 0; y < imp.getHeight(); y++){
					if(selection.getPolygon().contains(x,y)){
						intensitySum += imp.getStack().getVoxel(x, y, imp.getZ()-1);
						COM [0] += imp.getStack().getVoxel(x, y, imp.getZ()-1) * x * imp.getCalibration().pixelWidth;
						COM [1] += imp.getStack().getVoxel(x, y, imp.getZ()-1) * y * imp.getCalibration().pixelHeight;						
					}
				}
			}
		}
		COM [0] /= intensitySum;
		COM [1] /= intensitySum;
		return COM;
	}
	
	private static void thresholdImage (ImagePlus imp, String algorithm){
		//threshold image
		IJ.setAutoThreshold(imp, (algorithm + " dark"));
		double minThreshold = imp.getProcessor().getMinThreshold();
		double imageMax = Math.pow(2.0,imp.getBitDepth())-1.0;
					
		for(int x = 0; x < imp.getWidth(); x++){
			for(int y = 0; y < imp.getHeight(); y++){
				if(imp.getStack().getVoxel(x, y, 0) >= minThreshold){
					imp.getStack().setVoxel(x, y, 0, imageMax);
				}else{
					imp.getStack().setVoxel(x, y, 0, 0.0);
				}
			}
		}
//		IJ.log("bin ");userCheck(imp);
	}
	
	private static void optimal8BitConversion (ImagePlus imp){
		//set displayrange from minimum to maximum and then convert
		
		double min = Double.POSITIVE_INFINITY, max = 0.0;
		for(int x = 0; x < imp.getWidth(); x++){
			for(int y = 0; y < imp.getHeight(); y++){
				if(imp.getStack().getVoxel(x, y, 0)>max){
					max = imp.getStack().getVoxel(x, y, 0);
				}
				if(imp.getStack().getVoxel(x, y, 0)<min){
					min = imp.getStack().getVoxel(x, y, 0);
				}
			}
		}
		imp.setDisplayRange(min, max);
		if(imp.getBitDepth() != 8){
			ImageConverter impConv = new ImageConverter(imp);
			impConv.convertToGray8();
		}		
//		IJ.log("8bit " + max);userCheck(impMax);	
	}
	
	/**
	 * For tethered sperm. Reverse Traces whose last point is closer to the median start point than the first point
	 * Method improved on 29.06.2020 (build with version 0.3.2, copied from SpermQ v0.2.1): use distance instead of sum of x and y difference to decide for closer point
	 * @author Jan N. Hansen (github: hansenjn)
	 * */
	public static void reverseReversedTraces(ArrayList<trace> traces, ProgressDialog progress){
		//remove 0 traces
		for(int i = traces.size()-1; i >= 0; i--){
			if(traces.get(i).getTracePoints().size()==0){
				traces.remove(i);
			}
		}	
		traces.trimToSize();
		
		double [] stX = new double [traces.size()],
				stY = new double [traces.size()];
		
		for(int i = 0; i < traces.size(); i++){
			stX [i] = traces.get(i).getTracePoints().get(0).getX();
			stY [i] = traces.get(i).getTracePoints().get(0).getY();
		}
		
		double mediStX = tools.getMedian(stX);
		double mediStY = tools.getMedian(stY);
		
		for(int i = 0; i < traces.size(); i++){
			if(Math.sqrt(Math.pow(mediStX - stX [i], 2.0) + Math.pow(mediStY - stY [i],2.0))
				< Math.sqrt(Math.pow(mediStX - traces.get(i).getTracePoints().get(traces.get(i).getTracePoints().size()-1).getX(), 2.0)
					+ Math.pow(mediStY - traces.get(i).getTracePoints().get(traces.get(i).getTracePoints().size()-1).getY(), 2.0))){
			}else{
				//invert tracepoints
				traces.get(i).reverseTracePoints();
			}
			progress.updateBarText("Checking for reversed traces: trace nr " + i);
		}		
	}
		
	/**
	 * For freely swimming sperm. Method introduced on 29.06.2020 (build with version 0.3.2, copied from SpermQ v0.2.1, 
	 * originally developed on 14.08.2019 for SpermQ version v1.0.9)
	 * @author Jan N. Hansen (github: hansenjn)
	 * @description The method determines the median head position in the last traces (number of traces defined by 
	 *  frameDistanceForComparison / 2) and checks whether the head is still closer than the tail; 
	 *  if this not the case, the trace is inverted. Next the method checks whether the first point
	 *  moves towards non-ciliary regions considering the last traces (number of traces defined by 
	 *  frameDistanceForComparison) 
	 * */
	public static void reverseReversedTracesOfFree(ArrayList<trace> traces, ProgressDialog progress,
			int frameDistanceForComparison){
		//remove 0 traces
		for(int i = traces.size()-1; i >= 0; i--){
			if(traces.get(i).getTracePoints().size()==0){
				traces.remove(i);
				progress.notifyMessage("trace at t = " + i + " has been removed (no points)", ProgressDialog.LOG);
			}
		}	
		traces.trimToSize();
		
		boolean [] correctFirst = new boolean [traces.size()];
		correctFirst [0] = false;
		double stX, stY, eX, eY, ostX, ostY;
		double [] ostXs = new double [(int)(frameDistanceForComparison/2.0)];
		double [] ostYs = new double [(int)(frameDistanceForComparison/2.0)];
		int counter;
		for(int i = 1; i < traces.size(); i++){
			stX = traces.get(i).getTracePoints().get(0).getX();
			stY = traces.get(i).getTracePoints().get(0).getY();
			
			eX = traces.get(i).getTracePoints().get(traces.get(i).getTracePoints().size()-1).getX();
			eY = traces.get(i).getTracePoints().get(traces.get(i).getTracePoints().size()-1).getY();
			
			Arrays.fill(ostXs, Double.POSITIVE_INFINITY);
			Arrays.fill(ostYs, Double.POSITIVE_INFINITY);
			counter = 0;
			for(int j = 1; j < (int)(frameDistanceForComparison/2.0)+1 && i-j >=0; j++){
				ostXs [j-1] = traces.get(i-j).getTracePoints().get(0).getX();
				ostYs [j-1]= traces.get(i-j).getTracePoints().get(0).getY();
				counter ++;
			}
			
			ostX = tools.getMedianOfRange(ostXs, 0, counter-1);
			ostY = tools.getMedianOfRange(ostYs, 0, counter-1);
			
			if(Math.sqrt(Math.pow(stX-ostX, 2.0) + Math.pow(stY-ostY, 2.0))
					<= Math.sqrt(Math.pow(eX-ostX, 2.0) + Math.pow(eY-ostY, 2.0))){
				//all is ok, no inversion needed
			}else{
				//invert trace's points
				traces.get(i).reverseTracePoints();
			}
			progress.updateBarText("Checking for reversed traces (free sperm): trace nr " + i);
		}
		
		//check whether the first point moves towards non-ciliary regions
		if(frameDistanceForComparison<traces.size()){
			int counterWrong = 0;
			double xC, yC, xCAll, yCAll;
			int ct;
			for(int i = frameDistanceForComparison; i < traces.size(); i++){
				xCAll = 0.0; yCAll = 0.0; ct = 0;
				for(int j = -1; j <= 1 && i + j < traces.size(); j++){
					ct ++;
					xC = traces.get(i+j).getTracePoints().get(0).getX() 
							- traces.get(i+j).getTracePoints().get(traces.get(i+j).getTracePoints().size()-1).getX();
					yC = traces.get(i+j).getTracePoints().get(0).getY() 
							- traces.get(i+j).getTracePoints().get(traces.get(i+j).getTracePoints().size()-1).getY();
					xCAll += traces.get(i+j).getTracePoints().get(0).getX() - xC/2.0;
					yCAll += traces.get(i+j).getTracePoints().get(0).getY() - yC/2.0;
				}
				xCAll /= (double) ct;
				yCAll /= (double) ct;
						
				stX = traces.get(i-frameDistanceForComparison).getTracePoints().get(0).getX();
				stY = traces.get(i-frameDistanceForComparison).getTracePoints().get(0).getY();
				
				eX = traces.get(i-frameDistanceForComparison).getTracePoints().get(
						traces.get(i-frameDistanceForComparison).getTracePoints().size()-1).getX();
				eY = traces.get(i-frameDistanceForComparison).getTracePoints().get(
						traces.get(i-frameDistanceForComparison).getTracePoints().size()-1).getY();
				
				if(Math.sqrt(Math.pow(xCAll-stX, 2.0) + Math.pow(yCAll-stY, 2.0))
						> Math.sqrt(Math.pow(xCAll-eX, 2.0) + Math.pow(yCAll-eY, 2.0))){
					counterWrong ++;
				}	
				progress.updateBarText("Checking general reversement: trace nr " + i);
			}
			if((double) counterWrong / (double) (traces.size() - frameDistanceForComparison) > 0.5){
				for(int i = 0; i < traces.size(); i++){
					//invert tracepoints
					traces.get(i).reverseTracePoints();
				}
			}
		}		
	}
	
	/**
	 * Find and save a common trace starting point (recommended for tethered sperm)
	 * */
	public static void unifyStartPoints(ArrayList<trace> traces, ProgressDialog progress){
		//remove 0 traces
		for(int i = traces.size()-1; i >= 0; i--){
			if(traces.get(i).getTracePoints().size()==0){
				traces.remove(i);
			}
		}	
		traces.trimToSize();
		
		double [] stX = new double [traces.size()],
				stY = new double [traces.size()];
		
		for(int i = 0; i < traces.size(); i++){
			stX [i] = traces.get(i).getTracePoints().get(0).getX();
			stY [i] = traces.get(i).getTracePoints().get(0).getY();
		}
		
		double mediStX = tools.getMedian(stX);
		double mediStY = tools.getMedian(stY);
		
		//save new start point
		for(int i = 0; i < traces.size(); i++){
			traces.get(i).getTracePoints().add(0, new trackPoint(mediStX, mediStY));
			progress.updateBarText("Saving merged start point: trace nr " + i);
		}
	}
	
	public static void adjustPointsViaNormalVector(ImagePlus imp, ArrayList<trace> traces, ProgressDialog progress, boolean saveNormalVectorRoiSet,
			String saveDirectory, int vectorSize, double normalRad, boolean smoothNormal, boolean preventHeadFromCorr, int preventPoints){
		if(imp.getNSlices()>4)	IJ.error("to many slices for correction via normal vector...");
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		int calculationRadius = (int)Math.round(normalRad/imp.getCalibration().pixelWidth);
		
		RoiManager rm = RoiManager.getInstance();
		if (rm==null) rm = new RoiManager();
		
		int channel = 1;	//no multichannel implemented so far!
				
		//initialize variables
		final int halfSize = (vectorSize-(vectorSize%2))/2;
		
		double x1, y1, x2, y2;
		int t;
		ArrayList<trackPoint> points;
		int nrOfPoints;
		double [] radii = new double [4];
		int radius;
		
		int gaussMismatch = 0, undefined = 0;
		double resultsParameters [][] = new double [5][4];
		int slicesIncluded = 0;
		double normalPointsX [][];
		double normalPointsY [][];
	
		CurveFitter cf;
		double [] parameters;	
		double rSquare;
		
	//calculate
		
		String undefinedPoints;
		for(int i = 0; i < traces.size(); i++){
			traces.get(i).xyCorrected = true;
			undefinedPoints = "";
			if(traces.get(i).getTracePoints().size()>0){
				if(progress.isStopped())	return;
				progress.updateBarText("Improving x,y position " + (i+1) + "/" + traces.size());
				
				t = traces.get(i).getFrame();
				points = traces.get(i).getTracePoints();
				nrOfPoints = points.size();			
				
				if(saveNormalVectorRoiSet){
					rm.reset();
				}
				
				
				
				//generate vectors
				x1 = (double)points.get(0).getX();
				y1 = (double)points.get(0).getY();
				x2 = (double)points.get(0).getX();
				y2 = (double)points.get(0).getY();
				for(int j = nrOfPoints-1; j >= 0; j--){
					//calculate tangential vector					
					searchUpstream: for(int vu = halfSize; vu >= 0; vu--){
						if(j-vu >= 0){
							x1 = (double)points.get(j-vu).getX();
							y1 = (double)points.get(j-vu).getY();
							break searchUpstream;
						}
					}
					
					searchDownstream: for(int vd = halfSize; vd >= 0; vd--){
						if(j+vd < nrOfPoints){
							x2 = (double)points.get(j+vd).getX();
							y2 = (double)points.get(j+vd).getY();
							break searchDownstream;
						}
					}					
					points.get(j).setVectors(x2-x1, y2-y1);
				}
				
				//fit along normal vector
				gaussMismatch = 0;
				undefined = 0;				
				for(int j = nrOfPoints-1; j >= 0; j--){
					//adjust calculation radius if touching borders
					radii [0] = Math.sqrt(Math.pow(((double)(w-2)*imp.getCalibration().pixelWidth - points.get(j).getX()) / (points.get(j).getNormalVector()[0]*imp.getCalibration().pixelWidth),2.0));
//					if(radii [0] < 0)	radii [0] = Double.MAX_VALUE;
					radii [1] = Math.sqrt(Math.pow(((double)(h-2)*imp.getCalibration().pixelHeight - points.get(j).getY()) / (points.get(j).getNormalVector()[1]*imp.getCalibration().pixelHeight),2.0));
//					if(radii [1] < 0)	radii [1] = Double.MAX_VALUE;
					radii [2] = Math.sqrt(Math.pow((points.get(j).getX()) / (points.get(j).getNormalVector()[0]*imp.getCalibration().pixelWidth),2.0));
//					if(radii [2] < 0)	radii [2] = Double.MAX_VALUE;
					radii [3] = Math.sqrt(Math.pow((points.get(j).getY()) / (points.get(j).getNormalVector()[1]*imp.getCalibration().pixelHeight),2.0));
//					if(radii [3] < 0)	radii [3] = Double.MAX_VALUE;
					
					Arrays.sort(radii);
					radius = (int)Math.round(radii [0]);
					if(radius > calculationRadius){
//						IJ.log("r1" + radius);
						radius = calculationRadius;
//						IJ.log("r2" + radius);
					}else{
//						IJ.log(" x " + points.get(j).getX() + " y " + points.get(j).getY() + " w " + w + " h " + h);
//						IJ.log(" ri " + radii [0]);
//						IJ.log(" ri " + radii [1]);
//						IJ.log(" ri " + radii [2]);
//						IJ.log(" ri " + radii [3]);
//						IJ.log(" nx " + points.get(j).getNormalVector()[0]);
//						IJ.log(" ny " + points.get(j).getNormalVector()[1]);
					}
					
					//save vector Roi
					if(saveNormalVectorRoiSet){
						Line vectorRoi = points.get(j).getVectorLine(imp);
						vectorRoi.setName("vector " + j);
						rm.addRoi(vectorRoi);	
						
						Line normalVectorRoi = points.get(j).getNormalVectorLine(imp,calculationRadius);
						normalVectorRoi.setName("normal vector " + j);
						rm.addRoi(normalVectorRoi);	
					}
					
					//Fit along normal vector
					slicesIncluded = 0;
					normalPointsX = new double [imp.getNSlices()][radius*2+1];
					normalPointsY = new double [imp.getNSlices()][radius*2+1];
					
					for(int s = 0; s < imp.getNSlices(); s++){
						//obtain data points for fitting		
						try{
							normalPointsX [s][radius] = 0.0;
							normalPointsY [s][radius] = impProcessing.getInterpolatedIntensity2D(imp, points.get(j).getX(),	//x 
							points.get(j).getY(),	//y
							imp.getStackIndex(channel, s+1, t+1)-1);	//intensity		
						}catch(Exception e){
							IJ.log("Error s " + s + " x " + points.get(j).getX() + "y " + points.get(j).getY() + " stacki " + (imp.getStackIndex(channel, s+1, t+1)-1));
							IJ.log("w" + imp.getWidth() + " h" + imp.getHeight() + " ss" + imp.getStackSize());
						}
												
							
						for(int p = 0; p < radius; p++){
							normalPointsX [s][p + 1 + radius] = (p+1) * (points.get(j).getNormalVectorLength() * imp.getCalibration().pixelWidth);
							normalPointsY [s][p + 1 + radius] = impProcessing.getInterpolatedIntensity2D(imp,
									(double)points.get(j).getX() + (p+1) * (points.get(j).getNormalVector()[0] * imp.getCalibration().pixelWidth),
									(double)points.get(j).getY() + (p+1) * (points.get(j).getNormalVector()[1] * imp.getCalibration().pixelHeight),
									imp.getStackIndex(channel, s+1, t+1)-1);			
							
							normalPointsX [s][radius - 1 - p] = (-1) * (p+1) * (points.get(j).getNormalVectorLength() * imp.getCalibration().pixelWidth);
							normalPointsY [s][radius - 1 - p] = impProcessing.getInterpolatedIntensity2D(imp,
									(double)points.get(j).getX() - (p+1) * (points.get(j).getNormalVector()[0] * imp.getCalibration().pixelWidth),
									(double)points.get(j).getY() - (p+1) * (points.get(j).getNormalVector()[1] * imp.getCalibration().pixelHeight),
									imp.getStackIndex(channel, s+1, t+1)-1);
						}
						
						//eventually smooth normal
							if(smoothNormal){
								for(int n = 0; n < normalPointsY [s].length; n++){
									double newIntensity = normalPointsY [s][n];								
									if(n==0){
										if(normalPointsY [s].length>1){
											newIntensity += normalPointsY [s][n+1];
											newIntensity /= 2.0;
										}									
									}else if(n==normalPointsY [s].length-1){
										newIntensity += normalPointsY [s][n-1];
										newIntensity /= 2.0;
									}else{
										newIntensity += normalPointsY [s][n+1];
										newIntensity += normalPointsY [s][n-1];
										newIntensity /= 3.0;
									}
									normalPointsY [s][n] = newIntensity;	
								}
							}
						
						//gauss fit						
//							if(s == 0 && t == 0 && j == 30)	showAsPlot(normalPointsX[s],normalPointsY[s]);
							
							cf = new CurveFitter(normalPointsX [s], normalPointsY [s]);
							cf.setMaxIterations(1000);	
							cf.doFit(CurveFitter.GAUSSIAN);
													
							parameters = cf.getParams();	
							// case GAUSSIAN: p[0]+(p[1]-p[0])*Math.exp(-(x-p[2])*(x-p[2])/(2.0*p[3]*p[3]))
							rSquare= cf.getRSquared();
					        
						//filter fits
							//1. offset >= 0.0, 2. width < total calculated normal width, 3. rSquare > 0.8, 4. maximumShift<radius
							
							if(rSquare>=0.8
								&& parameters [0] >= 0.0
								&& parameters [3] <= imp.getCalibration().pixelWidth * radius * 2.0 * points.get(j).getNormalVectorLength()
								&& Math.abs(parameters [2]) <= imp.getCalibration().pixelWidth * radius * points.get(j).getNormalVectorLength()){
								
								for(int r = 0; r < 4; r++){
									resultsParameters [r][slicesIncluded] = parameters [r];
								}							
								resultsParameters [4][slicesIncluded] = s;
								slicesIncluded++;
							}else{
//								showAsPlot(normalPointsX[s],normalPointsY[s]);
								gaussMismatch++;
							}								
					}
					
					//check if one slice was successfull 
					if(slicesIncluded!=0){
						if(preventHeadFromCorr && j < preventPoints){
							points.get(j).widthGaussFit(radius, normalPointsX, normalPointsY, resultsParameters, slicesIncluded, false);	
						}else{
							points.get(j).widthGaussFit(radius, normalPointsX, normalPointsY, resultsParameters, slicesIncluded, true);	
						}						
					}else{
						if(preventHeadFromCorr && j < preventPoints){
							//Do not remove head points, if prevented
						}else{
							//remove point
							undefinedPoints += " " + j + ",";
							points.remove(j);
							points.trimToSize();
							undefined++;
						}
					}			
				} 			
				
				if(saveNormalVectorRoiSet){
					rm.runCommand("Save", saveDirectory + "NormVecRois_" + t + ".zip");	
				}	
				
				progress.addToBar(0.1*(1.0/traces.size()));
				if(gaussMismatch>0)	progress.notifyMessage("t=" + i + ": " + gaussMismatch + " fits with r < 0.8 in xy correction!",ProgressDialog.LOG);
				if(undefined>0)	progress.notifyMessage("t=" + i + ": " + undefined + " points removed as being undefined:" + undefinedPoints,ProgressDialog.LOG);			
			}
			if(i%50 == 0)	System.gc();
		}
	}
	
	public static void updateWidthFitData(ImagePlus imp, ArrayList<trace> traces, ProgressDialog progress, boolean saveNormalVectorRoiSet,
			String saveDirectory, int vectorSize, double normalRad, boolean smoothNormal, boolean preventHeadFromCorr, int preventPoints){
		if(imp.getNSlices()>4)	IJ.error("too many slices for correction via normal vector...");
		int w = imp.getWidth();
		int h = imp.getHeight();
		int calculationRadius = (int)Math.round(normalRad / imp.getCalibration().pixelWidth);
		
		RoiManager rm = RoiManager.getInstance();
		if (rm==null) rm = new RoiManager();
		
		int channel = 1;	//no multichannel implemented so far!
			
		//initialize variables
			final int halfSize = (vectorSize-(vectorSize%2))/2;
			
			String undefinedPoints;
			double x1, y1, x2, y2;
			int t;
			ArrayList<trackPoint> points;
			int nrOfPoints;
			double [] radii = new double [4];
			int radius;
			
			int gaussMismatch = 0, undefined = 0;
			double resultsParameters [][] = new double [5][4];
			int slicesIncluded = 0;
			double normalPointsX [][];
			double normalPointsY [][];
		
			CurveFitter cf;
			double [] parameters;	
			double rSquare;
			
		//calculate
		for(int i = 0; i < traces.size(); i++){
			if(traces.get(i).getTracePoints().size()>0){
				undefinedPoints = "";
				if(progress.isStopped())	return;
				progress.updateBarText("Obtaining width fit data " + (i+1) + "/" + traces.size());
				
				t = traces.get(i).getFrame();
				points = traces.get(i).getTracePoints();
				nrOfPoints = points.size();			
					
				if(saveNormalVectorRoiSet){
					rm.reset();
				}
				
				//generate vectors				
				x1 = points.get(0).getX();
				y1 = points.get(0).getY();
				x2 = points.get(0).getX();
				y2 = points.get(0).getY();
				for(int j = points.size()-1; j >= 0; j--){
					//generate vectors
					searchUpstream: for(int vu = halfSize; vu >= 0; vu--){
//						if(j-vu >= 0 && get2DDistance(points.get(j), points.get(j-vu)) < maxPointDist){ 	//TODO max DIst?
						if(j-vu >= 0){
							x1 = points.get(j-vu).getX();
							y1 = points.get(j-vu).getY();
							break searchUpstream;
						}
					}
					
					searchDownstream: for(int vd = halfSize; vd >= 0; vd--){
//						if(j+vd < nrOfPoints && get2DDistance(points.get(j), points.get(j+vd)) < maxPointDist){ 	//TODO max DIst?
						if(j+vd < nrOfPoints){
							x2 = points.get(j+vd).getX();
							y2 = points.get(j+vd).getY();
							break searchDownstream;
						}
					}					
					//TODO test for x2-x1==0 and eventually remove point
					points.get(j).setVectors(x2-x1, y2-y1);
				}
				
				gaussMismatch = 0;
				undefined = 0;				
				for(int j = nrOfPoints-1; j >= 0; j--){
					//adjust calculation radius if touching borders
					radii [0] = Math.sqrt(Math.pow(((double)(w-2)*imp.getCalibration().pixelWidth - points.get(j).getX()) / (points.get(j).getNormalVector()[0]*imp.getCalibration().pixelWidth),2.0));
//					if(radii [0] < 0)	radii [0] = Double.MAX_VALUE;
					radii [1] = Math.sqrt(Math.pow(((double)(h-2)*imp.getCalibration().pixelHeight - points.get(j).getY()) / (points.get(j).getNormalVector()[1]*imp.getCalibration().pixelHeight),2.0));
//					if(radii [1] < 0)	radii [1] = Double.MAX_VALUE;
					radii [2] = Math.sqrt(Math.pow((points.get(j).getX()) / (points.get(j).getNormalVector()[0]*imp.getCalibration().pixelWidth),2.0));
//					if(radii [2] < 0)	radii [2] = Double.MAX_VALUE;
					radii [3] = Math.sqrt(Math.pow((points.get(j).getY()) / (points.get(j).getNormalVector()[1]*imp.getCalibration().pixelHeight),2.0));
//					if(radii [3] < 0)	radii [3] = Double.MAX_VALUE;
					
					Arrays.sort(radii);
					radius = (int)Math.round(radii [0]);					
					if(radius > calculationRadius){
//						IJ.log("r1" + radius);
						radius = calculationRadius;
//						IJ.log("r2" + radius);
					}else{
//						IJ.log(" x " + points.get(j).getX() + " y " + points.get(j).getY() + " w " + w + " h " + h);
//						IJ.log(" ri " + radii [0]);
//						IJ.log(" ri " + radii [1]);
//						IJ.log(" ri " + radii [2]);
//						IJ.log(" ri " + radii [3]);
//						IJ.log(" nx " + points.get(j).getNormalVector()[0]);
//						IJ.log(" ny " + points.get(j).getNormalVector()[1]);
					}
					
					
					//save vector Roi-set
					if(saveNormalVectorRoiSet){
						Line vectorRoi = points.get(j).getVectorLine(imp);
						vectorRoi.setName("vector " + j);
						rm.addRoi(vectorRoi);	
						
						Line normalVectorRoi = points.get(j).getNormalVectorLine(imp,calculationRadius);
						normalVectorRoi.setName("normal vector " + j);
						rm.addRoi(normalVectorRoi);	
					}
					
					
					slicesIncluded = 0;
					normalPointsX = new double [imp.getNSlices()][radius*2+1];
					normalPointsY = new double [imp.getNSlices()][radius*2+1];					
					for(int s = 0; s < imp.getNSlices(); s++){
						//obtain data points for fitting						
							normalPointsX [s][radius] = 0.0;
							normalPointsY [s][radius] = impProcessing.getInterpolatedIntensity2D(imp, points.get(j).getX(),	//x 
									points.get(j).getY(),	//y
									imp.getStackIndex(channel, s+1, t+1)-1);	//intensity						
							
							for(int p = 0; p < radius; p++){
								normalPointsX [s][p + 1 + radius] = (p+1) * (points.get(j).getNormalVectorLength() * imp.getCalibration().pixelWidth);
								normalPointsY [s][p + 1 + radius] = impProcessing.getInterpolatedIntensity2D(imp,
										(double)points.get(j).getX() + (p+1) * (points.get(j).getNormalVector()[0] * imp.getCalibration().pixelWidth),
										(double)points.get(j).getY() + (p+1) * (points.get(j).getNormalVector()[1] * imp.getCalibration().pixelHeight),
										imp.getStackIndex(channel, s+1, t+1)-1);			
								
								normalPointsX [s][radius - 1 - p] = (-1) * (p+1) * (points.get(j).getNormalVectorLength() * imp.getCalibration().pixelWidth);
								normalPointsY [s][radius - 1 - p] = impProcessing.getInterpolatedIntensity2D(imp,
										(double)points.get(j).getX() - (p+1) * (points.get(j).getNormalVector()[0] * imp.getCalibration().pixelWidth),
										(double)points.get(j).getY() - (p+1) * (points.get(j).getNormalVector()[1] * imp.getCalibration().pixelHeight),
										imp.getStackIndex(channel, s+1, t+1)-1);
							}
						
						//eventually smooth normal
							if(smoothNormal){
								for(int n = 0; n < normalPointsY [s].length; n++){
									if(n==0){
										if(normalPointsY [s].length>1){
											normalPointsY [s][n] += normalPointsY [s][n+1];
											normalPointsY [s][n] /= 2.0;
										}									
									}else if(n==normalPointsY [s].length-1){
										normalPointsY [s][n] += normalPointsY [s][n-1];
										normalPointsY [s][n] /= 2.0;
									}else{
										normalPointsY [s][n] += normalPointsY [s][n+1];
										normalPointsY [s][n] += normalPointsY [s][n-1];
										normalPointsY [s][n] /= 3.0;
									}
								}
							}
						
						//gauss fit						
//							if(s == 0 && t == 0 && j == 30)	showAsPlot(normalPointsX[s],normalPointsY[s]);
							
//							if(t == 0){
//								IJ.log("s" + s + " t" + t + " j" + j + "");
//								showAsPlot(normalPointsX[s],normalPointsY[s]);
//							}
							
							cf = new CurveFitter(normalPointsX [s], normalPointsY [s]);
							cf.setMaxIterations(1000);	
							cf.doFit(CurveFitter.GAUSSIAN);
													
							parameters = cf.getParams();	
							// case GAUSSIAN: p[0]+(p[1]-p[0])*Math.exp(-(x-p[2])*(x-p[2])/(2.0*p[3]*p[3]))
							rSquare= cf.getRSquared();
					        
						//filter fits
							//1. offset >= 0.0, 2. width < total calculated normal width, 3. rSquare > 0.8, 4. maximumShift<radius
							
							if(rSquare>=0.8
								&& parameters [0] >= 0.0
								&& parameters [3] <= imp.getCalibration().pixelWidth * radius * 2.0 * points.get(j).getNormalVectorLength()
								&& Math.abs(parameters [2]) <= imp.getCalibration().pixelWidth * radius * points.get(j).getNormalVectorLength()){
								
								for(int r = 0; r < 4; r++){
									resultsParameters [r][slicesIncluded] = parameters [r];
								}							
								resultsParameters [4][slicesIncluded] = s;
								slicesIncluded++;
							}else{
	//							showAsPlot(normalPointsX[s],normalPointsY[s]);
								gaussMismatch++;
							}	
					}
					
					if(slicesIncluded!=0){
						//no adjustment!
						points.get(j).widthGaussFit(radius, normalPointsX, normalPointsY, resultsParameters, slicesIncluded, false);	
					}else{
						if(preventHeadFromCorr && j < preventPoints){
							//Do not remove head points, if prevented
						}else{
							//remove point
							undefinedPoints += " " + j + ",";
							points.remove(j);
							points.trimToSize();
							undefined++;
						}						
					}			
				} 			
				
				if(saveNormalVectorRoiSet){
					rm.runCommand("Save", saveDirectory + "NormVecRoisIter2_" + t + ".zip");	
				}
				
				progress.addToBar(0.1*(1.0/traces.size()));
				if(gaussMismatch>0)	progress.notifyMessage("t=" + i + ": " + gaussMismatch + " width fits with r < 0.8!",ProgressDialog.LOG);
				if(undefined>0)	progress.notifyMessage("t=" + i + ": " + undefined + " points removed as being undefined in width fit:" + undefinedPoints,ProgressDialog.LOG);
			}
		}	
	}

//	public static double [] projectPointToLine3D (double px, double py, double pz, double lx1, double ly1, double lz1, double lx2, double ly2, double lz2){
//		//calculate 
//		double m = (ly2-ly1)/(lx2-lx1);			
//		
//		//find dislocation of normal line (m*x+b=y)
//		double bNormal = py - ((1.0/m) * px);		
//		
//		//calculate new x
//		double x = ((ly1 - (m * lx1))-bNormal)/((1.0/m)-m);
//		
//		return new double [] {x,(1.0/m) * x + bNormal};
//		
//		//LONG WAY CALCULATION
////		//calculate 
////		double m = (ly2-ly1)/(lx2-lx1);
////		//find dislocation of line (m*x+b=y)
////		double b = ly1 - (m * lx1);
////				
////		//calculate slope of normal line
////		double mNormal = 1.0/m;
////		//find dislocation of normal line (m*x+b=y)
////		double bNormal = py - (mNormal * px);
////		
////		//calculate new x
////		double x = (b-bNormal)/(mNormal-m);
////		double y = mNormal * x + bNormal;
//	}
		
	public static void determineZInfo (ImagePlus imp, ArrayList<trace> traces, double [] slicePositions, ProgressDialog progress, boolean useNormalMax){
		if(imp.getNSlices()!=slicePositions.length)	IJ.error("wrong slice number compared to indicated slice positions for z-gauss fit...");
		if(imp.getNChannels()>1)	IJ.error("gauss fit z not implemented for > 1 channel...");
		
		double [] slicesArray = slicePositions.clone();
		Arrays.sort(slicesArray);
		final double sumOfInterplaneDistances = slicesArray [slicesArray.length-1] - slicesArray [0];
//		IJ.log("stackdepth " + slicesArray [slicesArray.length-1] + "-" + slicesArray [0]);
		
		//initialize
		double [] pX = new double [imp.getNSlices()];
		for(int s = 0; s < imp.getNSlices(); s++){
			pX [s] = slicePositions [s];
		}
		double [] pY = new double [imp.getNSlices()];
		double [] parameters;
		double rSquare;
		CurveFitter cf;
		
		for(int i = 0; i < traces.size(); i++){
			int excluded = 0;
			if(progress.isStopped()) return;			
			progress.updateBarText("Analyzing z position " + (i+1) + "/" + traces.size());
						
			for(int j = traces.get(i).getTracePoints().size()-1; j >= 0; j--){				
				//create point array
					for(int s = 0; s < imp.getNSlices(); s++){
						//get interpolated intensity
						if(useNormalMax){
							pY [s] = traces.get(i).getTracePoints().get(j).getNormalMaxIntensity(s);
						}else{
							pY [s] = traces.get(i).getTracePoints().get(j).getNormalPointIntensity(s);
						}						
					}
				
				//gauss fit in z
					cf = new CurveFitter(pX, pY);
					cf.setMaxIterations(1000);	
					cf.doFit(CurveFitter.GAUSSIAN);
					
					parameters = cf.getParams();	
					// case GAUSSIAN: p[0]+(p[1]-p[0])*Math.exp(-(x-p[2])*(x-p[2])/(2.0*p[3]*p[3]))
					rSquare = cf.getRSquared();
			       
        		//filter parameters by acceptance criteria
        			// 1. maximum within the 4 slices
        			// 2. width smaller than sum of interplane-distances
	        		if(parameters [2] <= slicesArray[slicesArray.length-1]
	        				&& parameters [2] >= slicesArray[0]
	        				&& parameters [3] <= sumOfInterplaneDistances
	        				&& rSquare > 0.8){
	        			traces.get(i).getTracePoints().get(j).setZGaussFit4P(parameters);
//	        			traces.get(i).addZGauss4PData(parameters, j);
		        			
	        		}else{
	        			excluded++;
//	        			progress.notifyMessage("t=" + i + " - p " + j + ": z gauss fit excluded"
//	        			+ " rSquare=" + rSquare
//	        			+ "	offset=" + parameters[0]
//		        		+ "	height=" + parameters[1]
//	        			+ "	center=" + parameters[2]
//						+ "	width=" + parameters[3]
//	        			+"",ProgressDialog.LOG);	 
	        			traces.get(i).getTracePoints().remove(j);
	        			traces.get(i).getTracePoints().trimToSize();
//	        			traces.get(i).removePoint(j);
	        		}	        		
			}	
			
			traces.get(i).g4pz = true;
			progress.addToBar(0.3*(1.0/traces.size()));
			if(excluded>0)	progress.notifyMessage("t=" + i + ": " + excluded + "z gauss fits excluded",ProgressDialog.LOG);
		}
	}
	
	public static void determineZInfoInWidthsOnly(ImagePlus imp, ArrayList<trace> traces, double [] slicePositions, ProgressDialog progress, String LUTPath, double widthStep){
		if(imp.getNSlices()!=slicePositions.length)	IJ.error("wrong slice number compared to indicated slice positions for z-gauss fit...");
		if(imp.getNChannels()>1)	IJ.error("gauss fit z not implemented for > 1 channel...");
		
		//read information about width LUT
		ImagePlus widthLUT = IJ.openImage(LUTPath);
		final int maxArcLength = widthLUT.getWidth();
		final int nrWidthSteps = widthLUT.getHeight();
		int widthCenter = (int)(nrWidthSteps/2);
		double minWidth = 0.0, maxWidth = 0.0, resolution = 0.0;
		
		try {
			FileReader fr = new FileReader(new File(LUTPath.substring(0, LUTPath.lastIndexOf(".")) + "_info.txt"));
			BufferedReader br = new BufferedReader(fr);
			String line = "";
			for(int i = 0; i < 9; i++){
				line = br.readLine();
			}
						
			if(line.contains("#")){
				minWidth = Double.parseDouble(line.substring(line.indexOf("#")+1, line.indexOf("-")));
				maxWidth = Double.parseDouble(line.substring(line.indexOf("-")+1, line.indexOf("+")));
				resolution = Double.parseDouble(line.substring(line.indexOf("+")+1, line.lastIndexOf("c")));
				widthCenter = Integer.parseInt(line.substring(line.lastIndexOf("c")+1));
			}
			
			br.close();
			fr.close();
		}catch (IOException e) {
			progress.notifyMessage("problem with text loading", ProgressDialog.ERROR);
			e.printStackTrace();
		}	
		
		//initialize
		double [][] widthParams;
		int excluded;		
		int counter;
		double [] distanceToWidthSteps = new double [nrWidthSteps];
		double minimumDifferencePlus; 
		int zPosOfMinimumPlus;
		double minimumDifferenceMinus; 
		int zPosOfMinimumMinus;
		
		//analysis
		for(int i = 0; i < traces.size(); i++){
			if(progress.isStopped()) return;		
			progress.updateBarText("Analyzing z position " + (i+1) + "/" + traces.size());
			
			excluded = 0;
			for(int j = traces.get(i).getTracePoints().size()-1; j >= 0; j--){
				if(progress.isStopped()) return;
				//get width fit z
	        		int arcLPosLUT = (int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(NOZ) * resolution);
	        		if(arcLPosLUT < maxArcLength){
	        			widthParams = traces.get(i).getTracePoints().get(j).getXYGaussWidth();
	        			//calculate fit to z position
	        			
	        			for(int widths = 0; widths < widthParams.length; widths++){
	        				//initialize
	        				for(int zPos = 0; zPos < nrWidthSteps; zPos++){
	        					distanceToWidthSteps [zPos] = -1.0;
	        				}
	        				
	        				//compare to width from table +-2 arclength, above center 
	        				minimumDifferencePlus = Double.MAX_VALUE;
	        				zPosOfMinimumPlus = Integer.MIN_VALUE;
	        				{
	        					for(int zPos = widthCenter; zPos < nrWidthSteps; zPos++){
		        					counter = 0;		        					
		        					for(int add = 0; add < 3 && arcLPosLUT + add < maxArcLength; add++){
		        						if(0.0 != widthLUT.getStack().getVoxel(arcLPosLUT + add,zPos,(int)Math.round(widthParams[widths][1]))){
		        							distanceToWidthSteps [zPos] += Math.sqrt(Math.pow(widthParams[widths][0] 
			        								- tools.getValueFromEncodedIntensity8bit(minWidth,maxWidth,
			        										widthLUT.getStack().getVoxel(arcLPosLUT + add,zPos,(int)Math.round(widthParams[widths][1])))
			        							,2.0));
			        						counter++;
//			        						if(j == 4 && widths == 0){
//			        							IJ.log("x " + arcLPosLUT + add 
//			        									+ "y" + zPos + " z" + (int)Math.round(widthParams[widths][1]));			
//			        							IJ.log("p " + j + " widthFrom lut " + getWidthFromIntensity(minWidth,maxWidth,
//			        										widthLUT.getStack().getVoxel(arcLPosLUT + add,zPos,(int)Math.round(widthParams[widths][1]))));
//			        							IJ.log("difference " + Math.sqrt(Math.pow(widthParams[widths][0] 
//				        								- getWidthFromIntensity(minWidth,maxWidth,
//				        										widthLUT.getStack().getVoxel(arcLPosLUT + add,zPos,(int)Math.round(widthParams[widths][1])))
//				        							,2.0)));
//			        						}
		        						}else{
//		        							if(j == 4){
//		        								IJ.log("intensity " + widthLUT.getStack().getVoxel(arcLPosLUT + add,
//				        								zPos,
//				        								(int)Math.round(widthParams[widths][1])));
//			        							IJ.log("x " + arcLPosLUT + add 
//			        									+ "y" + zPos + " z" + (int)Math.round(widthParams[widths][1]));			        							
//		        							}		        							
		        						}
		        					}
		        					for(int add = 1; add < 3 && arcLPosLUT - add >= 0; add++){
		        						if(0.0 != widthLUT.getStack().getVoxel(arcLPosLUT - add,
		        								zPos,
		        								(int)Math.round(widthParams[widths][1]))){
		        							distanceToWidthSteps [zPos] += Math.sqrt(Math.pow(widthParams[widths][0] 
			        								- tools.getValueFromEncodedIntensity8bit(minWidth,maxWidth,
			        										widthLUT.getStack().getVoxel(arcLPosLUT - add,zPos,(int)Math.round(widthParams[widths][1])))
			        								,2.0));
			        						counter++;
		        						}
		        						
		        					}   
		        					if(counter!=0){
//		        						if(j==4)	IJ.log("distToWS " + distanceToWidthSteps [zPos]);
		        						distanceToWidthSteps [zPos] += 1.0;
		        						distanceToWidthSteps [zPos] /= counter;
		        					}	        					
		        				}
		        				
		        				//find LUT field with smallest difference (above center)
		        				
		        				for(int zPos = nrWidthSteps-1; zPos >= widthCenter; zPos--){
		        					if(distanceToWidthSteps [zPos] != -1.0 && distanceToWidthSteps [zPos] < 1.0 
		        							&& distanceToWidthSteps [zPos] < minimumDifferencePlus){
		        						zPosOfMinimumPlus = zPos;
		        						minimumDifferencePlus = distanceToWidthSteps [zPos];
		        					}
		        				}
	        				}
	        				
	        				
	        				//compare to width from table +-2 arclength, above center
	        				minimumDifferenceMinus = Double.MAX_VALUE;
	        				zPosOfMinimumMinus = Integer.MIN_VALUE;
	        				{
	        					for(int zPos = widthCenter; zPos > 0; zPos--){
		        					counter = 0;
		        					for(int add = 0; add < 3 && arcLPosLUT + add < maxArcLength; add++){
//		        						IJ.log("run" + add);
		        						if(0.0 != widthLUT.getStack().getVoxel(arcLPosLUT + add,
		        								zPos,
		        								(int)Math.round(widthParams[widths][1]))){
		        							distanceToWidthSteps [zPos] += Math.sqrt(Math.pow(widthParams[widths][0] 
			        								- tools.getValueFromEncodedIntensity8bit(minWidth,maxWidth,widthLUT.getStack().getVoxel(arcLPosLUT + add,
			        									zPos,
			        									(int)Math.round(widthParams[widths][1])))
			        								,2.0));
			        						counter++;
//		        						}else{
//		        							IJ.log("problem");
		        						}      						
		        					}
		        					for(int add = 1; add < 3 && arcLPosLUT - add >= 0; add++){
		        						if(0.0 != widthLUT.getStack().getVoxel(arcLPosLUT - add,
		        								zPos,
		        								(int)Math.round(widthParams[widths][1]))){
		        							distanceToWidthSteps [zPos] += Math.sqrt(Math.pow(widthParams[widths][0] 
			        								- tools.getValueFromEncodedIntensity8bit(minWidth,maxWidth,widthLUT.getStack().getVoxel(arcLPosLUT - add,
			        										zPos,
			        									(int)Math.round(widthParams[widths][1])))
			        								,2.0));
			        						counter++;
//		        						}else{
//		        							IJ.log("problem");
		        						}
		        						
		        					}   
		        					if(counter!=0){
//		        						if(j==4)	IJ.log("distToWS " + distanceToWidthSteps [zPos]);
		        						distanceToWidthSteps [zPos] += 1.0;
		        						distanceToWidthSteps [zPos] /= counter;
		        					}else{
//		        						IJ.log("counter 0");
		        					}
		        				}
		        				
		        				//find LUT field with smallest difference (above center)
		        				for(int zPos = 0; zPos <= widthCenter; zPos++){
		        					if(distanceToWidthSteps [zPos] != -1.0 && distanceToWidthSteps [zPos] < 1.0 
		        							&& distanceToWidthSteps [zPos] < minimumDifferenceMinus){
		        						zPosOfMinimumMinus = zPos;
		        						minimumDifferenceMinus = distanceToWidthSteps [zPos];
		        					}
		        				}
		        				
		        				
	        				}
	        				
	        				//save z based on LUT field with smallest difference (above center)
	        				if(zPosOfMinimumMinus != Integer.MIN_VALUE && zPosOfMinimumPlus != Integer.MIN_VALUE){
//	        					IJ.log(widthCenter + " r " + (widthCenter - zPosOfMinimumMinus));
	        					traces.get(i).getTracePoints().get(j).addZWidthFit(slicePositions[(int)Math.round(widthParams[widths][1])] 
	        								+ ((widthCenter - zPosOfMinimumMinus) * widthStep),
	        							slicePositions[(int)Math.round(widthParams[widths][1])] 
	        								- ((zPosOfMinimumPlus - widthCenter) * widthStep));
//	        				}else if(zPosOfMinimumMinus == Integer.MIN_VALUE && zPosOfMinimumPlus != Integer.MIN_VALUE){
//	        					traces.get(i).getTracePoints().get(j).addZWidthFit(slicePositions[(int)Math.round(widthParams[widths][1])] 
//        								+ ((zPosOfMinimumPlus - widthCenter) * widthStep),
//        							slicePositions[(int)Math.round(widthParams[widths][1])] 
//        								- ((zPosOfMinimumPlus - widthCenter) * widthStep));
//	        					keepPoint = true;
//	        				}else if(zPosOfMinimumMinus != Integer.MIN_VALUE && zPosOfMinimumPlus == Integer.MIN_VALUE){
//	        					traces.get(i).getTracePoints().get(j).addZWidthFit(slicePositions[(int)Math.round(widthParams[widths][1])] 
//        								+ ((widthCenter - zPosOfMinimumMinus) * widthStep),
//    								slicePositions[(int)Math.round(widthParams[widths][1])] 
//            								+ ((widthCenter - zPosOfMinimumMinus) * widthStep));
//	        					keepPoint = true;
//	        				}else{
//	        					IJ.log("no width match...");
	        				}
	        			}        			
	        		}
	        		
	        	//eventually exclude point
//	        		if(i == 0){
//						IJ.log("p" + j);
//						IJ.log("   z " + traces.get(i).getTracePoints().get(j).getZ(G4PZ));
//					}
	        		if(traces.get(i).getTracePoints().get(j).zDeterminable() == false){
//	        				|| traces.get(i).getTracePoints().get(j).getZ(G4PZ) < -10000.0			//TODO find cause for these values
//	        				|| traces.get(i).getTracePoints().get(j).getZ(G4PZ) > 10000.0){
//	        			if(traces.get(i).getTracePoints().get(j).getZ(G4PZ) < -10000.0)	IJ.log("t" + traces.get(i).getFrame() + " p" + j + " <-10000");
//	        			if(traces.get(i).getTracePoints().get(j).getZ(G4PZ) > 10000.0)	IJ.log("t" + traces.get(i).getFrame() + " p" + j + " >10000");
	        			excluded++;
	        			traces.get(i).getTracePoints().remove(j);
	        			traces.get(i).getTracePoints().trimToSize();
//	        			IJ.log("excluded! " + j);
	        		}	  
			}	
			
			traces.get(i).g4pz = true;
			if(i%50 == 0)	System.gc();
			progress.addToBar(0.3*(1.0/traces.size()));
			if(excluded>0)	progress.notifyMessage("t=" + i + ": " + excluded + " z fits excluded",ProgressDialog.LOG);
		}		
		widthLUT.close();
	
	}
	
	public static void determineZInfoIncWidth (ImagePlus imp, ArrayList<trace> traces, double [] slicePositions, ProgressDialog progress, boolean useNormalMax, String LUTPath, double widthStep){
		if(imp.getNSlices()!=slicePositions.length)	IJ.error("wrong slice number compared to indicated slice positions for z-gauss fit...");
		if(imp.getNChannels()>1)	IJ.error("gauss fit z not implemented for > 1 channel...");
		
		double [] slicesArray = Arrays.copyOf(slicePositions, slicePositions.length);
		Arrays.sort(slicesArray);
		double sumOfInterplaneDistances = slicesArray [slicesArray.length-1] - slicesArray [0];
//		IJ.log("stackdepth " + slicesArray [slicesArray.length-1] + "-" + slicesArray [0]);
		//read information about width LUT
		ImagePlus widthLUT = IJ.openImage(LUTPath);
		final int maxArcLength = widthLUT.getWidth();
		final int nrWidthSteps = widthLUT.getHeight();
		int widthCenter = (int)(nrWidthSteps/2);
		double minWidth = 0.0, maxWidth = 0.0, resolution = 0.0;
		
		try {
			FileReader fr = new FileReader(new File(LUTPath.substring(0, LUTPath.lastIndexOf(".")) + "_info.txt"));
			BufferedReader br = new BufferedReader(fr);
			String line = "";
			for(int i = 0; i < 9; i++){
				line = br.readLine();
			}
						
			if(line.contains("#")){
				minWidth = Double.parseDouble(line.substring(line.indexOf("#")+1, line.indexOf("-")));
				maxWidth = Double.parseDouble(line.substring(line.indexOf("-")+1, line.indexOf("+")));
				resolution = Double.parseDouble(line.substring(line.indexOf("+")+1, line.lastIndexOf("c")));
				widthCenter = Integer.parseInt(line.substring(line.lastIndexOf("c")+1));
			}
			
			br.close();
			fr.close();
		}catch (IOException e) {
			progress.notifyMessage("problem with text loading", ProgressDialog.ERROR);
			e.printStackTrace();
		}	
		
		//initialize
		double [] pX = new double [imp.getNSlices()];
		for(int s = 0; s < imp.getNSlices(); s++){
			pX [s] = slicePositions [s];
		}
		double [] pY = new double [imp.getNSlices()];
		CurveFitter cf;
		double [] parameters;
		double rSquare;
		double [][] widthParams;
		int excluded;		
		int counter;
		double [] distanceToWidthSteps = new double [nrWidthSteps];
		double minimumDifferencePlus; 
		int zPosOfMinimumPlus;
		double minimumDifferenceMinus; 
		int zPosOfMinimumMinus;
		
		//analysis
		for(int i = 0; i < traces.size(); i++){
			if(progress.isStopped()) return;		
			progress.updateBarText("Analyzing z position " + (i+1) + "/" + traces.size());
			
			excluded = 0;
			for(int j = traces.get(i).getTracePoints().size()-1; j >= 0; j--){
				if(progress.isStopped()) return;
				//create point array
					for(int s = 0; s < imp.getNSlices(); s++){						
						//get interpolated intensity
						if(useNormalMax){
							pY [s] = traces.get(i).getTracePoints().get(j).getNormalMaxIntensity(s);
						}else{
							pY [s] = traces.get(i).getTracePoints().get(j).getNormalPointIntensity(s);
						}						
					}
				
				//gauss fit in z
					cf = new CurveFitter(pX, pY);
					cf.setMaxIterations(1000);	
					cf.doFit(CurveFitter.GAUSSIAN);
					
					parameters = cf.getParams();	
					// case GAUSSIAN: p[0]+(p[1]-p[0])*Math.exp(-(x-p[2])*(x-p[2])/(2.0*p[3]*p[3]))
					rSquare= cf.getRSquared();
			       
        		//fiter fits
					// 1. maximum within the 4 slices
        			// 2. width smaller than sum of interplane-distances
	        		if(parameters [2] <= slicesArray[slicesArray.length-1]
	        				&& parameters [2] >= slicesArray[0]
	        				&& parameters [3] <= sumOfInterplaneDistances
	        				&& rSquare > 0.8){
	        			traces.get(i).getTracePoints().get(j).setZGaussFit4P(parameters);
//	        		}else{
//	        			progress.notifyMessage("t=" + i + " - p " + j + ": z gauss fit excluded"
//	    	        			+ " rSquare=" + rSquare
//	    	        			+ "	offset=" + parameters[0]
//	    		        		+ "	height=" + parameters[1]
//	    	        			+ "	center=" + parameters[2]
//	    						+ "	width=" + parameters[3]
//	    	        			+"",ProgressDialog.LOG);	 
	        		}
	        		
	        		
	        	//get width fit z
	        		int arcLPosLUT = (int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(NOZ) * resolution);
	        		if(arcLPosLUT < maxArcLength){
	        			widthParams = traces.get(i).getTracePoints().get(j).getXYGaussWidth();
	        			//calculate fit to z position
	        			
	        			for(int widths = 0; widths < widthParams.length; widths++){
	        				//initialize
	        				for(int zPos = 0; zPos < nrWidthSteps; zPos++){
	        					distanceToWidthSteps [zPos] = -1.0;
	        				}
	        				
	        				//compare to width from table +-2 arclength, above center 
	        				minimumDifferencePlus = Double.MAX_VALUE;
	        				zPosOfMinimumPlus = Integer.MIN_VALUE;
	        				{
	        					for(int zPos = widthCenter; zPos < nrWidthSteps; zPos++){
		        					counter = 0;		        					
		        					for(int add = 0; add < 3 && arcLPosLUT + add < maxArcLength; add++){
		        						if(0.0 != widthLUT.getStack().getVoxel(arcLPosLUT + add,zPos,(int)Math.round(widthParams[widths][1]))){
		        							distanceToWidthSteps [zPos] += Math.sqrt(Math.pow(widthParams[widths][0] 
			        								- tools.getValueFromEncodedIntensity8bit(minWidth,maxWidth,
			        										widthLUT.getStack().getVoxel(arcLPosLUT + add,zPos,(int)Math.round(widthParams[widths][1])))
			        							,2.0));
			        						counter++;
//			        						if(j == 4 && widths == 0){
//			        							IJ.log("x " + arcLPosLUT + add 
//			        									+ "y" + zPos + " z" + (int)Math.round(widthParams[widths][1]));			
//			        							IJ.log("p " + j + " widthFrom lut " + getWidthFromIntensity(minWidth,maxWidth,
//			        										widthLUT.getStack().getVoxel(arcLPosLUT + add,zPos,(int)Math.round(widthParams[widths][1]))));
//			        							IJ.log("difference " + Math.sqrt(Math.pow(widthParams[widths][0] 
//				        								- getWidthFromIntensity(minWidth,maxWidth,
//				        										widthLUT.getStack().getVoxel(arcLPosLUT + add,zPos,(int)Math.round(widthParams[widths][1])))
//				        							,2.0)));
//			        						}
		        						}else{
//		        							if(j == 4){
//		        								IJ.log("intensity " + widthLUT.getStack().getVoxel(arcLPosLUT + add,
//				        								zPos,
//				        								(int)Math.round(widthParams[widths][1])));
//			        							IJ.log("x " + arcLPosLUT + add 
//			        									+ "y" + zPos + " z" + (int)Math.round(widthParams[widths][1]));			        							
//		        							}		        							
		        						}
		        					}
		        					for(int add = 1; add < 3 && arcLPosLUT - add >= 0; add++){
		        						if(0.0 != widthLUT.getStack().getVoxel(arcLPosLUT - add,
		        								zPos,
		        								(int)Math.round(widthParams[widths][1]))){
		        							distanceToWidthSteps [zPos] += Math.sqrt(Math.pow(widthParams[widths][0] 
			        								- tools.getValueFromEncodedIntensity8bit(minWidth,maxWidth,
			        										widthLUT.getStack().getVoxel(arcLPosLUT - add,zPos,(int)Math.round(widthParams[widths][1])))
			        								,2.0));
			        						counter++;
		        						}
		        						
		        					}   
		        					if(counter!=0){
//		        						if(j==4)	IJ.log("distToWS " + distanceToWidthSteps [zPos]);
		        						distanceToWidthSteps [zPos] += 1.0;
		        						distanceToWidthSteps [zPos] /= counter;
		        					}	        					
		        				}
		        				
		        				//find LUT field with smallest difference (above center)
		        				
		        				for(int zPos = nrWidthSteps-1; zPos >= widthCenter; zPos--){
		        					if(distanceToWidthSteps [zPos] != -1.0 && distanceToWidthSteps [zPos] < 1.0 
		        							&& distanceToWidthSteps [zPos] < minimumDifferencePlus){
		        						zPosOfMinimumPlus = zPos;
		        						minimumDifferencePlus = distanceToWidthSteps [zPos];
		        					}
		        				}
	        				}
	        				
	        				
	        				//compare to width from table +-2 arclength, above center
	        				minimumDifferenceMinus = Double.MAX_VALUE;
	        				zPosOfMinimumMinus = Integer.MIN_VALUE;
	        				{
	        					for(int zPos = widthCenter; zPos > 0; zPos--){
		        					counter = 0;
		        					for(int add = 0; add < 3 && arcLPosLUT + add < maxArcLength; add++){
//		        						IJ.log("run" + add);
		        						if(0.0 != widthLUT.getStack().getVoxel(arcLPosLUT + add,
		        								zPos,
		        								(int)Math.round(widthParams[widths][1]))){
		        							distanceToWidthSteps [zPos] += Math.sqrt(Math.pow(widthParams[widths][0] 
			        								- tools.getValueFromEncodedIntensity8bit(minWidth,maxWidth,widthLUT.getStack().getVoxel(arcLPosLUT + add,
			        									zPos,
			        									(int)Math.round(widthParams[widths][1])))
			        								,2.0));
			        						counter++;
//		        						}else{
//		        							IJ.log("problem");
		        						}      						
		        					}
		        					for(int add = 1; add < 3 && arcLPosLUT - add >= 0; add++){
		        						if(0.0 != widthLUT.getStack().getVoxel(arcLPosLUT - add,
		        								zPos,
		        								(int)Math.round(widthParams[widths][1]))){
		        							distanceToWidthSteps [zPos] += Math.sqrt(Math.pow(widthParams[widths][0] 
			        								- tools.getValueFromEncodedIntensity8bit(minWidth,maxWidth,widthLUT.getStack().getVoxel(arcLPosLUT - add,
			        										zPos,
			        									(int)Math.round(widthParams[widths][1])))
			        								,2.0));
			        						counter++;
//		        						}else{
//		        							IJ.log("problem");
		        						}
		        						
		        					}   
		        					if(counter!=0){
//		        						if(j==4)	IJ.log("distToWS " + distanceToWidthSteps [zPos]);
		        						distanceToWidthSteps [zPos] += 1.0;
		        						distanceToWidthSteps [zPos] /= counter;
		        					}else{
//		        						IJ.log("counter 0");
		        					}
		        				}
		        				
		        				//find LUT field with smallest difference (above center)
		        				for(int zPos = 0; zPos <= widthCenter; zPos++){
		        					if(distanceToWidthSteps [zPos] != -1.0 && distanceToWidthSteps [zPos] < 1.0 
		        							&& distanceToWidthSteps [zPos] < minimumDifferenceMinus){
		        						zPosOfMinimumMinus = zPos;
		        						minimumDifferenceMinus = distanceToWidthSteps [zPos];
		        					}
		        				}
		        				
		        				
	        				}
	        				
	        				//save z based on LUT field with smallest difference (above center)
	        				if(zPosOfMinimumMinus != Integer.MIN_VALUE && zPosOfMinimumPlus != Integer.MIN_VALUE){
//	        					IJ.log(widthCenter + " r " + (widthCenter - zPosOfMinimumMinus));
	        					traces.get(i).getTracePoints().get(j).addZWidthFit(slicePositions[(int)Math.round(widthParams[widths][1])] 
	        								+ ((widthCenter - zPosOfMinimumMinus) * widthStep),
	        							slicePositions[(int)Math.round(widthParams[widths][1])] 
	        								- ((zPosOfMinimumPlus - widthCenter) * widthStep));
//	        				}else if(zPosOfMinimumMinus == Integer.MIN_VALUE && zPosOfMinimumPlus != Integer.MIN_VALUE){
//	        					traces.get(i).getTracePoints().get(j).addZWidthFit(slicePositions[(int)Math.round(widthParams[widths][1])] 
//        								+ ((zPosOfMinimumPlus - widthCenter) * widthStep),
//        							slicePositions[(int)Math.round(widthParams[widths][1])] 
//        								- ((zPosOfMinimumPlus - widthCenter) * widthStep));
//	        					keepPoint = true;
//	        				}else if(zPosOfMinimumMinus != Integer.MIN_VALUE && zPosOfMinimumPlus == Integer.MIN_VALUE){
//	        					traces.get(i).getTracePoints().get(j).addZWidthFit(slicePositions[(int)Math.round(widthParams[widths][1])] 
//        								+ ((widthCenter - zPosOfMinimumMinus) * widthStep),
//    								slicePositions[(int)Math.round(widthParams[widths][1])] 
//            								+ ((widthCenter - zPosOfMinimumMinus) * widthStep));
//	        					keepPoint = true;
//	        				}else{
//	        					IJ.log("no width match...");
	        				}
	        			}        			
	        		}
	        		
	        	//eventually exclude point
//	        		if(i == 0){
//						IJ.log("p" + j);
//						IJ.log("   z " + traces.get(i).getTracePoints().get(j).getZ(G4PZ));
//					}
	        		if(traces.get(i).getTracePoints().get(j).zDeterminable() == false){
//	        				|| traces.get(i).getTracePoints().get(j).getZ(G4PZ) < -10000.0			//TODO find cause for these values
//	        				|| traces.get(i).getTracePoints().get(j).getZ(G4PZ) > 10000.0){
//	        			if(traces.get(i).getTracePoints().get(j).getZ(G4PZ) < -10000.0)	IJ.log("t" + traces.get(i).getFrame() + " p" + j + " <-10000");
//	        			if(traces.get(i).getTracePoints().get(j).getZ(G4PZ) > 10000.0)	IJ.log("t" + traces.get(i).getFrame() + " p" + j + " >10000");
	        			excluded++;
	        			traces.get(i).getTracePoints().remove(j);
	        			traces.get(i).getTracePoints().trimToSize();
//	        			IJ.log("excluded! " + j);
	        		}	  
			}	
			
			traces.get(i).g4pz = true;
			if(i%50 == 0)	System.gc();
			progress.addToBar(0.3*(1.0/traces.size()));
			if(excluded>0)	progress.notifyMessage("t=" + i + ": " + excluded + " z fits excluded",ProgressDialog.LOG);
		}		
		widthLUT.close();
	}
	
	public static double [] getWidthZCorrectionFactor (ArrayList<trace> traces, ProgressDialog progress) {
		/**
		 * This method compares the z value obtained by gaussian fit over z-axis and the z value obtained using at width-to-Z LUT to calculate
		 * a correction factor and apply this to further calculations of the z-value
		 * */
		
		int [] nrOfCorrValues = new int [10];
		Arrays.fill(nrOfCorrValues, 0);
		int index;
		for(int i = 0; i < traces.size(); i++){
			for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
				index = (int)(10.0*((double)j/(double)traces.get(i).getTracePoints().size()));
				nrOfCorrValues [index] ++;
			}			
		}
		Arrays.sort(nrOfCorrValues);
		
		double [][] corrValues = new double [10][nrOfCorrValues[nrOfCorrValues.length-1]];
		int [] corrValueCounter = new int [10];
		Arrays.fill(corrValueCounter, 0);		
		
		double corr;
		for(int i = 0; i < traces.size(); i++){
			if(progress.isStopped()) return new double [] {0.0, 0.0};
			progress.updateBarText("Calibrating width-z-calculation " + (i+1) + "/" + traces.size());
			
//			traces.get(i).getTracePoints().trimToSize();
			for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
				traces.get(i).getTracePoints().get(j).getZ(PUREZ);
				corr = traces.get(i).getTracePoints().get(j).calculateWidthCorrFactor();
				if(corr != 0.0){
					index = (int)(10.0*((double)j/(double)traces.get(i).getTracePoints().size()));
					corrValues [index][corrValueCounter[index]] = corr;
					corrValueCounter [index]++;
				}				 
			}
		}
		
		//trim array
//		IJ.log("array " + corrValues.length + " counter " + corrValueCounter);
		double [] corrFactor = new double [10]; 
		Arrays.fill(corrFactor, 1.0);
		double [] medianArray;
		for(int j = 0; j < 10; j++){
//			IJ.log(j + " lo:" + corrValues [j].length);
			medianArray = Arrays.copyOfRange(corrValues [j], 0, corrValueCounter [j]);
//			IJ.log(j + " ln:" + medianArray.length + " ct " + corrValueCounter [j]);
			if(medianArray.length > 1){
				corrFactor [j] = tools.getMedian(medianArray);
			}else{
				corrFactor [j] = 0.0;
			}
		}
		
		//fill undefined ranges using calculations for neighbored ranges
		for(int j = 0; j < 10; j++){
			if(corrFactor [j] == 0.0){
				corr = 0.0;
				fillingHoles: for(int l = 1; l < 10; l++){
					if(j+l<10){
						corr = corrFactor [j+l];
					}
					if(j-l>0){
						if(corr != 0.0){
							corr += corrFactor [j-l];
							corr /= 2.0;
						}else{
							corr = corrFactor [j-l];
						}
					}
					if(corr != 0.0){
						corrFactor [j] = corr;
						break fillingHoles;
					}
				}
			}
		}
		
		//write 1.0 if 0.0
		for(int j = 0; j < 10; j++){
			if(corrFactor [j] == 0.0){
				corrFactor [j] = 1.0;
			}
		}
		
		for(int i = 0; i < traces.size(); i++){
			for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
				traces.get(i).getTracePoints().get(j).setWidhtZCorrectionFactor(corrFactor [(int)(10.0*((double)j/(double)traces.get(i).getTracePoints().size()))]);
			}
		}
		
		return corrFactor;
	}
	
	public static void interpolateZLinear(ArrayList<trace> traces, ProgressDialog progress, double acceptedDist, int plusMinusPoints){		
		for(int i = 0; i < traces.size(); i++){
			if(traces.get(i).getTracePoints().size()>0){
				if(progress.isStopped()) return;			
				progress.updateBarText("Interpolate z position " + (i+1) + "/" + traces.size());				
				traces.get(i).interpolateZLinear(acceptedDist, plusMinusPoints);				
				progress.addToBar(0.2*(1.0/traces.size()));
			}			
		}
	}
	
	/**
	 * Filters the ArrayList of traces and ecludes traces which are not oriented or do not contain points
	 * */
	public static void removeProblematicTraces(ProgressDialog progress, ArrayList<trace> traces){
		for(int i = 0; i < traces.size(); i++){  			
  			if(traces.get(i).getTracePoints().size() > 0 && traces.get(i).oriented){				  	
			}else{
				progress.notifyMessage("Trace " + traces.get(i).getFrame() + " needs to be removed (no points present)", ProgressDialog.LOG);
				traces.remove(i);
  				traces.trimToSize();
			}
  		}	
	}
	
	public static void saveTraceImage(ImagePlus imp, ArrayList<trace> traces, int encoding, double zMin, double zMax, String path, int preventPoints){
		if(imp.getNChannels()>1)	IJ.error("generation of trace image not implemented for multi-channel stacks...");
		int channel = 0;
		String ending = "_ti";
		
		double [][][][] saveImage = new double [imp.getWidth()][imp.getHeight()][imp.getNFrames()][2];
		for(int x = 0; x < imp.getWidth(); x++){
			for(int y = 0; y < imp.getHeight(); y++){
				for(int t = 0; t < imp.getNFrames(); t++){
					saveImage [x][y][t][0] = 0.0;
					saveImage [x][y][t][1] = 0.0;
				}
			}
		}
		
		ImagePlus traceImp = IJ.createHyperStack("trace image", imp.getWidth(), imp.getHeight(), imp.getNChannels(), 1, imp.getNFrames(), 8);
		traceImp.setCalibration(imp.getCalibration());
		final double maxIntensity = 255.0;
		
		//save positions
		if (encoding == NOZ){
			if(!traces.get(0).xyCorrected){
				ending = "_ti_rawskl";	
			}		
			
			for(int i = 0; i < traces.size(); i++){
				int t = traces.get(i).getFrame();
				for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
					int z = traceImp.getStackIndex(channel+1, 1, t+1)-1;
					saveImage[(int)Math.round(traces.get(i).getTracePoints().get(j).getX()/imp.getCalibration().pixelWidth)]
							[(int)Math.round(traces.get(i).getTracePoints().get(j).getY()/imp.getCalibration().pixelHeight)]
							[z][0] = maxIntensity;
					saveImage[(int)Math.round(traces.get(i).getTracePoints().get(j).getX()/imp.getCalibration().pixelWidth)]
							[(int)Math.round(traces.get(i).getTracePoints().get(j).getY()/imp.getCalibration().pixelHeight)]
							[z][1] = 1.0;
					
				}			
			}
		}else{
			//get z min / max
			double  min = Double.MAX_VALUE, max = 0.0;
			
			for(int i = 0; i < traces.size(); i++){
				if(traces.get(i).getTracePoints().size() > 0 && traces.get(i).oriented){	//Double.isNaN(traces.get(i).getThetaDegree(encoding))==false
					for(int j = preventPoints; j < traces.get(i).getTracePoints().size(); j++){
						if(traces.get(i).getTracePoints().get(j).getZ(encoding) > max){
							max = traces.get(i).getTracePoints().get(j).getZ(encoding);
						}
						if(traces.get(i).getTracePoints().get(j).getZ(encoding) < min){
							min = traces.get(i).getTracePoints().get(j).getZ(encoding);
						}
					}
				}
			}
						
			ending = "_ti_zC";
			if(encoding == MEDIANZ){
				ending += "median";
			}else if(encoding == MEANZ){
				ending += "mean";
			}
			
			for(int i = 0; i < traces.size(); i++){
				int t = traces.get(i).getFrame();				
				for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
					int z = traceImp.getStackIndex(channel+1, 1, t+1)-1;
					saveImage[(int)Math.round(traces.get(i).getTracePoints().get(j).getX()/imp.getCalibration().pixelWidth)]
							[(int)Math.round(traces.get(i).getTracePoints().get(j).getY()/imp.getCalibration().pixelHeight)]
							[z][0] += tools.getEncodedIntensity8bit(traces.get(i).getTracePoints().get(j).getZ(encoding),min,max);
					saveImage[(int)Math.round(traces.get(i).getTracePoints().get(j).getX()/imp.getCalibration().pixelWidth)]
							[(int)Math.round(traces.get(i).getTracePoints().get(j).getY()/imp.getCalibration().pixelHeight)]
							[z][1] += 1.0;
				}			
			}
			
			//save image metadata
			TextPanel tp = new TextPanel("Metadata");
			tp.append("Image information for image:	" + path.substring(path.lastIndexOf(System.getProperty("file.separator"))) + ending + ".tif");
			tp.append("");
			tp.append("dimension		minimum	maximum");
			tp.append("gray value (1.0-254.0):	z information [um]	" + constants.df6US.format(min) + "	" + constants.df6US.format(max));
			tp.append("");
			tp.append("code#"+min+"-"+max+"");
			addFooter(tp);
		  	tp.saveAs(path + ending + "_info.txt");
		}		
		
		//write to image
		for(int x = 0; x < imp.getWidth(); x++){
			for(int y = 0; y < imp.getHeight(); y++){
				for(int t = 0; t < imp.getNFrames(); t++){
					int z = traceImp.getStackIndex(channel+1, 1, t+1)-1;
					if(saveImage[x][y][z][1]>0.0){						
						traceImp.getStack().setVoxel(x, y, z, saveImage[x][y][z][0]/saveImage[x][y][z][1]); 
					}
				}
			}
		}
		
		IJ.run(traceImp, "Fire", "");		
		IJ.saveAsTiff(traceImp, path + ending + ".tif");
		traceImp.changes = false;
		traceImp.close();
	}

	public static void saveXYWidthGraph(ArrayList<trace> traces, double resolution, String path){	//resolution = scaling factor
		int slices = 4;
		
		//get parameters
			double length = 0.0;
			double min = Double.POSITIVE_INFINITY;
			double max = 0.0;
			int tMax = 0;
			for(int i = 0; i < traces.size(); i++){
				if(traces.get(i).getFrame() > tMax)	tMax = traces.get(i).getFrame();
				if(traces.get(i).getTracePoints().size() > 0){	//Double.isNaN(traces.get(i).getThetaDegree(encoding))==false
					for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
						if(traces.get(i).getTracePoints().get(j).getArcLengthPos(NOZ)> length){
							length = traces.get(i).getTracePoints().get(j).getArcLengthPos(NOZ);
						}
					}
				}				
				if(traces.get(i).getMinWidth() < min)	min = traces.get(i).getMinWidth();
				if(traces.get(i).getMaxWidth() > max)	max = traces.get(i).getMaxWidth();
			}
						
		//create Image
			int width = (int)Math.round(length * resolution)+1;
			ImagePlus codeImp = IJ.createImage("Width graph", width, tMax+1, slices, 8);
//			ImagePlus codeImp2 = IJ.createImage("Width graph 2", width, (tMax+1)*5, slices, 8);
			
			double values [][][] = new double [width][slices][2];
			for(int i = 0; i < traces.size(); i++){
				for(int a = 0; a < width; a++){
					for(int s = 0; s < slices; s++){
						values [a][s][0] = 0.0;
						values [a][s][1] = 0.0;
					}
				}
								
				double [][] widths;
				for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
					widths = traces.get(i).getTracePoints().get(j).getXYGaussWidth();
					for(int s = 0; s < widths.length; s++){
//						try{
							values [(int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(NOZ) * resolution)]
									[(int)Math.round(widths [s][1])][0] += tools.getEncodedIntensity8bit (widths [s][0], min, max);
							values [(int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(NOZ) * resolution)]
									[(int)Math.round(widths [s][1])][1] += 1.0;
//						}catch(Exception e){
//							IJ.log("alp " + (int)Math.round(points.get(j).getArcLengthPos(NOZ) * resolution));
//							IJ.log("maxalp " + length + " *r " + width);
//							IJ.log("s" + (int)Math.round(widths [s][1]));
//						}						
					}				
				}
				for(int a = 0; a < width; a++){
					for(int s = 0; s < slices; s++){
						if(values[a][s][1]!=0.0){
//							if(values[a][s][1]>1.0)	IJ.log(">1.0");
							codeImp.getStack().setVoxel(a, traces.get(i).getFrame(), s, (values[a][s][0] / values[a][s][1]));
//							for(int m = 0; m < 5; m++){
//								codeImp2.getStack().setVoxel(a, t*5+m, s, (values[a][s][0] / values[a][s][1]));
//							}	
						}														
					}
				}	
			}
			codeImp.getCalibration().pixelWidth = (1.0/resolution);
			IJ.run(codeImp, "Fire", "");
//			IJ.run(codeImp2, "Fire", "");
			
			IJ.saveAsTiff(codeImp, path + "_w.tif");
			codeImp.changes = false;
			codeImp.close();
			
//			IJ.saveAsTiff(codeImp2, path + "_w_5fT.tif");
//			codeImp2.changes = false;
//			codeImp2.close();
			
		//save image metadata
			TextPanel tp = new TextPanel("Metadata");
			tp.append("Image information for image:	" + path.substring(path.lastIndexOf(System.getProperty("file.separator"))) + "_w.tif");
			tp.append("");
			tp.append("dimension		minimum	maximum");
			tp.append("x axis:	xy arc length	0	" + constants.df0.format((int)Math.round(length*resolution)+1) + "	scaling: " + constants.df6US.format(resolution) + "x");
			tp.append("y axis:	time [frame]	0	" + constants.df0.format(tMax));
			tp.append("z axis:	slice	0	" + constants.df0.format(slices-1));
			tp.append("gray value (1.0-254.0):	xy gauss fit width	" + constants.df6US.format(min) + "	" + constants.df6US.format(max));
			tp.append("");
			tp.append("code#"+min+"-"+max+"+"+resolution);
			addFooter(tp);
		  	tp.saveAs(path + "_w_info.txt");
	}
	
	public static void saveGaussMaxGraphs(ArrayList<trace> traces, double calibration, String path){	//resolution = scaling factor
		int slices = 4;
		double resolution = (1.0/calibration);
		//get parameters
			double length = 0.0;			
			int tMax = 0;
			for(int i = 0; i < traces.size(); i++){
				if(traces.get(i).getFrame() > tMax)	tMax = traces.get(i).getFrame();
				if(traces.get(i).getTracePoints().size() > 0){	//Double.isNaN(traces.get(i).getThetaDegree(encoding))==false
					for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
						if(traces.get(i).getTracePoints().get(j).getArcLengthPos(NOZ)> length){
							length = traces.get(i).getTracePoints().get(j).getArcLengthPos(NOZ);
						}
					}
				}
			}
						
		//create Image
			int width = (int)Math.round(length * resolution)+1;
			ImagePlus codeImp = IJ.createImage("Height graph", width, tMax+1, slices, 16);
			double table [][][] = new double [slices][width][tMax+1];
			for(int a = 0; a < width; a++){
				for(int s = 0; s < slices; s++){
					Arrays.fill(table[s][a], Double.NaN);
				}
			}
			
			double values [][][] = new double [width][slices][2];
			for(int i = 0; i < traces.size(); i++){
				for(int a = 0; a < width; a++){
					for(int s = 0; s < slices; s++){
						values [a][s][0] = 0.0;
						values [a][s][1] = 0.0;
					}
				}
								
				double [][] heights;
				for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
					heights = traces.get(i).getTracePoints().get(j).getXYGaussHeight();
					for(int s = 0; s < heights.length; s++){
//						try{
							values [(int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(NOZ) * resolution)]
									[(int)Math.round(heights [s][1])][0] += heights [s][0];
							values [(int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(NOZ) * resolution)]
									[(int)Math.round(heights [s][1])][1] += 1.0;
//						}catch(Exception e){
//							IJ.log("alp " + (int)Math.round(points.get(j).getArcLengthPos(NOZ) * resolution));
//							IJ.log("maxalp " + length + " *r " + width);
//							IJ.log("s" + (int)Math.round(widths [s][1]));
//						}						
					}				
				}
				for(int a = 0; a < width; a++){
					for(int s = 0; s < slices; s++){
						if(values[a][s][1]!=0.0){
//							if(values[a][s][1]>1.0)	IJ.log(">1.0");
							codeImp.getStack().setVoxel(a, traces.get(i).getFrame(), s, (values[a][s][0] / values[a][s][1]));
							table [s][a][traces.get(i).getFrame()] = (values[a][s][0] / values[a][s][1]);
//							for(int m = 0; m < 5; m++){
//								codeImp2.getStack().setVoxel(a, t*5+m, s, (values[a][s][0] / values[a][s][1]));
//							}	
						}														
					}
				}	
			}
			codeImp.getCalibration().pixelWidth = (1.0/resolution);
			IJ.run(codeImp, "Fire", "");
			
			IJ.saveAsTiff(codeImp, path + "_GaussH.tif");
			codeImp.changes = false;
			codeImp.close();
						
		//save image metadata
			TextPanel tp = new TextPanel("Metadata");
			tp.append("Image information for image:	" + path.substring(path.lastIndexOf(System.getProperty("file.separator"))) + "_w.tif");
			tp.append("");
			tp.append("dimension		minimum	maximum");
			tp.append("x axis:	xy arc length	0	" + constants.df0.format((int)Math.round(length*resolution)+1) + "	scaling: " + constants.df6US.format(resolution) + "x");
			tp.append("y axis:	time [frame]	0	" + constants.df0.format(tMax));
			tp.append("z axis:	slice	0	" + constants.df0.format(slices-1));
			tp.append("gray value:	xy gauss fit height");
			addFooter(tp);
		  	tp.saveAs(path + "_GaussH_info.txt");
		  	
		  	//save also as Table	TODO
		  	for(int s = 0; s < slices; s++){
		  		tp = new TextPanel("Kymograph for the parameter ");
				tp.append("Kymograph");
				tp.append("");
				
				String appText = "frame";
				for(int a = 0; a < width; a++){
					appText += "	" + constants.df6US.format((double)a * calibration);
				}
				tp.append(appText);
				
				for(int i = 0; i <= tMax; i++){
					appText = constants.df0.format(i);
					for(int a = 0; a < width; a++){
						appText += "	";
						if(!Double.isNaN(table [s][a][i])){
							appText += constants.df6US.format(table [s][a][i]);
						}														
					}
					tp.append(appText);
				}
				addFooter(tp);
			  	tp.saveAs(path + "_GaussH.txt");
		  	} 
		  	
	}
	
	public static void saveOrientedTraceImage(ImagePlus imp, ArrayList<trace> traces, int encoding, String path, double calibration){
		String ending = "";
		
		//get min / max values
		int	tMax = 0;
		for(int i = 0; i < traces.size(); i++){			
			if(!(traces.get(i).getTracePoints().size() > 0)) continue;
			if(!traces.get(i).oriented)	continue;			
			if(traces.get(i).getFrame() > tMax)	tMax = traces.get(i).getFrame();
		}
		
		//initialize image and array
		int width = (int)Math.round((KYMAX_X - KYMIN_X)/calibration) + 1,
			height = (int)Math.round((KYMAX_Y - KYMIN_Y)/calibration) + 1;
		
		ImagePlus traceImp = IJ.createImage("Ori Image", width, height, tMax+1, 8);
		traceImp.setCalibration(imp.getCalibration());
		
		//save positions		
		ending = "_ori";
		if(encoding == MEDIANZ){
			ending += "ZCmedian";
		}else if(encoding == MEANZ){
			ending += "ZCmean";
		}else if(encoding == PUREZ){
			ending += "ZC";
		}
		
		try{
			double [][][] saveImage = new double [width][height][2];
			int cx, cy;
			for(int i = 0; i < traces.size(); i++){
				if(!(traces.get(i).getTracePoints().size() > 0)) continue;
				if(!traces.get(i).oriented)	continue;
				
				for(int x = 0; x < width; x++){
					for(int y = 0; y < height; y++){
							saveImage [x][y][0] = 0.0;
							saveImage [x][y][1] = 0.0;
					}
				}
				
				for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
					cx = (int)Math.round((traces.get(i).getTracePoints().get(j).get2DOrientedVector()[0]-KYMIN_X)/calibration);
					cy = (int)Math.round((traces.get(i).getTracePoints().get(j).get2DOrientedVector()[1]-KYMIN_Y)/calibration);
					if(cx >= 0 && cx < saveImage.length && cy >= 0 && cy < saveImage[0].length){
						saveImage [cx][cy][0] += tools.getEncodedIntensity8bit(traces.get(i).getTracePoints().get(j).getZ(encoding),KYMIN_Z,KYMAX_Z);
						saveImage [cx][cy][1] += 1.0;
					}					
				}
				
				//write to image			
				for(int x = 0; x < width; x++){
					for(int y = 0; y < height; y++){
						if(saveImage[x][y][1] > 0.0){
							traceImp.getStack().setVoxel(x, y, traces.get(i).getFrame(), saveImage[x][y][0] / saveImage[x][y][1]); 
						}
					}
				
				}
			}
		}catch(Exception e){
//			IJ.log("error w" + width + " h" + height);
		}
						
		IJ.run(traceImp, "Fire", "");		
		IJ.saveAsTiff(traceImp, path + ending + ".tif");
		traceImp.changes = false;
		traceImp.close();
		
		//save image metadata
		TextPanel tp = new TextPanel("Metadata");
		tp.append("Image information for image:	" + path.substring(path.lastIndexOf(System.getProperty("file.separator"))) + ending + ".tif");
		tp.append("");
		tp.append("dimension		minimum	maximum");
		tp.append("gray value (1.0-254.0):	z information [um]	" + constants.df6US.format(KYMIN_Z) + "	" + constants.df6US.format(KYMAX_Z));
		tp.append("");
		tp.append("code#"+KYMIN_Z+"-"+KYMAX_Z+"");
		addFooter(tp);
	  	tp.saveAs(path + ending + "_info.txt");
	}
	
	public static void saveKymograph(ArrayList<trace> traces, double calibration, String path, int encoding, int kymoType, int excludeHeadPoints, boolean oriented3D){
		double min = getKymographMin(kymoType), max = getKymographMax(kymoType);		
		String xyzText = getKymographTxtLabel(kymoType);		
		int arcLengthEncoding = NOZ;
		if(oriented3D){
			arcLengthEncoding = encoding;
		}
		//TODO MAYBE INTEGRATE WHEN Z WORKS BETTER?
//		switch(kymoType){
//			case KYMOCANGLE2D:
//			case KYMOCURV2D:
//			case KYMODZ:
//			case KYMOX:
//			case KYMOY:
//			case KYMOZ:
//				arcLengthEncoding = NOZ;
//				break;
//			case KYMOCANGLE3D:
//			case KYMOCURV3D:
//				arcLengthEncoding = encoding;
//				break;
//		}
		
		//test function
		if(min == 0.0 && max == 0.0)	IJ.log("minmax0 at kymoType" + kymoType);
		
		//get arc min max and tmax
			double arcLengthMax = Double.NEGATIVE_INFINITY;
			int	tMax = 0;
			for(int i = 0; i < traces.size(); i++){
				if(!traces.get(i).oriented)	continue;
				if(traces.get(i).getTracePoints().size() <= excludeHeadPoints)	continue;	
				if(traces.get(i).getFrame()>tMax)	tMax = traces.get(i).getFrame();
				
				for(int j = excludeHeadPoints; j < traces.get(i).getTracePoints().size(); j++){
					if(traces.get(i).getTracePoints().get(j).getArcLengthPos(arcLengthEncoding) > arcLengthMax){
						arcLengthMax = traces.get(i).getTracePoints().get(j).getArcLengthPos(arcLengthEncoding);
					}				
				}
			}			
			
		//generate image
			ImagePlus codeImp = IJ.createImage("Kymograph", (int)Math.round(arcLengthMax / calibration)+1, tMax+1, 1, 16);
			double values [][] = new double [(int)Math.round(arcLengthMax / calibration)+1][2];			
			
			for(int i = 0; i < traces.size(); i++){
				if(!traces.get(i).oriented)	continue;
				if(traces.get(i).getTracePoints().size() <= excludeHeadPoints)	continue;	
				// && Double.isNaN(traces.get(i).getThetaDegree(encoding))==false
				
				for(int a = 0; a < (int)Math.round(arcLengthMax / calibration)+1; a++){
					values [a][0] = 0.0;
					values [a][1] = 0.0;
				}
				
				for(int j = excludeHeadPoints; j < traces.get(i).getTracePoints().size(); j++){
					try{
						values[(int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(arcLengthEncoding) / calibration)][0]
								+= tools.getEncodedIntensity16bit (getKymoTypeParameter(traces.get(i).getTracePoints().get(j), kymoType, encoding, oriented3D), min, max);
						values[(int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(arcLengthEncoding) / calibration)][1] += 1.0;											
					}catch (Exception e){
						IJ.log("T" + i + ": problem in kymograph generation..." + " j" + j + " - " + (int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(arcLengthEncoding) / calibration) 
						+ " (allowed length " + values.length + ")");
					}
				}	
				
				for(int a = 0; a < (int)Math.round(arcLengthMax / calibration)+1; a++){
//							IJ.log("j " + j + " i " + i + "xyz" + xyz + " enc " + encoding);
//							IJ.log("al " + points.get(j).getArcLengthPos(encoding));
//							IJ.log("al " + (int)Math.round(points.get(j).getArcLengthPos(encoding) * resolution));
//							IJ.log("corrI " + getCorrectedIntensity8bit (points.get(j).getOrientedVector(encoding)[xyz], min, max));					
					if(values[a][1]>0.0){
						codeImp.getStack().setVoxel(a, traces.get(i).getFrame(), 0, values[a][0]/values[a][1]);
					}
				}								
			}
			
			IJ.run(codeImp, "Fire", "");
			
			codeImp.getCalibration().pixelWidth = calibration;
			impProcessing.setOptimalDisplayRange(codeImp, false);
			
			String ending = "_k" + xyzText;
			if((kymoType == KYMOZ || kymoType == KYMODZ)  && encoding==MEDIANZ){
				 ending += "medi";
			}else if((kymoType == KYMOZ || kymoType == KYMODZ) && encoding==MEANZ){
				ending += "mean";
			}
			
			IJ.saveAsTiff(codeImp, path + ending + ".tif");
			codeImp.changes = false;
			codeImp.close();
			
		//save image metadata
			TextPanel tp = new TextPanel("Metadata");
			tp.append("Image information for image:	" + path.substring(path.lastIndexOf(System.getProperty("file.separator"))) + ending + ".tif");
			tp.append("");
			tp.append("dimension		minimum	maximum");
			if(kymoType == KYMOCANGLE3D || kymoType == KYMOCURV3D){
				tp.append("x axis:	3D arc length	0	" + constants.df0.format((int)Math.round(arcLengthMax/calibration)+1) + "	calibration [um/px]: " + constants.df6US.format(calibration) + "");
			}else{
				tp.append("x axis:	2D arc length	0	" + constants.df0.format((int)Math.round(arcLengthMax/calibration)+1) + "	calibration [um/px]: " + constants.df6US.format(calibration) + "");
			}			
			tp.append("y axis:	time [frame]	0	" + constants.df0.format(tMax));
			tp.append("gray value (1.0-65534.0):	oriented " + xyzText + " position	" + constants.df6US.format(min) + "	" + constants.df6US.format(max));
			tp.append("");
			tp.append("code#"+min+"-"+max+"+"+(1.0/calibration));
			addFooter(tp);
		  	tp.saveAs(path + ending + "_info.txt");
	}
	
	/**
	 * @return kymograph minimum displayed value
	 * */
	private static double getKymographMin(int kymoType){
		switch(kymoType){
		case KYMOX:
			return KYMIN_X;
		case KYMOY:
			return KYMIN_Y;
		case KYMOZ:
			return KYMIN_Z;
		case KYMOMAXINTENSITY:
			return KYMIN_MAXINTENSITY;
		case KYMOCURV2D:
			return KYMIN_CURV;
		case KYMOCURV3D:
			return KYMIN_CURV;
		case KYMOCANGLE2D:
			return KYMIN_CANGLE2D;
		case KYMOCANGLE3D:
			return KYMIN_CANGLE3D;
		case KYMODZ:
			return KYMIN_DZ;
		case KYMOTANGENTANGLE:
			return KYMIN_CANGLE2D;
		default:
			return -1000.0;
		}
	}
	
	/**
	 * @return kymograph maximum displayed value
	 * */
	private static double getKymographMax(int kymoType){
		switch(kymoType){
		case KYMOX:
			return KYMAX_X;
		case KYMOY:
			return KYMAX_Y;
		case KYMOZ:
			return KYMAX_Z;
		case KYMOMAXINTENSITY:
			return KYMAX_MAXINTENSITY;
		case KYMOCURV2D:
			return KYMAX_CURV;
		case KYMOCURV3D:
			return KYMAX_CURV;
		case KYMOCANGLE2D:
			return KYMAX_CANGLE2D;
		case KYMOCANGLE3D:
			return KYMAX_CANGLE3D;
		case KYMODZ:
			return KYMAX_DZ;
		case KYMOTANGENTANGLE:
			return KYMAX_CANGLE2D;
		default:
			return -1000.0;
		}
	}
	
	/**
	 * @return label of displayed values
	 * */
	private static String getKymographTxtLabel(int kymoType){
		switch(kymoType){
		case KYMOX:
			return "X";
		case KYMOY:
			return "Y";
		case KYMOZ:
			return "Z";
		case KYMOMAXINTENSITY:
			return "MaxI";
		case KYMOCURV2D:
			return "Curv2D";
		case KYMOCURV3D:
			return "Curv3D";
		case KYMOCANGLE2D:
			return "cA2D";
		case KYMOCANGLE3D:
			return "cA3D";
		case KYMODZ:
			return "dZ";
		case KYMOTANGENTANGLE:
			return "tAng";
		default:
			return "";
		}
	}
	
	public static double getKymoTypeParameter (trackPoint p, int type, int encoding, boolean oriented3D){
		switch(type){
		case KYMOX: 
			if(oriented3D){
				return p.get3DOrientedVector(encoding)[0]; 
			}else{
				return p.get2DOrientedVector()[0];
			}
		case KYMOY: 
			if(oriented3D){
				return p.get3DOrientedVector(encoding)[1]; 
			}else{
				return p.get2DOrientedVector()[1];
			} 
		case KYMOZ: 
			if(oriented3D){
				return p.get3DOrientedVector(encoding)[2]; 
			}else{
				return p.getZ(encoding);
			}
		case KYMOMAXINTENSITY: return p.getGaussFitMax();
		case KYMOCURV2D: return p.getCurvature2D();
		case KYMOCURV3D: return p.getCurvature3D();
		case KYMOCANGLE2D: return p.getCAngle2D();
		case KYMOCANGLE3D: return p.getCAngle3D();
//		case KYMODZ: return p.getDZ(encoding);
		default: return 0.0;
		}
	}
	
	public static void saveXYZCoordinates(ArrayList<trace> traces, double calibration, String path, int encoding, ProgressDialog progress){
		//get arc min max and tmax
//			double arcLengthMin = Double.POSITIVE_INFINITY;
			double arcLengthMax = Double.NEGATIVE_INFINITY;
			int	tMax = 0;
			progress.notifyMessage("Create XYZ Coordinate lists: nr of timepoints to save: " + traces.size(), ProgressDialog.LOG);
			for(int i = 0; i < traces.size(); i++){
//				if(!traces.get(i).oriented)	continue;
				if(traces.get(i).getTracePoints().size() <= 0){
					progress.notifyMessage("Create XYZ Coordinate lists: trace " + i + " harbors 0 points", ProgressDialog.LOG);
					continue;	
				}
				if(traces.get(i).getFrame() > tMax){
					tMax = traces.get(i).getFrame();
				}
				
				for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
					if(traces.get(i).getTracePoints().get(j).getArcLengthPos(encoding) > arcLengthMax){
						arcLengthMax = traces.get(i).getTracePoints().get(j).getArcLengthPos(encoding);
					}
//					if(traces.get(i).getTracePoints().get(j).getArcLengthPos() < arcLengthMin){
//						arcLengthMin = traces.get(i).getTracePoints().get(j).getArcLengthPos();
//					}					
				}
			}			
			
			progress.notifyMessage("Create XYZ Coordinate lists: arc length maximum: " + constants.dfdialog.format(arcLengthMax), ProgressDialog.LOG);
			
		//save all
			double values [][] = new double [(int)Math.round(arcLengthMax / calibration)+1][2];
			TextPanel tp = new TextPanel("Results");
			String appText = "";
			//save X
			{
				tp.clear();
				tp.append("time	arc length (micron)");
			  	appText = "";
			  	for(int a = 0; a < (int)Math.round(arcLengthMax / calibration)+1; a++){
			  		appText += "	" + constants.df3US.format(a*calibration);
			  	}
			  	tp.append(appText);
			  	
			  	for(int i = 0; i < traces.size(); i++){
			  		appText = "" + i;
					if(traces.get(i).getTracePoints().size() <= 0){		//!traces.get(i).oriented || 
						tp.append(appText);
						continue;	
					}
					// && Double.isNaN(traces.get(i).getThetaDegree(encoding))==false
					
					for(int a = 0; a < (int)Math.round(arcLengthMax / calibration)+1; a++){
						values [a][0] = 0.0;
						values [a][1] = 0.0;
					}
					
					for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
						try{
							values[(int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(encoding) / calibration)][0]
									+= traces.get(i).getTracePoints().get(j).getX();
							values[(int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(encoding) / calibration)][1] += 1.0;											
						}catch (Exception e){
							IJ.log("T" + i + ": problem in kymograph generation..." + " j" + j + " - " 
									+ (int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(encoding) / calibration) + " length " + values.length);
						}
					}	
					
					for(int a = 0; a < (int)Math.round(arcLengthMax / calibration)+1; a++){
//								IJ.log("j " + j + " i " + i + "xyz" + xyz + " enc " + encoding);
//								IJ.log("al " + points.get(j).getArcLengthPos());
//								IJ.log("al " + (int)Math.round(points.get(j).getArcLengthPos() * resolution));
//								IJ.log("corrI " + getCorrectedIntensity8bit (points.get(j).getOrientedVector(encoding)[xyz], min, max));					
						appText += "	";
						if(values[a][1]>0.0){
							appText += constants.df6US.format(values[a][0]/values[a][1]);
						}						
					}
					tp.append(appText);
				}
			  	multi_focal_tools.addFooter(tp);		
			  	tp.saveAs(path + "_coordX.txt");
			}
			
			//save Y
			{
				tp.clear();
				tp.append("time	arc length (micron)");
			  	appText = "";
			  	for(int a = 0; a < (int)Math.round(arcLengthMax / calibration)+1; a++){
			  		appText += "	" + constants.df3US.format(a*calibration);
			  	}
			  	tp.append(appText);
			  				  	
			  	for(int i = 0; i < traces.size(); i++){
			  		appText = "" + i;
			  		if(traces.get(i).getTracePoints().size() <= 0){		//!traces.get(i).oriented || 
						tp.append(appText);
						continue;	
					}
					// && Double.isNaN(traces.get(i).getThetaDegree(encoding))==false
					
					for(int a = 0; a < (int)Math.round(arcLengthMax / calibration)+1; a++){
						values [a][0] = 0.0;
						values [a][1] = 0.0;
					}
					
					for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
						try{
							values[(int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(encoding) / calibration)][0]
									+= traces.get(i).getTracePoints().get(j).getY();
							values[(int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(encoding) / calibration)][1] += 1.0;											
						}catch (Exception e){
							IJ.log("T" + i + ": problem in kymograph generation..." + " j" + j + " - " 
									+ (int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(encoding) / calibration) + " length " + values.length);
						}
					}	
					
					for(int a = 0; a < (int)Math.round(arcLengthMax / calibration)+1; a++){					
						appText += "	";
						if(values[a][1]>0.0){
							appText += constants.df6US.format(values[a][0]/values[a][1]);
						}
					}
					tp.append(appText);							
				}
			  	multi_focal_tools.addFooter(tp);		
			  	tp.saveAs(path + "_coordY.txt");
			}
			
			//save Z
			{
				tp.clear();
				tp.append("time	arc length (micron)");
			  	appText = "";
			  	for(int a = 0; a < (int)Math.round(arcLengthMax / calibration)+1; a++){
			  		appText += "	" + constants.df3US.format(a*calibration);
			  	}
			  	tp.append(appText);
			  				  	
			  	for(int i = 0; i < traces.size(); i++){
			  		appText = "" + i;
			  		if(traces.get(i).getTracePoints().size() <= 0){		//!traces.get(i).oriented || 
						tp.append(appText);
						continue;	
					}
					
					for(int a = 0; a < (int)Math.round(arcLengthMax / calibration)+1; a++){
						values [a][0] = 0.0;
						values [a][1] = 0.0;
					}
					
					for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
						try{
							values[(int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(encoding) / calibration)][0]
									+= traces.get(i).getTracePoints().get(j).getZ(encoding);
							values[(int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(encoding) / calibration)][1] += 1.0;											
						}catch (Exception e){
							IJ.log("T" + i + ": problem in kymograph generation..." + " j" + j + " - " 
									+ (int)Math.round(traces.get(i).getTracePoints().get(j).getArcLengthPos(encoding) / calibration) + " length " + values.length);
						}
					}	
					
					for(int a = 0; a < (int)Math.round(arcLengthMax / calibration)+1; a++){					
						appText += "	";
						if(values[a][1]>0.0){
							appText += constants.df6US.format(values[a][0]/values[a][1]);
						}
					}
					tp.append(appText);							
				}
			  	multi_focal_tools.addFooter(tp);		
			  	tp.saveAs(path + "_coordZ.txt");
			}
					  	
	}
	
	public static void saveHeadRotationMatrixImage(ArrayList<trace> traces, double [] slicePositions, double calibration, String path){
		double [] positionsCopy = new double [slicePositions.length]; 
		for(int i = 0; i < slicePositions.length; i++){
			positionsCopy [i] = slicePositions [i];
		}
		Arrays.sort(positionsCopy);
		
		//find maximum
		double min = Double.MAX_VALUE;
		double max = 0.0;
		int tMax = 0;
		double [][] matrix = new double [41][4];
		for(int i = 0; i < traces.size(); i++){
			if(!traces.get(i).oriented)	continue;
			if(!traces.get(i).hrVset)	continue;
			
			if(tMax < traces.get(i).getFrame())	tMax = traces.get(i).getFrame();
			
			matrix = traces.get(i).getHeadRotationMatrix();
			for(int x = 0; x < matrix.length; x++){
				for(int y = 0; y < matrix [0].length; y++){
					if(matrix [x][y] > max){
						max = matrix [x][y];
					}
					if(matrix [x][y] < min){
						min = matrix [x][y];
					}
				}
			}			
		}
		
		ImagePlus codeImp = IJ.createImage("HRimage", matrix [0].length, matrix.length, tMax+1, 8);
		int s;
		for(int i = 0; i < traces.size(); i++){
			if(!traces.get(i).oriented)	continue;
			if(!traces.get(i).hrVset)	continue;
			
			matrix = traces.get(i).getHeadRotationMatrix();
			for(int x = 0; x < matrix.length; x++){
				for(int y = 0; y < matrix [0].length; y++){
					s = tools.getIndexOfClosestValue(positionsCopy, slicePositions[x]);
					codeImp.getStack().setVoxel(y, s, traces.get(i).getFrame(), tools.getEncodedIntensity8bit(matrix[x][y], min, max));
				}
			}			
		}
		
		IJ.run(codeImp, "Fire", "");
		IJ.saveAsTiff(codeImp, path + "HRI.tif");
		codeImp.changes = false;
		codeImp.close();
		
		//save image metadata
		TextPanel tp = new TextPanel("Metadata");
		tp.append("Image information for image:	" + path.substring(path.lastIndexOf(System.getProperty("file.separator"))) + "HRI.tif");
		tp.append("");
		tp.append("dimension		minimum	maximum");
		tp.append("x axis:	points normal to orientation vector	0	9	calibration [um/px]: " + constants.df6US.format(calibration) + "");
		tp.append("y axis:	slice	0	4");
		tp.append("gray value (1.0-254.0):	interpolated intensity	" + constants.df6US.format(min) + "	" + constants.df6US.format(max));
		tp.append("");
		tp.append("code#"+min+"-"+max+"+"+(1.0/calibration));
		addFooter(tp);
	  	tp.saveAs(path + "HRI_info.txt");
		
	}
	
	public static double [][] getAndSaveFrequencies(ArrayList<trace> traces, String path, double calibration, int groupedTimesteps, int kymoType, int encoding,
			double sampleRate, double neglectedInitialArclength, boolean oriented3D, double speciesLength){
		int arcLengthEncoding = encoding;
		
		//get maximum arcLength and time-step
			double length = 0.0;
			int tMax = 0;
			for(int i = 0; i < traces.size(); i++){				
				if(traces.get(i).getFrame() > tMax)	tMax = traces.get(i).getFrame();
				if(traces.get(i).getTracePoints().size() > 0 && traces.get(i).oriented){
					for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
						if(traces.get(i).getTracePoints().get(j).getArcLengthPos(arcLengthEncoding) > length){
							length = traces.get(i).getTracePoints().get(j).getArcLengthPos(arcLengthEncoding);
						}
					}					
				}
			}
						
		//create Image
			int width = (int)Math.round(length / calibration)+1;
			int neglectedWidth = (int)Math.round(neglectedInitialArclength / calibration);
			if(neglectedInitialArclength <= 0.0)	neglectedWidth = -1;
			
		//initialize
			double [] medianVs = new double [] {0.0,0.0,0.0};
			int counter;
			
			double values [][][] = new double [width][2][groupedTimesteps];
			double valueColumn [] = new double [groupedTimesteps];
			double frequencies [][][] = new double [width][tMax + 1 - groupedTimesteps + 1][5];
			int calArcLPos = 0;
			
			for(int orT = 0; orT < traces.size() && traces.get(orT).getFrame() <= tMax + 1 - groupedTimesteps; orT++){
				for(int t = 0; t < groupedTimesteps; t++){
					for(int a = neglectedWidth + 1; a < width; a++){
						values [a][0][t] = 0.0;
						values [a][1][t] = 0.0;
					}					
				}	
				
			
				//align
				for(int i = orT; i < traces.size() && traces.get(i).getFrame() < traces.get(orT).getFrame() + groupedTimesteps; i++){
					for(int j = 0; j < traces.get(i).getTracePoints().size(); j++){
						calArcLPos = (int)(traces.get(i).getTracePoints().get(j).getArcLengthPos(arcLengthEncoding)/calibration);
						values [calArcLPos][0][traces.get(i).getFrame()-traces.get(orT).getFrame()] 
								+= getKymoTypeParameter(traces.get(i).getTracePoints().get(j), kymoType, encoding, oriented3D);
						values [calArcLPos][1][traces.get(i).getFrame()-traces.get(orT).getFrame()] += 1.0;
						if(calArcLPos+1 < width){
							values [calArcLPos+1][0][traces.get(i).getFrame()-traces.get(orT).getFrame()] 
									+= getKymoTypeParameter(traces.get(i).getTracePoints().get(j), kymoType, encoding, oriented3D);
							values [calArcLPos+1][1][traces.get(i).getFrame()-traces.get(orT).getFrame()] += 1.0;
						}
						if(calArcLPos-1 >= 0){
							values [calArcLPos-1][0][traces.get(i).getFrame()-traces.get(orT).getFrame()] 
									+= getKymoTypeParameter(traces.get(i).getTracePoints().get(j), kymoType, encoding, oriented3D);
							values [calArcLPos-1][1][traces.get(i).getFrame()-traces.get(orT).getFrame()] += 1.0;
						}
						
					}
				}
				
				//calculate values and determine frequency
				for(int t = 0; t < groupedTimesteps; t++){
					for(int a = neglectedWidth + 1; a < width; a++){
						if(values [a][1][t] != 0.0){
							values [a][0][t] /= values [a][1][t];
						}
					}
				}
								
				//get frequencies
				for(int a = neglectedWidth + 1; a < width; a++){
//					if(a == (int)(width*3.0/4.0)){
//						tools2D.showAsPlot(values [a][0]);
//						IJ.log("arclpos " + a + " = " + (a*calibration));
//					}
					{
						//get smoothed graph
						for(int t = 0; t < groupedTimesteps; t++){
							if(values [a][1][t] != 0.0){
								counter = 1;
								medianVs [0] = values [a][0][t];
								if(t + 1 < groupedTimesteps
										&& values [a][1][t+1] != 0.0){
									medianVs [counter] = values [a][0][t+1];
									counter++;
								}
								if(t != 0
										&& values [a][1][t-1] != 0.0){
									medianVs [counter] = values [a][0][t-1];
									counter++;
								}
								if(counter==1){
									valueColumn [t] = medianVs [0];
								}else if(counter == 2){
									valueColumn [t] = (medianVs [0] + medianVs [1]) / 2.0;
								}else{
									valueColumn [t] = tools.getMedian(medianVs);
								}
							}else{
								//fill holes in graph
								valueColumn [t] = 0.0;
								fillingHoles: for(int l = 1; l < 10; l++){
									if(t + l < values[0].length){
										valueColumn [t] = values [a][0][t+l];
									}
									if(t - l >= 0){
										if(valueColumn [t] != 0.0){
											valueColumn [t] += values [a][0][t-l];
											valueColumn [t] /= 2.0;
										}else{
											valueColumn [t] = values [a][0][t-l];
										}
									}
									if(valueColumn [t] != 0.0){
										break fillingHoles;
									}
								}
							}								
						}
					}					
					if(kymoType == multi_focal_tools.KYMOCURV2D || kymoType == multi_focal_tools.KYMOCURV3D || kymoType == multi_focal_tools.KYMOCANGLE3D 
							|| kymoType == multi_focal_tools.KYMOMAXINTENSITY || kymoType == multi_focal_tools.KYMOX){
						frequencies [a][traces.get(orT).getFrame()] = getFrequenciesAndAmplitude(valueColumn, sampleRate, false, true);
					}else{
						frequencies [a][traces.get(orT).getFrame()] = getFrequenciesAndAmplitude(valueColumn, sampleRate, false, false);
					}					
					//display frequency
//					if(a == (int)(width*3.0/4.0)){
//						tools2D.showAsPlot(valueColumn);
//						frequencies [a][traces.get(orT).getFrame()] = tools2D.getFrequency(valueColumn, sampleRate, true);
//					}
				}
			}
			
			//Median smoothing of Freqs and Amplitudes along flagellum
			medianVs = new double [5];			
			for(int t = 0; t < tMax + 1 - groupedTimesteps + 1; t++){
				for(int a = neglectedWidth+1; a < width; a++){
					for(int f = 0; f < 5; f++){
						if(frequencies [a][t][f] > 0.0){
							counter = 0;
							for(int add = a - 2; add <= a + 2; add++){
								if(add < width 
										&& add > neglectedWidth 
										&& frequencies [add][t][f] > 0.0){
									medianVs [counter] = frequencies [add][t][f];
									counter++;
								}
							}
							if(counter==1){
								frequencies [a][t][f] = medianVs [0];
							}else if(counter>1){
								frequencies [a][t][f] = tools.getMedianOfRange(medianVs, 0, counter-1);
							}						
//						}else{
//							IJ.log(getKymographTxtLabel(kymoType) + ": " + a + "-" + t + "-" + f + ": " + frequencies [a][t][f]);
						}
					}
					
					
				}
			}
			
			double [][] freqs = new double [3][3];
			int [][] freqsCt = new int [3][3];
		  	for(int f1 = 0; f1 < 3; f1++){
		  		freqs [0][f1] = 0.0;
		  		freqs [1][f1] = 0.0;
		  		freqs [2][f1] = 0.0;
		  		freqsCt [0][f1] = 0;
		  		freqsCt [1][f1] = 0;
		  		freqsCt [2][f1] = 0;
		  	}
		  	
			for(int a = neglectedWidth+1; a < speciesLength; a++){
				for(int t = 0; t < tMax + 1 - groupedTimesteps + 1; t++){
					if(frequencies [a][t][0] > 0.0){						
						freqs [0][(int)((a/(double)speciesLength)*3.0)] += frequencies [a][t][0];
						freqsCt [0][(int)((a/(double)speciesLength)*3.0)] ++;
						freqs [1][(int)((a/(double)speciesLength)*3.0)] += frequencies [a][t][2];
						freqsCt [1][(int)((a/(double)speciesLength)*3.0)] ++;
						freqs [2][(int)((a/(double)speciesLength)*3.0)] += frequencies [a][t][4];
						freqsCt [2][(int)((a/(double)speciesLength)*3.0)] ++;
					}					
				}
			}

			ImagePlus codeImpF = IJ.createImage("Frequency graph", width, tMax + 1 - groupedTimesteps + 1, 5, 16);

			for(int a = neglectedWidth+1; a < width; a++){
				for(int t = 0; t < tMax + 1 - groupedTimesteps + 1; t++){
					if(frequencies [a][t][0] > 0.0){
						codeImpF.getStack().setVoxel(a, t, 0, tools.getEncodedIntensity16bit(frequencies [a][t][0], KYMIN_FREQ, KYMAX_FREQ));
					}					
					if(frequencies [a][t][1] > 0.0){
						codeImpF.getStack().setVoxel(a, t, 1, tools.getEncodedIntensity16bit(frequencies [a][t][1], KYMIN_AMPLITUDE, KYMAX_AMPLITUDE));
					}
					if(frequencies [a][t][2] > 0.0){
						codeImpF.getStack().setVoxel(a, t, 2, tools.getEncodedIntensity16bit(frequencies [a][t][2], KYMIN_FREQ, KYMAX_FREQ));
					}
					if(frequencies [a][t][3] > 0.0){
						codeImpF.getStack().setVoxel(a, t, 3, tools.getEncodedIntensity16bit(frequencies [a][t][3], KYMIN_AMPLITUDE, KYMAX_AMPLITUDE));
					}
					if(frequencies [a][t][4] > 0.0){
						codeImpF.getStack().setVoxel(a, t, 4, tools.getEncodedIntensity16bit(frequencies [a][t][4], KYMIN_FREQ, KYMAX_FREQ));
					}
				}
			}
			codeImpF.getCalibration().pixelWidth = calibration;
			IJ.run(codeImpF, "Fire", "");
			
			IJ.saveAsTiff(codeImpF, path  + "_" + getKymographTxtLabel(kymoType) + "_f.tif");
			codeImpF.changes = false;
			codeImpF.close();
			
		//save image metadata
			TextPanel tp = new TextPanel("Metadata");
			tp.append("Image information for image:	" + path.substring(path.lastIndexOf(System.getProperty("file.separator"))) + "_f.tif");
			tp.append("");
			tp.append("dimension		minimum	maximum");
			tp.append("x axis:	xy arc length	0	" + constants.df0.format((int)Math.round(length/calibration)+1) + "	calibration: " + constants.df6US.format(calibration) + "");
			tp.append("y axis:	time [frame]	0 - " + constants.df0.format(groupedTimesteps) 
				+ "	" + constants.df0.format(tMax-1 - groupedTimesteps) + " - " + constants.df0.format(tMax-1));
			tp.append("z axis:	0,2 = 1st and 2nd largest frequency in FFT	1,3 = FFT amplitudes of 1st and 2nd largest frequency");
			tp.append("gray value of z=0,2 (1.0-65534.0):	frequency	" + constants.df6US.format(KYMIN_FREQ) + "	" + constants.df6US.format(KYMAX_FREQ));
			tp.append("gray value of z=1,3 (1.0-65534.0):	amplitude	" + constants.df6US.format(KYMIN_AMPLITUDE) + "	" + constants.df6US.format(KYMAX_AMPLITUDE));
			tp.append("");
			tp.append("code#"+KYMIN_FREQ+"&"+KYMAX_FREQ+"*"+KYMIN_AMPLITUDE+"!"+KYMAX_AMPLITUDE+"+"+calibration);
			addFooter(tp);
		  	tp.saveAs(path  + "_" + getKymographTxtLabel(kymoType) + "_f_info.txt");
		  	
		//calculate frequencies by group
		  	for(int f1 = 0; f1 < 3; f1++){
		  		freqs [0][f1] /= (double) freqsCt [0][f1];
		  		freqs [1][f1] /= (double) freqsCt [1][f1];
		  		freqs [2][f1] /= (double) freqsCt [2][f1];
		  	}
		  	return freqs;
		  	
	}
	
	/**
	 * @deprecated
	 * */
	public static double getFrequency (double [] values, double sampleRate, boolean showPlot){
		DoubleFFT_1D fftDo = new DoubleFFT_1D(values.length);	//Package unter: http://incanter.org/docs/parallelcolt/api/edu/emory/mathcs/jtransforms/fft/DoubleFFT_1D.html
        double [] fft = new double[values.length * 2];
        double [] magnitude = new double[values.length];
        System.arraycopy(values, 0, fft, 0, values.length); 
                       fftDo.realForwardFull(fft);
                       
        double real, imaginary;
        for(int j = 0; j < values.length; j++){
        	real = fft[2*j];
        	imaginary = fft[2*j+1];
        	magnitude [j] = Math.sqrt(real*real+imaginary*imaginary);
        }
        
        //display as plot
        if(showPlot)	tools.showAsPlot(magnitude);
        
        //output maximum frequency found (from index 2 on)    	
        return (tools.getMaximumIndexWithinRange(magnitude, 2, (magnitude.length/2)) * sampleRate / magnitude.length) ;
	}
	
	/**
	 * @performs a 1D frequency analysis of the array @param values
	 * @param showPlot: if true displays a plot of the found FFT
	 * @returns an array containing {highest freq. peak, amplitude of highest freq. peak,
	 * 								 2nd highest freq. peak, amplitude of 2nd highest freq. peak}
	 * */
	public static double [] getFrequenciesAndAmplitude (double [] values, double sampleRate, boolean showPlot, boolean normalizePlusMinus){
		//DoubleFFT package from: http://incanter.org/docs/parallelcolt/api/edu/emory/mathcs/jtransforms/fft/DoubleFFT_1D.html
				DoubleFFT_1D fftDo = new DoubleFFT_1D(values.length);	
		        double [] fft = new double[values.length * 2];
		        double [] magnitude = new double[values.length];
		        System.arraycopy(values, 0, fft, 0, values.length); 
		        
		        //normalization of values to +/- range
		      		if(normalizePlusMinus && tools.getMinimumWithinRange(fft, 0, fft.length-1) >= 0.0){
		      			double max = tools.getMaximum(fft);
		      			for(int i = 0; i < fft.length; i++){
		      				fft [i] -= max;
		      			}
		      		}
		        
		        fftDo.realForwardFull(fft);
		                       
		        double real, imaginary;
		        for(int j = 0; j < values.length; j++){
		        	real = fft[2*j];
		        	imaginary = fft[2*j+1];
		        	magnitude [j] = Math.sqrt(real*real+imaginary*imaginary);
		        }
		        
		        //display as plot
		        if(showPlot)	tools.showAsPlot(magnitude);
		        
		        //output maximum frequencies and amplitudes found (from index 2 on)
		        int [] freqs = tools.get2HighestMaximaIndicesWithinRange(magnitude, 2, (magnitude.length/2));
		        double [] ampl = new double [2];
		        double com = tools.getCenterOfMassOfRange(magnitude, 0, (magnitude.length/2)-1);
		        if(freqs[0] >= 0 && freqs [0] < magnitude.length){
		        	ampl [0] = magnitude [freqs[0]];
		        }else{
		        	ampl [0] = 0.0;
		        }
		        if(freqs[1] >= 0 && freqs [1] < magnitude.length){
		        	ampl [1] = magnitude [freqs[1]];
		        }else{
		        	ampl [1] = 0.0;
		        }
		        return new double [] {freqs [0] * sampleRate / magnitude.length, ampl [0],
			    		freqs [1] * sampleRate / magnitude.length, ampl [1], com * sampleRate / magnitude.length};
//		//DoubleFFT package from: http://incanter.org/docs/parallelcolt/api/edu/emory/mathcs/jtransforms/fft/DoubleFFT_1D.html
//		DoubleFFT_1D fftDo = new DoubleFFT_1D(values.length);	
//        double [] fft = new double[values.length * 2];
//        double [] magnitude = new double[values.length];
//        System.arraycopy(values, 0, fft, 0, values.length); 
//                       fftDo.realForwardFull(fft);
//                       
//        double real, imaginary;
//        for(int j = 0; j < values.length; j++){
//        	real = fft[2*j];
//        	imaginary = fft[2*j+1];
//        	magnitude [j] = Math.sqrt(real*real+imaginary*imaginary);
//        }
//        
//        //display as plot
//        if(showPlot)	tools.showAsPlot(magnitude);
//        
//        //output maximum frequencies and amplitudes found (from index 2 on)
//        int [] freqs = tools.get2HighestMaximaIndicesWithinRange(magnitude, 2, (magnitude.length/2));
//        double [] ampl = new double [2];
//        if(freqs[0] >= 0 && freqs [0] < magnitude.length){
//        	ampl [0] = magnitude [freqs[0]];
//        }else{
//        	ampl [0] = 0.0;
//        }
//        if(freqs[1] >= 0 && freqs [1] < magnitude.length){
//        	ampl [1] = magnitude [freqs[1]];
//        }else{
//        	ampl [1] = 0.0;
//        }
//        return new double [] {freqs [0] * sampleRate / magnitude.length, ampl [0],
//	    		freqs [1] * sampleRate / magnitude.length, ampl [1]};
	}
	
	public static double getTraceTypeParameter (trace t, int type, int encoding,
			double [] minIntensity, double [] maxIntensity, 
			double [] minPosition, double [] maxPosition){
		switch(type){
			case TRACE_THETA: 
				return t.getThetaDegree(encoding);
			case TRACE_HRMAXPOSITION: 
				return tools.getAverage(t.getHeadRotationMaximumPositions());
			case TRACE_HRMAXINTENSITY: 
				return tools.getAverage(t.getHeadRotationMaximumIntensities());
			case TRACE_HRANGLE: 
				return tools.getAverage(t.getHeadRotationAngle(minIntensity, maxIntensity, minPosition, maxPosition));
			default: return 0.0;
		}
	}
	
	public static String getTraceTypeParameterText (int type, int encoding){
		switch(type){
			case TRACE_THETA:
				switch(encoding){
					case NOZ: return "Th2D";
					case PUREZ: return "Th3D";
					case MEANZ: return "Th3DMean";
					case MEDIANZ: return "Th3DMedi";
				}
				break;
			case TRACE_HRMAXPOSITION: 
				return "HRMaxPos";
			case TRACE_HRMAXINTENSITY: 
				return "HRMaxInt";
			case TRACE_HRANGLE: 
				return "HRAng";
			default: return null;
		}
		return null;
	}
	
	public static double [] getAndSaveGroupedFrequenciesTraceParam (ArrayList<trace> traces, int traceParamType, int encoding, 
			int groupedTimesteps, double sampleRate,
			String savePath,
			double [] minIntensity, double [] maxIntensity, 
			double [] minPosition, double [] maxPosition){
		
		//get max time-step
			int tMax = 0;
			for(int i = 0; i < traces.size(); i++){				
				if(traces.get(i).getFrame() > tMax)	tMax = traces.get(i).getFrame();
			}
			
		//initialize
			double [] medianVs = new double [] {0.0,0.0,0.0};
			int counter;
			
			double values [] = new double [groupedTimesteps];
			double valueColumn [] = new double [groupedTimesteps];
			double frequencies [][] = new double [tMax + 1 - groupedTimesteps + 1][5];
			for(int i = 0; i < frequencies.length; i++){
				Arrays.fill(frequencies [i], 0.0);
			}
			
			for(int orT = 0; orT < traces.size() && traces.get(orT).getFrame() <= tMax + 1 - groupedTimesteps; orT++){
				for(int t = 0; t < groupedTimesteps; t++){
					values [t] = Double.NEGATIVE_INFINITY;										
				}	
				
				//get values
				for(int i = orT; i < traces.size() && traces.get(i).getFrame() < traces.get(orT).getFrame() + groupedTimesteps; i++){
					values [traces.get(i).getFrame()-traces.get(orT).getFrame()] 
							= getTraceTypeParameter(traces.get(i), traceParamType, encoding, 
							minIntensity, maxIntensity, minPosition, maxPosition);
				}
								
				//get frequencies
					{
						//get smoothed graph
						for(int t = 0; t < groupedTimesteps; t++){
							if(values [t] != Double.NEGATIVE_INFINITY){
								counter = 1;
								medianVs [0] = values [t];
								if(t + 1 < groupedTimesteps 
										&& values [t+1] != Double.NEGATIVE_INFINITY){
									medianVs [counter] = values [t+1];
									counter++;
								}
								if(t != 0
										&& values [t-1] != Double.NEGATIVE_INFINITY){
									medianVs [counter] = values [t-1];
									counter++;
								}
								if(counter==1){
									valueColumn [t] = medianVs [0];
								}else if(counter == 2){
									valueColumn [t] = (medianVs [0] + medianVs [1]) / 2.0;
								}else{
									valueColumn [t] = tools.getMedian(medianVs);
								}
							}else{
								//fill holes in graph
								valueColumn [t] = 0.0;
								fillingHoles: for(int l = 1; l < 10; l++){
									if(t + l < values.length){
										valueColumn [t] = values [t+l];
									}
									if(t - l >= 0){
										if(valueColumn [t] != 0.0){
											valueColumn [t] += values [t-l];
											valueColumn [t] /= 2.0;
										}else{
											valueColumn [t] = values [t-l];
										}
									}
									if(valueColumn [t] != 0.0){
										break fillingHoles;
									}
								}
							}								
						}
					}						
					frequencies [traces.get(orT).getFrame()] = getFrequenciesAndAmplitude(valueColumn, sampleRate, false, true);
				}
			
			double [] freqs = new double [3];
			int [] freqsCt = new int [3];
	  		freqs [0] = 0.0;
	  		freqs [1] = 0.0;
	  		freqs [2] = 0.0;
	  		freqsCt [0] = 0;
	  		freqsCt [1] = 0;
	  		freqsCt [2] = 0;
		  	
			for(int t = 0; t < tMax + 1 - groupedTimesteps + 1; t++){
				if(frequencies [t][0] > 0.0){						
					freqs [0] += frequencies [t][0];
					freqsCt [0] ++;
					freqs [1] += frequencies [t][2];
					freqsCt [1] ++;
					freqs [2] += frequencies [t][4];
					freqsCt [2] ++;
				}					
			}
			
			TextPanel tp = new TextPanel("Freq Results");
			tp.append("Frequencies for trace parameter:	" + getTraceTypeParameterText(traceParamType, encoding));
			tp.append("Analyzed image:	" + savePath.substring(savePath.lastIndexOf(System.getProperty("file.separator"))) + "_f.tif");
			tp.append("");
			tp.append("t [frame]" + "	" + "Primary freq." + "	" + "Amplitude of primary freq. peak in FFT"
					+ "	" + "Secondary freq." + "	" + "Amplitude of secondary freq. peak in FFT" + "	"
					+ "	" + "COM freq.	");
			
			String txt;
			for(int t = 0; t < tMax + 1 - groupedTimesteps + 1; t++){
				txt = "" + constants.df0.format(t) + "-" + constants.df0.format(t+groupedTimesteps);
				txt += "	";
				if(frequencies [t][0] > 0.0){
					txt += "" + frequencies [t][0];
				}					
				txt += "	";
				if(frequencies [t][1] > 0.0){
					txt += "" + frequencies [t][1];
				}
				txt += "	";
				if(frequencies [t][2] > 0.0){
					txt += "" + frequencies [t][2];
				}
				txt += "	";
				if(frequencies [t][3] > 0.0){
					txt += "" + frequencies [t][3];
				}
				txt += "	";
				if(frequencies [t][4] > 0.0){
					txt += "" + frequencies [t][4];
				}
				tp.append(txt);
			}
			addFooter(tp);
		  	tp.saveAs(savePath + "_" + getTraceTypeParameterText(traceParamType, encoding) + "_f.txt");
		  	
		//calculate frequencies by group
	  		freqs [0] /= (double) freqsCt [0];
	  		freqs [1] /= (double) freqsCt [1];
	  		freqs [2] /= (double) freqsCt [2];
		  	return freqs;
	}
	
	static void addFooter (TextPanel tp){
		tp.append("");
		tp.append("Datafile was generated by '"+multi_focal.PLUGINNAME+"', (\u00a9 2013-" + constants.dateY.format(new Date()) + ": Jan N Hansen (jan.hansen@caesar.de) \u0026 Jan F Jikeli (jan.jikeli@caesar.de))");
		tp.append("The plug-in '"+multi_focal.PLUGINNAME+"' is distributed in the hope that it will be useful,"
				+ " but WITHOUT ANY WARRANTY; without even the implied warranty of"
				+ " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.");
		tp.append("Skeleton results are generated using the ImageJ plug-ins 'skeletonize' and 'analyze skeleton'"
		+ " by: Ignacio Arganda-Carreras, Rodrigo Fernandez-Gonzalez, Arrate Munoz-Barrutia, Carlos Ortiz-De-Solorzano,"
		+ " '3D reconstruction of histological sections: Application to mammary gland tissue', Microscopy Research and Technique,"
		+ " Volume 73, Issue 11, pages 1019-1029, October 2010.");
		tp.append("Plug-in version:	"+multi_focal.PLUGINVERSION);	  	
	}
	
	/**
	 * @returns the distance between two trackPoints, @param p and @param q
	 * @param encoding specifies which kind of determined Z is used to calculate the distance
	 * */
	public static double getDistance (trackPoint p, trackPoint q, int encoding){
		if(p == null)	IJ.log("p=0");
		if(q == null)	IJ.log("q=0");
		return Math.sqrt(Math.pow(p.getX()-q.getX(),2.0)+Math.pow(p.getY()-q.getY(),2.0)+Math.pow(p.getZ(encoding)-q.getZ(encoding),2.0));
	}
	
	public static double [] get2DVectorFromPoints (double [] values){
		double slope = 0.0;
		int counter = 0;
		for(int i = 0; i < values.length; i++){
			for(int j = i+1; j < values.length; j++){
				slope += (values [j] - values [i])/(j-i);
				counter ++;
			}
		}
		return new double [] {1.0, (slope/counter)};
	}	
	
	public static double getCurvatureFactor(trackPoint p, trackPoint q, int encoding){
		double tf = tools.mathAbs(p.getArcLengthPos(encoding)-q.getArcLengthPos(encoding));
		if(tf == 0.0)	return 0.0;
		return tf / getDistance(p, q, encoding);
	}
}