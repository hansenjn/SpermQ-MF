package spermQ_mf;
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import ij.IJ;
import ij.ImagePlus;
import spermQ_mf.jnhsupport.*;

public class trace {
	private ArrayList<trackPoint> points;
	private int frame;	// 0 <= frame <= imp.getNFrames()-1
	private double [][] orientationVector;
	private double [][] origin;
	private double [] theta;
	private double hrCOMAngle;
	private double hrMaximaLineAngle;
	private double [] hrAngle;
	private double [] hrCOMVector;
	private double [] hrMaximumPosition;
	private double [] hrMaximumIntensity;
	private double [] hrMaximaLineVector;
	private double [] hrImageCOM;
	private double [][] hrMatrix;
	public boolean xyCorrected, g4pz, linearInterpolatedZ, hrVset, oriented;	
	
	
	public trace(ArrayList<trackPoint> pointList, int timeStep, double minDistance, double maxDistance){
		pointList.trimToSize();
		points = new ArrayList<trackPoint>(pointList);
		this.frame = timeStep;
		
		orientationVector = new double [4][3];
		origin = new double [4][3];
		theta = new double [4];
		hrCOMAngle = Double.NEGATIVE_INFINITY;
		hrMaximaLineAngle = Double.NEGATIVE_INFINITY;
		hrAngle = new double [4];
		hrCOMVector = new double [3];
		hrMaximumPosition = new double [4];
		hrMaximumIntensity = new double [4];
		hrImageCOM = new double [2];
		hrMatrix = new double [41][4];
		hrVset = false;
		oriented = false;
		xyCorrected = false;
		g4pz = false;
		linearInterpolatedZ = false;	
		
		this.updateAllOrientationVectors(minDistance, maxDistance);
	}
	
	public void reverseTracePoints(){
		Collections.reverse(points);		
	}
	
	//copy function
//	public trace (trace another) {
//		points = new ArrayList<trackPoint>(another.getTracePoints());
//		this.frame = another.frame;
//		
//		orientationVector = new double [4][3];
//		origin = new double [4][3];
//		theta = new double [4];
//		
//		g4pz = false;
//		linearInterpolated = false;	
//		
//		this.updateOrientationVector(15);
//	  }
	
	public void setFrame(int t){
		this.frame = t;
	}
	
	public int getFrame(){
		return frame;
	}
	
	public ArrayList<trackPoint> getTracePoints(){
		return points;
	}
	
	public void upscalePoints(int factor){
		/**
		 * Add points between all adjacent points to increase preciseness of the trace
		 * factor: multiplication factor -> between two adjacent points (factor-1) points will be added
		 * */
		
		points.trimToSize();
		ArrayList<trackPoint> pointList = new ArrayList<trackPoint>(points.size()*factor); 
		for(int i = 0; i < points.size()-1; i++){
			trackPoint p = new trackPoint(points.get(i).getX(), points.get(i).getY());
			pointList.add(p);
			for(int j = 1; j < factor; j++){
				double fraction = (j)/(double)(factor);
				pointList.add(new trackPoint(tools.getInterpolatedValue1D(fraction, 0, 1, p.getX(), points.get(i+1).getX()),
						tools.getInterpolatedValue1D(fraction, 0, 1, p.getY(), points.get(i+1).getY())));
			}			
		}
		pointList.add(new trackPoint(points.get(points.size()-1).getX(), points.get(points.size()-1).getY()));
		
		//replace points array with new array
		points.clear();
		points.ensureCapacity(pointList.size());
		for(int i = 0; i < pointList.size(); i++){
			points.add(pointList.get(i));
		}
		points.trimToSize();
	}
	
	public void sortPointsAndRemoveOutliers(ProgressDialog progress, int type){
		//create sortedList (list)
		ArrayList<trackPoint> list = new ArrayList<trackPoint>(points.size());
		{
			list.add(points.get(0));
			points.remove(0);
			points.trimToSize();
			trackPoint end = points.get(points.size()-1);
			
			int index = 0;
			double distance;
			trackPoint p;
			int pIndex;
			ordering: while(!points.isEmpty()){
				distance = Double.POSITIVE_INFINITY;
				p = null;
				pIndex = -1;
				for(int i = points.size()-1; i >= 0; i--){
					if(multi_focal_tools.getDistance(points.get(i),list.get(index),type) < distance){
						p = points.get(i);
						pIndex = i;
						distance = multi_focal_tools.getDistance(points.get(i),list.get(index),type);
					}
				}
				points.remove(pIndex);
				list.add(p);
				index++;
				if(p.equals(end)){
					if(points.size()>0){
						progress.notifyMessage("frame " + frame + ": points resorted, " + points.size() + " outliers removed...",  ProgressDialog.LOG);
					}					
//					IJ.log("Outliers removed in " + frame + ": "+ points.size());
					points.clear();					
					break ordering;
				}
			}
		}
		list.trimToSize();
		points = new ArrayList<trackPoint> (list);
		System.gc();
	}
	
	public void interpolateXYLinear(int plusMinusRange, int smoothingType){
		int nrOfPoints = points.size();
		
		//create point coordinate list
		double coords [][] = new double [nrOfPoints][2];
		for(int j = 0; j < nrOfPoints; j++){
			coords [j][0] = points.get(j).getX();
			coords [j][1] = points.get(j).getY();
		}
			
		int ctReqSize = (plusMinusRange)*(plusMinusRange)+1;
		
		double data [][] = new double [ctReqSize][3];	//x,y,distanceToPrevious
		double distances [] = new double [ctReqSize];
		int interPolCt;
		int index1, index2;
		for(int j = 1; j < nrOfPoints; j++){
			for(int m = 0; m < ctReqSize; m++){
				data [m][0] = Double.MAX_VALUE;
				data [m][1] = Double.MAX_VALUE;
				data [m][2] = Double.MAX_VALUE;
				distances [m] = Double.MAX_VALUE;
			}
			
			interPolCt = 0;
			data [interPolCt][0] = coords[j][0];
			data [interPolCt][1] = coords[j][1];
			data [interPolCt][2] = -1.0;
			distances [interPolCt] = data [interPolCt][2];
			interPolCt++;
			
			for(int d1 = -1; d1 >= -plusMinusRange; d1--){
				for(int d2 = 1; d2 <= plusMinusRange; d2++){
//					if(d1!=0&&d2!=0){
						if(j+d1>=0&&j+d1<nrOfPoints&&j+d2>=0&&j+d2<nrOfPoints){
							double [] newP = tools.projectPointToLine2D(coords[j][0],coords[j][1],
									coords[j+d1][0],coords[j+d1][1],
									coords[j+d2][0],coords[j+d2][1]);
							data [interPolCt][0] = newP [0];
							data [interPolCt][1] = newP [1];
							//determine shift distance
							data [interPolCt][2] = multi_focal_tools.getDistance(new trackPoint(coords[j][0],coords[j][1]), new trackPoint(data[interPolCt][0],data[interPolCt][1]), multi_focal_tools.NOZ);				
							distances [interPolCt] = data [interPolCt][2];
							interPolCt++;
						}
//					}
				}
			}			
			
			if(interPolCt != 1){
//				if((interPolCt+1)!=ctReqSize)	IJ.log("p" + j + ": " + (interPolCt+1) + "-" + ctReqSize);
				

				Arrays.sort(distances);
				if(smoothingType == multi_focal_tools.MINDIST){
					
					index1 = tools.getIndexOfClosestValue(tools.getArrayColumn(data, 2), distances[1]);
					points.get(j).resetXCoordinate(data[index1][0]);
					points.get(j).resetYCoordinate(data[index1][1]);
					
				}else if(smoothingType == multi_focal_tools.MAXDIST){
					index1 = tools.getIndexOfClosestValue(tools.getArrayColumn(data, 2), distances[interPolCt-1]);
					points.get(j).resetXCoordinate(data[index1][0]);
					points.get(j).resetYCoordinate(data[index1][1]);
					
				}else if(smoothingType == multi_focal_tools.MEAN){
					double x = 0.0, y = 0.0;
					for(int i = 0; i < interPolCt; i++){
						x += data [i][0];
						y += data [i][1];
					}
					
					points.get(j).resetXCoordinate(x/(double)(interPolCt));
					points.get(j).resetYCoordinate(y/(double)(interPolCt));
					
				}else if(smoothingType == multi_focal_tools.MEDIAN){	
//					double [] x = new double [interPolCt+1], y  = new double [interPolCt+1];
//					for(int i = 0; i <= interPolCt; i++){
////						if(data[i][0] == Double.MAX_VALUE) IJ.log("problem " + i);
////						if(data[i][1] == Double.MAX_VALUE) IJ.log("problem " + i);
//						x [i] = data [i][0];
//						y [i] = data [i][1];
//					}
//					Arrays.sort(x);
//					Arrays.sort(y);
//					
//					if((interPolCt+1)%2==0){
//						points.get(j).resetXCoordinate((x[(interPolCt)/2-1]+x[(interPolCt)/2])/2.0);
//						points.get(j).resetYCoordinate((y[(interPolCt)/2-1]+y[(interPolCt)/2])/2.0);
//					}else{
//						points.get(j).resetXCoordinate(x[(interPolCt)/2]);
//						points.get(j).resetYCoordinate(y[(interPolCt)/2]);
//					}
					
					if(interPolCt%2 == 0){
						index1 = tools.getIndexOfClosestValue(tools.getArrayColumn(data, 2), distances[(interPolCt)/2-1]);
						index2 = tools.getIndexOfClosestValue(tools.getArrayColumn(data, 2), distances[(interPolCt)/2]);
						
						points.get(j).resetXCoordinate((data[index1][0]+data[index2][0])/2.0);
						points.get(j).resetYCoordinate((data[index1][1]+data[index2][1])/2.0);
					}else{
						index1 = tools.getIndexOfClosestValue(tools.getArrayColumn(data, 2), distances[(interPolCt)/2]);
						points.get(j).resetXCoordinate(data[index1][0]);
						points.get(j).resetYCoordinate(data[index1][1]);
					}
					
				}else{
					IJ.error("wrong type of xy smoothing");
				}
			}
			
			
		}
	}
	
	public void updateAllOrientationVectors(double minDistance, double maxDistance){
		/**
		 * updates all origin and orientation vectors for which z-information is available
		 * maxDistance defines the maximum distance of points from the first point that are included into calculation
		 * */
		
		if(points.size()==0)	return;
				
		if(g4pz==false){
			this.updateOrientationVector(multi_focal_tools.NOZ, minDistance, maxDistance);
		}else{
			this.updateOrientationVector(multi_focal_tools.PUREZ, minDistance, maxDistance);
						
			if(linearInterpolatedZ){
				this.updateOrientationVector(multi_focal_tools.MEDIANZ, minDistance, maxDistance);
				this.updateOrientationVector(multi_focal_tools.MEANZ, minDistance, maxDistance);							
			}
		}				
	}
	
	/**
	 * updates origin and orientation vector of the trace according to the z-information selected in type
	 * Possible types are multi_focal_tools.NOZ, .G4PZ, .G4PZMEAN, and .G4PZMEDIAN
	 * maxDistance defines the maximum distance of points from the first point that are included into calculation
	 * */
	public void updateOrientationVector (int type, double minDistance, double maxDistance){
		origin [type][0] = points.get(0).getX();
		origin [type][1] = points.get(0).getY();
		origin [type][2] = points.get(0).getZ(type);
		
		//find startPoint
		double [] newVector = {0.0,0.0,0.0};
		int counter = 0, skippedPoints = 0;
		searchStart: for(int i = 1; i < points.size(); i++){
			newVector [0] = points.get(i).getX() - origin [type][0];
			newVector [1] = points.get(i).getY() - origin [type][1];
			newVector [2] = 0;
			if(tools.getVectorLength(newVector) >= minDistance){
				origin [type][0] = points.get(i).getX();
				origin [type][1] = points.get(i).getY();
				origin [type][2] = points.get(i).getZ(type);	
				counter = i;
				skippedPoints = i;
				break searchStart;
			}
		}
		
		
		orientationVector [type][0] = 0.0;
		orientationVector [type][1] = 0.0;
		orientationVector [type][2] = 0.0;
		
		
		
		searching: while(true){
			if(points.size()==counter+1)	break searching;
			
			newVector [0] = points.get(counter+1).getX() - origin [type][0];
			newVector [1] = points.get(counter+1).getY() - origin [type][1];
			newVector [2] = 0; //points.get(counter+1).getZ(type) - origin [type][2];
//			if(Math.sqrt(Math.pow(xComponent, 2.0)+Math.pow(yComponent, 2.0)+Math.pow(zComponent, 2.0)) > maxDistance){
//			if(tools.getVectorLength(newVector) < minDistance){
//				counter++; 
//				skippedPoints++;
//				continue;
//			}
//			
			if(tools.getVectorLength(newVector) > maxDistance){
				break searching;
			}
					
			newVector = tools.getNormalizedVector(newVector);
//			if(Double.isNaN(newVector [0]))	IJ.log("T" + frame + " newV0 NaN:! type " + type + " ct " + counter);
//			if(Double.isNaN(newVector [1]))	IJ.log("T" + frame + " newV1 NaN:! type " + type + " ct " + counter);
//			if(Double.isNaN(newVector [2]))	IJ.log("T" + frame + " newV2 NaN:! type " + type + " ct " + counter);
			
			if(!Double.isNaN(newVector [0]) && !Double.isNaN(newVector [1]) && !Double.isNaN(newVector [2])){
				orientationVector [type][0] += newVector [0];
				orientationVector [type][1] += newVector [1];
				orientationVector [type][2] += newVector [2];
			}else{
				skippedPoints++;
			}
			counter++;	
		}
		counter -= skippedPoints;
		if(counter!=0){
			orientationVector [type][0] /= (double) counter;
			orientationVector [type][1] /= (double) counter;
			orientationVector [type][2] /= (double) counter;
//			if(Double.isNaN(orientationVector [type][0]))	IJ.log("T" + frame + " ov0 NaN: origin " + origin [type][0] + " - " + origin [type][1] + " - " + origin [type][2] + "! type " + type + " ct " + counter);
//			if(Double.isNaN(orientationVector [type][1]))	IJ.log("T" + frame + " ov1 NaN: origin " + origin [type][0] + " - " + origin [type][1] + " - " + origin [type][2] + "! type " + type + " ct " + counter);
//			if(Double.isNaN(orientationVector [type][2]))	IJ.log("T" + frame + " ov2 NaN: origin " + origin [type][0] + " - " + origin [type][1] + " - " + origin [type][2] + "! type " + type + " ct " + counter);
//			if(Double.isNaN(orientationVector [type][0]) || Double.isNaN(orientationVector [type][1]) || Double.isNaN(orientationVector [type][2])){
//				IJ.log("problem " + origin [type][0] + "-" + origin [type][1] + "-" + origin [type][2] + "! " + frame + " type" + type);
//			}
			if(orientationVector [type][0]==0.0
					&&orientationVector [type][1]==0.0
					&&orientationVector [type][2]==0.0){
				
			}
					
					
			theta [type] = tools.getAbsoluteAngle(orientationVector[type], constants.X_AXIS);
			oriented = true;
		}else{
			IJ.log("T"  + frame + " no orientation possible (0 points in proximity)");
		}
			
//		if(Double.isNaN(orientationVector [type][0]))	IJ.log("T" + frame + "ov0 NaN" + origin [type][0] + "-" + origin [type][1] + "-" + origin [type][2] + "! type" + type + " ct" + counter);
//		if(Double.isNaN(orientationVector [type][1]))	IJ.log("T" + frame + "ov1 NaN" + origin [type][0] + "-" + origin [type][1] + "-" + origin [type][2] + "! type" + type + " ct" + counter);
//		if(Double.isNaN(orientationVector [type][2]))	IJ.log("T" + frame + "ov2 NaN" + origin [type][0] + "-" + origin [type][1] + "-" + origin [type][2] + "! type" + type + " ct" + counter);
//		if(orientationVector [type][0]==0.0
//				&&orientationVector [type][1]==0.0
//				&&orientationVector [type][2]==0.0){
//			IJ.log("problem " + origin [type][0] + "-" + origin [type][1] + "-" + origin [type][2] + "! " + frame + " type" + type);
//		}
//		if(counter==0)IJ.log("problem ct=0 " + origin [type][0] + "-" + origin [type][1] + "-" + origin [type][2] + "! " + frame + " type" + type);
//		if(Double.isNaN(theta [type])){
//			IJ.log("theta NAN ");
//			IJ.log("  origin " + origin [type][0] + "-" + origin [type][1] + "-" + origin [type][2] + "! " + frame + " type" + type + " ct" + counter);
//			IJ.log("  ovec " + orientationVector [type][0] + "-" + orientationVector [type][1] + "-" + orientationVector [type][2] + "! " + frame + " type" + type + " ct" + counter);
//		}
	}
	
	
	
	public double [][] getOrientationVector (){
		return orientationVector;
	}
	
	public void calculateOrientedPoints3D(){		
		for(int i = 0; i < points.size(); i++){
			points.get(i).setOrientedVector3D (orientationVector, origin, constants.X_AXIS, this.theta, multi_focal_tools.NOZ);
		}
		if(g4pz){
			for(int i = 0; i < points.size(); i++){
				points.get(i).setOrientedVector3D (orientationVector, origin, constants.X_AXIS, this.theta, multi_focal_tools.PUREZ);
			}
			if(linearInterpolatedZ){
				for(int i = 0; i < points.size(); i++){
					points.get(i).setOrientedVector3D (orientationVector, origin, constants.X_AXIS, this.theta, multi_focal_tools.MEDIANZ);
					points.get(i).setOrientedVector3D (orientationVector, origin, constants.X_AXIS, this.theta, multi_focal_tools.MEANZ);
				}
			}
		}
	}
		
	public void calculateOrientedPoints2D(){		
		for(int i = 0; i < points.size(); i++){
			points.get(i).setOrientedVector2D(origin [multi_focal_tools.NOZ], theta[multi_focal_tools.NOZ]);
		}
	}
	
	public void orientToHeadPosition(){		
		if(g4pz){
			for(int i = 0; i < points.size(); i++){
				points.get(i).setOrientedVector3D (points.get(i).get3DOrientedVector(multi_focal_tools.PUREZ), origin,
						constants.Z_AXIS, this.hrCOMAngle, multi_focal_tools.PUREZ);
			}
			if(linearInterpolatedZ){
				for(int i = 0; i < points.size(); i++){
					points.get(i).setOrientedVector3D (points.get(i).get3DOrientedVector(multi_focal_tools.MEDIANZ), origin,
							constants.Z_AXIS, this.hrCOMAngle, multi_focal_tools.MEDIANZ);
					points.get(i).setOrientedVector3D (points.get(i).get3DOrientedVector(multi_focal_tools.MEANZ), origin,
							constants.Z_AXIS, this.hrCOMAngle, multi_focal_tools.MEANZ);
				}
			}
		}
	}
	
	public void interpolateZLinear(double acceptedDist, int plusMinusPoints){	
//		IJ.log("smoothing started");
		int nrOfPoints = points.size();
		int plusMinusPointsThird = (int)((double)plusMinusPoints/3.0);
		
		//count possible variants
		int ipVariants = 1;
		for(int d1 = -plusMinusPoints; d1 <= plusMinusPoints; d1++){
			for(int d2 = d1+1; d2 <= plusMinusPoints; d2++){
				ipVariants++;
			}
		}
		
		double meanZ;
		int interPolCounter;
		for(int j = 0; j < nrOfPoints; j++){	
			double medianZ [] = new double [ipVariants];
			medianZ [0] = points.get(j).getZ(multi_focal_tools.PUREZ);
			for(int m = 1; m < ipVariants; m++){	
				medianZ [m] = Double.MAX_VALUE;
			}
			
			interPolCounter = 1;
			meanZ = medianZ [0];
//			if(j == 10) IJ.log("medi0 " + median [0]);
			
			double distance1 = 0.0;
			double distance2 = 0.0;
			double value1 = 0.0;
			double value2 = 0.0;
			
			for(int d1 = -plusMinusPoints; d1 <= plusMinusPoints; d1++){
				for(int d2 = d1+1; d2 <= plusMinusPoints; d2++){
					if(d1!=0&&d2!=0){
						if(j+d1>0&&j+d1<nrOfPoints&&j+d2>0&&j+d2<nrOfPoints){
							if(tools.mathAbs(d2 - d1) > plusMinusPointsThird){
								distance1 = (tools.mathAbs(d1)/d1) * multi_focal_tools.getDistance(points.get(j+d1), points.get(j), multi_focal_tools.NOZ);
								distance2 = (tools.mathAbs(d2)/d2) * multi_focal_tools.getDistance(points.get(j+d2), points.get(j), multi_focal_tools.NOZ);

								if(Math.abs(distance1) <= acceptedDist 
										&& Math.abs(distance2) <= acceptedDist
										&& tools.mathAbs(distance2-distance1) > 0.2){
									
									value1 = points.get(j+d1).getZ(multi_focal_tools.PUREZ);								
									if(j+d1-1>=0){
										value1 += points.get(j+d1-1).getZ(multi_focal_tools.PUREZ);
										if(j+d1+1 < points.size()){
											value1 += points.get(j+d1+1).getZ(multi_focal_tools.PUREZ);
											value1 /= 3.0;
										}else{
											value1 /= 2.0;
										}
									}else if(j+d1+1 < points.size()){
										value1 += points.get(j+d1+1).getZ(multi_focal_tools.PUREZ);
										value1 /= 2.0;
									}
									
									value2 = points.get(j+d2).getZ(multi_focal_tools.PUREZ);
									if(j+d2-1>=0){
										value2 += points.get(j+d2-1).getZ(multi_focal_tools.PUREZ);
										if(j+d2+1 < points.size()){
											value2 += points.get(j+d2+1).getZ(multi_focal_tools.PUREZ);
											value2 /= 3.0;
										}else{
											value2 /= 2.0;
										}
									}else if(j+d2+1 < points.size()){
										value2 += points.get(j+d2+1).getZ(multi_focal_tools.PUREZ);
										value2 /= 2.0;
									}
									
									medianZ [interPolCounter] = tools.getInterpolatedValue1D(0.0,
											distance1, distance2, value1, value2);
																		
									if(medianZ [interPolCounter] < -100.0 || medianZ [interPolCounter] > 100.0 ){	
										//TODO eventually needs adjustment
										IJ.log("t" + this.getFrame() + " median at " + j + "=" + medianZ [interPolCounter] + " - "
												+ " d1=" + d1 + " d2=" + d2);
										IJ.log("zd1 " + points.get(j+d1).getZ(multi_focal_tools.PUREZ) + " " + "zd2 " + points.get(j+d2).getZ(multi_focal_tools.PUREZ));
										IJ.log("distance1 " + distance1 + " - distance2 " + distance2);
										
										medianZ [interPolCounter] = 0.0;
									}else{
										meanZ += medianZ [interPolCounter];
										interPolCounter++;
									}							
								}	
							}														
						}
					}
				}
			}
			meanZ /= (double) interPolCounter;	
			points.get(j).setInterpolatedZ(meanZ, false);
			if(interPolCounter!=1){
				Arrays.sort(medianZ);
				try{
					if((interPolCounter)%2==0){
						points.get(j).setInterpolatedZ((medianZ[(interPolCounter)/2-1]+medianZ[(interPolCounter)/2])/2.0, true);
					}else{
						points.get(j).setInterpolatedZ(medianZ[(int)((double)(interPolCounter)/2.0)], true);
					}
				}catch(Exception e){
					IJ.log("problem at t" +frame +" - ct="+ interPolCounter);
					IJ.log("j" + j);
				}
			}
		}				
		linearInterpolatedZ = true;
	}
	
	public void removePoint(int index){
		points.remove(index);
		points.trimToSize();
	}
	
	public void trimList(){
		points.trimToSize();
	}
	
	public void calculateArcLengthPositions(){
		if(points.size()==0)	return;
		double arcLength = 0.0;
		points.get(0).setArcLengthPos(0.0);
		for(int i = 1; i < points.size(); i++){
			arcLength += multi_focal_tools.getDistance(points.get(i),points.get(i-1),multi_focal_tools.NOZ);
			points.get(i).setArcLengthPos(arcLength);
		}	
		double [] length;
		if(g4pz && linearInterpolatedZ == false){
			length = new double [1];
			length [0] = 0.0;
			points.get(0).setArcLength3DPos(length);
			for(int i = 1; i < points.size(); i++){
				length [0] += multi_focal_tools.getDistance(points.get(i),points.get(i-1),multi_focal_tools.PUREZ);
				points.get(i).setArcLength3DPos(length);
			}
		}else if(g4pz && linearInterpolatedZ){
			length = new double [] {0.0, 0.0, 0.0};
			points.get(0).setArcLength3DPos(length);
			for(int i = 1; i < points.size(); i++){
				length [multi_focal_tools.PUREZ-1] += multi_focal_tools.getDistance(points.get(i),points.get(i-1),multi_focal_tools.PUREZ);
				length [multi_focal_tools.MEDIANZ-1] += multi_focal_tools.getDistance(points.get(i),points.get(i-1),multi_focal_tools.MEDIANZ);
				length [multi_focal_tools.MEANZ-1] += multi_focal_tools.getDistance(points.get(i),points.get(i-1),multi_focal_tools.MEANZ);
				points.get(i).setArcLength3DPos(length);
			}
		}
	}
	
	public void calculateDAnglesAndCurvatures(double upstreamDist, int preventPoints, int encoding){
		double halfUpDist = upstreamDist / 2.0;
		trackPoint p = points.get(0), q = points.get(0);
		
		//calculate Tensions
		for(int j = 0; j < points.size(); j++){	
			searchUpstream: for(int vu = 0; vu <= j; vu++){
				p = points.get(j-vu);
				if(tools.mathAbs(p.getArcLengthPos(encoding)-points.get(j).getArcLengthPos(encoding))/halfUpDist>=1){
					break searchUpstream;
				}
			}			
			searchUpstream: for(int vu = 0; vu + j < points.size(); vu++){
				q = points.get(j+vu);
				if(tools.mathAbs(q.getArcLengthPos(encoding)-points.get(j).getArcLengthPos(encoding))/halfUpDist>=1){
					break searchUpstream;
				}
			}			
			if(encoding == multi_focal_tools.NOZ){
				points.get(j).setCurvature2D(multi_focal_tools.getCurvatureFactor(p, q, encoding));
			}else{
				points.get(j).setCurvature3D(multi_focal_tools.getCurvatureFactor(p, q, encoding));
			}
		}
		
		//calculate dAngle
		p = points.get(0);
		q = points.get(0);
		double angle;
		for(int j = 0; j < points.size(); j++){
			searchUpstream: for(int vu = 0; vu <= j; vu++){
				p = points.get(j-vu);
				if(j-vu >= preventPoints) q = points.get(j-vu);
				
				if(tools.mathAbs(p.getArcLengthPos(encoding)-points.get(j).getArcLengthPos(encoding))/upstreamDist>=1){
					break searchUpstream;
				}
			}	
			if(encoding == multi_focal_tools.NOZ){
				angle = tools.get2DAngle(points.get(j).getVector(), p.getVector());
				points.get(j).setCAngle2D(Math.toDegrees(angle));
			}else{
				angle = tools.getAbsoluteAngle(new double [] {points.get(j).getVector()[0],points.get(j).getVector()[1],points.get(j).getZ(encoding)},
						new double [] {p.getVector()[0],p.getVector()[1],p.getZ(encoding)});
				points.get(j).setCAngle3D(Math.toDegrees(angle));
			}			
		}
	}
	
	public void determineHeadRotation2D (ImagePlus imp, double [] slicePositions, int plusMinusRange){
		//TODO establish for all encoding strategies!
		this.hrMatrix = new double [4][plusMinusRange*2+1];
		
		Arrays.fill(this.hrCOMVector, 0.0);
		{
			double [] normalOV = tools.getNormalVectorXY(this.orientationVector[multi_focal_tools.NOZ]);
			normalOV = tools.getNormalizedVector(normalOV);
			double [] OV = tools.getNormalizedVector(this.orientationVector[multi_focal_tools.NOZ]);
			
			double centerZPos = 0.0;
			for(int s = 0; s < imp.getNSlices(); s++){
				centerZPos += slicePositions [s];
			}					
			centerZPos /= 4.0;
			
			double intensitySum = 0.0;
			double intensity = 0.0;
			double point [] = new double [] {0.0,0.0};
			double normalCorrection [] = {0.0,0.0,0.0,0.0};
			double normalVector [] = {0.0,0.0};
			{
				int normalCorrCt [] = {0,0,0,0};
				int normalVectorCt = 0;
				double widthParams [][];
				for(int j = 0; j <= 4 && j < points.size(); j++){
					widthParams = points.get(j).getXYGaussCenter();
					for(int i = 0; i < widthParams.length; i++){
						if(widthParams [i][0] > 0.0){
							normalCorrection [(int)widthParams[i][1]] += widthParams [i][0];
							normalCorrCt [(int)widthParams[i][1]] ++;							
						}						
					}	
					normalVector [0] += points.get(j).getNormalVector() [0];
					normalVector [1] += points.get(j).getNormalVector() [1];
					normalVectorCt ++;
				}
				normalVector [0] /= normalVectorCt;
				normalVector [1] /= normalVectorCt;
				normalVector = tools.getNormalizedVector(normalVector);
				for(int j = 0; j < 4; j++){
					if(normalCorrCt [j] != 0 && normalCorrection [j] != 0.0){
						normalCorrection [j] /= normalCorrCt [j];
					}					
				}
			}
			
			for(int s = 0; s < imp.getNSlices(); s++){
				for(int i = -plusMinusRange; i <= plusMinusRange; i++){
					for(int j = -4; j <= 4; j++){
						point [0] = normalVector [0] * normalCorrection [s] + this.origin [multi_focal_tools.NOZ][0] + normalOV [0] * i * imp.getCalibration().pixelWidth;
						point [1] = normalVector [1] * normalCorrection [s] + this.origin [multi_focal_tools.NOZ][1] + normalOV [1] * i * imp.getCalibration().pixelHeight;
						
						intensity = impProcessing.getInterpolatedIntensity2D(imp,
								point [0] + OV [0] * j * imp.getCalibration().pixelWidth,
								point [1] + OV [1] * j * imp.getCalibration().pixelHeight,
								imp.getStackIndex(1, s+1, frame+1)-1);
						
						this.hrCOMVector [0] += intensity * (point [0] - this.origin [multi_focal_tools.NOZ][0]);
						this.hrCOMVector [1] += intensity * (point [1] - this.origin [multi_focal_tools.NOZ][1]);
						this.hrCOMVector [2] += intensity * (slicePositions [s] - centerZPos);
						//TODO establish for all encoding strategies!
						intensitySum += intensity;
						
						hrMatrix [s][plusMinusRange+i] += intensity;
					}
					hrMatrix [s][plusMinusRange+i] /= 9;
				}
				
				hrMaximumPosition [s] = (tools.getMaximumIndex(hrMatrix [s]) - plusMinusRange) * imp.getCalibration().pixelWidth;
				hrMaximumIntensity [s] = tools.getMaximum(hrMatrix [s]);
			}
			
			hrCOMVector [0] /= intensitySum;
			hrCOMVector [1] /= intensitySum;
			hrCOMVector [2] /= intensitySum;		
			hrImageCOM [0] = hrCOMVector [0] + plusMinusRange * imp.getCalibration().pixelWidth;
			hrImageCOM [1] = hrCOMVector [1] + plusMinusRange * imp.getCalibration().pixelWidth;
			hrCOMVector = tools.getNormalizedVector(hrCOMVector);
			hrCOMAngle = tools.getAbsoluteAngle(constants.Z_AXIS, hrCOMVector);
			hrMaximaLineVector = multi_focal_tools.get2DVectorFromPoints(this.hrMaximumPosition);
			hrMaximaLineVector = tools.getNormalizedVector(hrMaximaLineVector);
			hrMaximaLineAngle = tools.get2DAngle(hrMaximaLineVector, constants.X_AXIS_2D);
		}
		System.gc();
		hrVset = true;
	}
		
	public double [] getHeadRotationAngle (double [] minIntensity, double [] maxIntensity, double [] minDisplacement, double [] maxDisplacement){
		double [] vector = {0.0,0.0};
		for(int i = 0; i < 4; i++){
			vector [0] = tools.getNormalizedValue(hrMaximumPosition [i], minDisplacement [i], maxDisplacement [i]) - 0.5;
			vector [1] = tools.getNormalizedValue(hrMaximumIntensity [i], minIntensity [i], maxIntensity [i]) - 0.5;
			hrAngle [i] = Math.toDegrees(tools.get2DAngle(vector, constants.X_AXIS_2D));
		}
		return hrAngle;
	}
	
	public double getHeadRotationCOMAngleInDegree(){
		if(!hrVset){
			return 0.0;
		}		
		return Math.toDegrees(hrCOMAngle);
	}
	
	public double getHeadRotationMaximaAngleInDegree(){
		if(!hrVset){
			return 0.0;
		}		
		return Math.toDegrees(hrMaximaLineAngle);
	}
	
	public double [] getHeadRotationCOMVector(){
		if(!hrVset){
			return new double [] {0.0,0.0,0.0};
		}		
		return hrCOMVector;
	}
	
	public double [] getHeadRotationMaximaLineVector(){
		if(!hrVset){
			return new double [] {0.0,0.0};
		}		
		return hrMaximaLineVector;
	}
	
	public double [] getHeadRotationMaximumPositions(){
		if(!hrVset){
			return new double [] {0.0,0.0,0.0,0.0};
		}		
		return hrMaximumPosition;
	}
	
	public double [] getHeadRotationMaximumIntensities(){
		if(!hrVset){
			return new double [] {0.0,0.0,0.0,0.0};
		}		
		return hrMaximumIntensity;
	}
	
	public double [] getHeadRotationMatrixCOM(){
		if(!hrVset){
			return new double [] {0.0,0.0};
		}
		
		return hrImageCOM;
	}
	
	public double [][] getHeadRotationMatrix(){
		if(!hrVset){
			IJ.error("head rotation vector set not initialized");
		}		
		return hrMatrix;
	}
	
	public double getArcLengthXY(){
		double length = 0.0;
		for(int i = 1; i < points.size(); i++){
			length += multi_focal_tools.getDistance(points.get(i), points.get(i-1),multi_focal_tools.NOZ);
		}
		return length;
	}
	
	public double getArcLengthXZ(int type){
		double length = 0.0;
		for(int i = 1; i < points.size(); i++){
			length += Math.sqrt(Math.pow(points.get(i).getX()-points.get(i-1).getX(), 2.0)
					+ Math.pow(points.get(i).getZ(type)-points.get(i-1).getZ(type), 2.0));
		}	
		return length;
	}
	
	public double getArcLengthYZ(int type){
		double length = 0.0;
		for(int i = 1; i < points.size(); i++){
			length += Math.sqrt(Math.pow(points.get(i).getY()-points.get(i-1).getY(), 2.0)
					+ Math.pow(points.get(i).getZ(type)-points.get(i-1).getZ(type), 2.0));
		}
		return length;
	}
	
	public double getArcLength3D(int type){
		double length = 0.0;
		for(int i = 1; i < points.size(); i++){
			length += multi_focal_tools.getDistance(points.get(i),points.get(i-1),type);
		}
		return length;
	}	
	
	public double getMinWidth(){
		double minWidth = Double.POSITIVE_INFINITY;
		double [][] widths;
		for(int i = 0; i < points.size(); i++){
			widths = points.get(i).getXYGaussWidth().clone();
			for(int j = 0; j < widths.length; j++){
				if(widths [j][0] < minWidth)	minWidth = widths [j][0];
			}
		}
		return minWidth;
	}
	
	public double getMaxWidth(){
		double maxWidth = 0.0;
		double [][] widths;
		for(int i = 0; i < points.size(); i++){
			widths = points.get(i).getXYGaussWidth().clone();
			for(int j = 0; j < widths.length; j++){
				if(widths [j][0] > maxWidth)	maxWidth = widths [j][0];
			}
		}
		return maxWidth;
	}
	
	public double getThetaDegree(int type){
		if(Double.isNaN(Math.toDegrees(theta[type]))){
			IJ.log("Time " + frame + ": in degree = " + Math.toDegrees(theta[type]) + " - orth " + theta[type]);
		}
		return Math.toDegrees(theta[type]);
//		return theta[type];
	}
}

