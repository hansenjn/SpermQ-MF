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

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Line;
import spermQ_mf.jnhsupport.*;

public class trackPoint {
	// properties
	private double cXp, cYp;
	private double [] pointIntensity;
	
	private double [] vector;
	private double vectorLength;
	private double [] normalVector;
	private double normalVectorLength;
	private double [] normalMaxIntensity;
	private double [][] widthGaussParams;	//[][] gaussresults+4=sliceindex,resultnr
	private double correction;
	private double [] gauss4PParams; 
	
	private double gaussZ4P;
	private double gaussHeight4P;
	private double [][] widthZ;
	private double zByWidth;
	private double widthToZCorrFactor = 1.0;
	private double finalZ = Double.NEGATIVE_INFINITY;	
	private int finalZIncl = 0;
	private String selectedZCombo = "";
	private double selectedZComboSD = 0.0;
	private double curvature2D, curvature3D;	
	private double cAngle2D, cAngle3D;	
	
	private double meanInterpolatedZ, medianInterpolatedZ;
	
	private double arcLengthXYPos;
	private double [] arcLengthXYZPos;
	private double [] orientedVectorBy2D;
	private double [][] orientedVectorBy3D;
	
	public trackPoint(double px, double py){
		cXp = px;
		cYp = py;	
		
		vector = new double [] {0.0,0.0};
		vectorLength = 0.0;
		normalVector = new double [] {0.0,0.0};
		normalVectorLength = 0.0;
		correction = 0.0;
		
		gaussZ4P = Double.NEGATIVE_INFINITY;	
		gaussHeight4P = Double.NEGATIVE_INFINITY;	
		widthZ = new double [1][2];
		widthZ [0][0] = Double.NEGATIVE_INFINITY;
		widthZ [0][1] = Double.NEGATIVE_INFINITY;
		zByWidth = Double.NEGATIVE_INFINITY;
		orientedVectorBy2D = new double [2];
		orientedVectorBy3D = new double [4][3];
				
		arcLengthXYZPos = new double [3];
	}	
	
	//*********************
	//SETS
	//*********************	
	public void setVectors(double vx, double vy){
		//set normalized vector
			vectorLength = Math.sqrt(Math.pow(vx,2.0)+Math.pow(vy,2.0));
			vector [0] = vx / vectorLength;
			vector [1] = vy / vectorLength;
			vectorLength = Math.sqrt(Math.pow(vector[0],2.0)+Math.pow(vector[1],2.0));
			
		//define normal vector
			normalVector = tools.getNormalVectorXY(vector);
			normalVectorLength = Math.sqrt(Math.pow(normalVector[0],2.0)+Math.pow(normalVector[1],2.0));
	}
	
	public void widthGaussFit (int calculationRadius, double [][] xNormalPoints, double [][] yNormalPoints, double [][] parameters, int nrOfFits, boolean correctCoords){
//		normalLineRadius = calculationRadius;
		widthGaussParams = new double [nrOfFits][5];
		if(vectorLength!=0.0){
			//y offset, height, center, width, slicePos
			double shift = 0.0;
			for(int i = 0; i < nrOfFits; i++){
				shift += parameters [2][i];
				
				//copy parameters and determine z
				for(int j = 0; j < 5; j++){
					widthGaussParams [i][j] = parameters [j][i];
				}				
			}
			shift /= nrOfFits;
			if(correctCoords){
				this.correctCoordinates((shift/normalVectorLength) * normalVector[0], (shift/normalVectorLength) * normalVector[1]);
			}else{
				shift = 0.0;
			}
			correction = shift;
			
			//save interpolated point intensity and normal max intensity
			pointIntensity = new double [xNormalPoints.length];
			normalMaxIntensity = new double [xNormalPoints.length];
								
			for(int i = 0; i < xNormalPoints.length; i++){
				if(Math.abs((int)(shift/normalVectorLength) + (int)(shift/Math.abs(shift))) <= calculationRadius){
					pointIntensity [i] = tools.getInterpolatedValue1D(shift, xNormalPoints[i][calculationRadius + 1 + (int)(shift/normalVectorLength)],
							xNormalPoints[i][calculationRadius + (int)(shift/normalVectorLength) + (int)(shift/Math.abs(shift))],
							yNormalPoints[i][calculationRadius + (int)(shift/normalVectorLength)],
							yNormalPoints[i][calculationRadius + (int)(shift/normalVectorLength) + (int)(shift/Math.abs(shift))]);
				}else{
					pointIntensity [i] = tools.getInterpolatedValue1D(shift, xNormalPoints[i][calculationRadius + 1 + (int)(shift/normalVectorLength)],
							xNormalPoints[i][calculationRadius + (int)(shift/normalVectorLength) - (int)(shift/Math.abs(shift))],
							yNormalPoints[i][calculationRadius + (int)(shift/normalVectorLength)],
							yNormalPoints[i][calculationRadius + (int)(shift/normalVectorLength) - (int)(shift/Math.abs(shift))]);
				}
				
				normalMaxIntensity [i] = 0.0;
				for(int j = 0; j < yNormalPoints[i].length; j++){
					if(yNormalPoints[i][j]>normalMaxIntensity[i]){
						normalMaxIntensity[i] = yNormalPoints[i][j];
					}					
				}
			}		
		}else{
			IJ.error("trackPoint not initialized!");
		}		
	}
	
	private void correctCoordinates(double xShift, double yShift){
		cXp = cXp + xShift;
		cYp = cYp + yShift;
	}
	
	public void resetXCoordinate(double x){
		cXp = x;
	}
	
	public void resetYCoordinate(double y){
		cYp = y;
	}
	
	public void setZGaussFit4P (double [] parameters){
		gauss4PParams = new double [parameters.length];
		for(int i = 0; i < parameters.length; i++){
			gauss4PParams [i] = parameters [i] + 0.0;
		}		
		//y offset, height, center, width
		gaussZ4P = gauss4PParams [2];
		gaussHeight4P = gauss4PParams [0] + gauss4PParams [1];  
	}
	
	public void addZWidthFit (double wZminus, double wZplus){
		if(widthZ [0][0] == Double.NEGATIVE_INFINITY){
			widthZ [0][0] = wZminus;
			widthZ [0][1] = wZplus;			
		}else{
			double [][] arrayCopy = new double [widthZ.length][2];
			for(int i = 0; i < widthZ.length; i++){
				arrayCopy [i][0] = widthZ [i][0];
				arrayCopy [i][1] = widthZ [i][1];
			}	
			
			widthZ = new double [widthZ.length + 1][2];
			widthZ [0][0] = wZminus;
			widthZ [0][1] = wZplus;			
			for(int i = 0; i < arrayCopy.length; i++){
				widthZ [i+1][0] = arrayCopy [i][0];
				widthZ [i+1][1] = arrayCopy [i][1];
			}
		}	
	}
		
	public void setInterpolatedZ (double interpolated, boolean median){
		if(median){
			medianInterpolatedZ = interpolated;
		}else{
			meanInterpolatedZ = interpolated;
		}
	}
	
	public void setArcLengthPos(double pos){
		arcLengthXYPos = pos;
	}
	
	public void setArcLength3DPos(double [] pos){
		for(int i = 0; i < pos.length; i++){
			arcLengthXYZPos [i] = pos [i];
		}
	}
	
	public void setOrientedVector3D(double [] referenceVector, double [][] origin, double [] targetVector, double angle, int smoothing){				
		//point to rotate x,y,z
		double x = this.getX();
		double y = this.getY();
		double z = 0.0;		if(smoothing!=multi_focal_tools.NOZ) z = this.getZ(smoothing);
		
		//base point of the line
		double a = origin [smoothing][0];
		double b = origin [smoothing][1];
		double c = origin [smoothing][2];
		
		//normalized vector (u^2+v^2+w^2 = 1) of the line to be rotated around 
		//vector to be rotated around equals the vector of the cross-product of the reference vector and the vector of the target line
		double [] areaNormalVector = tools.crossProduct(referenceVector, targetVector);
		double [] normalizedNormalVector = tools.getNormalizedVector(areaNormalVector);
		double u = normalizedNormalVector [0];
		double v = normalizedNormalVector [1];
		double w = normalizedNormalVector [2];
		
		//rotation
		orientedVectorBy3D [smoothing][0] = (a*(v*v+w*w)-u*(b*v+c*w-u*x-v*y-w*z))*(1-Math.cos(angle))
				+ x * Math.cos(angle)
				+ ((-1)*c*v+b*w-w*y+v*z) * Math.sin(angle) - origin [smoothing][0];
		
		orientedVectorBy3D [smoothing][1] = (b*(u*u+w*w)-v*(a*u+c*w-u*x-v*y-w*z))*(1-Math.cos(angle))
				+ y * Math.cos(angle)
				+ (c*u-a*w+w*x-u*z) * Math.sin(angle) - origin [smoothing][1];
		
		orientedVectorBy3D [smoothing][2] = (c*(u*u+v*v)-w*(a*u+b*v-u*x-v*y-w*z))*(1-Math.cos(angle))
				+ z * Math.cos(angle)
				+ ((-1)*b*u+a*v-v*x+u*y) * Math.sin(angle) - origin [smoothing][2];	
	}
	
	public void setOrientedVector3D(double [][] referenceVector, double [][] origin, double [] targetVector, double [] angle, int smoothing){
		this.setOrientedVector3D(referenceVector [smoothing], origin, targetVector, angle [smoothing], smoothing);
	}
	
	public void setOrientedVector2D(double [] origin, double angle){
		orientedVectorBy2D [0] = (this.getX() - origin [0]) * Math.cos(angle) + (this.getY() - origin [1]) * Math.sin(angle);
		orientedVectorBy2D [1] = -1 * (this.getX() - origin [0]) * Math.sin(angle) + (this.getY() - origin [1]) * Math.cos(angle);	
	}
			
	public void setWidhtZCorrectionFactor(double factor){
		widthToZCorrFactor = factor;
		finalZ = Double.NEGATIVE_INFINITY;
	}
	
	/**
	 * @param curvature3d the curvature3D to set
	 */
	public void setCurvature3D(double curvature3d) {
		curvature3D = curvature3d;
	}
	
	/**
	 * @param curvature2d the curvature2D to set
	 */
	public void setCurvature2D(double curvature2d) {
		curvature2D = curvature2d;
	}
	
	/**
	 * @param cAngle3D the cAngle3D to set
	 */
	public void setCAngle3D(double cAngle3D) {
		this.cAngle3D = cAngle3D;
	}

	/**
	 * @param cAngle2D the cAngle2D to set
	 */
	public void setCAngle2D(double cAngle2D) {
		this.cAngle2D = cAngle2D;
	}
	
	//*********************
	//GETS
	//*********************
	/**
	 * @deprecated
	 * */
	public double getOrX(){
		return cXp;
	}
	
	/**
	 * @deprecated
	 * */
	public double getOrY(){
		return cYp;
	}
	
	public double getX(){
		return cXp;	
	}
	
	public double getY(){
		return cYp;		
	}
	

	public double getNormalVectorLength(){
		return normalVectorLength;
	}
	
	public double [] getNormalVector(){
		return normalVector;
	}
	
	public Line getVectorLine(ImagePlus imp){
		if(vectorLength==0.0) return null;
		double x = cXp / imp.getCalibration().pixelWidth,
				y = cYp / imp.getCalibration().pixelHeight;
		return new Line(x, y, x + vector [0], y + vector [1]);
	}
	
	public Line getNormalVectorLine(ImagePlus imp, double radius){
		if(normalVectorLength==0.0) return null;
		double x = cXp / imp.getCalibration().pixelWidth,
			y = cYp / imp.getCalibration().pixelHeight;
		return new Line(x + normalVector [0]*(radius), 		//x1
				y + normalVector [1]*(radius),				//y1
				x + (-1.0) * normalVector [0] * (radius), 	//x2
				y + (-1.0) * normalVector [1] * (radius));	//y2
	}
		
	public double getNormalMaxIntensity(int slice){
		return normalMaxIntensity [slice];
	}
		
	
	public double getNormalPointIntensity(int slice){
		return pointIntensity [slice];
	}
	
	public double [][] getXYGaussWidth (){
		//widthGaussParams: y offset, height, center, width, slicePos
		double [][] parameters = new double [widthGaussParams.length][2];	//0 = width, 1 = slicePos
		for(int i = 0; i < widthGaussParams.length; i++){
			parameters [i][0] = widthGaussParams [i][3];
			parameters [i][1] = widthGaussParams [i][4];
		}
		return parameters;
	}
	
	public double [][] getXYGaussCenter (){
		//widthGaussParams: y offset, height, center, width, slicePos
		double [][] parameters = new double [widthGaussParams.length][2];	//0 = center, 1 = slicePos
		for(int i = 0; i < widthGaussParams.length; i++){
			if(correction != 0.0){
				parameters [i][0] = widthGaussParams [i][2] - correction;
			}else{
				parameters [i][0] = 0.0;
			}
			parameters [i][1] = widthGaussParams [i][4];
		}
		return parameters;
	}
	
	public double [][] getXYGaussHeight (){
		//widthGaussParams: y offset, height, center, width, slicePos
		double [][] parameters = new double [widthGaussParams.length][2];	//0 = maximum, 1 = slicePos
		for(int i = 0; i < widthGaussParams.length; i++){
			parameters [i][0] = widthGaussParams [i][1];
			parameters [i][1] = widthGaussParams [i][4];
		}
		return parameters;
	}
	
	public double getArcLengthPos(int type){
		if(type == multi_focal_tools.NOZ){
			return arcLengthXYPos;
		}else{
			return arcLengthXYZPos [type-1];
		}		
	}
	
	public double getGaussFitMax(){
		if(gaussHeight4P != Double.NEGATIVE_INFINITY){
			return gaussHeight4P;
		}
		return 0.0;
	}
	
	public double getZ(int smoothing){
		if(smoothing == multi_focal_tools.NOZ){
			return 0.0;
		}else if(smoothing == multi_focal_tools.PUREZ){
			if(gaussZ4P == Double.NEGATIVE_INFINITY && widthZ[0][0] == Double.NEGATIVE_INFINITY){
				IJ.error("z not initialized");
				return 0.0;
			}else if(widthZ[0][0] == Double.NEGATIVE_INFINITY && gaussZ4P != Double.NEGATIVE_INFINITY){
//				IJ.log("no width fit");
				return gaussZ4P;			
			}else if(finalZ != Double.NEGATIVE_INFINITY && finalZIncl == widthZ.length){
//				IJ.log("finalZ returned");
				return finalZ;				
			}else{//if width fit data are available calculate best variant by standard deviation
				//generate array
				int arraySize = 0;
				if(gaussZ4P != Double.NEGATIVE_INFINITY){
					arraySize++;
				}
				arraySize += widthZ.length;
				finalZIncl = widthZ.length;
				
				//calculate nr of combinations
				int nrOfCombinations = (int) Math.pow(2, widthZ.length);
//				for(int i = arraySize; i > 0; i--){
//					nrOfCombinations *= i;
//				}
				double [][] zValues = new double [nrOfCombinations][arraySize];
//				String [] combinations = new String [nrOfCombinations];
//				for(int nrC = 0; nrC < nrOfCombinations; nrC++){
//					combinations [nrC] = "";
//				}
				//create all combinations
				int comboCounter = 0;
				for(int wZPs = 0; wZPs < widthZ.length; wZPs++){
					comboCounter = 0;
					for(int nrC = 0; nrC < nrOfCombinations/Math.pow(2,(wZPs+1)); nrC++){					
						for(int i = 0; i < Math.pow(2,wZPs); i++){
//							combinations [comboCounter] += "0";
							zValues [comboCounter][wZPs] = widthZ [wZPs][0] * widthToZCorrFactor;	//minus
							comboCounter++;
						}
						for(int i = 0; i < Math.pow(2,wZPs); i++){
//							combinations [comboCounter] += "1";
							zValues [comboCounter][wZPs] = widthZ [wZPs][1] * widthToZCorrFactor;	//plus
							comboCounter++;
						}
					}
				}
				
//				{
//					IJ.log("comboCounter " + comboCounter + " arraylength " + nrOfCombinations);
//					String output = "";
//					for(int i = 0; i < nrOfCombinations; i++){
//						output = "z" + i + "";
//						for(int j = 0; j < widthZ.length; j++){
//							output += " " + zValues [i][j];
//						}		
//						IJ.log(output);
//					}
//				}
								
				double minSD = Double.MAX_VALUE;
				int selCombo = -1;
//				IJ.log("combinations");
				for(int nrC = 0; nrC < nrOfCombinations; nrC++){
//					IJ.log("combination " + combinations [nrC]);
//					IJ.log("zvPre " + zValues [nrC][arraySize-1]);
					if(gaussZ4P != Double.NEGATIVE_INFINITY)	zValues [nrC][arraySize-1] = gaussZ4P;
//					IJ.log("zvPost " + zValues [nrC][arraySize-1]);
					
					//calculate SD
					double SD = tools.getSD(zValues[nrC]);
					if(SD < minSD){
						minSD = SD;
						selCombo = nrC;						
					}
				}
				
				try{
					finalZ = tools.getMedian(zValues[selCombo]);
					if(gaussZ4P != Double.NEGATIVE_INFINITY){
						zByWidth = tools.getAverageOfRange(zValues[selCombo], 0, arraySize-2);
					}else{
						zByWidth = tools.getAverage(zValues[selCombo]);
					}
//					finalZ = multi_focal_tools.getAverage(zValues[nrC]);
				}catch(Exception e){
					IJ.log("nrOfC " + nrOfCombinations + " - selCombo " + selCombo + " - GaussZ4P " + gaussZ4P
							+ " - minSD " + minSD + " - widthZ[0][0] " + (widthZ[0][0] * widthToZCorrFactor) + " - widthZ[0][1]" + (widthZ[0][1] * widthToZCorrFactor) 
							+ " - comboCounter " + comboCounter + " - widthZLength " + widthZ.length
							+ " - arraysize " + arraySize);
					IJ.error("index error");
				}				
				
				for(int i = 0; i < zValues[selCombo].length; i++){
					this.selectedZCombo += (", " + zValues[selCombo][i]);
				}
				if(selectedZCombo.length()>0)	selectedZCombo = selectedZCombo.substring(2);
				
//				if(finalZ == Double.NEGATIVE_INFINITY)	IJ.log("finalZ is negative infinity");
//				if(finalZ == Double.POSITIVE_INFINITY)	IJ.log("finalZ is positive infinity");
				
				return finalZ;
				
//				IJ.log("fZ " + finalZ);
//				for(int i = 0; i < zValues[selCombo].length; i++){
//					IJ.log("   " + i + " - " + zValues[selCombo][i]);
//				}
//				IJ.log("allZ ");
//				for(int i = 0; i < widthZ.length; i++){
//					IJ.log("   m " + i + " - " + widthZ[i][0]);
//					IJ.log("   p " + i + " - " + widthZ[i][1]);
//				}
//				return finalZ;
			}			
		}else if(smoothing == multi_focal_tools.MEDIANZ){
//			if(medianInterpolatedZGauss4P==0.0)	IJ.log("returned medi " + medianInterpolatedZGauss4P);
			return medianInterpolatedZ;
		}else if(smoothing == multi_focal_tools.MEANZ){
			return meanInterpolatedZ;
		}
		return 0.0;
	}

	public double calculateWidthCorrFactor(){
		if(gaussZ4P == Double.NEGATIVE_INFINITY || zByWidth == Double.NEGATIVE_INFINITY){
			return 0.0;
		}
//		try{
			return gaussZ4P / zByWidth;
//		}catch (Exception e){			
//			return 0.0;
//		}
	}
	
	public boolean zDeterminable(){
		if(gaussZ4P != Double.NEGATIVE_INFINITY)	return true;
		if(widthZ.length>1)	return true;
		return false;
	}
	
	public double getZComboSD(){
		return this.selectedZComboSD;
	}
	
	public String getZCombination(){
		return this.selectedZCombo;
	}
	
	
	public double [] get3DOrientedVector(int type){
		return orientedVectorBy3D [type];
	}
	
	public double [] get2DOrientedVector(){
		return orientedVectorBy2D;
	}

	/**
	 * @return the curvature3D
	 */
	public double getCurvature3D() {
		return curvature3D;
	}
	
	/**
	 * @return the curvature2D
	 */
	public double getCurvature2D() {
		return curvature2D;
	}

	/**
	 * @return the cAngle3D
	 */
	public double getCAngle3D() {
		return cAngle3D;
	}

	/**
	 * @return the cAngle2D
	 */
	public double getCAngle2D() {
		return cAngle2D;
	}

	/**
	 * @return the tangential vector of the point
	 */
	public double [] getVector() {
		return vector;
	}

}
