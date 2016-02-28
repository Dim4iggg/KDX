import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.math3.stat.regression.*;

import jxl.Cell;
import jxl.CellType; 
import jxl.Sheet;
import jxl.Workbook;
import jxl.read.biff.BiffException;
import jxl.write.WritableCell;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;
import jxl.write.Number;

public class Main {

	private static ArrayList<DataPoint> points;
	private static String dataPath = "Eingabewerte.xls";
	
	private static int NUM_PSEUDOPOINTS = 2;
	private static int NUM_FITTINGPOS = 5;
	private static int O = 2;
			
	public static void main(String[] args) {
	
	    points  = new ArrayList<>();
	    
	  //read in data
	    try {
			ReadData();
		} catch (BiffException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	    //Visualizer.AddPointSet(points, "Sample Points");
	    
	    //generate equally distributed pseudo points (eg 5 points)
	    ArrayList<DataPoint> Mu = GeneratePseudoPoints();
	    //Visualizer.AddPointSet(Mu, "Pseudo Points");
	    
	    //chose F from S  (eg. every second point) -> here every point! F=S
	    //combine the locations of a subset of the historical instances in S with a set of different time points
	    ArrayList<DataPoint> F = GenerateFittingPositions();
	    //Visualizer.AddPointSet(F, "Fitting Positions");
	    
	    //visualize data
	    try {
			//AnalysisLauncher.open(new Visualizer());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
	    
	    FitKDX(points, F, Mu, new double[10][10], new double[3][3]);
	    
	    Visualizer app = new Visualizer(
                "KDX", Visualizer.PLOTTYPE.SCATTER);
        app.pack();
        app.setVisible(true);
        
	}


	private static void ReadData() throws BiffException, IOException
	{
		Workbook workbook = Workbook.getWorkbook(new File(dataPath)); 
		
		Sheet sheet =  workbook.getSheet(0);
		
		//read in data
		int row = 0;
		String s = "points: ";
		try
		{
		Cell cell = sheet.getCell(0,0);
		while(!cell.getContents().isEmpty())  //each row is a DataPoint
		{
			DataPoint point = new DataPoint();
			for(int column = 0; column<DataPoint.DIMENSIONS+1; column++)  //each row is a parameter
			{
				cell = sheet.getCell(column, row);
				CellType type = cell.getType();
				double zahl = Double.parseDouble(cell.getContents());
				
				if(column != DataPoint.TIME_COLUMN)
				{
					//set next attribute value
					point.AddData(zahl);
					s += zahl;
				}
				else
				{ 
					//set time
					point.SetTime(zahl);
					s += ", t=" + zahl + "; ";
				}
				
			}
			points.add(point);
			row++;
		}
		}
		catch(java.lang.ArrayIndexOutOfBoundsException e){};
		
		System.out.println(s);
	}
	
	private static ArrayList<DataPoint> GeneratePseudoPoints()
	{
		ArrayList<DataPoint> mu = new ArrayList<>();
		String s = "pseudo points: ";
		for(int i=0; i<NUM_PSEUDOPOINTS; i++)
		{
			DataPoint pPoint = new DataPoint();
			pPoint.SetTime(1);  //TODO: change?
			
			
			//distribute equally in each dimension
			for(int d=0; d<DataPoint.DIMENSIONS; d++)
			{
				double diff = DataPoint.maxValues[d] - DataPoint.minValues[d];  //eg 100
				double step = diff/(double)(NUM_PSEUDOPOINTS-1);  //if n=5  step = 25
				double pseudoValue = DataPoint.minValues[d] + step*i;  //points at 0, 25, 50, 75, 100
				pPoint.AddData(pseudoValue); 
				s += pseudoValue + ", t=" + pPoint.time + "; ";
			}
			mu.add(pPoint);
		}
		System.out.println(s);
		return mu;
	}
	
	private static ArrayList<DataPoint> GenerateFittingPositions()
	{
		//a set F of N !historical! (!IN THE PAST!), spatio-temporal fitting positions
		//combine the locations of a subset of the historical instances in S with a set of different time points.
		ArrayList<DataPoint> F = new ArrayList<>();
		
		if(NUM_FITTINGPOS > points.size())
		{
			NUM_FITTINGPOS = points.size(); //max as many as real points
		}
		String s = "fitting positions: ";
		for(int i=0; i<NUM_FITTINGPOS; i++)
		{
			DataPoint pPoint = new DataPoint();

			//distribute equally in time  //TODO: how to select time for fitting positions?
			double diff = DataPoint.maxTime - DataPoint.minTime;  
			double step = diff/(double)(NUM_FITTINGPOS-1);  
			double time = DataPoint.minTime + step*i;  
			pPoint.SetTime(time); 
			pPoint.values = points.get(i).values;  //just copy the position of a real point
			
			s += pPoint.values[0] + ", t=" + time + "; ";
			
			F.add(pPoint);
		}
		System.out.println(s);
	    return F;
	}
	
	private static double[] FitKDX(ArrayList<DataPoint> S, ArrayList<DataPoint> F, ArrayList<DataPoint> Mu, double[][] spatialBWidth, double[][] temporalBWidth)
	{
		int n = S.size();
		int N = F.size();
		int m = Mu.size();
		int M = (m-1)*(O+1);
		double[] k = new double[N];
		double[][] K = new double[N][m*(O+1)];
		
		int e=1;
		for(DataPoint p: S)
		{
			Visualizer.AddDensityLineSet(p,GaussianSpatialDensityKernel(p.values, p.values, spatialBWidth), "", spatialBWidth, 0);
			Visualizer.AddDensityLinePoint(p, GaussianSpatialDensityKernel(p.values, p.values, spatialBWidth), "S"+e, spatialBWidth,0);
			e++;
		}
		
		//get historic density estimates
		for(int j = 0; j<N; j++)
		{
			//k is a (N × 1)-vector that is obtained by spatiotemporal
			//density estimation for the N fitting positions !using the sample S!
			//as reference instances and pre-tuned spatial and temporal bandwidths.
			k[j] =  KDE(F.get(j).values, spatialBWidth, temporalBWidth , S);  
			
			System.out.println( "k[" + j + "] = " + k[j]);
			//add density line to visualization
			//Visualizer.AddDensityLineSet(F.get(j),k[j], j+"s point", spatialBWidth, 1);
			Visualizer.AddDensityLinePoint(F.get(j), k[j], "F"+j, spatialBWidth,1);
			
			
			//for each pseudo point compute its density with a kernel 
			// as the sum of gaussian densities of neighboring sample points
			for(int i=0; i<m; i++)
			{
				double kernel = GaussianSpatialDensityKernel(F.get(j).values, Mu.get(i).values, spatialBWidth);
				
				for(int o=0; o<O; o++)
				{
					K[j][i+o*m] = kernel * Math.pow(F.get(j).time, o);  //Kj,i+o·m = Ki(xj) · to j
				}
			}
		}
		
		double[] p = new double[N];
		double[][] P = new double[N][(m-1) * (O+1)];
		
		//create design matrix P and corresponding vector p
		for(int j=0; j<N; j++)
		{
			p[j] = k[j] - K[j][m];
			for(int i=0; i< m; i++)  //TODO: check im pseudocode steht m-1!
			{
				for(int o=0; o<O; o++)
				{
					P[j][i+o*(m-1)] = (K[j][i+o*m] - K[j][m+o*m])* Math.pow(F.get(j).time, o);
					if(P[j][i+o*(m-1)] == 0)
					{
						int a=0;
					}
				}
			}
		}
		
		//add regularisation terms to p and P
		for(int i=0; i<(m-1); i++)
		{
			for(int o=0; o<O; o++)
			{
				//p[N+i+o*(m-1)] = Cio; //TODO: Cio was ist das?
				//P[N+i+o*(m-1)][i+o*(m-1)] = 0; //TODO: da fehlt etwas!
				//P[N+i+o*(m-1)][i+o*(m-1)] = 1;
			}
		}
		
		//TODO: REMOVE! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//	if(true)
		//	return k;
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		
		
		//fit regression coefficients
		//double[] beta = new double[(m-1)*O+1];
		double[] beta = Regress(P, p, m, O);

		double[] alpha = new double[m*(O+1)]; //TODO: check size of alpha
		//reconstruct regression coefficient of alpha
		for(int o=0; o<O; o++)
		{
			double sum =0;
			for(int i=0; i<m; i++)
			{
				if(i<m-1)
				{
					alpha[m*o+i] = beta[i] + o * (m-1);
					sum += alpha[m*o+i];
				}
				else if( o==0)
				{
					alpha[m*o+i] = 1 - sum;
				}
				else  //i = m-1 && o>0
				{
					alpha[m*o+i] = -sum;
				}
			}
		}
		
		//TODO: draw area using alpha, and mu 
		
		return alpha; 
	}
	
	private static double SpatioTemporalDensity(int N, double[] pointValues, double pointTime, double[][] spatialBWidth, double[][] temporalBWidth , ArrayList<DataPoint> S)
	{
		//TODO: k berechnen wie?
		double k = 0;		
		return k;
	}
	
	public static double KDE(double[] pointValues, double[][] spatialBWidth, double[][] temporalBWidth , ArrayList<DataPoint> S)
	{
		double k = 0;
		
		for(int s=0; s<S.size(); s++)
		{
			k += GaussianSpatialDensityKernel(pointValues, S.get(s).values, spatialBWidth);
		}
		
		k = k/S.size();   //kernel density estimation, see aggarwal
		
		return k;
	}		
			
	public static double GaussianSpatialDensityKernel(double[] xValues, double[] muValues, double[][] spatialBWidth)
	{
		//TODO: wie wird die bandbreite angegeben - SUMi - was genau bedeutet sie?
		double summand = Math.pow((2*Math.PI), - (double)DataPoint.DIMENSIONS/2)* Math.pow(spatialBWidth.length, -0.5)* Math.exp(-0.5 * XSquare(xValues, muValues)); 
		
		return summand;
	}
	
	private static double XSquare(double[] xValues, double[] muValues)
	{
		//(x - mui)*(x - mui)
		double[] dif = new double[xValues.length];
		double sum = 0;
		
		for(int i=0; i<xValues.length; i++)
		{
			dif[i] = xValues[i] - muValues[i];
			sum += dif[i]*dif[i];
		}
		
		return sum;		
	}
	
	private static double[] Regress(double[][] P, double[] p, int m, int O)
	{
		
		//double[] beta = new double[(m-1)*(O+1)];
		
		OLSMultipleLinearRegression regression2 = new OLSMultipleLinearRegression();
	    regression2.setNoIntercept(true);
	    
	 /*   double[] a = {
	            2.6,
	            1.6,
	            4.0,
	            3.0,
	            4.9
	    };
	    double[][] b =
	            {
	                    { 1, 1.2 },
	                    { 1, 3.0  },
	                    { 1, 4.5  },
	                    { 1, 5.8  },
	                    { 1, 7.2  },
	            };
	    */
	    
	    regression2.newSampleData(p, P);

	    double[] regressionParameters = regression2.estimateRegressionParameters();

	    for (int i = 0; i < regressionParameters.length; i++) {
	        double regressionParameter = regressionParameters[i];
	        System.out.println("beta[" + i  + "] = " + regressionParameter);
	    }
	    
	    if(regressionParameters.length != (m-1)*O+1)
	    {
	    	System.out.println("wrong length of beta");
	    }
	    
		return regressionParameters;
	}

	
}
