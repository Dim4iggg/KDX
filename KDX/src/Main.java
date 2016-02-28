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
	    Visualizer.AddPointSet(points, "Sample Points");
	    
	    //generate equally distributed pseudo points (eg 5 points)
	    ArrayList<DataPoint> Mu = GeneratePseudoPoints(5);
	    Visualizer.AddPointSet(Mu, "Pseudo Points");
	    
	    //chose F from S  (eg. every second point) -> here every point! F=S
	    //combine the locations of a subset of the historical instances in S with a set of different time points
	    ArrayList<DataPoint> F = GenerateFittingPositions(5);
	    Visualizer.AddPointSet(F, "Fitting Positions");
	    
	    //visualize data
	    try {
			//AnalysisLauncher.open(new Visualizer());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
	    
	    FitKDX(points, F, Mu, new double[10][10], new double[3][3] , 2);
	    
	    Visualizer app = new Visualizer(
                "KDX", Visualizer.PLOTTYPE.ALL);
        app.pack();
        app.setVisible(true);
        
	}


	private static void ReadData() throws BiffException, IOException
	{
		Workbook workbook = Workbook.getWorkbook(new File(dataPath)); 
		
		Sheet sheet =  workbook.getSheet(0);
		
		//read in data
		int row = 0;
		try
		{
		Cell cell = sheet.getCell(0,0);
		while(!cell.getContents().isEmpty())  //each row is a DataPoint
		{
			DataPoint point = new DataPoint();
			System.out.println("new data point");
			for(int column = 0; column<DataPoint.DIMENSIONS+1; column++)  //each row is a parameter
			{
				cell = sheet.getCell(column, row);
				CellType type = cell.getType();
				double zahl = Double.parseDouble(cell.getContents());
				
				if(column == DataPoint.TIME_COLUMN)
				{
					//set time
					point.SetTime(zahl);
					System.out.println("adding time: " + zahl);
				}
				else
				{ 
					//set next attribute value
					point.AddData(zahl);
					System.out.println("adding value:" + zahl);
				}
				
			}
			points.add(point);
			row++;
		}
		}
		catch(java.lang.ArrayIndexOutOfBoundsException e){};
		
		System.out.println(points.size() + " points found in " + dataPath);
		
	}
	
	private static ArrayList<DataPoint> GeneratePseudoPoints(int num)
	{
		ArrayList<DataPoint> mu = new ArrayList<>();
		for(int i=0; i<num; i++)
		{
			DataPoint pPoint = new DataPoint();
			pPoint.SetTime(0);  //TODO: change?
			
			String s = "new pseudo point: ";
			//distribute equally in each dimension
			for(int d=0; d<DataPoint.DIMENSIONS; d++)
			{
				double diff = DataPoint.maxValues[d] - DataPoint.minValues[d];  //eg 100
				double step = diff/(double)(num-1);  //if n=5  step = 25
				double pseudoValue = DataPoint.minValues[d] + step*i;  //points at 0, 25, 50, 75, 100
				pPoint.AddData(pseudoValue); 
				s += pseudoValue + "; ";
			}
			mu.add(pPoint);
			System.out.println(s);
			
		}
		return mu;
	}
	
	private static ArrayList<DataPoint> GenerateFittingPositions(int num)
	{
		//a set F of N !historical! (!IN THE PAST!), spatio-temporal fitting positions
		//combine the locations of a subset of the historical instances in S with a set of different time points.
		ArrayList<DataPoint> F = new ArrayList<>();
		
		if(num > points.size())
		{
			num = points.size(); //max as many as real points
		}
		for(int i=0; i<num; i++)
		{
			DataPoint pPoint = new DataPoint();
		
			String s = "new fitting position at time: ";
			//distribute equally in time
			double diff = DataPoint.maxTime - DataPoint.minTime;  //eg 100
			double step = diff/(double)(num-1);  //if n=5  step = 25
			double time = DataPoint.minTime + step*i;  //points at 0, 25, 50, 75, 100
			pPoint.SetTime(time); 
			s += time + "; ";
			
			pPoint.values = points.get(i).values;  //just copy the position of a real point
			
			F.add(pPoint);
			System.out.println(s);
			
		}
	    return F;
	}
	
	private static double[] FitKDX(ArrayList<DataPoint> S, ArrayList<DataPoint> F, ArrayList<DataPoint> Mu, double[][] spatialBWidth, double[][] temporalBWidth , int O)
	{
		int n = S.size();
		int N = F.size();
		int m = Mu.size();
		int M = (m-1)*(O+1);
		double[] k = new double[N];
		double[][] K = new double[N][m*(O+1)];
		
		//get historic density estimates
		for(int j = 0; j<N; j++)
		{
			 //k[j] =  KDE(F.get(j).values, spatialBWidth, temporalBWidth , S);
			//apply gaussian density function for each point
			k[j] = GaussianSpatialDensityKernel(F.get(j).values, F.get(j).values, spatialBWidth);
			System.out.println( "k[" + j + "] = " + k[j]);
			//add density line to visualization
			Visualizer.AddDensityLineSet(F.get(j),k[j], j+"s point", spatialBWidth);
			
			//for each pseudo point compute its density with a kernel 
			// as the sum of gaussian densities of neighboring sample points
			for(int i=0; i<m; i++)
			{
				double kernel = GaussianSpatialDensityKernel(F.get(j).values, Mu.get(i).values, spatialBWidth);
				
				
				for(int o=0; o<O; o++)
				{
					K[j][i+o*m] = kernel * Math.pow(F.get(j).time, o);  //Kj,i+o�m = Ki(xj) � to j
				}
			}
		}
		
		double[] p = new double[N];
		double[][] P = new double[N][(m-1) * (O+1)];
		
		//create design matrix P and corresponding vector p
		for(int j=0; j<N; j++)
		{
			p[j] = k[j] - K[j][m];
			for(int i=0; i< (m-1); i++)
			{
				for(int o=0; o<O; o++)
				{
					P[j][i+o*(m-1)] = (K[j][i+o*m] - K[j][m+o*m])* Math.pow(F.get(j).time, o);
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
		if(true)
			return k;
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
		//TODO: k berechnen wie?
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
		
		//double[] beta = new double[(m-1)*O+1];
		
		OLSMultipleLinearRegression regression2 = new OLSMultipleLinearRegression();
	    regression2.setNoIntercept(true);
	    regression2.newSampleData(p, P);

	    double[] regressionParameters = regression2.estimateRegressionParameters();

	    for (int i = 0; i < regressionParameters.length; i++) {
	        double regressionParameter = regressionParameters[i];
	        System.out.println(i + " " + regressionParameter);
	    }
	    
	    if(regressionParameters.length != (m-1)*O+1)
	    {
	    	System.out.println("wrong length of beta");
	    }
	    
		return regressionParameters;
	}

	
}
