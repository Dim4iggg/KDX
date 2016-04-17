import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.math3.linear.QRDecomposition;
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

public class KDX {

	private static ArrayList<DataPoint> points;
	private static ArrayList<DataPoint> Mu;
	private static ArrayList<DataPoint> F;
	
	private static String dataPath = "Eingabewerte2.xls";
	
	//this many pseudo points will be added at FIXED positions in EVERY time -> num points * num times
	private static int NUM_PSEUDOPOINTS = 10;
	private static int NUM_FITTINGPOS = 6;
	private static int O = 5;
			
	private static double regWeight = 1;
	
	private static double[] alpha;
	private static double[][] spatialBWidth;
	
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
	    
	  //  DO(0,0);
	    
	    //identity matrix
	    spatialBWidth = new double[DataPoint.DIMENSIONS][DataPoint.DIMENSIONS];//TODO: what to put in this matrix?
	    for(int i=0; i< spatialBWidth.length; i++)
	    {
	    	 for(int j=0; j< spatialBWidth[0].length; j++)
	    	 {
	    		 if(i==j)
	    			 spatialBWidth[i][j] = 1;
	    	 }
	    }
	    
	    
	    //generate equally distributed pseudo points (eg 5 points)
	    Mu = GeneratePseudoPoints();
	    //Visualizer.AddPointSet(Mu, "Pseudo Points");
	    
	    //chose F from S  (eg. every second point) -> here every point! F=S
	    //combine the locations of a subset of the historical instances in S with a set of different time points
	    F = GenerateFittingPositions();
	    F = points;
	    
	    //visualize data
	    try {
			//AnalysisLauncher.open(new Visualizer());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
	    
	    
	    FitKDX(F, Mu, spatialBWidth, new double[10][10]);
	    
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
					point.AddData(zahl, true);
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
		//copy the points at the same position for each time step (in the past)
		//for(int t=0; t<DataPoint.timePoints.size(); t++)
		{
			for(int i=0; i<NUM_PSEUDOPOINTS; i++)    
			{
				DataPoint pPoint = new DataPoint();
				//pPoint.SetTime(DataPoint.timePoints.get(1)); 
				pPoint.SetTime(DataPoint.timePoints.get((int)(Math.random()*DataPoint.timePoints.size()))); 
				
				
				//distribute equally in each dimension
				for(int d=0; d<DataPoint.DIMENSIONS; d++)
				{
					double diff = DataPoint.maxValues[d] - DataPoint.minValues[d] ;  //eg 100
					double step = diff/(double)(NUM_PSEUDOPOINTS-1);  //if n=5  step = 25
					double pseudoValue = DataPoint.minValues[d] + step*i ;  //points at 0, 25, 50, 75, 100
					pPoint.AddData(pseudoValue,false ); 
					s += pseudoValue + ", t=" + pPoint.time + "; ";
				}
				mu.add(pPoint);
				
				//add pseudo point Mu(i) to visualization
				Visualizer.AddDensityLinePoint(pPoint, GaussianSpatialDensityKernel(pPoint.values, pPoint.values, spatialBWidth), "Mu"+i, spatialBWidth,1);
			}
		}
		System.out.println(s);
		return mu;
	}
	
	private static ArrayList<DataPoint> GenerateFittingPositions()
	{
		//a set F of N !historical! (!IN THE PAST!), spatio-temporal fitting positions
		//combine the locations of a subset of the historical instances in S with a set of different time points.
		ArrayList<DataPoint> F = new ArrayList<>();
		
		String s = "fitting positions: ";
		for(int t = 0; t< DataPoint.timePoints.size(); t++)
		{
			for(int i=0; i<NUM_FITTINGPOS; i++)
			{ 
				double time = DataPoint.timePoints.get(t);  
				
				DataPoint pPoint = new DataPoint();
				pPoint.SetTime(time); 
				
				//distribute equally in each dimension
				for(int d=0; d<DataPoint.DIMENSIONS; d++)
				{
					double diff = DataPoint.maxValues[d] - DataPoint.minValues[d] +2;  //eg 100
					double step = diff/(double)(NUM_FITTINGPOS-1);  //if n=5  step = 25
					double pseudoValue = DataPoint.minValues[d]-1 + step*i;  //points at 0, 25, 50, 75, 100
					pPoint.AddData(pseudoValue, false); 
					s += pseudoValue + ", t=" + pPoint.time + "; ";
				}				
				F.add(pPoint);
			}
		}
		System.out.println(s);
	    return F;
	}
	
	private static double[] FitKDX( ArrayList<DataPoint> F, ArrayList<DataPoint> Mu, double[][] spatialBWidth, double[][] temporalBWidth)
	{
		//int n = S.size();
		//System.out.println("\n|S| = " + n);
		
		int N = F.size();
		System.out.println("|F| = N = " + N);
		
		int m = Mu.size();
		System.out.println("|Mu| = m = " + m + "\n");
		
		double[] k = new double[N];
		double[][] K = new double[N][m];   
		
		
		String kMatrix = "\n K: \n";
		//get historic density estimates 
		for(int j = 0; j<N; j++)  //for all fitting positions F
		{
			//k is a (N * 1)-vector 
			k[j] =  GaussianSpatialDensityKernel( F.get(j).values, F.get(j).values, spatialBWidth);
			Visualizer.AddDensityLineSet(F.get(j),GaussianSpatialDensityKernel(F.get(j).values, F.get(j).values, spatialBWidth), "", spatialBWidth, 0);
			Visualizer.AddDensityLinePoint(F.get(j), GaussianSpatialDensityKernel(F.get(j).values, F.get(j).values, spatialBWidth), "F"+j, spatialBWidth,0);
			F.get(j).SetDensity(k[j]);

			System.out.println( "k[" + j + "] = " + k[j]);
		
			//for each pseudo point compute its density with a kernel 
			// as the sum of gaussian densities of neighboring sample points
			
			for(int i=0; i<m; i++)
			{
				//double kernel = GaussianSpatialDensityKernel(F.get(j).values, Mu.get(i).values, spatialBWidth);
				double kernel = GaussianSpatialDensityKernel( F.get(j).values, Mu.get(i).values, spatialBWidth); //>< if time same
				//spatial density of m-th pseudo-point's kernel at f-th fitting position
				K[j][i] = kernel * k[j];

				kMatrix = kMatrix + round(K[j][i], 3) + " | ";
				
				
			}
			
			kMatrix = kMatrix + "\n";
		}
		
		for(int i=0; i<m; i++)
		{
			double muDen = 0;
			for(int j = 0; j<N; j++)  //for all fitting positions F
			{
				muDen += K[j][i];
			}
			Mu.get(i).SetDensity(muDen);
		}
		
		String h ="";
		for(int z=0; z<K.length; z++)
		{
			for(int w=0; w< K[0].length; w++)
			{
				h += round(K[z][w], 3) + " | ";
			}
			h+="\n";
		}
		System.out.println("K = \n" + h);	
			
		double[] p = new double[N];
		double[][] P = new double[N][(m-1) * (O+1)];
		
		String pMatrix = "\n P: \n";
		//create design matrix P and corresponding vector p
		for(int j=0; j<N; j++)
		{
			p[j] = k[j] - K[j][m-1]; 
			System.out.println("p[j] = " + p[j]);
			
			for(int o=0; o<=O; o++)
			{
				for(int i=0; i< m-1; i++)  
				{
					double q = (K[j][i] - K[j][m-1])* Math.pow(F.get(j).time, o);
					P[j][i+o*(m-1)] = q;
					//if(P[j][i+o*(m-1)] == 0)
					//	P[j][i+o*(m-1)] += 0.0000001; 

					pMatrix = pMatrix + round(P[j][i+o*(m-1)], 3) + " | ";
				}
			}
			pMatrix = pMatrix + "\n";
		}
		//System.out.println(pMatrix);
		
		//add regularisation terms to p and P
		/*for(int i=0; i<(m-1); i++)
		{
			//for(int j=0; j<N; j++)
			{
				for(int o=0; o<=O; o++)
				{
					//p[(N/DataPoint.timePoints.size()-1)+i+o*(m-1)] = 3; //TODO: Cio was ist das?  p oder P?
					if( (i+o*(m-1))% ((N/DataPoint.timePoints.size() -1)+i+o*(m-1)) == 0) //diagonal
					//if( j!=0 && (j/(N)) == i )
					{
						//P[(N/DataPoint.timePoints.size() -1)+i+o*(m-1)][i+o*(m-1)] = 1; //set diagonal element to 1
						//P[j][i+o*(m-1)] = 5;
					}
					else
					{
						//P[(N/DataPoint.timePoints.size()-1)+i+o*(m-1)][i+o*(m-1)] = 0; //set row to zero
						//P[j][i+o*(m-1)] = 8;
					}
				}	
			}
		} */
		
		String g ="";
		for(int z=0; z<P.length; z++)
		{
			for(int w=0; w< P[0].length; w++)
			{
				g += round(P[z][w], 3) + " | ";
			}
			g+="\n";
		}
		System.out.println("P = \n" + g);
		
		for(int z=0; z<p.length; z++)
		{
			System.out.println("p[]" + p[z]);
		}
		//TODO: REMOVE! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		//if(true)
		//	return k;
//		<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		
		
		//fit regression coefficients
		double[] beta = Regress(P, p, O);

		alpha = new double[m*(O+1)]; 
		
		//reconstruct regression coefficient of alpha
		for(int o=0; o<=O; o++)
		{
			double sum =0;
			for(int i=0; i<m; i++)
			{
				if(i == m-1)
				{
					if(o==0 )
					{
						alpha[i] = 1 - sum; 
						sum = 0;
					}
					else //o>0
					{
						alpha[m*o+i] = - sum;
						sum = 0;
					}
				}
				else 
				{
					alpha[m*o+i] = beta[i+o*(m-1)];
					sum += alpha[m*o+i];
				}
				
				System.out.print("\n alpha[m*o+i] = " + alpha[m*o+i]);
			}
		}
		return alpha; 
	}

	public static double KDE(DataPoint f, double[][] spatialBWidth, ArrayList<DataPoint> S)
	{
		double[] pointValues = f.values;
		double time = f.time;
		double k = 0;
		
		for(int s=0; s<S.size(); s++)
		{
			if(Math.abs(time - S.get(s).time)< 0.5)
			{
				double value = GaussianSpatialDensityKernel(S.get(s).values, pointValues, spatialBWidth);
				k = k+value;
			}
		}
		
		//k = k/S.size();   //kernel density estimation, see aggarwal
		
		return k;
	}		
			
	public static double GaussianSpatialDensityKernel(double[] xValues, double[] muValues, double[][] spatialBWidth)
	{
		//TODO: wie wird die bandbreite angegeben - SUMi - was genau bedeutet sie?
		double summand = Math.pow((2*Math.PI), - (double)DataPoint.DIMENSIONS/2)* Math.pow(spatialBWidth.length, -0.5)* Math.exp(-0.5 * XSquare(xValues, muValues, spatialBWidth)); 
		
		return summand;
	}
	
	private static double XSquare(double[] xValues, double[] muValues, double[][] spatialBWidth)
	{
		//(x - mui)*(x - mui)
		double[] dif = new double[xValues.length];
		double sum = 0;

	
        
		for(int i=0; i<xValues.length; i++)
		{
			double s = 0;
			
			dif[i] = xValues[i] - muValues[i];
			for(int j=0; j < spatialBWidth[0].length; j++ )
			{
				s += spatialBWidth[i][j];
			}
			sum += dif[i]* s * dif[i];
		}
		
		return sum;		
	}
	
	private static double[] Regress(double[][] P, double[] p,  int O)
	{
		
		//double[] beta = new double[(m-1)*(O+1)];
		
		OLSMultipleLinearRegression regression2 = new OLSMultipleLinearRegression();
	    regression2.setNoIntercept(false);
	   
	    MLRegression regression = new MLRegression(P, p);
	    
	   MillerUpdatingRegression regression3 = new MillerUpdatingRegression(P[0].length, false, 0.01);
	    regression3.addObservations(P, p);
	    
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
	    //regression3.addObservations(P, p);
	    
	    double[] regressionParameters3 = new double[P[0].length];
	    for(int i =0 ; i<P[0].length; i++)
	    {
	    	regressionParameters3[i] = regression.beta.get(i,0);
	    }
	    
	    double[] regressionParameters = regression2.estimateRegressionParameters();

	    for (int i = 0; i < regressionParameters3.length; i++) {
	        double regressionParameter = regressionParameters3[i];
	        System.out.println("beta[" + i  + "] = " + regressionParameter);
	    }
    
	    check(P, p, regressionParameters3);
	    
		return regressionParameters3;
	}

	
	private static void check(double[][] P, double[] p, double[] b)
	{
		for(int i=0; i< b.length; i++)
		{
			double summand = 0;
			for(int j=0; j< P[0].length; j++)
			{
				summand = summand + P[i][j]*b[i];
			}
			
			double abw = p[i] - summand;
			System.out.println("abweichung ["+i+"] = " +abw);
		}
	}
	
	public static double Extrapolation(double[] xValues, double time)
	{
		O = alpha.length/Mu.size() -1;
		
		double f = 0;		
		for(int i=0; i<Mu.size(); i++)
		{
			double mult = 0;
			//if( Math.abs((Mu.get(i).time-time))<= 1)
			//{				
				for(int l=0; l<O; l++)
				{					
					double summand = alpha[Mu.size()*l+i] * Math.pow(time, l);
					//wegen zweiten summenzeichen
					mult = mult+summand;
				}			
			//wegen ersten summenzeichen
			f = f + (mult * GaussianSpatialDensityKernel(xValues, Mu.get(i).values, spatialBWidth));	
			//}
		}

		return f;
	}

	public static double round(double value, int places) {
	    if (places < 0) throw new IllegalArgumentException();

	    long factor = (long) Math.pow(10, places);
	    value = value * factor;
	    long tmp = Math.round(value);
	    return (double) tmp / factor;
	}

	public static double DO(double x, double z)
	{
		 double[] k = {
		            1,
		           2,
		           3
		    };
		    double[][] K =
		            {
		                    { 0.33,   0.66  },
		                    { 0.5,  0.5  },
		                    { 0.66,  0.33  }
		            };
		    
		    double[] a = Regress(K, k, 3);

		 
		check(K, k, a);   
		    
		return 1;
	}

	public static double wtemporal(double delta_t)
	{
		double ht = 3;	
		return Math.max(1/ht - Math.abs(delta_t)/Math.pow(ht,2), 0);
	}
	
	public static double Extra(double[] x, double z)
	{
		O = alpha.length/Mu.size() -1;
		
		double ht = 1;
		double delta_t = 3;
		//double wtemporal = (delta_t)* Math.max(1/ht - Math.abs(delta_t)/Math.pow(ht,2),0); 
		
		// Q *M Matrix (only 1 querry point at once)
		double[][] wt2 = new double[1][Mu.size()];
		
		 for( int i = 0; i<wt2.length; i++)
		 {	
			for(int j=0; j< wt2[0].length; j++)
			{
				//wt = wtemporal(repmat(RefT',length(QueryT),1)-repmat(QueryT,1,length(RefT))); 
				double w = Mu.get(j).time - z;
				double wFunc = wtemporal(w);
				wt2[i][j] = wFunc;
				for(int o = 0; o<O; o++)
				{
									
					wt2[i][j] = wt2[i][j] + Math.pow(z, o)*alpha[Mu.size()*o+j];	
				}
				if(wt2[i][j] < 0)
				{
					wt2[i][j] = 0; //non-negativity correction
				}
			}
		 }
	 
		 double sum = 0;
		//row-sum equals one correction
		 for(int j=0; j< wt2[0].length; j++)
		 {
			//containing the sum of each row.
			 
			sum += wt2[0][j];	 
		 }
		 
		 for(int j=0; j< wt2[0].length; j++)
		 {
			 wt2[0][j] = wt2[0][j]/sum;
		 }
		 //calculate spatial densities
		 //TODO: take Mus density into account?
		 double[][] ds = new double [1][Mu.size()];
		 for(int i = 0; i < Mu.size(); i++)
		 {
			 double kernel = GaussianSpatialDensityKernel(Mu.get(i).values, x, spatialBWidth); //>< if time same
			 ds [0][i] = kernel * Mu.get(i).density;
		 }
		 
		 //combine into spatio-temporal extrapolated densities:

		 //.* is elementwise multiplication
		 //ds.*wtnorm
		 //(Q x 1) 
		 double dst = 0;
		 for(int i = 0; i < Mu.size(); i++)
		 {
			// dst += ds[0][i]*Mu.get(i).values[0];  //dst = sum(ds.*wtnorm,2);
			 dst += ds[0][i]* wt2[0][i];
					 
		 }
		 
		 System.out.println(dst);
		 int g;
		 if(2-x[0] < 0.05)
			 g =0;
		 return dst;
	}
	
}











