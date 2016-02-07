import java.util.ArrayList;
import org.apache.commons.math3.stat.regression.*;

public class Main {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		OLSMultipleLinearRegression regression2 = new OLSMultipleLinearRegression();


		 double[] y = {
		            2.6,
		            1.6,
		            4.0,
		            3.0,
		            4.9
		    };
		    double[][] x2 =
		            {
		                    { 1, 1.2 },
		                    { 1, 3.0  },
		                    { 1, 4.5  },
		                    { 1, 5.8  },
		                    { 1, 7.2  },
		            };

	    regression2.setNoIntercept(true);
	    regression2.newSampleData(y, x2);

	    double[] regressionParameters = regression2.estimateRegressionParameters();

	    for (int i = 0; i < regressionParameters.length; i++) {
	        double regressionParameter = regressionParameters[i];
	        System.out.println(i + " " + regressionParameter);
	    }

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
			 k[j] =  SpatioTemporalDensity(N, F.get(j).values, F.get(j).time, spatialBWidth, temporalBWidth , S);
			
			for(int i=0; i<m; i++)
			{
				 GaussianSpatialDensityKernel(F.get(j).values, Mu.get(i).values, spatialBWidth, S);
				 
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
					P[j][i+o*(m-1)] = (K[j][i] - K[j][m])* Math.pow(F.get(j).time, o);
				}
			}
		}
		
		//add regularisation terms to p and P
		for(int i=0; i<(m-1); i++)
		{
			for(int o=0; o<O; o++)
			{
				p[N+i+o*(m-1)] = o; //TODO: Cio was ist das?
				P[N+i+o*(m-1)][i+o*(m-1)] = 0; //TODO: change! siehe paper
				P[N+i+o*(m-1)][i+o*(m-1)] = 1;
			}
		}
		
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
		return 0;
	}
	
	private static double GaussianSpatialDensityKernel(double[] xValues, double[] muValues, double[][] spatialBWidth, ArrayList<DataPoint> S)
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
		//TODO: beta regression berechnen - wie?
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
