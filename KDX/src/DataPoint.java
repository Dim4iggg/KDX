import java.util.ArrayList;



public class DataPoint {
	
	//how many dimensions does the data have
	public static int DIMENSIONS = 1;
	
	//in which column in the data is the time 
	public static int TIME_COLUMN = 1;
	
	
	
	//min and max values for each dimension (ALL data)
	public static double[] minValues;
	public static double[] maxValues;
	 
	//first and last time points (of ALL data)
	public static double minTime = Double.MAX_VALUE;
	public static double maxTime = Double.MIN_VALUE;
	
	//time points found through ALL data
	public static ArrayList<Double> timePoints;
	
	//time of THIS DataPoint
	public double time;
	
	//values of THIS DataPoint
	public double[] values;
	
	//for how many dimensions the values were read in 
	int indexFilled = 0;
	
	//computed density
	public double density;
	
	
	public DataPoint()
	{
		this.values = new double[DIMENSIONS];
		if(minValues == null)
		{
			//min and max values are saved for each dimension in order to compute the bandwidths
		    minValues = new double[DataPoint.DIMENSIONS];
		    maxValues = new double[DataPoint.DIMENSIONS];
		    for(int i = 0; i<DIMENSIONS; i++)
		    {
		    	minValues[i] = Double.MAX_VALUE;
		    	maxValues[i] = Double.MIN_VALUE;
		    }
		    //all time points are saved to contribute the pseudo points on them
		    timePoints = new ArrayList<Double>();
		}
		
	}
	
	public boolean AddData(double data, boolean updateValBoundaries)
	{
		if(indexFilled >= DIMENSIONS)
			return false;
		
		if(updateValBoundaries)
		{
			//new min value
			if(data < minValues[indexFilled])
			{
				minValues[indexFilled] = data;
			}
			//new max value
			if(data > maxValues[indexFilled])
			{
				maxValues[indexFilled] = data;
			}
		}
		
		this.values[indexFilled] = data;
		indexFilled++;
		return true;
	}
	
	public boolean SetTime(double time)
	{
		if(time<0)
			return false;
		
		if(time < DataPoint.minTime)
		{
			DataPoint.minTime = time;
		}
		if(time > DataPoint.maxTime)
		{
			DataPoint.maxTime = time;
		}
		
		//is this time point already saved?
		boolean contains = false;
		for(int i=0; i<timePoints.size(); i++)
		{
			if(timePoints.get(i) == time)
			{
				contains = true;
				break;
			}
		}
		//if not, save it
		if(!contains)
		{
			timePoints.add(time);
		}
		
		this.time = time;
		return true;
	}
	
	public void SetDensity(double den)
	{
		this.density = den;
	}
}
