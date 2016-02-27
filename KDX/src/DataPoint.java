
public class DataPoint {
	
	public static int DIMENSIONS = 1;
	public static int TIME_COLUMN = 1;
	
	public static double[] minValues;
	public static double[] maxValues;
	
	
	public double time;
	public double[] values;
	int indexFilled = 0;
	
	
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
		    
		}
		
	}
	
	public boolean AddData(double data)
	{
		if(indexFilled >= DIMENSIONS)
			return false;
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
		
		this.values[indexFilled] = data;
		indexFilled++;
		return true;
	}
	
	public boolean SetTime(double time)
	{
		if(time<0)
			return false;
		
		this.time = time;
		return true;
	}
}
