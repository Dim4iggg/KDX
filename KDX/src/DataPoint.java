
public class DataPoint {
	
	public static int DIMENSIONS = 1;
	
	public double time;
	public double[] values;
	
	public DataPoint(double[] newValues, double timePoint)
	{
		this.time = timePoint;
		
		if(newValues.length != DIMENSIONS)
		{
			//create an empty array
			this.values = new double[DIMENSIONS];
			System.out.println("Wrong number of values!");
			return;
		}
		//set values
		this.values = newValues;
	}
}
