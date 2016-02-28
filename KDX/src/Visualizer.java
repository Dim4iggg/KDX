import java.awt.BorderLayout;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Random;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;

import com.orsoncharts.Chart3D;
import com.orsoncharts.Chart3DFactory;
import com.orsoncharts.Chart3DPanel;
import com.orsoncharts.Colors;
import com.orsoncharts.axis.NumberAxis3D;
import com.orsoncharts.axis.NumberTickSelector;
import com.orsoncharts.data.Dataset3D;
import com.orsoncharts.data.DefaultKeyedValues;
import com.orsoncharts.data.KeyedValues;
import com.orsoncharts.data.category.CategoryDataset3D;
import com.orsoncharts.data.category.StandardCategoryDataset3D;
import com.orsoncharts.data.xyz.XYZDataset;
import com.orsoncharts.data.xyz.XYZSeries;
import com.orsoncharts.data.xyz.XYZSeriesCollection;
import com.orsoncharts.demo.LineChart3D1;
import com.orsoncharts.demo.ScatterPlot3D1;
import com.orsoncharts.demo.swing.DemoPanel;
import com.orsoncharts.demo.swing.ExitOnClose;
import com.orsoncharts.demo.swing.OrsonChartsDemo;
import com.orsoncharts.demo.swing.ScatterPlot3DDemo1;
import com.orsoncharts.demo.swing.StackedBarChart3DDemo3;
import com.orsoncharts.graphics3d.Dimension3D;
import com.orsoncharts.graphics3d.ViewPoint3D;
import com.orsoncharts.graphics3d.swing.DisplayPanel3D;
import com.orsoncharts.label.StandardXYZLabelGenerator;
import com.orsoncharts.plot.CategoryPlot3D;
import com.orsoncharts.plot.XYZPlot;
import com.orsoncharts.renderer.category.CategoryRenderer3D;
import com.orsoncharts.renderer.category.StackedBarRenderer3D;
import com.orsoncharts.renderer.category.StandardCategoryColorSource;
import com.orsoncharts.renderer.xyz.ScatterXYZRenderer;




public class Visualizer  extends JFrame
{
	public enum PLOTTYPE
	{
		LINE, AREA, SCATTER, ALL;
	}
	
	static XYZSeriesCollection scatterdataset = null;
	static CategoryDataset3D linedataset = null;
	
    public Visualizer(String title, PLOTTYPE ptype) {
        super(title);
        addWindowListener(new ExitOnClose());
        getContentPane().add(createDemoPanel(ptype));
    }

    
    /**
     * Returns a panel containing the content for the demo.  This method is
     * used across all the individual demo applications to allow aggregation 
     * into a single "umbrella" demo (OrsonChartsDemo).
     * 
     * @return A panel containing the content for the demo.
     */
    public static JPanel createDemoPanel(PLOTTYPE ptype) {
    	DemoPanel mainPanel = new DemoPanel(new BorderLayout());
    	mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));
    	mainPanel.setPreferredSize(OrsonChartsDemo.DEFAULT_CONTENT_SIZE);
        
    	//prepare scatter plot
        DemoPanel scatterPanel = new DemoPanel(new BorderLayout());
        scatterPanel.setPreferredSize(OrsonChartsDemo.DEFAULT_CONTENT_SIZE);
        Chart3D scatterChart = Visualizer.createScatterChart(scatterdataset);
    	Chart3DPanel scatterChartPanel = new Chart3DPanel(scatterChart);
    	scatterPanel.setChartPanel(scatterChartPanel);
    	scatterChartPanel.zoomToFit(OrsonChartsDemo.DEFAULT_CONTENT_SIZE);           
        scatterPanel.add(new DisplayPanel3D(scatterChartPanel ,false, false));
               
        //prepare line plot
        DemoPanel linePanel = new DemoPanel(new BorderLayout());
        linePanel.setPreferredSize(OrsonChartsDemo.DEFAULT_CONTENT_SIZE);
        Chart3D lineChart = Visualizer.createLineChart(linedataset); 
        Chart3DPanel lineChartPanel = new Chart3DPanel(lineChart);
        linePanel.setChartPanel(lineChartPanel);
        lineChartPanel.zoomToFit(OrsonChartsDemo.DEFAULT_CONTENT_SIZE);
        linePanel.add(new DisplayPanel3D(lineChartPanel ,false, false));
        
        //add panels 
        switch (ptype) 
        {
        
        case LINE:   	
        	mainPanel.add(linePanel);
        	break;
           
        case AREA:
        	break;
        	
        case SCATTER:
        	mainPanel.add(scatterPanel);
            break;
            
        case ALL:

        	mainPanel.add(scatterPanel);
            mainPanel.add(linePanel);
            break;   
        }

        return mainPanel;
    }
    
    public static void AddPointSet(ArrayList<DataPoint> points, String title)
    {
    	if(scatterdataset == null)
    	{
    		scatterdataset = new XYZSeriesCollection();
    	}
    	
    	 XYZSeries s = new XYZSeries(title);
    	 DataPoint p;
         for (int i = 0; i < points.size(); i++) 
         {
        	 p = points.get(i);
             s.add(p.values[0], 1, p.time);
         }
         scatterdataset.add(s);
    }

    public static void AddDensityLineSet(DataPoint point, double density, String title, double[][] spatialBWidth,int numRefs )
    {
    	if(linedataset == null)
    	{
    		linedataset = new StandardCategoryDataset3D();
    	}
    	
    	
    	DefaultKeyedValues<Double> series = new DefaultKeyedValues<Double>();
    	
    	if(numRefs <=0)
    		numRefs = 10;
    	else if(numRefs == 1) //point
    	{
        	series.put(point.values[0], density);
    		((StandardCategoryDataset3D) linedataset).addSeriesAsRow(title, point.time, series);   
    		return;
    	}

        // 5 points "left" from the sample point
        for(int i= numRefs-1;  i>0; i--)
        {
        	//copy all values
        	double[] lpoint =  new double[point.values.length];
        	for(int j =0; j<lpoint.length; j++)
        	{
        		lpoint[j] = point.values[j];
        	}
        	
        	//move in 1st dimension left //TODO: change to the dimension that is shown
        	double xPos = lpoint[0] - 0.5*(i+1);
        	lpoint[0] = xPos;
        	double val = Main.GaussianSpatialDensityKernel(point.values, lpoint, spatialBWidth);
        	series.put(xPos, val);
        } 
        
        series.put(point.values[0], density);
        
        // 5 points "right" from the sample point
        for(int i= 0; i<numRefs; i++)
        {
        	//copy all values
        	double[] rpoint =  new double[point.values.length];
        	for(int j =0; j<rpoint.length; j++)
        	{
        		rpoint[j] = point.values[j];
        	}
        	
        	//move in 1st dimension right //TODO: change to the dimension that is shown
        	double xPos = rpoint[0] + 0.5*(i+1);
        	rpoint[0] = xPos;
        	double val = Main.GaussianSpatialDensityKernel(point.values, rpoint, spatialBWidth);
        	series.put(xPos, val);
        }
        
        ((StandardCategoryDataset3D) linedataset).addSeriesAsRow(title, point.time, series);   
    }

    public static void AddDensityLinePoint(DataPoint point, double density, String title, double[][] spatialBWidth,int numRefs )
    {
    	if(scatterdataset == null)
    	{
    		scatterdataset = new XYZSeriesCollection();
    	}
    	 XYZSeries s = new XYZSeries(title);
    	 
    	 if(numRefs <=0)
     		numRefs = 50;
     	else if(numRefs == 1) //point
     	{
     		s.add(point.values[0], density, point.time);
     		 scatterdataset.add(s);   
     		return;
     	}
    	 
    	// 5 points "left" from the sample point
         for(int i= numRefs-1;  i>0; i--)
         {
         	//copy all values
         	double[] lpoint =  new double[point.values.length];
         	for(int j =0; j<lpoint.length; j++)
         	{
         		lpoint[j] = point.values[j];
         	}
         	
         	//move in 1st dimension left //TODO: change to the dimension that is shown
         	double xPos = lpoint[0] - 0.05*(i+1);
         	lpoint[0] = xPos;
         	double val = Main.GaussianSpatialDensityKernel(point.values, lpoint, spatialBWidth);
         	s.add(xPos, val, point.time);
         } 
         
         s.add(point.values[0], density, point.time);
         
         // 5 points "right" from the sample point
         for(int i= 0; i<numRefs; i++)
         {
         	//copy all values
         	double[] rpoint =  new double[point.values.length];
         	for(int j =0; j<rpoint.length; j++)
         	{
         		rpoint[j] = point.values[j];
         	}
         	
         	//move in 1st dimension right //TODO: change to the dimension that is shown
         	double xPos = rpoint[0] + 0.05*(i+1);
         	rpoint[0] = xPos;
         	double val = Main.GaussianSpatialDensityKernel(point.values, rpoint, spatialBWidth);
         	s.add(xPos, val, point.time);
         }
         scatterdataset.add(s);
    }
    
    
    /**
     * Creates a scatter chart based on the supplied dataset.
     * 
     * @param dataset  the dataset.
     * 
     * @return A scatter chart. 
     */
    public static Chart3D createScatterChart(XYZDataset dataset) {
        Chart3D chart = Chart3DFactory.createScatterChart("Data Points", 
                null, dataset, "X", "Density", "Time");
        XYZPlot plot = (XYZPlot) chart.getPlot();
        //plot.setDimensions(new Dimension3D(10.0, 4.0, 4.0));
        plot.setLegendLabelGenerator(new StandardXYZLabelGenerator(
                StandardXYZLabelGenerator.COUNT_TEMPLATE));
        plot.getZAxis().setInverted(true);
        ScatterXYZRenderer renderer = (ScatterXYZRenderer) plot.getRenderer();
        renderer.setSize(0.15);
        renderer.setColors(Colors.createIntenseColors());
        chart.setViewPoint(ViewPoint3D.createAboveLeftViewPoint(30));
        return chart;
    }
    
    /**
     * Creates a line chart with the supplied dataset.
     * 
     * @param dataset  the dataset.
     * 
     * @return A line chart.
     */
    public static Chart3D createLineChart(CategoryDataset3D dataset) {
        Chart3D chart = Chart3DFactory.createLineChart(
                "Density functions of Fitting Positions",null, dataset, "Time", null,  
                "Density");
        CategoryPlot3D plot = (CategoryPlot3D) chart.getPlot();
       // plot.setDimensions(new Dimension3D(18, 8, 4));
        plot.getRowAxis().setVisible(true);
        plot.getColumnAxis().setVisible(true);

        
        NumberAxis3D valueAxis = (NumberAxis3D) plot.getValueAxis();
        valueAxis.setTickSelector(new NumberTickSelector(true));
        CategoryRenderer3D renderer = plot.getRenderer();
        renderer.setColors(Colors.createFancyDarkColors());
       
        chart.setViewPoint(ViewPoint3D.createAboveLeftViewPoint(40));
        return chart;    
    }
}
