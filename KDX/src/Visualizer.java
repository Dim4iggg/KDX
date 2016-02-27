import java.awt.BorderLayout;
import java.awt.Color;
import java.util.Random;

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
import com.orsoncharts.renderer.category.StackedBarRenderer3D;
import com.orsoncharts.renderer.category.StandardCategoryColorSource;
import com.orsoncharts.renderer.xyz.ScatterXYZRenderer;




public class Visualizer  extends JFrame
{
	public enum PLOTTYPE
	{
		LINE, AREA, SCATTER;
	}
	
    public Visualizer(String title, PLOTTYPE ptype) {
        super(title);
        addWindowListener(new ExitOnClose());
        getContentPane().add(createDemoPanel(ptype));
    }

    
    public static void main(String[] args) {
    	Visualizer app = new Visualizer(
                "MY demo", PLOTTYPE.LINE);
        app.pack();
        app.setVisible(true);
    }
    
    /**
     * Returns a panel containing the content for the demo.  This method is
     * used across all the individual demo applications to allow aggregation 
     * into a single "umbrella" demo (OrsonChartsDemo).
     * 
     * @return A panel containing the content for the demo.
     */
    public static JPanel createDemoPanel(PLOTTYPE ptype) {
        DemoPanel content = new DemoPanel(new BorderLayout());
        content.setPreferredSize(OrsonChartsDemo.DEFAULT_CONTENT_SIZE);
        
        Chart3D chart = null;
        
        switch (ptype) 
        {
        
        case LINE:   
        	CategoryDataset3D dataset = Visualizer.createLineDataset();
            chart = Visualizer.createLineChart(dataset);
            break;
        case AREA:
        	
        case SCATTER:
        	XYZDataset dataset2 = Visualizer.createScatterDataset();
        	chart = Visualizer.createScatterChart(dataset2);
        	break;
  	
        }

        Chart3DPanel chartPanel = new Chart3DPanel(chart);
        content.setChartPanel(chartPanel);
        chartPanel.zoomToFit(OrsonChartsDemo.DEFAULT_CONTENT_SIZE);
        content.add(new DisplayPanel3D(chartPanel ,false, false));
        return content;
    }

    public static XYZDataset createScatterDataset() {
        XYZSeries s1 = createRandomSeries("S1", 15);
        XYZSeries s2 = createRandomSeries("S2", 50);
        XYZSeries s3 = createRandomSeries("S3", 150);
        XYZSeriesCollection dataset = new XYZSeriesCollection();
        dataset.add(s1);
        dataset.add(s2);
        dataset.add(s3);
        return dataset;
    }
    
    private static XYZSeries createRandomSeries(String name, int count) {
        XYZSeries s = new XYZSeries(name);
        for (int i = 0; i < count; i++) {
            s.add(Math.random() * 100, Math.random() / 100, Math.random() * 100);
        }
        return s;
    }
    
    public static CategoryDataset3D createLineDataset() {
        StandardCategoryDataset3D dataset = new StandardCategoryDataset3D();
        dataset.addSeriesAsRow("Safari", createLineData());
        dataset.addSeriesAsRow("Firefox", createLineData());
        dataset.addSeriesAsRow("Internet Explorer", createLineData());
        dataset.addSeriesAsRow("Chrome", createLineData());
        return dataset;
    }
    
    private static KeyedValues<Double> createLineData() {
        DefaultKeyedValues<Double> series = new DefaultKeyedValues<Double>();
        series.put("Nov-12", Math.random());
        series.put("Dec-12", Math.random());
        series.put("Jan-13", Math.random());
        return series;
    }
    
    /**
     * Creates a scatter chart based on the supplied dataset.
     * 
     * @param dataset  the dataset.
     * 
     * @return A scatter chart. 
     */
    public static Chart3D createScatterChart(XYZDataset dataset) {
        Chart3D chart = Chart3DFactory.createScatterChart("MY plot demo", 
                "Here could be your ad!", dataset, "X", "Y", "Z");
        XYZPlot plot = (XYZPlot) chart.getPlot();
        plot.setDimensions(new Dimension3D(10.0, 4.0, 4.0));
        plot.setLegendLabelGenerator(new StandardXYZLabelGenerator(
                StandardXYZLabelGenerator.COUNT_TEMPLATE));
        ScatterXYZRenderer renderer = (ScatterXYZRenderer) plot.getRenderer();
        renderer.setSize(0.15);
        renderer.setColors(Colors.createIntenseColors());
        chart.setViewPoint(ViewPoint3D.createAboveLeftViewPoint(40));
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
                "Here are some lines", 
                "made randomly", dataset, null, null, 
                "description");
        CategoryPlot3D plot = (CategoryPlot3D) chart.getPlot();
        plot.setDimensions(new Dimension3D(18, 8, 4));
        plot.getRowAxis().setVisible(false);
        NumberAxis3D valueAxis = (NumberAxis3D) plot.getValueAxis();
        valueAxis.setTickSelector(new NumberTickSelector(true));
        plot.getRenderer().setColors(Colors.createFancyDarkColors());
        chart.setViewPoint(ViewPoint3D.createAboveViewPoint(30));
        return chart;    
    }
}
