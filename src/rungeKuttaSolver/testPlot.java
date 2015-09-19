package rungeKuttaSolver;


import java.io.IOException;

import com.xeiam.xchart.BitmapEncoder;
import com.xeiam.xchart.BitmapEncoder.BitmapFormat;
import com.xeiam.xchart.Chart;
import com.xeiam.xchart.QuickChart;
import com.xeiam.xchart.SwingWrapper;

public class testPlot {
	 public static void main(String[] args) throws IOException
	 {
		 double t_0 = 0;
		 double lambda = -21;
		 double y_0 = 1;
		 double stepSize = 0.1;
		 double t_end = 10;
		 
		 DifferentialEquation equation = new Exp_DE(lambda);
		 forwardEuler solver = new forwardEuler(equation);
		 solver.run(y_0, t_0, t_end, stepSize);		 
		 
		 double[] xData = solver.t_values;
		 double[] yData = solver.y_values;

		    // Create Chart
		    Chart chart = QuickChart.getChart("Sample Chart", "X", "Y", "y(x)", xData, yData);

		    // Show it
		    new SwingWrapper(chart).displayChart();

		    // Save it
		    BitmapEncoder.saveBitmap(chart, "./Sample_Chart1", BitmapFormat.PNG);

		    // or save it in high-res
		    BitmapEncoder.saveBitmapWithDPI(chart, "./Sample_Chart_300_DPI1", BitmapFormat.PNG, 300);
	 }
}
