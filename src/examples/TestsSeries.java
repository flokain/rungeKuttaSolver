package examples;


import java.io.IOException;

import rungeKutta.EmbeddedRungeKutta;
import rungeKutta.RungeKutta;
import rungeKutta.Solver;

import com.xeiam.xchart.BitmapEncoder;
import com.xeiam.xchart.BitmapEncoder.BitmapFormat;
import com.xeiam.xchart.Chart;
import com.xeiam.xchart.QuickChart;
import com.xeiam.xchart.SwingWrapper;

import differentialEquations.DifferentialEquation;
import differentialEquations.Exp_DE;
import differentialEquations.Sin_DE;

public class TestsSeries {
	 public static void main(String[] args) throws IOException
	 {
		 double t_0 = 0;
		 double lambda = -27;
		 double[] y_0 = new double[] {1, 0};
		 double stepSize = 0.001;
		 double t_end = 1;
		 
		 DifferentialEquation equation = new Sin_DE(1);
		 
		 //Dorman_Prince
		 Solver solver = new EmbeddedRungeKutta("ode45",equation);
		 solver.run(y_0, t_0, t_end);	
		 
		 //forward euler
		 double[][] A = new double[][]{{0}};
		 double[] b = new double[]{1};
		 double[] c = new double[]{1};
		 RungeKutta solver2 = new RungeKutta(A, b, c);
		 solver2.setEquation(equation);
		 solver2.run(y_0, t_0, t_end, stepSize);
		 
		 //Heun
		 A = new double[][]{{0		,0	},
				 			{1		,0	}};
		 b = new double[]{0.5,0.5};
		 c = new double[]{0,1};
		 RungeKutta solver3 = new RungeKutta(A, b, c);
		 
		 solver3.setEquation(equation);
		 solver3.run(y_0, t_0, t_end, stepSize);

		 //RK4 (Simpson)
		 A = new double[][]{{0,		0,		0,		0	},
				 			{0.5,	0,		0,		0	},
				 			{0,		0.5,	0,		0	},
				 			{0,		0,		1,		0	}};
		 
		 b = new double[]	{1./6	,2./6,	2./6,	1./6};
		 c = new double[]	{0		,0.5,	0.5,	1};
		 RungeKutta solver4 = new RungeKutta(A, b, c);
		 
		 solver4.setEquation(equation);
		 solver4.run(y_0, t_0, t_end, stepSize);
		 
		 // get plot data
		 new Plot(solver);
		 
		 double[] xData = solver.getT_values();
		 double[] yData = solver.getY_values(0);
		 double[] xData2 = solver2.getT_values();
		 double[] yData2 = solver2.getY_values(0);
		 double[] xData3 = solver3.getT_values();
		 double[] yData3 = solver3.getY_values(0);
		 double[] xData4 = solver4.getT_values();
		 double[] yData4 = solver4.getY_values(0);
		 
		    // Create Chart
		    Chart chart = QuickChart.getChart("Dormand_Prince", "X", "Y", "y' = -27 * y relative tolerance= 1e-3", xData, yData);
		    Chart chart2 = QuickChart.getChart("forward euler", "X", "Y", "y' = -27 * y Stepsize= 0.1", xData2, yData2);
		    Chart chart3 = QuickChart.getChart("Heun", "X", "Y", "y' = -27 * y Stepsize= 0.1", xData3, yData3);
		    Chart chart4 = QuickChart.getChart("RK4 (Simpson)", "X", "Y", "y' = -27 * y Stepsize= 0.1", xData4, yData4); 
		    // Show it
		    new SwingWrapper(chart).displayChart();
		    new SwingWrapper(chart2).displayChart();
		    new SwingWrapper(chart3).displayChart();
		    new SwingWrapper(chart4).displayChart();
		    // Save it
		    BitmapEncoder.saveBitmap(chart, "./Sample_Chart1", BitmapFormat.PNG);

		    // or save it in high-res
		    BitmapEncoder.saveBitmapWithDPI(chart, "./Sample_Chart_300_DPI1", BitmapFormat.PNG, 300);
	 }
}
