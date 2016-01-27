package examples;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import rungeKutta.Solver;

import com.xeiam.xchart.Chart;
import com.xeiam.xchart.QuickChart;
import com.xeiam.xchart.SwingWrapper;

import differentialEquations.DifferentialEquation;
import differentialEquations.Exp_DE;
import differentialEquations.InitialValueProblem;
import differentialEquations.Sin_DE;


public class Plot{
	
	public List<Chart> charts;
	public SwingWrapper sw;
	
	public Plot(double[] xData, double[] yData){
		Chart chart = QuickChart.getChart("", "X", "Y", "", xData, yData);	
		this.show();
	}
	
	public Plot(Solver solver, List<InitialValueProblem> ivps) throws IOException{
		
		charts = new ArrayList<Chart>();
		for (InitialValueProblem ivp : ivps)
		{	
			solver.setEquation(ivp.getEquation());
			solver.run(ivp.getY_0(),ivp.getT_0(), 1);
			double[] xData = solver.getT_values();
			double[] yData = solver.getY_values(0);
			charts.add(QuickChart.getChart("edited with plotter", "X", "Y", "bla", xData, yData));
		}
		this.show();
	}
	
	public Plot(Solver solver) throws IOException
	{
		this(solver, standartProblems());
	}
	
	private void show(){
		new SwingWrapper(charts).displayChartMatrix("Standard Problems");
	}
	
	private static List<InitialValueProblem> standartProblems() throws IOException
	{
		List<InitialValueProblem> list = new ArrayList<InitialValueProblem>();
		
		list.add(new InitialValueProblem( new Sin_DE(1),new double[]{1,0},0));
		//list.add(new InitialValueProblem( new Exp_DE(1), 1, 0));
		
		return list;
		
	}
}