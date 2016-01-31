package examples;

import java.awt.print.Printable;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import rungeKutta.Solver;

import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.BitmapEncoder.BitmapFormat;
import org.knowm.xchart.internal.chartpart.Chart;
import org.knowm.xchart.QuickChart;
import org.knowm.xchart.SwingWrapper;

import differentialEquations.DifferentialEquation;
import differentialEquations.InitialValueProblem;
import differentialEquations.Sin_DE;


public class Plot{
	
	public List<Chart> charts;
	public SwingWrapper sw;
	
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
	public Plot(Solver solver,DifferentialEquation equation, double[]y_0, double t_0) throws IOException
	{
		charts = new ArrayList<Chart>();
		InitialValueProblem ivp = new InitialValueProblem(equation, y_0, t_0);
		solver.setEquation(ivp.getEquation());
		solver.run(ivp.getY_0(),ivp.getT_0(), 1);
		double[] xData = solver.getT_values();
		for (int i = 0; i < y_0.length; i++) 
		{
			double[] yData = solver.getY_values(i);
			charts.add(QuickChart.getChart("edited with plotter", "X", "Y" + i , "bla", xData, yData));
		}
		this.show();
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
