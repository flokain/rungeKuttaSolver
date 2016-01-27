package differentialEquations;

import java.io.IOException;

public class InitialValueProblem {
	
	private DifferentialEquation equation;
	private double[] y_0;
	private double t_0;
	
	public InitialValueProblem(DifferentialEquation equation,double[]y_0,double t_0) throws IOException 
	{
		equation.calculate(t_0, y_0); //Test for correct length of y
		
		this.equation = equation;
		this.t_0 = t_0;
		this.y_0 = y_0;
	}
	public InitialValueProblem(DifferentialEquation equation,double y_0,double t_0) throws IOException 
	{
		
		equation.calculate(t_0, new double[]{y_0}); //Test for correct length of y
		
		this.equation = equation;
		this.t_0 = t_0;
		this.y_0 = new double[]{y_0};
	}
	
	public DifferentialEquation getEquation() {
		return equation;
	}

	public void setEquation(DifferentialEquation equation) throws IOException {
		equation.calculate(this.t_0, this.y_0); //Test for correct length of y
		this.equation = equation;
	}

	public double[] getY_0() {
		return y_0;
	}

	public void setY_0(double[] y_0) throws IOException {
		
		equation.calculate(this.t_0, y_0); //Test for correct length of y
		this.y_0 = y_0;
	}
	
	public void setY_0(double y_0) throws IOException {
		
		equation.calculate(this.t_0, new double[]{y_0}); //Test for correct length of y
		this.y_0 = new double[]{y_0};
	}

	public double getT_0() {
		return t_0;
	}

	public void setT_0(double t_0) {
		this.t_0 = t_0;
	}

	
}

