package rungeKutta;

import java.io.IOException;

import differentialEquations.DifferentialEquation;

public abstract class Solver {

	protected double[][] y_values;
	protected double[] t_values;

	public double[][] getY_values() {
		return y_values;
	}
	public double[] getY_values(int i) throws IOException {
		if( i >= y_values[0].length || i < 0)
			throw new IOException("Solver GetY_value error: index exceeds array");
		
		double[] y = new double[y_values.length];
		for( int j = 0; j < y.length; j++)
			y[j]= y_values[j][i];
		
		return y;
	}
	public double[] getT_values() {
		return t_values;
	}
	
	// --------should also have public--------------- 
	 public void run(double[] y_0, double t_0, double t_end) throws IOException
	 {}
	 
	 public void run(double y_0, double t_0, double t_end) throws IOException
	 {
		 this.run(new double[]{y_0}, t_0 ,t_end);
	 }
	 
	 public void setEquation(DifferentialEquation equation) { }
	 public DifferentialEquation getEquation() {return null;}
}
