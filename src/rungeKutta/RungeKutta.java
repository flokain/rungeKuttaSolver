package rungeKutta;

import java.io.IOException;

import differentialEquations.DifferentialEquation;

public class RungeKutta extends Solver
{
	private double[] c;
	private double[][] A;
	private double[] b;
	private DifferentialEquation equation;

	
	public RungeKutta(double[][] A, double[] b, double[] c) throws IOException
	{
		this(A,b,c,null);
	}

	public RungeKutta(double[][] A, double[] b, double[] c , DifferentialEquation equation) throws IOException
	{
		//check dimensions
		if(c.length != A[0].length || b.length != A.length || c.length != b.length)
			throw new IOException("RungeKutta Intialization Error: Dimensions of A,b,c incompatible");

		//check conditions for b1 and b2 : sum(b_i)=1 
		double sum = 0;
		
		for (int i = 0; i < b.length;i++) 
			sum += b[i];
		
		if (Math.abs(sum - 1) > 0.000001) 
			throw new IOException("RungeKutta Intialization Error: sum of elements of b is not equal to 1");
		
		this.A = A;
		this.b = b;
		this.c = c;
		this.equation = equation;
	}
	
	public void run(double y_0, double t_0, double t_end, double stepSize) throws IOException
	{
		this.run(new double[]{y_0}, t_0 ,t_end ,stepSize);
	}
	
	public void run(double[] y_0, double t_0, double t_end, double stepSize) throws IOException
	{ 
		int stepsCount = (int) Math.floor( Math.abs( (t_end-t_0)/stepSize) );
		if (Math.abs(stepsCount) < 1) 
		{
			stepsCount = 1;
			stepSize = Math.abs(t_end-t_0);
		}
		int s = b.length;

		y_values = new double[stepsCount+1][y_0.length];
		t_values = new double[stepsCount+1];
		t_values[0] = t_0;
		y_values[0] = y_0;
		double h = Math.signum(t_end-t_0)*stepSize;
		
		for (int i = 0; i< stepsCount; i++)
		{
			double t = t_values[i]; 
			double y[] = y_values[i];
			double[][]k = new double[s][y_0.length];
			
			for (int j = 0; j < s; j ++)
			{
				double[] temp = new double[y.length];
				
				for (int m = 0; m < j; m ++)
				{
					temp = sum(temp,mul(A[j][m],k[m]));
				}
				
				k[j] = equation.calculate(t+c[j]*h, sum(y,mul(h,temp))) ;
				y_values[i+1] = sum(y_values[i+1],mul(b[j],k[j]));	
			}
			t_values[i+1] =  t_values[i]+h;
			y_values[i+1]= sum(y_values[i],mul(h,y_values[i+1]));
		}
	}	

	public void run(double[] y_0, double t_0, double t_end) throws IOException
	{
		run(y_0,  t_0,  t_end, (t_end-t_0)/1000);
	}
	
	public void setEquation(DifferentialEquation equation)
	{
		this.equation = equation;
	}
	
	public double[] getC() {
		return c;
	}
	
	public double[][] getA() {
		return A;
	}
	public double[] getB() {
		return b;
	}
	public DifferentialEquation getEquation() throws IOException {
		
		if (equation == null) throw new IOException("Equation not set");
		return equation;
	}
	
	private double[] sum(double[]a,double[]b)
	{
		if (a.length != b.length)
			new IOException("double[]-sum exception: diffrent length");
		double[] c = new double[a.length];
		for (int i = 0; i < c.length; i++)
			c[i]= a[i] + b[i];
		return c;
	}
	
	private double[] mul(double a,double[]b)
	{
		double[] c = new double[b.length];
		for (int i = 0; i < c.length; i++)
			c[i]= a * b[i];
		return c;
	}
	
}
