package rungeKuttaSolver;

public class Exp_DE implements DifferentialEquation 
{
	private double lambda;

	public Exp_DE(double lambda)
	{
		this.lambda = lambda;
	}

	@Override
	public double[] calculate(double t,double[] y)
	{
		return new double[]{lambda*y[0]};	//exp(lambda*t) solves y'= y*t				
	}
}