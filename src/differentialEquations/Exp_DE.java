package differentialEquations;


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
		return new double[]{lambda*y[0]};	// a * exp(lambda*t) solves y'= y*lambda				
	}
}