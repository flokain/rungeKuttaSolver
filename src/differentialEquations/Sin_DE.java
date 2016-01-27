package differentialEquations;

import java.io.IOException;


public class Sin_DE implements DifferentialEquation 
{
	private double omega;

	public Sin_DE(double omega)
	{
		this.omega = omega;
	}

	@Override
	public double[] calculate(double t,double[] y) throws IOException
	{
		if (y.length != 2) 
			throw new IOException("This is a "+Integer.toString(y.length)+" dimensional Differential Equation y must have the length "+Integer.toString(y.length));
		
		double dy  = y[1]; 			            //     y' = (y')
		double ddy = -omega*omega*y[0]; 	    //  (y')' = -omega^2*(y)
		return new double[]{dy,ddy};			// solution y = a*sin( omega *t) + b*cos( omega *t)
	}
}