package differentialEquations;

import java.io.IOException;

public interface DifferentialEquation
{
	public double[] calculate(double t, double[] y) throws IOException;
}
