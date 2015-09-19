package rungeKuttaSolver;

import java.util.ArrayList;

public class forwardEuler 
{
	public double[] y_values;
	public double[] t_values;
	DifferentialEquation function;
	public forwardEuler(DifferentialEquation function)
	{
		this.function = function;
	}
	public void run(double y_0, double t_0, double t_end, double stepSize)
	{
		int stepsCount = (int)((t_end-t_0)/stepSize + 1);
		
		y_values = new double[stepsCount];
		t_values = new double[stepsCount];
		
		t_values[0] = t_0;
		y_values[0] = y_0;
		for(int i = 1; i< stepsCount ; i++)
		{
			double time = t_0+i*stepSize;
			double dy = function.calculate(time,new double[]{y_values[i-1]})[0];
			t_values[i] = time;
			y_values[i]= y_values[i-1] + dy*stepSize;
		}	
	}

}
