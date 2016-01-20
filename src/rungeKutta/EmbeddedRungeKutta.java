package rungeKutta;

import java.io.IOException;
import differentialEquations.DifferentialEquation;

public class EmbeddedRungeKutta {
	
	RungeKutta leadRK;
	RungeKutta controlRK;
	int order;
	
	private double[][] y;
	private double[] t;
	
	public EmbeddedRungeKutta(double[][] A, double[] b_lead, double[] b_control, double[] c, DifferentialEquation equation) throws IOException
	{
		//check dimensions
		if(c.length != A[0].length || b_lead.length != A.length || b_control.length != A.length || c.length != b_lead.length)
			throw new IOException("RungeKutta Intialization Error: Dimensions of A,b,c incompatible");

		//check conditions for b1 and b2 : sum(b_i)=1 
		double sum1 = 0;
		double sum2 = 0;
		
		for (int i = 0; i < b_lead.length;i++) 
		{
			sum1 += b_lead[i];
			sum2 += b_control[i];
		}
		
		if (Math.abs(sum1 - 1) > 0.000001 || Math.abs(sum2 - 1) > 0.000001) 
			throw new IOException("RungeKutta Intialization Error: sum of elements of b is not equal to 1");
		
		leadRK = new RungeKutta(A, b_lead, c, equation);
		controlRK = new RungeKutta(A, b_control, c, equation);
	}
		
	public EmbeddedRungeKutta(String method, DifferentialEquation equation) throws IOException
	{
		double[] c = null;
		double[] b_lead = null;
		double[] b_control = null;
		double[][] A = null;
		switch (method)
		{
			case "Dormand–Prince":
			case "ode45":
			{
				this.order = 5;
				
				c = 	  	new double[]{0, 1./5, 3./10, 4./5, 8./9, 1., 1.};
				b_lead = 	new double[]{35./384, 		0,  500./1113, 	  125./192,   -2187./6784,     11./84,	   0    };
				b_control =new double[]{5179./57600,	0,	7571./16695,  393./640,   -92097./339200,  187./2100,  1./40};
				 
				A = new double[][] {
						 {0,			0,				0,				0,			0,				0,		0},
						 {1./5,       	0,				0,     			0,      	0,				0,		0},
					     {3./40,        9./40,   		0,  			0,    		0,         		0,		0},
					     {44./45,       -56./15,        32./9,    		0,     		0,     			0,		0},
					     {19372./6561,  -25360./2187,	64448./6561,    -212./729,  0,         		0,		0},
					     {9017./3168,   -355./33,       46732./5247,    49./176,    -5103./18656,   0,  	0},
					     {35./384,      0,       		500./1113,     	125./192,   -2187./6784, 	11./84,	0}};
				 
				break;
			}
		}
		/*-------------standard constructor code-------------*/
		// (because otherwise it has to be used in the first line, stupid :(
		
		// use standard parameters from ode45 from Matlab (20.1.2015):
		
			
		//check dimensions
		if(c.length != A[0].length || b_lead.length != A.length || b_control.length != A.length || c.length != b_lead.length)
			throw new IOException("RungeKutta Intialization Error: Dimensions of A,b,c incompatible");
	
		//check conditions for b1 and b2 : sum(b_i)=1 
		double sum1 = 0;
		double sum2 = 0;
		
		for (int i = 0; i < b_lead.length;i++) 
		{
			sum1 += b_lead[i];
			sum2 += b_control[i];
		}
		
		if (Math.abs(sum1 - 1) > 0.000001 || Math.abs(sum2 - 1) > 0.000001) 
			throw new IOException("RungeKutta Intialization Error: sum of elements of b is not equal to 1");
		
		leadRK = new RungeKutta(A, b_lead, c, equation);
		controlRK = new RungeKutta(A, b_control, c, equation);
		
		
		
		/*-------------END standard constructor END-------------*/
	}

	public EmbeddedRungeKutta(String method) throws IOException
	{
		this(method,null);
	}
	
	//solver
	public void run(double[] y_0, double t_0, double t_end, double stepSize_min, double tolerance ,double propability, double deltaStepsize_min, double stepSize_start) throws IOException
	{
		if ( Math.abs(stepSize_start) > Math.abs(t_end-t_0) )
			throw new IOException("Stepsize for the start is to big, choose it smaller than t_end-t_0");
		
		int y_size = y_0.length;
		
		t = new double[1];
		y = new double[1][y_size];
		
		t[0] = t_0;
		y[0] = y_0;
		
		double h = stepSize_start;
		int order = leadRK.getB().length;
		
		int s = 0; //last index <=> y.length-1
		while( t[s] < t_end)
		{
			leadRK.run(y[s], t[s], t[s]+h, h);
			double[] lead_y = leadRK.getY_values()[1];
			
			controlRK.run(y[s], t[s], t[s]+h, h);
			double[] control_y = controlRK.getY_values()[1];

			double estimatedError = h * this.calcRelError(y[s],lead_y,control_y,1e-3); 
			
			if( estimatedError / h <= tolerance || h <= stepSize_min) // tolerance or minimal stepSize reached => accept solution, make a step
			{	
				s = s+1;
				double[][] y_tmp = new double[s+1][y_size];
				double[] t_tmp = new double[s+1];
				
				for( int i = 0; i<s;i++)
				{
					t_tmp[i] = t[i];
					y_tmp[i] = y[i];
				}
				y = y_tmp;
				t = t_tmp;
				
				y[s] = lead_y;
				t[s] = t[s-1]+h;
				
			  /* //Code from script Melenk
				 h = Math.max( stepSize_min, 
							  Math.min( deltaStepsize_max * h, 
									    propability * Math.pow( (tolerance/estimatedError * Math.pow(h,order) ), order-1) ) );
			  */
				//code from matlab ode45
				h = Math.max( stepSize_min, 
						  			h * Math.max( deltaStepsize_min, 
						  							propability * Math.pow(tolerance/estimatedError, 1./order) 
						  						) 
						  	);
				
				if( t[s]+h > t_end)
					h= t_end-t[s];
			}
			else // estimated Error was to big, try again at same y,t with smaller stepSize
			{
				h = h/2;
			}
			
		}
	}

	public void run(double[] y_0, double t_0, double t_end, double stepSize_min, double tolerance ,double propability, double deltaStepsize_min) throws IOException
	{
		run(y_0, t_0, t_end, stepSize_min, tolerance, propability, deltaStepsize_min, guessFirstStep(y_0,t_0, t_end, stepSize_min, tolerance, propability));
	}
	public void run(double[] y_0, double t_0, double t_end) throws IOException
	{
		// stepSize_min = 1e-16;    //precision of double
		// deltaStepsize_max = 0.1; // matlab ode45
		// tolerance = 1e-3; 		// matlab ode45
		// propability = 0.8;		// matlab ode45
		run(y_0, t_0, t_end, 1e-16, 1e-3, 0.8, 0.1, guessFirstStep(y_0,t_0, t_end, 1e-16, 1e-3, 0.8));
	}

	public void run(double y_0, double t_0, double t_end) throws IOException
	{
		run(new double[]{y_0}, t_0, t_end);
	}
	
	
	// Setter
	public void setEquation(DifferentialEquation equation)
	{
		leadRK.setEquation(equation);
		controlRK.setEquation(equation);
	}
	
	//Geter
	public double[][] getY_values() {
		return y;
	}
	public double[] getY_values(int i) throws IOException {
		if( i >= y[0].length || i < 0)
			throw new IOException("Solver GetY_value error: index exceeds array");
		
		double[] y_ = new double[y.length];
		for( int j = 0; j < y_.length; j++)
			y_[j]= y[j][i];
		
		return y_;
	}
	public double[] getT_values() {
		return t;
	}
	
	// Internal functions
	private double[] diff(double[] a,double[] b)
	{
		if (a.length != b.length)
			new IOException("double[]-sum exception: diffrent length");
		double[] c = new double[a.length];
		for (int i = 0; i < c.length; i++)
			c[i]= a[i] - b[i];
		return c;
	}
	
	private double calcRelError(double[] y, double[]y_lead, double[] y_control, double absTolerance)
	{
		
		if (y_lead.length != y_control.length || y.length != y_control.length )
			new IOException("double[]-calcError exception: diffrent length");
		
		double c = 0;
		double[] yDiff = this.diff(y_lead, y_control);
		
		for (int i = 0; i < y.length; i++)
			c = Math.max(c, yDiff[i] / Math.max( Math.max( Math.abs(y_lead[i]) , Math.abs(y[i]) ),absTolerance) );
		
		return c;
	}

	private double guessFirstStep(double[] y_0, double t_0,double t_end, double stepSize_min, double tolerance, double propability)
	{
		double[] dy = leadRK.getEquation().calculate(t_0, y_0);
		double rh = 0;
		double h = 0.1 * (t_end-t_0);
		
		for (int i = 0; i < dy.length; i++)
			rh = Math.max(rh, Math.abs(dy[i] / Math.max( Math.abs( y_0[i] ), 1e-3) ) ); //  1e-3 is actually absolute tolerance / tolerance
		
		rh = rh /  ( propability * Math.pow(tolerance,1./order) );
		
		if(h*rh > 1) h = 1/rh; 
		
		h = Math.max(stepSize_min,h);	
		return h;
	}
}
