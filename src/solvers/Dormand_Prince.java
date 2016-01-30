package solvers;

import java.io.IOException;

import differentialEquations.DifferentialEquation;
import rungeKutta.EmbeddedRungeKutta;

public class Dormand_Prince extends EmbeddedRungeKutta{

	public Dormand_Prince() throws IOException {
		super(Methods.DORMAND_PRINCE);
	}
	public Dormand_Prince(DifferentialEquation equation) throws IOException {
		super(Methods.DORMAND_PRINCE,equation);
	}
}
