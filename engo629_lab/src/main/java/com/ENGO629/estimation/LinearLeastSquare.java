package com.ENGO629.estimation;

import java.util.ArrayList;
import java.util.stream.IntStream;

import org.ejml.simple.SimpleMatrix;

import com.ENGO629.models.Satellite;

public class LinearLeastSquare {
	private final static double SpeedofLight = 299792458;

	public static double[] process(ArrayList<Satellite> satList, boolean isWLS) throws Exception {

		// Satellite count
		int n = satList.size();
		// Weight matrix
		double[][] weight = new double[n][n];
		/*
		 * If 'isWLS' flag is true, the estimation method is WLS and weight matrix will
		 * be based on elevation angle otherwise identity matrix will assigned for LS
		 */
		if (isWLS) {
			for (int i = 0; i < n; i++) {
				double elevAngle = Math.toRadians(satList.get(i).getElevation());
				double var = 1 / Math.pow(Math.sin(elevAngle), 2);
				weight[i][i] = 1 / var;
			}

		} else {
			IntStream.range(0, n).forEach(i -> weight[i][i] = 1);
		}
		// variable to store estimated Rx position and clk offset
		double[] estEcefClk = new double[4];
		/*
		 * Error variable based on norm value deltaX vector, intially assigned a big
		 * value
		 */
		double error = Double.MAX_VALUE;
		// Threshold to stop iteration or regression
		double threshold = 1e-3;
		// Minimum 4 satellite are required to proceed
		if (n >= 4) {

			while (error >= threshold) {

				// Misclosure vector
				double[][] deltaPR = new double[n][1];
				// Jacobian or Design Matrix
				double[][] h = new double[n][4];
				// Iterate through each satellite, to compute LOS vector and Approx pseudorange
				for (int i = 0; i < n; i++) {

					Satellite sat = satList.get(i);
					double[] satECEF = sat.getECEF();
					double PR = sat.getPseduorange();
					// Approx Geometric Range
					double approxGR = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> satECEF[j] - estEcefClk[j])
							.map(j -> Math.pow(j, 2)).reduce((j, k) -> j + k).getAsDouble());
					// Approx Pseudorange Range
					double approxPR = approxGR + (SpeedofLight * estEcefClk[3]);
					deltaPR[i][0] = approxPR - PR;
					int index = i;
					IntStream.range(0, 3).forEach(j -> h[index][j] = (satECEF[j] - estEcefClk[j]) / approxGR);
					h[i][3] = 1;
				}
				// Least Squares implementation
				SimpleMatrix H = new SimpleMatrix(h);
				SimpleMatrix Ht = H.transpose();
				SimpleMatrix W = new SimpleMatrix(weight);
				SimpleMatrix HtWHinv = (Ht.mult(W).mult(H)).invert();
				SimpleMatrix DeltaPR = new SimpleMatrix(deltaPR);
				SimpleMatrix DeltaX = HtWHinv.mult(Ht).mult(W).mult(DeltaPR);
				// updating Rx state vector, by adding deltaX vector
				IntStream.range(0, 3).forEach(i -> estEcefClk[i] = estEcefClk[i] + DeltaX.get(i, 0));
				estEcefClk[3] += (-DeltaX.get(3, 0)) / SpeedofLight;
				// Recomputing error - norm of deltaX vector
				error = Math.sqrt(IntStream.range(0, 3).mapToDouble(i -> Math.pow(DeltaX.get(i, 0), 2)).reduce(0,
						(i, j) -> i + j));

			}
			/*
			 * Regression is completed, error is below threshold, successfully estimated Rx
			 * Position and Clk Offset
			 */

			return estEcefClk;
		}

		throw new Exception("Satellite count is less than 4, can't compute user position");

	}

}
