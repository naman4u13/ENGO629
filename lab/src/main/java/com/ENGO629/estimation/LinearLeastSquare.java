package com.ENGO629.estimation;

import java.util.ArrayList;
import java.util.stream.IntStream;

import org.ejml.simple.SimpleMatrix;

import com.ENGO629.models.Satellite;

public class LinearLeastSquare {
	private final static double SpeedofLight = 299792458;

	public static double[] process(ArrayList<Satellite> satList, boolean isWLS) throws Exception {

		int n = satList.size();
		double[][] weight = new double[n][n];
		if (isWLS) {
			for (int i = 0; i < n; i++) {
				double elevAngle = Math.toRadians(satList.get(i).getElevation());
				double var = 1 / Math.pow(Math.sin(elevAngle), 2);
				weight[i][i] = 1 / var;
			}
		} else {
			IntStream.range(0, n).forEach(i -> weight[i][i] = 1);
		}
		double[] estEcefClk = new double[4];
		double error = Double.MAX_VALUE;
		// Get Millimeter Accuracy, actually it takes atleast 5 iterations to converge
		// to the result which is accurate more than micrometer scale
		double threshold = 1e-3;
		if (n >= 4) {

			while (error >= threshold) {

				double[][] deltaPR = new double[n][1];

				double[][] h = new double[n][4];
				for (int i = 0; i < n; i++) {

					Satellite sat = satList.get(i);
					double[] satECEF = sat.getECEF();
					double PR = sat.getPseduorange();

					double approxGR = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> satECEF[j] - estEcefClk[j])
							.map(j -> Math.pow(j, 2)).reduce((j, k) -> j + k).getAsDouble());

					double approxPR = approxGR + (SpeedofLight * estEcefClk[3]);
					deltaPR[i][0] = approxPR - PR;
					int index = i;
					IntStream.range(0, 3).forEach(j -> h[index][j] = (satECEF[j] - estEcefClk[j]) / approxGR);
					h[i][3] = 1;
				}
				SimpleMatrix H = new SimpleMatrix(h);
				SimpleMatrix Ht = H.transpose();
				SimpleMatrix W = new SimpleMatrix(weight);

				SimpleMatrix HtWHinv = (Ht.mult(W).mult(H)).invert();

				SimpleMatrix DeltaPR = new SimpleMatrix(deltaPR);
				SimpleMatrix DeltaX = HtWHinv.mult(Ht).mult(W).mult(DeltaPR);
				IntStream.range(0, 3).forEach(i -> estEcefClk[i] = estEcefClk[i] + DeltaX.get(i, 0));
				estEcefClk[3] += (-DeltaX.get(3, 0)) / SpeedofLight;
				error = Math.sqrt(IntStream.range(0, 3).mapToDouble(i -> Math.pow(DeltaX.get(i, 0), 2)).reduce(0,
						(i, j) -> i + j));

			}
			return estEcefClk;
		}

		throw new Exception("Satellite count is less than 4, can't compute user position");

	}

}
