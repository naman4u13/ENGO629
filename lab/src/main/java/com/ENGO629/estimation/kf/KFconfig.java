package com.ENGO629.estimation.kf;

import java.util.stream.IntStream;

public class KFconfig extends KF {

	private final double SpeedofLight = 299792458;
	private final double c2 = SpeedofLight * SpeedofLight;
	// Typical Allan Variance Coefficients for TCXO (low quality)
	private final double h0 = 2E-19;
	private final double h_2 = 2E-20;
	private final double sf = h0 / 2;
	private final double sg = 2 * Math.PI * Math.PI * h_2;

	public void configSPP(double deltaT) {

		double[][] F = new double[8][8];
		double[][] Q = new double[8][8];
		IntStream.range(0, 8).forEach(x -> F[x][x] = 1);
		IntStream.range(0, 4).forEach(x -> F[x][x + 4] = deltaT);
		double[] q = new double[] { 9, 9, 4, sg * c2 };
		for (int i = 0; i < 4; i++) {
			Q[i][i] = q[i] * Math.pow(deltaT, 3) / 3;
			Q[i][i + 4] = q[i] * Math.pow(deltaT, 2) / 2;
			Q[i + 4][i] = q[i] * Math.pow(deltaT, 2) / 2;
			Q[i + 4][i + 4] = q[i] * deltaT;
		}
		Q[3][3] += (sf * c2 * deltaT);
		super.configure(F, Q);
	}

}
