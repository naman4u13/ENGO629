package com.ENGO629.estimation.kf;

import java.util.stream.IntStream;

import org.ejml.simple.SimpleMatrix;

import com.ENGO629.util.LatLonUtil;

public class KFconfig extends KF {

	private final double SpeedofLight = 299792458;
	private final double c2 = SpeedofLight * SpeedofLight;
	// Typical Allan Variance Coefficients for TCXO (low quality)
	private final double h0 = 2E-19;
	private final double h_2 = 2E-20;
	private final double sf = h0 / 2;
	private final double sg = 2 * Math.PI * Math.PI * h_2;

	public void configSPP(double deltaT, Flag flag) {

		double[] ecef = new double[] { getState().get(0), getState().get(1), getState().get(2) };
		SimpleMatrix R = getR(ecef);
		if (flag == Flag.POSITION) {

			double[][] F = new double[5][5];
			double[][] Q = new double[5][5];
			IntStream.range(0, 5).forEach(i -> F[i][i] = 1);
			F[3][4] = deltaT;
			double[] qENU = new double[] { 100, 100, 0 };
			double[] qECEF = enuToEcef(R, qENU);
			double[] q = qECEF;
			// double[] q = new double[] { 1000, 1000, 1000 };
			IntStream.range(0, 3).forEach(i -> Q[i][i] = q[i]);
			Q[3][3] = ((sf * deltaT) + ((sg * Math.pow(deltaT, 3)) / 3)) * c2;
			Q[3][4] = (sg * Math.pow(deltaT, 2)) * c2 / 2;
			Q[4][3] = (sg * Math.pow(deltaT, 2)) * c2 / 2;
			Q[4][4] = sg * deltaT * c2;
			super.configure(F, Q);

		} else if (flag == Flag.VELOCITY) {
			double[][] F = new double[8][8];
			double[][] Q = new double[8][8];
			IntStream.range(0, 8).forEach(i -> F[i][i] = 1);
			IntStream.range(0, 4).forEach(i -> F[i][i + 4] = deltaT);
			double[] qENU = new double[] { 4, 4, 0 };
			double[] qECEF = enuToEcef(R, qENU);
			// double[] q = new double[] { 25, 25, 25, sg * c2 };
			double[] q = new double[] { qECEF[0], qECEF[1], qECEF[2], sg * c2 };
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

	private SimpleMatrix getR(double[] ecef) {
		double[] llh = LatLonUtil.ecef2lla(ecef);
		double lat = Math.toRadians(llh[0]);
		double lon = Math.toRadians(llh[1]);
		double[][] r = new double[][] {
				{ -Math.sin(lon), -Math.sin(lat) * Math.cos(lon), Math.cos(lat) * Math.cos(lon) },
				{ Math.cos(lon), -Math.sin(lat) * Math.sin(lon), Math.cos(lat) * Math.sin(lon) },
				{ 0, Math.cos(lat), Math.sin(lat) } };
		SimpleMatrix R = new SimpleMatrix(r);
		return R;
	}

	private double[] enuToEcef(SimpleMatrix R, double[] _enuParam) {

		SimpleMatrix enuParam = new SimpleMatrix(3, 1, false, _enuParam);
		SimpleMatrix _ecefParam = R.mult(enuParam);
		double[] ecefParam = new double[] { _ecefParam.get(0), _ecefParam.get(1), _ecefParam.get(2) };
		return ecefParam;
	}

}
