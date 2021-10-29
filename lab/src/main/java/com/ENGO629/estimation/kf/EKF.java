package com.ENGO629.estimation.kf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.stream.IntStream;

import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.ENGO629.models.Satellite;

public class EKF {

	private KFconfig kfObj;
	private double prObsNoiseVar;

	public EKF() {
		kfObj = new KFconfig();
	}

	public ArrayList<double[]> process(HashMap<Integer, ArrayList<Satellite>> satMap, ArrayList<Integer> timeList,
			double[] intialECEF) throws Exception {

		double[][] x = new double[][] { { intialECEF[0] }, { intialECEF[1] }, { intialECEF[2] }, { 0 }, { 0 }, { 0 },
				{ 0 }, { 0 } };
		double[][] P = new double[8][8];
		IntStream.range(0, 3).forEach(i -> P[i][i] = 25);
		IntStream.range(3, 8).forEach(i -> P[i][i] = 1e13);

		prObsNoiseVar = 100;
		kfObj.setState_ProcessCov(x, P);
		return iterateSPP(satMap, timeList);

	}

	private ArrayList<double[]> iterateSPP(HashMap<Integer, ArrayList<Satellite>> satMap, ArrayList<Integer> timeList)
			throws Exception {
		ArrayList<double[]> ecefList = new ArrayList<double[]>();
		int time = timeList.get(0);
		for (int i = 1; i < timeList.size(); i++) {
			int currentTime = timeList.get(i);
			ArrayList<Satellite> satList = satMap.get(currentTime);
			double deltaT = currentTime - time;
			runFilter(deltaT, satList);
			SimpleMatrix x = kfObj.getState();
			SimpleMatrix P = kfObj.getCovariance();
			double[] estECEF = new double[] { x.get(0), x.get(1), x.get(2) };
			ecefList.add(estECEF);
			if (!MatrixFeatures_DDRM.isPositiveDefinite(P.getMatrix())) {

				throw new Exception("PositiveDefinite test Failed");
			}
			time = currentTime;

		}

		return ecefList;
	}

	private void runFilter(double deltaT, ArrayList<Satellite> satList) {

		int n = satList.size();

		kfObj.configSPP(deltaT);
		kfObj.predict();

		SimpleMatrix x = kfObj.getState();
		double[] estECEF = new double[] { x.get(0), x.get(1), x.get(2) };
		double rxClkOff = x.get(3);// in meters

		// H is the Jacobian matrix of partial derivatives Observation StateModel(h) of
		// with
		// respect to x
		double[][] H = getJacobian(satList, estECEF);
		double[][] z = new double[n][1];
		IntStream.range(0, n).forEach(i -> z[i][0] = satList.get(i).getPseduorange());
		double[][] ze = new double[n][1];
		IntStream.range(0, n)
				.forEach(i -> ze[i][0] = Math
						.sqrt(IntStream.range(0, 3).mapToDouble(j -> estECEF[j] - satList.get(i).getECEF()[j])
								.map(j -> j * j).reduce(0, (j, k) -> j + k))
						+ rxClkOff);
		double[][] R = new double[n][n];
		IntStream.range(0, n).forEach(i -> R[i][i] = prObsNoiseVar);
		kfObj.update(z, R, ze, H);

	}

	private double[][] getJacobian(ArrayList<Satellite> satList, double[] estECEF) {
		int n = satList.size();
		double[][] H = new double[n][8];

		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			// Line of Sight vector
			double[] LOS = IntStream.range(0, 3).mapToDouble(j -> sat.getECEF()[j] - estECEF[j]).toArray();
			// Geometric Range
			double GR = Math.sqrt(Arrays.stream(LOS).map(j -> j * j).reduce(0.0, (j, k) -> j + k));
			// Converting LOS to unit vector
			final int _i = i;
			IntStream.range(0, 3).forEach(j -> H[_i][j] = -LOS[j] / GR);
			H[i][3] = 1;
		}

		return H;

	}
}
