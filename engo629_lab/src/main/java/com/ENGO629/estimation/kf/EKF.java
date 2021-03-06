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
			double[] intialECEF, Flag flag) throws Exception {

		int n = 0;
		/* constant position model - state vector(n=5) -> (x,y,z,cdt,cdt_dot) */
		/*
		 * constant velocity model - state vector(n=8) ->
		 * (x,y,z,cdt,x_dot,y_dot,z_dot,cdt_dot)
		 */
		if (flag == Flag.POSITION) {
			n = 5;
		} else if (flag == Flag.VELOCITY) {
			n = 8;
		}
		double[][] x = new double[n][1];
		double[][] P = new double[n][n];
		/*
		 * state XYZ is intialized WLS generated Rx position estimated using first epoch
		 * data, a-priori estimate error covariance matrix for state XYZ is therefore
		 * assigned 25 m^2 value. Other state variables are assigned infinite(big)
		 * variance
		 */
		IntStream.range(0, 3).forEach(i -> x[i][0] = intialECEF[i]);
		IntStream.range(0, 3).forEach(i -> P[i][i] = 25);
		IntStream.range(3, n).forEach(i -> P[i][i] = 1e13);

		kfObj.setState_ProcessCov(x, P);
		// Begin iteration or recursion
		return iterate(satMap, timeList, flag);

	}

	private ArrayList<double[]> iterate(HashMap<Integer, ArrayList<Satellite>> satMap, ArrayList<Integer> timeList,
			Flag flag) throws Exception {
		ArrayList<double[]> ecefList = new ArrayList<double[]>();
		int time = timeList.get(0);
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {
			int currentTime = timeList.get(i);
			ArrayList<Satellite> satList = satMap.get(currentTime);
			double deltaT = currentTime - time;
			// Perform Predict and Update
			runFilter(deltaT, satList, flag);
			// Fetch Posteriori state estimate and estimate error covariance matrix
			SimpleMatrix x = kfObj.getState();
			SimpleMatrix P = kfObj.getCovariance();
			double[] estECEF = new double[] { x.get(0), x.get(1), x.get(2) };
			// Add position estimate to the list
			ecefList.add(estECEF);
			/*
			 * Check whether estimate error covariance matrix is positive semidefinite
			 * before further proceeding
			 */
			if (!MatrixFeatures_DDRM.isPositiveDefinite(P.getMatrix())) {

				throw new Exception("PositiveDefinite test Failed");
			}
			time = currentTime;

		}

		return ecefList;
	}

	private void runFilter(double deltaT, ArrayList<Satellite> satList, Flag flag) {

		// Satellite count
		int n = satList.size();

		// Assign Q and F matrix
		kfObj.config(deltaT, flag);
		kfObj.predict();

		SimpleMatrix x = kfObj.getState();
		double[] estECEF = new double[] { x.get(0), x.get(1), x.get(2) };
		double rxClkOff = x.get(3);// in meters

		/*
		 * H is the Jacobian matrix of partial derivatives Observation StateModel(h) of
		 * with respect to x
		 */
		double[][] H = getJacobian(satList, estECEF, x.numRows());
		// Measurement vector
		double[][] z = new double[n][1];
		IntStream.range(0, n).forEach(i -> z[i][0] = satList.get(i).getPseduorange());
		// Estimated Measurement vector
		double[][] ze = new double[n][1];
		IntStream.range(0, n)
				.forEach(i -> ze[i][0] = Math
						.sqrt(IntStream.range(0, 3).mapToDouble(j -> estECEF[j] - satList.get(i).getECEF()[j])
								.map(j -> j * j).reduce(0, (j, k) -> j + k))
						+ rxClkOff);
		// Measurement Noise
		double[][] R = new double[n][n];
		for (int i = 0; i < n; i++) {
			double elevAngle = Math.toRadians(satList.get(i).getElevation());
			double var = 1 / Math.pow(Math.sin(elevAngle), 2);
			R[i][i] = var;
		}
		// Perform Update Step
		kfObj.update(z, R, ze, H);

	}

	private double[][] getJacobian(ArrayList<Satellite> satList, double[] estECEF, int stateN) {
		int n = satList.size();
		double[][] H = new double[n][stateN];

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
