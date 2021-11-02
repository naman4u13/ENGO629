package com.ENGO629;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.stream.IntStream;

import com.ENGO629.estimation.LinearLeastSquare;
import com.ENGO629.estimation.kf.EKF;
import com.ENGO629.estimation.kf.Flag;
import com.ENGO629.models.Satellite;
import com.ENGO629.util.GraphPlotter;
import com.ENGO629.util.LatLonUtil;

public class MainApp {

	public static void main(String[] args) throws IOException {
		try {

			// Path to store output file
			String path = "D:\\projects\\eclipse_projects\\UCalgary\\ENGO629\\results\\output";
			File output = new File(path + ".txt");
			PrintStream stream;
			stream = new PrintStream(output);
			System.setOut(stream);
			ArrayList<Integer> timeList = new ArrayList<Integer>();
			ArrayList<double[]> rxEcefList = new ArrayList<double[]>();
			// Map based variable to store enu errors for each algorithm
			HashMap<String, ArrayList<double[]>> EnuMap = new HashMap<String, ArrayList<double[]>>();
			String fileName = "satpos_meas.txt";
			InputStream is = MainApp.class.getClassLoader().getResourceAsStream(fileName);
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String line = null;
			// Map based variable to store satellite data, (k,v) -> (epoch,satellite list)
			HashMap<Integer, ArrayList<Satellite>> satMap = new HashMap<Integer, ArrayList<Satellite>>();
			// Parsing "satpos_meas.txt" file
			while ((line = br.readLine()) != null) {
				double[] data = Arrays.stream(line.split("\\s+")).mapToDouble(i -> Double.parseDouble(i)).toArray();
				satMap.computeIfAbsent((int) data[0], k -> new ArrayList<Satellite>()).add(new Satellite(data));
			}
			fileName = "refpos.txt";
			is = MainApp.class.getClassLoader().getResourceAsStream(fileName);
			br = new BufferedReader(new InputStreamReader(is));
			// Choose which type of PVT method to implement - LS, WLS, EKF
			int estimatorType = 3;
			// Parsing 'refpos.txt' line by line
			while ((line = br.readLine()) != null) {
				double[] data = Arrays.stream(line.split("\\s+")).mapToDouble(i -> Double.parseDouble(i)).toArray();
				// epoch GPS time
				int t = (int) data[0];
				// True receiver(rx) ECEF position
				double[] rxECEF = Arrays.copyOfRange(data, 1, 4);
				// Get satellite list for that epoch
				ArrayList<Satellite> satList = satMap.get(t);
				// variable to store estimated rx position and clk offset
				double[] estEcefClk = null;
				switch (estimatorType) {
				case 1:
					// Implement LS method
					estEcefClk = LinearLeastSquare.process(satList, false);
					EnuMap.computeIfAbsent("LS", k -> new ArrayList<double[]>()).add(estimateENU(estEcefClk, rxECEF));
					break;
				case 2:
					// Implement WLS method
					estEcefClk = LinearLeastSquare.process(satList, true);
					EnuMap.computeIfAbsent("WLS", k -> new ArrayList<double[]>()).add(estimateENU(estEcefClk, rxECEF));

					break;
				case 3:
					// Implement all PVT methods
					estEcefClk = LinearLeastSquare.process(satList, false);
					EnuMap.computeIfAbsent("LS", k -> new ArrayList<double[]>()).add(estimateENU(estEcefClk, rxECEF));

					estEcefClk = LinearLeastSquare.process(satList, true);
					EnuMap.computeIfAbsent("WLS", k -> new ArrayList<double[]>()).add(estimateENU(estEcefClk, rxECEF));

					break;

				}

				rxEcefList.add(rxECEF);
				timeList.add(t);

			}

			if (estimatorType == 4 || estimatorType == 3) {
				EKF ekf = new EKF();
				double[] initialECEF = LinearLeastSquare.process(satMap.get(timeList.get(0)), false);
				// Implement EKF based on receiver’s position and clock offset errors as a
				// random walk process
				ArrayList<double[]> estEcefList_pos = ekf.process(satMap, timeList, initialECEF, Flag.POSITION);
				// Implement EKF based on receiver’s velocity and clock drift errors as a random
				// walk process
				ArrayList<double[]> estEcefList_vel = ekf.process(satMap, timeList, initialECEF, Flag.VELOCITY);
				for (int i = 1; i < timeList.size(); i++) {
					EnuMap.computeIfAbsent("EKF - pos. random walk", k -> new ArrayList<double[]>())
							.add(estimateENU(estEcefList_pos.get(i - 1), rxEcefList.get(i)));

					EnuMap.computeIfAbsent("EKF - vel. random walk", k -> new ArrayList<double[]>())
							.add(estimateENU(estEcefList_vel.get(i - 1), rxEcefList.get(i)));

				}
			}

			// Calculate Accuracy Metrics
			HashMap<String, ArrayList<double[]>> GraphEnuMap = new HashMap<String, ArrayList<double[]>>();
			for (String key : EnuMap.keySet()) {
				ArrayList<Double>[] errList = new ArrayList[5];
				IntStream.range(0, 5).forEach(i -> errList[i] = new ArrayList<Double>());
				ArrayList<double[]> enuList = EnuMap.get(key);
				int n = enuList.size();

				for (int i = 0; i < n; i++) {
					double[] enu = enuList.get(i);
					// error in East direction
					errList[0].add(Math.sqrt(enu[0] * enu[0]));
					// error in North direction
					errList[1].add(Math.sqrt(enu[1] * enu[1]));
					// error in Up direction
					errList[2].add(Math.sqrt(enu[2] * enu[2]));
					// 3d error
					errList[3].add(Math.sqrt(Arrays.stream(enu).map(j -> j * j).sum()));
					// 2d error
					errList[4].add(Math.sqrt((enu[0] * enu[0]) + (enu[1] * enu[1])));

				}

				GraphEnuMap.put(key, enuList);

				// RMSE
				System.out.println("\n" + key + " RMSE");
				System.out.println("RMS - ");
				System.out.println(" E - " + RMS(errList[0]));
				System.out.println(" N - " + RMS(errList[1]));
				System.out.println(" U - " + RMS(errList[2]));
				System.out.println(" 3d Error - " + RMS(errList[3]));
				System.out.println(" 2d Error - " + RMS(errList[4]));

				// 95th Percentile
				IntStream.range(0, 5).forEach(i -> Collections.sort(errList[i]));
				int q95 = (int) (n * 0.95);

				System.out.println("\n" + key + " 95%");
				System.out.println("RMS - ");
				System.out.println(" E - " + errList[0].get(q95));
				System.out.println(" N - " + errList[1].get(q95));
				System.out.println(" U - " + errList[2].get(q95));
				System.out.println(" 3d Error - " + errList[3].get(q95));
				System.out.println(" 2d Error - " + errList[4].get(q95));

			}

			// Plot Error Graphs
			GraphPlotter.graphENU(GraphEnuMap, timeList);

		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}

	}

	public static double[] estimateError(double[] estEcefClk, double[] rxECEF) {

		double[] rxLL = LatLonUtil.ecef2lla(rxECEF);

		double ecefErr = Math.sqrt(IntStream.range(0, 3).mapToDouble(x -> rxECEF[x] - estEcefClk[x])
				.map(x -> Math.pow(x, 2)).reduce(0, (a, b) -> a + b));
		double[] estLL = LatLonUtil.ecef2lla(estEcefClk);

		// Great Circle Distance
		double gcErr = LatLonUtil.getHaversineDistance(estLL, rxLL);
		double[] err = new double[] { ecefErr, gcErr };
		return err;
	}

	// Transform from ECEF to ENU
	public static double[] estimateENU(double[] estEcefClk, double[] rxECEF) {
		double[] enu = LatLonUtil.ecef2enu(estEcefClk, rxECEF);
		return enu;
	}

	public static double RMS(ArrayList<Double> list) {
		return Math.sqrt(list.stream().mapToDouble(x -> x * x).average().orElse(Double.NaN));
	}

	public static double MAE(ArrayList<Double> list) {
		return list.stream().mapToDouble(x -> x).average().orElse(Double.NaN);
	}

}
