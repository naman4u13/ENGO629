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

import org.jfree.ui.RefineryUtilities;

import com.ENGO629.estimation.LinearLeastSquare;
import com.ENGO629.estimation.kf.EKF;
import com.ENGO629.estimation.kf.Flag;
import com.ENGO629.models.Satellite;
import com.ENGO629.util.GraphPlotter;
import com.ENGO629.util.LatLonUtil;

public class MainApp {

	public static void main(String[] args) throws IOException {
		try {

			String path = "D:\\projects\\eclipse_projects\\UCalgary\\ENGO629\\results\\output2";
			File output = new File(path + ".txt");
			PrintStream stream;
			stream = new PrintStream(output);
			System.setOut(stream);
			ArrayList<Integer> timeList = new ArrayList<Integer>();
			ArrayList<double[]> rxEcefList = new ArrayList<double[]>();
			HashMap<String, ArrayList<double[]>> ErrMap = new HashMap<String, ArrayList<double[]>>();
			HashMap<String, ArrayList<double[]>> EnuMap = new HashMap<String, ArrayList<double[]>>();
			String fileName = "satpos_meas.txt";
			InputStream is = MainApp.class.getClassLoader().getResourceAsStream(fileName);
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String line = null;
			HashMap<Integer, ArrayList<Satellite>> satMap = new HashMap<Integer, ArrayList<Satellite>>();
			while ((line = br.readLine()) != null) {
				double[] data = Arrays.stream(line.split("\\s+")).mapToDouble(i -> Double.parseDouble(i)).toArray();
				satMap.computeIfAbsent((int) data[0], k -> new ArrayList<Satellite>()).add(new Satellite(data));
			}
			fileName = "refpos.txt";
			is = MainApp.class.getClassLoader().getResourceAsStream(fileName);
			br = new BufferedReader(new InputStreamReader(is));
			int estimatorType = 3;
			while ((line = br.readLine()) != null) {
				double[] data = Arrays.stream(line.split("\\s+")).mapToDouble(i -> Double.parseDouble(i)).toArray();
				int t = (int) data[0];
				double[] rxECEF = Arrays.copyOfRange(data, 1, 4);
				ArrayList<Satellite> satList = satMap.get(t);
				double[] estEcefClk = null;
				switch (estimatorType) {
				case 1:
					estEcefClk = LinearLeastSquare.process(satList, false);
					ErrMap.computeIfAbsent("LS", k -> new ArrayList<double[]>()).add(estimateError(estEcefClk, rxECEF));
					EnuMap.computeIfAbsent("LS", k -> new ArrayList<double[]>()).add(estimateENU(estEcefClk, rxECEF));
					break;
				case 2:
					estEcefClk = LinearLeastSquare.process(satList, true);
					ErrMap.computeIfAbsent("WLS", k -> new ArrayList<double[]>())
							.add(estimateError(estEcefClk, rxECEF));
					EnuMap.computeIfAbsent("WLS", k -> new ArrayList<double[]>()).add(estimateENU(estEcefClk, rxECEF));

					break;
				case 3:
					estEcefClk = LinearLeastSquare.process(satList, false);
					ErrMap.computeIfAbsent("LS", k -> new ArrayList<double[]>()).add(estimateError(estEcefClk, rxECEF));
					EnuMap.computeIfAbsent("LS", k -> new ArrayList<double[]>()).add(estimateENU(estEcefClk, rxECEF));

					estEcefClk = LinearLeastSquare.process(satList, true);
					ErrMap.computeIfAbsent("WLS", k -> new ArrayList<double[]>())
							.add(estimateError(estEcefClk, rxECEF));
					EnuMap.computeIfAbsent("WLS", k -> new ArrayList<double[]>()).add(estimateENU(estEcefClk, rxECEF));

					break;

				}
				rxEcefList.add(rxECEF);
				timeList.add(t);

			}
			if (estimatorType == 4 || estimatorType == 3) {
				EKF ekf = new EKF();
				double[] initialECEF = LinearLeastSquare.process(satMap.get(timeList.get(0)), false);
				ArrayList<double[]> estEcefList_pos = ekf.process(satMap, timeList, initialECEF, Flag.POSITION);
				ArrayList<double[]> estEcefList_vel = ekf.process(satMap, timeList, initialECEF, Flag.VELOCITY);
				for (int i = 1; i < timeList.size(); i++) {
					ErrMap.computeIfAbsent("EKF - pos. random walk", k -> new ArrayList<double[]>())
							.add(estimateError(estEcefList_pos.get(i - 1), rxEcefList.get(i)));
					EnuMap.computeIfAbsent("EKF - pos. random walk", k -> new ArrayList<double[]>())
							.add(estimateENU(estEcefList_pos.get(i - 1), rxEcefList.get(i)));

					ErrMap.computeIfAbsent("EKF - vel. random walk", k -> new ArrayList<double[]>())
							.add(estimateError(estEcefList_vel.get(i - 1), rxEcefList.get(i)));
					EnuMap.computeIfAbsent("EKF - vel. random walk", k -> new ArrayList<double[]>())
							.add(estimateENU(estEcefList_vel.get(i - 1), rxEcefList.get(i)));
				}
			}
			HashMap<String, ArrayList<Double>> GraphErrMap = new HashMap<String, ArrayList<Double>>();
			for (String key : ErrMap.keySet()) {

				ArrayList<double[]> errList = ErrMap.get(key);
				ArrayList<Double> ecefErrList = new ArrayList<Double>();
				ArrayList<Double> llErrList = new ArrayList<Double>();
				for (int i = 0; i < errList.size(); i++) {
					ecefErrList.add(errList.get(i)[0]);
					llErrList.add(errList.get(i)[1]);
				}
				String header = "order - ";
				String minECEF = "ECEF - ";
				String minLL = "LL - ";
				String rmsECEF = "ECEF - ";
				String maeECEF = "ECEF - ";
				String rmsLL = "LL - ";
				String maeLL = "LL - ";

				minECEF += Collections.min(ecefErrList) + " ";
				minLL += Collections.min(llErrList) + " ";
				rmsECEF += RMS(ecefErrList) + " ";
				rmsLL += RMS(llErrList) + " ";
				maeECEF += MAE(ecefErrList) + " ";
				maeLL += MAE(llErrList) + " ";
				GraphErrMap.put(key + " ECEF off", ecefErrList);
				GraphErrMap.put(key + " LL off", llErrList);

				System.out.println("\n" + key);
				System.out.println(header);
				System.out.println("MIN - ");
				System.out.println(minECEF);
				System.out.println(minLL);
				System.out.println("RMS - ");
				System.out.println(rmsECEF);
				System.out.println(rmsLL);
				System.out.println("MAE - ");
				System.out.println(maeECEF);
				System.out.println(maeLL);

			}
			HashMap<String, ArrayList<double[]>> GraphEnuMap = new HashMap<String, ArrayList<double[]>>();
			for (String key : EnuMap.keySet()) {

				ArrayList<double[]> enuList = EnuMap.get(key);
				int n = enuList.size();
				double err = 0;
				double horizontalErr = 0;
				double eErr = 0;
				double nErr = 0;
				double uErr = 0;

				for (int i = 0; i < n; i++) {
					double[] enu = enuList.get(i);
					eErr += enu[0] * enu[0];
					nErr += enu[1] * enu[1];
					uErr += enu[2] * enu[2];
					err += Arrays.stream(enu).map(j -> j * j).sum();
					horizontalErr += (enu[0] * enu[0]) + (enu[1] * enu[1]);
				}
				eErr = Math.sqrt(eErr / n);
				nErr = Math.sqrt(nErr / n);
				uErr = Math.sqrt(uErr / n);
				err = Math.sqrt(err / n);
				horizontalErr = Math.sqrt(horizontalErr / n);

				GraphEnuMap.put(key, enuList);

				System.out.println("\n" + key);
				System.out.println("RMS - ");
				System.out.println(" E - " + eErr);
				System.out.println(" N - " + nErr);
				System.out.println(" U - " + uErr);
				System.out.println(" 3d Error - " + err);
				System.out.println(" 2d Error - " + horizontalErr);

			}

			GraphPlotter chart = new GraphPlotter("GPS PVT Error - ", "Error Estimate(m)", GraphErrMap, timeList);
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);

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
