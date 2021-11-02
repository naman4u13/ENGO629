package com.ENGO629.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ui.ApplicationFrame;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

public class GraphPlotter extends ApplicationFrame {
	String path = "D:\\projects\\eclipse_projects\\UCalgary\\ENGO629\\results\\";

	public GraphPlotter(String applicationTitle, String chartTitle, HashMap<String, ArrayList<double[]>> dataMap,
			ArrayList<Integer> timeList) throws IOException {
		super(applicationTitle);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(chartTitle, "GPS-time", chartTitle,
				createDataset3dErr(dataMap, timeList));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String applicationTitle, String chartTitle, HashMap<String, double[]> dataMap,
			ArrayList<Integer> timeList, boolean flag) throws IOException {
		super(applicationTitle);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(chartTitle, "GPS-time", chartTitle,
				createDatasetENU(dataMap, timeList));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String applicationTitle, HashMap<String, ArrayList<double[]>> dataMap) throws IOException {
		super(applicationTitle);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createScatterPlot("2D Error", "East(m)", "North(m)",
				createDataset2dErr(dataMap));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + "2d Error" + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public static void graphENU(HashMap<String, ArrayList<double[]>> dataMap, ArrayList<Integer> timeList)
			throws IOException {

		String[] chartNames = new String[] { "E", "N", "U" };
		for (int i = 0; i < 3; i++) {
			final int index = i;
			HashMap<String, double[]> subDataMap = new HashMap<String, double[]>();
			for (String key : dataMap.keySet()) {
				ArrayList<double[]> data = dataMap.get(key);
				double[] arr = data.stream().mapToDouble(j -> j[index]).toArray();
				subDataMap.put(key, arr);
			}
			GraphPlotter chart = new GraphPlotter("GPS PVT Error - ", chartNames[i] + "(m)", subDataMap, timeList,
					true);
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
		}
		GraphPlotter chart = new GraphPlotter("2D-Error", dataMap);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter("3d-Error", "3d-Error", dataMap, timeList);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

	}

	private XYDataset createDataset3dErr(HashMap<String, ArrayList<double[]>> dataMap, ArrayList<Integer> timeList) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (String key : dataMap.keySet()) {
			final XYSeries series = new XYSeries(key);
			ArrayList<double[]> list = dataMap.get(key);
			for (int i = 0; i < list.size(); i++) {
				double[] data = list.get(i);
				double err = Math.sqrt(Arrays.stream(data).map(j -> j * j).sum());
				series.add(timeList.get(i), Double.valueOf(err));
			}
			dataset.addSeries(series);
		}

		return dataset;

	}

	private XYDataset createDatasetENU(HashMap<String, double[]> dataMap, ArrayList<Integer> timeList) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (String key : dataMap.keySet()) {
			final XYSeries series = new XYSeries(key);
			double[] data = dataMap.get(key);
			for (int i = 0; i < data.length; i++) {
				series.add(timeList.get(i), Double.valueOf(data[i]));
			}
			dataset.addSeries(series);
		}

		return dataset;

	}

	private XYDataset createDataset2dErr(HashMap<String, ArrayList<double[]>> dataMap) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (String key : dataMap.keySet()) {
			final XYSeries series = new XYSeries(key);
			ArrayList<double[]> list = dataMap.get(key);
			for (int i = 0; i < list.size(); i++) {
				double[] data = list.get(i);
				series.add(data[0], data[1]);
			}
			dataset.addSeries(series);
		}

		return dataset;

	}
}
