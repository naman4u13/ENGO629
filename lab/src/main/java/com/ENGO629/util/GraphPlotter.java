package com.ENGO629.util;

import java.util.ArrayList;
import java.util.HashMap;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ui.ApplicationFrame;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class GraphPlotter extends ApplicationFrame {

	public GraphPlotter(String applicationTitle, String chartTitle, HashMap<String, ArrayList<Double>> dataMap,
			ArrayList<Integer> timeList) {
		super(applicationTitle);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(chartTitle, "GPS-time", chartTitle,
				createDatasetCS(dataMap, timeList));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);

		setContentPane(chartPanel);

	}

	private XYDataset createDatasetCS(HashMap<String, ArrayList<Double>> dataMap, ArrayList<Integer> timeList) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (String key : dataMap.keySet()) {
			final XYSeries series = new XYSeries(key);
			ArrayList<Double> data = dataMap.get(key);
			for (int i = 0; i < data.size(); i++) {
				series.add(timeList.get(i), data.get(i));
			}
			dataset.addSeries(series);
		}

		return dataset;

	}
}
