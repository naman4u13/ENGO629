package com.ENGO629.models;

public class Satellite {

	private int t;
	private int prn;
	private double[] ecef;

	private double pseduorange;
	private double phase;
	private double elevation;

	public Satellite(int t, int prn, double x, double y, double z, double pseduorange, double phase, double elevation) {
		super();
		this.t = t;
		this.prn = prn;
		this.ecef = new double[] { x, y, z };
		this.pseduorange = pseduorange;
		this.phase = phase;
		this.elevation = elevation;
	}

	public Satellite(double[] data) {
		this.t = (int) data[0];
		this.prn = (int) data[1];
		this.ecef = new double[] { data[2], data[3], data[4] };
		this.pseduorange = data[5];
		this.phase = data[6];
		this.elevation = data[7];

	}

	public int getT() {
		return t;
	}

	public int getPrn() {
		return prn;
	}

	public double[] getECEF() {
		return ecef;
	}

	public double getPseduorange() {
		return pseduorange;
	}

	public double getPhase() {
		return phase;
	}

	public double getElevation() {
		return elevation;
	}

}
