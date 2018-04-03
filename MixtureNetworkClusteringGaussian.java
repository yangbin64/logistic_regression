package mixtureNetwork;

import java.io.*;
import java.util.*;

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.special.*;

public class MixtureNetworkClusteringGaussian<TID> {
	
	String filenameInput = "E:\\RR\\MixtureNetwork\\Correlation.txt";
	String pathnameOutput = "E:\\RR\\MixtureNetwork\\";
	
	int CLUSTER_NUM = 20;
	double TAU0 = 1;
	double TAU = 256;
	double GAMMA = 0.02;
	
	HashMap<TID, HashMap<TID, Double>> mapWeightedNetwork = new HashMap<TID, HashMap<TID, Double>>();

	HashMap<TID, Integer> mapLabel = new HashMap<TID, Integer>();
	int[] n = new int[CLUSTER_NUM];
	double[][] s = new double[CLUSTER_NUM][CLUSTER_NUM];
	double[][] s2 = new double[CLUSTER_NUM][CLUSTER_NUM];
	
	public MixtureNetworkClusteringGaussian() {
		
	}
	
	public MixtureNetworkClusteringGaussian(String fnInput, String pnOutput, int clusternum) {
		filenameInput = fnInput;
		pathnameOutput = pnOutput;
		CLUSTER_NUM = clusternum;
		n = new int[CLUSTER_NUM];
		s = new double[CLUSTER_NUM][CLUSTER_NUM];
		s2 = new double[CLUSTER_NUM][CLUSTER_NUM];
	}
	
	public void parse(String[] args) {
		if (args.length == 6) {
			System.out.println("Read parameters!");
			filenameInput = args[0];
			pathnameOutput = args[1];
			CLUSTER_NUM = Integer.valueOf(args[2]);
			TAU0 = Double.valueOf(args[3]);
			TAU = Double.valueOf(args[4]);
			GAMMA = Double.valueOf(args[5]);
		} else {
			System.out.println("Use default parameters!");
		}
		
		System.out.println("filenameInput=" + filenameInput);
		System.out.println("pathnameOutput=" + pathnameOutput);
		System.out.println("CLUSTER_NUM=" + CLUSTER_NUM);
		System.out.println("TAU0=" + TAU0);
		System.out.println("TAU=" + TAU);
		System.out.println("GAMMA=" + GAMMA);
	}
	
	public TID conv(String s) {
		return (TID)Integer.valueOf(s);
	}
	
	public double convertScore(double d) {
		//double s = 0.01;
		//double dd = (d+1+s)/(2+s*2);
		//return -Math.log(1.0/dd-1);
		return d;
	}
	
	public void read() {
		UniformIntegerDistribution uni = new UniformIntegerDistribution(0, CLUSTER_NUM-1);
		try {
			FileReader fr = new FileReader(filenameInput);
			BufferedReader buf = new BufferedReader(fr);
			String s;
			
			int count = 0;
			while ((s=buf.readLine())!=null) {
				count++;
				if (count%10000000==0) {
					System.out.println("" + count + " records are read.");
				}
				String[] ss=s.split("\t");
				TID id1 = conv(ss[0]);
				TID id2 = conv(ss[1]);
				double sc = convertScore(Double.valueOf(ss[2]));
				if (!mapWeightedNetwork.containsKey(id1)) {
					mapWeightedNetwork.put(id1, new HashMap<TID, Double>());
				}
				mapWeightedNetwork.get(id1).put(id2, sc);
				// init labels
				if (!mapLabel.containsKey(id1)) {
					int randomLabel = uni.sample();
					mapLabel.put(id1, randomLabel);
					n[randomLabel]++;
				}
			}
			
			fr.close();
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public void read2D() {
		UniformIntegerDistribution uni = new UniformIntegerDistribution(0, CLUSTER_NUM-1);
		try {
			FileReader fr = new FileReader(filenameInput);
			BufferedReader buf = new BufferedReader(fr);
			String s;
			
			ArrayList<TID> IDList = new ArrayList<TID>();
			
			boolean bFirstLine = true;
			while ((s=buf.readLine())!=null) {
				String[] ss=s.split("\t");
				
				if (bFirstLine) {
					for (int i=1; i<ss.length; i++) {
						TID id1 = conv(ss[i]);
						mapWeightedNetwork.put(id1, new HashMap<TID, Double>());
						IDList.add(id1);
					}
					bFirstLine = false;
				} else {
					TID id2 = conv(ss[0]);
					for (int i=1; i<ss.length; i++) {
						double sc = convertScore(Double.valueOf(ss[i]));
						TID id1 = IDList.get(i-1);
						mapWeightedNetwork.get(id1).put(id2, sc);
						// init labels
						if (!mapLabel.containsKey(id1)) {
							int randomLabel = uni.sample();
							mapLabel.put(id1, randomLabel);
							n[randomLabel]++;
						}
					}
				}
			}
			
			fr.close();
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public void initNetwork() {
		read2D();
		
		// init s and s2
		for (Iterator<Map.Entry<TID, HashMap<TID, Double>>> it1=mapWeightedNetwork.entrySet().iterator(); it1.hasNext();) {
			Map.Entry<TID, HashMap<TID, Double>> entry1 = it1.next();
			TID id1 = entry1.getKey();
			int label1 = mapLabel.get(id1);
			HashMap<TID, Double> scList = entry1.getValue();
			
			for (Iterator<Map.Entry<TID, Double>> it2=scList.entrySet().iterator(); it2.hasNext();) {
				Map.Entry<TID, Double> entry2 = it2.next();
				TID id2 = entry2.getKey();
				int label2 = mapLabel.get(id2);
				double sc = entry2.getValue();
				s[label1][label2] += sc;
				s2[label1][label2] += (sc*sc);
			}
		}
	}
	
	public double getLogGaussian(double snew, double s2new, int nnew, double s, double s2, int n) {
		double logGaussiannew = -TAU/2*(s2new-snew*snew/(nnew+TAU0/TAU));
		double logGaussian = -TAU/2*(s2-s*s/(n+TAU0/TAU));
		return logGaussiannew - logGaussian;
	}
	
	public void check() {
		int[] nc = new int[CLUSTER_NUM];
		int[][] sc = new int[CLUSTER_NUM][CLUSTER_NUM];
		int[][] s2c = new int[CLUSTER_NUM][CLUSTER_NUM];
		
		// compute nc
		for (Iterator<Map.Entry<TID, Integer>> it=mapLabel.entrySet().iterator(); it.hasNext();) {
			Map.Entry<TID, Integer> entry = it.next();
			int label = entry.getValue();
			nc[label]++;
		}
		
		// compute mc and Mc
		for (Iterator<Map.Entry<TID, HashMap<TID, Double>>> it1=mapWeightedNetwork.entrySet().iterator(); it1.hasNext();) {
			Map.Entry<TID, HashMap<TID, Double>> entry1 = it1.next();
			TID id1 = entry1.getKey();
			int label1 = mapLabel.get(id1);
			HashMap<TID, Double> setEntry2 = entry1.getValue();
			for (Iterator<Map.Entry<TID, Double>> it2=setEntry2.entrySet().iterator(); it2.hasNext();) {
				Map.Entry<TID, Double> entry2 = it2.next();
				TID id2 = entry2.getKey();
				int label2 = mapLabel.get(id2);
				double score = entry2.getValue();
				sc[label1][label2] += score;
				s2c[label1][label2] += (score*score);
			}
		}
		
		// compare
		for (int i=0; i<CLUSTER_NUM; i++) {
			if (nc[i]!=n[i]) {
				System.out.println("" + nc[i] + "!=" + n[i] + ", ");
				System.out.println("NOT match!");
			} else {
				System.out.println("" + nc[i] + "==" + n[i] + ", ");
			}
			
			for (int j=0; j<CLUSTER_NUM; j++) {
				if (sc[i][j] == s[i][j]) {
					System.out.print("[" + sc[i][j] + "==" + s[i][j] + "]");
				} else {
					System.out.print("[" + sc[i][j] + "!=" + s[i][j] + "]");
					System.out.println("NOT match!");
				}
			}
			System.out.println();
			
			for (int j=0; j<CLUSTER_NUM; j++) {
				if (s2c[i][j] == s2[i][j]) {
					System.out.print("[" + s2c[i][j] + "==" + s2[i][j] + "]");
				} else {
					System.out.print("[" + s2c[i][j] + "!=" + s2[i][j] + "]");
					System.out.println("NOT match!");
				}
			}
			System.out.println();
		}
	}
	
	public void clustering() {
		int changed = 0;
		int unchanged = 0;

		for (Iterator<Map.Entry<TID, Integer>> it=mapLabel.entrySet().iterator(); it.hasNext();) {
			
			Map.Entry<TID, Integer> entry = it.next();
			TID id = entry.getKey();
			int z = entry.getValue();
			
			double[] sum = new double[CLUSTER_NUM];
			double[] sum2 = new double[CLUSTER_NUM];
			for (Iterator<Map.Entry<TID, Double>> itn=mapWeightedNetwork.get(id).entrySet().iterator(); itn.hasNext();) {
				Map.Entry<TID, Double> entryn = itn.next();
				TID idn = entryn.getKey();
				int zz = mapLabel.get(idn);
				double score = entryn.getValue();
				sum[zz] += score;
				sum2[zz] += (score*score);
			}
			
			int[] singletons = new int[CLUSTER_NUM];
			double[] probabilities = new double[CLUSTER_NUM];
			double[] logProb = new double[CLUSTER_NUM];
			for (int z_new=0; z_new<CLUSTER_NUM; z_new++) {
				singletons[z_new] = z_new;
				if (z_new==z) {
					logProb[z_new] = Math.log(n[z_new]-1+GAMMA);
				} else {
					logProb[z_new] = Math.log(n[z_new]+GAMMA);
					for (int a=0; a<CLUSTER_NUM; a++) {
						if (a==z) {
							// [z][z]
							double w = mapWeightedNetwork.get(a).get(a);
							double snew = s[a][a]-sum[a]*2 + w;
							double s2new = s2[a][a]-sum2[a]*2 + w*w;
							logProb[z_new] += getLogGaussian(snew, s2new, (n[a]-1)*(n[a]-1), s[a][a], s2[a][a], n[a]*n[a]);
							if (Double.isNaN(logProb[z_new])) {
								System.out.println("NaN");
							}
						} else if (a==z_new) {
							// [z_new][z_new]
							double w = mapWeightedNetwork.get(a).get(a);
							double snew = s[a][a]+sum[a]*2 + w;
							double s2new = s2[a][a]+sum2[a]*2 + w*w;
							logProb[z_new] += getLogGaussian(snew, s2new, (n[a]+1)*(n[a]+1), s[a][a], s2[a][a], n[a]*n[a]);
							if (Double.isNaN(logProb[z_new])) {
								System.out.println("NaN");
							}
							// [z][z_new] and [z_new][z]
							double snew1 = s[z][z_new]-sum[z_new]+sum[z] - w;
							double s2new1 = s2[z][z_new]-sum2[z_new]+sum2[z] - w*w;
							logProb[z_new] += getLogGaussian(snew1, s2new1, (n[z_new]+1)*(n[z]-1), s[z][z_new], s2[z][z_new], n[z_new]*n[z])*2;
							if (Double.isNaN(logProb[z_new])) {
								System.out.println("NaN");
							}
						} else {
							// [z][a] and [a][z]
							double snew1 = s[z][a]-sum[a];
							double s2new1 = s2[z][a]-sum2[a];
							logProb[z_new] += getLogGaussian(snew1, s2new1, n[a]*(n[z]-1), s[z][a], s2[z][a], n[a]*n[z])*2;
							if (Double.isNaN(logProb[z_new])) {
								System.out.println("NaN");
							}
							// [z_new][a] and [a][z_new]
							double snew2 = s[z_new][a]+sum[a];
							double s2new2 = s2[z_new][a]+sum2[a];
							logProb[z_new] += getLogGaussian(snew2, s2new2, n[a]*(n[z_new]+1), s[z_new][a], s2[z_new][a], n[a]*n[z_new])*2;
							if (Double.isNaN(logProb[z_new])) {
								System.out.println("NaN");
							}
						}
					}
				}
			}
			double maxLogProb = 0;
			for (int z_new=0; z_new<CLUSTER_NUM; z_new++) {
				if (logProb[z_new]>maxLogProb) {
					maxLogProb = logProb[z_new];
				}
			}
			if (maxLogProb>200) {
				double d = maxLogProb-200;
				for (int z_new=0; z_new<CLUSTER_NUM; z_new++) {
					logProb[z_new] -= d;
				}
			}
			for (int z_new=0; z_new<CLUSTER_NUM; z_new++) {
				probabilities[z_new] = Math.exp(logProb[z_new]);
			}
			EnumeratedIntegerDistribution discrete = new EnumeratedIntegerDistribution(singletons, probabilities);
			int z_new = discrete.sample();
			
			// Update
			if (z_new!=z) {
				double[] sumn = new double[CLUSTER_NUM];
				double[] sum2n = new double[CLUSTER_NUM];
				for (Iterator<Map.Entry<TID, Double>> itn=mapWeightedNetwork.get(id).entrySet().iterator(); itn.hasNext();) {
					Map.Entry<TID, Double> entryn = itn.next();
					TID idn = entryn.getKey();
					int zz = mapLabel.get(idn);
					double score = entryn.getValue();
					sumn[zz] += score;
					sum2n[zz] += (score*score);
				}
				for (int a=0; a<CLUSTER_NUM; a++) {
					if (a==z) {
						// [z][z]
						double w = mapWeightedNetwork.get(a).get(a);
						s[a][a] = s[a][a]-sumn[a]*2 + w;
						s2[a][a] = s2[a][a]-sumn[a]*2 + w*w;
					} else if (a==z_new) {
						// [z_new][z_new]
						double w = mapWeightedNetwork.get(a).get(a);
						s[a][a] = s[a][a]+sumn[a]*2 + w;
						s2[a][a] = s2[a][a]+sum2n[a]*2 + w*w;
						// [z][z_new] and [z_new][z]
						s[z][z_new] = s[z][z_new]-sumn[z_new]+sumn[z] - w;
						s2[z][z_new] = s2[z][z_new]-sum2n[z_new]+sum2n[z] - w*w;
						s[z_new][z] = s[z_new][z]-sumn[z_new]+sumn[z] - w;
						s2[z_new][z] = s2[z_new][z]-sum2n[z_new]+sum2n[z] - w*w;
					} else {
						// [z][a] and [a][z]
						s[z][a] = s[z][a]-sumn[a];
						s2[z][a] = s2[z][a]-sum2n[a];
						s[a][z] = s[a][z]-sumn[a];
						s2[a][z] = s2[a][z]-sum2n[a];
						// [z_new][a] and [a][z_new]
						s[z_new][a] = s[z_new][a]+sumn[a];
						s2[z_new][a] = s2[z_new][a]+sum2n[a];
						s[a][z_new] = s[a][z_new]+sumn[a];
						s2[a][z_new] = s2[a][z_new]+sum2n[a];
					}
				}
				mapLabel.put(id, z_new);
				n[z]--;
				n[z_new]++;
				
				//System.out.println("" + z + "=>" + z_new);
				
				//check();
			}
			
			if (z_new==z) {
				unchanged++;
			} else {
				changed++;
			}
		}
		System.out.println("Unchanged:" + unchanged + "   Changed:" + changed);
	}
	
	public void output(int r) {
		String filenameZ = pathnameOutput + "z_" + r + ".out";
		String filenamem = pathnameOutput + "m_" + r + ".out";
		String filenameM = pathnameOutput + "mm_" + r + ".out";
		String filenamemCheck = pathnameOutput + "m_check_" + r + ".out";
		String filenameMCheck = pathnameOutput + "mm_check_" + r + ".out";
		
		try {
			FileWriter fwZ = new FileWriter(filenameZ);
			for (Iterator<Map.Entry<TID, Integer>> it=mapLabel.entrySet().iterator(); it.hasNext();) {
				Map.Entry<TID, Integer> entry = it.next();
				TID id = entry.getKey();
				int label = entry.getValue();
				//fwZ.write(id + "," + label + "\n");
				fwZ.write(id + "\t" + label + "\n");
			}
			fwZ.close();
			
			FileWriter fwm = new FileWriter(filenamem);
			//FileWriter fwM = new FileWriter(filenameM);
			for (int i=0; i<CLUSTER_NUM; i++) {
				fwm.write("\t" + i);
				//fwM.write("\t" + i);
			}
			fwm.write("\n");
			//fwM.write("\n");
			for (int i=0; i<CLUSTER_NUM; i++) {
				fwm.write("" + i);
				//fwM.write("" + i);
				for (int j=0; j<CLUSTER_NUM; j++) {
					fwm.write("\t" + (s[i][j]/n[i]/n[j]));
					//fwM.write("\t" + s2[i][j]);
				}
				fwm.write("\n");
				//fwM.write("\n");
			}
			fwm.close();
			//fwM.close();
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public void outputNetwork(int r) {
		String filenameNetwork = pathnameOutput + "network_" + r + ".out";
		
		try {
			FileWriter fwNetwork = new FileWriter(filenameNetwork);
			for (int i=0; i<CLUSTER_NUM; i++) {
				for (Iterator<Map.Entry<TID, Integer>> it1=mapLabel.entrySet().iterator(); it1.hasNext();) {
					Map.Entry<TID, Integer> entry1 = it1.next();
					//TID id1 = entry1.getKey();
					int label1 = entry1.getValue();
					if (label1==i) {
						for (int j=0; j<CLUSTER_NUM; j++) {
							for (Iterator<Map.Entry<TID, Integer>> it2=mapLabel.entrySet().iterator(); it2.hasNext();) {
								Map.Entry<TID, Integer> entry2 = it2.next();
								//TID id2 = entry2.getKey();
								int label2 = entry2.getValue();
								if (label2==j) {
									fwNetwork.write("" + mapWeightedNetwork.get(label1).get(label2) + ",");
								}
							}
						}
						fwNetwork.write("\n");
					}
				}
			}
			fwNetwork.close();
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		MixtureNetworkClusteringGaussian<String> mnc = new MixtureNetworkClusteringGaussian<String>();
		mnc.parse(args);
		mnc.initNetwork();
		mnc.output(0);
		mnc.outputNetwork(0);
		int r=1;
		for (r=1; r<1000; r++) {
			System.out.println("Starting round " + r + "...");
			mnc.clustering();
		}
		mnc.output(r);
		mnc.outputNetwork(r);
	}

}
