package mixtureNetwork;

import java.io.*;
import java.util.*;

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.special.*;

public class MixtureNetworkClustering<TID> {
	
	int CLUSTER_NUM = 10;
	final double ALPHA = 0.1;
	final double BETA = 0.1;
	final double GAMMA = 0.1;
	
	String filenameInput = "E:\\Brasil\\Genre\\03JaccardRank.out";
	String pathnameOutput = "E:\\Brasil\\ClusteringGenre\\";
	
	HashMap<TID, HashSet<TID>> mapBinaryNetwork = new HashMap<TID, HashSet<TID>>();

	HashMap<TID, Integer> mapLabel = new HashMap<TID, Integer>();
	int[] n = new int[CLUSTER_NUM];
	int[][] m = new int[CLUSTER_NUM][CLUSTER_NUM];
	int[][] M = new int[CLUSTER_NUM][CLUSTER_NUM];
	
	public MixtureNetworkClustering() {
		
	}
	
	public MixtureNetworkClustering(String fnInput, String pnOutput, int clusternum) {
		filenameInput = fnInput;
		pathnameOutput = pnOutput;
		CLUSTER_NUM = clusternum;
		n = new int[CLUSTER_NUM];
		m = new int[CLUSTER_NUM][CLUSTER_NUM];
		M = new int[CLUSTER_NUM][CLUSTER_NUM];
	}
	
	public TID conv(String s) {
		return (TID)Integer.valueOf(s);
	}
	
	public void initNetwork() {
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
				if (!mapBinaryNetwork.containsKey(id1)) {
					mapBinaryNetwork.put(id1, new HashSet<TID>());
				}
				mapBinaryNetwork.get(id1).add(id2);
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
		
		// init m and M
		for (Iterator<Map.Entry<TID, HashSet<TID>>> it1=mapBinaryNetwork.entrySet().iterator(); it1.hasNext();) {
			Map.Entry<TID, HashSet<TID>> entry1 = it1.next();
			TID id1 = entry1.getKey();
			int label1 = mapLabel.get(id1);
			HashSet<TID> setEntry2 = entry1.getValue();
			for (Iterator<TID> it2=setEntry2.iterator(); it2.hasNext();) {
				TID id2 = it2.next();
				int label2 = mapLabel.get(id2);
				m[label1][label2]++;
			}
		}
		
		for (int i=0; i<CLUSTER_NUM; i++) {
			for (int j=0; j<CLUSTER_NUM; j++) {
				M[i][j] = n[i]*n[j]-m[i][j];
			}
		}
	}
	
	public double getLogBeta(int mnew, int Mnew, int m, int M) {
		double logBeta = Beta.logBeta(mnew+ALPHA, Mnew+BETA);
		logBeta -= Beta.logBeta(m+ALPHA, M+BETA);
		return logBeta;
	}
	
	public void check() {
		int[] nc = new int[CLUSTER_NUM];
		int[][] mc = new int[CLUSTER_NUM][CLUSTER_NUM];
		int[][] Mc = new int[CLUSTER_NUM][CLUSTER_NUM];
		
		// compute nc
		for (Iterator<Map.Entry<TID, Integer>> it=mapLabel.entrySet().iterator(); it.hasNext();) {
			Map.Entry<TID, Integer> entry = it.next();
			int label = entry.getValue();
			nc[label]++;
		}
		
		// compute mc and Mc
		for (Iterator<Map.Entry<TID, HashSet<TID>>> it1=mapBinaryNetwork.entrySet().iterator(); it1.hasNext();) {
			Map.Entry<TID, HashSet<TID>> entry1 = it1.next();
			TID id1 = entry1.getKey();
			int label1 = mapLabel.get(id1);
			HashSet<TID> setEntry2 = entry1.getValue();
			for (Iterator<TID> it2=setEntry2.iterator(); it2.hasNext();) {
				TID id2 = it2.next();
				int label2 = mapLabel.get(id2);
				mc[label1][label2]++;
			}
		}
		
		for (int i=0; i<CLUSTER_NUM; i++) {
			for (int j=0; j<CLUSTER_NUM; j++) {
				Mc[i][j] = nc[i]*nc[j]-mc[i][j];
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
				if (mc[i][j] == m[i][j]) {
					System.out.print("[" + mc[i][j] + "==" + m[i][j] + "]");
				} else {
					System.out.print("[" + mc[i][j] + "!=" + m[i][j] + "]");
					System.out.println("NOT match!");
				}
			}
			System.out.println();
			
			for (int j=0; j<CLUSTER_NUM; j++) {
				if (Mc[i][j] == M[i][j]) {
					System.out.print("[" + Mc[i][j] + "==" + M[i][j] + "]");
				} else {
					System.out.print("[" + Mc[i][j] + "!=" + M[i][j] + "]");
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
			
			int[] counts = new int[CLUSTER_NUM];
			for (Iterator<TID> itn=mapBinaryNetwork.get(id).iterator(); itn.hasNext();) {
				TID idn = itn.next();
				int zz = mapLabel.get(idn);
				counts[zz]++;
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
							int mm = m[a][a]-counts[a]*2+1;
							int MM = M[a][a]-(n[a]-counts[a])*2;
							logProb[z_new] += getLogBeta(mm, MM, m[a][a], M[a][a]);
							if (Double.isNaN(logProb[z_new])) {
								System.out.println("NaN");
							}
						} else if (a==z_new) {
							// [z_new][z_new]
							int m1 = m[a][a]+counts[a]*2+1;
							int M1 = M[a][a]+(n[a]-counts[a])*2;
							logProb[z_new] += getLogBeta(m1, M1, m[a][a], M[a][a]);
							if (Double.isNaN(logProb[z_new])) {
								System.out.println("NaN");
							}
							// [z][z_new] and [z_new][z]
							int m2 = m[z][z_new]-counts[z_new]+counts[z]-1;
							int M2 = M[z][z_new]-(n[z_new]-counts[z_new])+(n[z]-counts[z]);
							logProb[z_new] += getLogBeta(m2, M2, m[z][z_new], M[z][z_new])*2;
							if (Double.isNaN(logProb[z_new])) {
								System.out.println("NaN");
							}
						} else {
							// [z][a] and [a][z]
							int m1 = m[z][a]-counts[a];
							int M1 = M[z][a]-(n[a]-counts[a]);
							logProb[z_new] += getLogBeta(m1, M1, m[z][a], M[z][a])*2;
							if (Double.isNaN(logProb[z_new])) {
								System.out.println("NaN");
							}
							// [z_new][a] and [a][z_new]
							int m2 = m[z_new][a]+counts[a];
							int M2 = M[z_new][a]+(n[a]-counts[a]);
							logProb[z_new] += getLogBeta(m2, M2, m[z_new][a], M[z_new][a])*2;
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
				int[] countsn = new int[CLUSTER_NUM];
				for (Iterator<TID> itn=mapBinaryNetwork.get(id).iterator(); itn.hasNext();) {
					TID idn = itn.next();
					int zz = mapLabel.get(idn);
					countsn[zz]++;
				}
				for (int a=0; a<CLUSTER_NUM; a++) {
					if (a==z) {
						// [z][z]
						m[a][a] = m[a][a]-countsn[a]*2+1;
						M[a][a] = M[a][a]-(n[a]-countsn[a])*2;
					} else if (a==z_new) {
						// [z_new][z_new]
						m[a][a] = m[a][a]+countsn[a]*2+1;
						M[a][a] = M[a][a]+(n[a]-countsn[a])*2;
						// [z][z_new] and [z_new][z]
						m[z][z_new] = m[z][z_new]-countsn[z_new]+countsn[z]-1;
						M[z][z_new] = M[z][z_new]-(n[z_new]-countsn[z_new])+(n[z]-countsn[z]);
						m[z_new][z] = m[z_new][z]-countsn[z_new]+countsn[z]-1;
						M[z_new][z] = M[z_new][z]-(n[z_new]-countsn[z_new])+(n[z]-countsn[z]);
					} else {
						// [z][a] and [a][z]
						m[z][a] = m[z][a]-countsn[a];
						M[z][a] = M[z][a]-(n[a]-countsn[a]);
						m[a][z] = m[a][z]-countsn[a];
						M[a][z] = M[a][z]-(n[a]-countsn[a]);
						// [z_new][a] and [a][z_new]
						m[z_new][a] = m[z_new][a]+countsn[a];
						M[z_new][a] = M[z_new][a]+(n[a]-countsn[a]);
						m[a][z_new] = m[a][z_new]+countsn[a];
						M[a][z_new] = M[a][z_new]+(n[a]-countsn[a]);
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
				fwZ.write(id + "," + label + "\n");
			}
			fwZ.close();
			
			FileWriter fwm = new FileWriter(filenamem);
			FileWriter fwM = new FileWriter(filenameM);
			for (int i=0; i<CLUSTER_NUM; i++) {
				for (int j=0; j<CLUSTER_NUM; j++) {
					fwm.write(m[i][j] + ",");
					fwM.write(M[i][j] + ",");
				}
				fwm.write("\n");
				fwM.write("\n");
			}
			fwm.close();
			fwM.close();
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
					TID id1 = entry1.getKey();
					int label1 = entry1.getValue();
					if (label1==i) {
						for (int j=0; j<CLUSTER_NUM; j++) {
							for (Iterator<Map.Entry<TID, Integer>> it2=mapLabel.entrySet().iterator(); it2.hasNext();) {
								Map.Entry<TID, Integer> entry2 = it2.next();
								TID id2 = entry2.getKey();
								int label2 = entry2.getValue();
								if (label2==j) {
									if (mapBinaryNetwork.get(id1).contains(id2)) {
										fwNetwork.write("" + label1 + ",");
									} else {
										fwNetwork.write(",");
									}
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
		MixtureNetworkClustering<String> mnc = new MixtureNetworkClustering<String>();
		mnc.initNetwork();
		mnc.output(0);
		mnc.outputNetwork(0);
		for (int r=1; r<=100; r++) {
			System.out.println("Starting round " + r + "...");
			mnc.clustering();
			mnc.output(r);
			mnc.outputNetwork(r);
		}
	}

}
