package em;

import java.io.*;
//import java.util.*;

import org.apache.commons.math3.distribution.*;
//import org.apache.commons.math3.linear.*;

public class GaussianMixtureStandard {

//	String pathname;
	String filename_training;
	String filename_result;
	String filename_testing;
	
	int FEATURE_NUM; // number of features
	int CLUSTER_NUM; // number of clusters
	int DATA_SIZE; // data size
	
	class Cluster{
		// Parameters
		double alpha = 0;
		double[] mu = new double[FEATURE_NUM];
		double[][] sigma = new double[FEATURE_NUM][FEATURE_NUM];
		MultivariateNormalDistribution normal = null;
		
		// Work area
		double Nk = 0;
		double[] sum_mu = new double[FEATURE_NUM];
		double[][] sum_sigma = new double[FEATURE_NUM][FEATURE_NUM];
		
		public Cluster() {
			UniformRealDistribution uniform = new UniformRealDistribution();
			
			for (int f=0; f<FEATURE_NUM; f++) {
				alpha = 1.0/CLUSTER_NUM;
				mu[f] = uniform.sample();
				sigma[f][f] = 10000000000.0;
			}
			updateDistribution();
		}
		
		public void initWorkArea() {
			Nk = 0;
			for (int f=0; f<FEATURE_NUM; f++) {
				sum_mu[f] = 0;
				for (int f2=0; f2<FEATURE_NUM; f2++) {
					sum_sigma[f][f2] = 0;
				}
			}
		}
		
		public double getDensity(double[] p) {
			return normal.density(p);
		}
		
		public void updateAlpha() {
			alpha = Nk/DATA_SIZE;
		}
		
		public void updateMu() {
			for (int f=0; f<FEATURE_NUM; f++) {
				mu[f] = sum_mu[f]/Nk;
			}
		}
		
		public void updateSigma() {
			for (int f=0; f<FEATURE_NUM; f++) {
				for (int f2=f; f2<FEATURE_NUM; f2++) {
					sigma[f][f2] = sum_sigma[f][f2]/Nk;
					sigma[f2][f] = sigma[f][f2];
				}
				sigma[f][f] += 0.2;
			}
		}
		
		public void updateDistribution() {
			try {
				normal = new MultivariateNormalDistribution(mu, sigma);
			} catch (Exception e) {
				System.out.println(e);
		        
				for (int f=0; f<FEATURE_NUM; f++) {
					for (int f2=f+1; f2<FEATURE_NUM; f2++) {
						if (sigma[f][f2]!=sigma[f2][f]) {
							System.out.println("" + sigma[f][f2] + "," + sigma[f2][f]);
						}
					}
				}
				
//				RealMatrix covarianceMatrix = new Array2DRowRealMatrix(sigma);

		        // Covariance matrix eigen decomposition.
//		        final EigenDecomposition covMatDec = new EigenDecomposition(covarianceMatrix);

		        // Compute and store the inverse.
//		        RealMatrix covarianceMatrixInverse = covMatDec.getSolver().getInverse();
			}
		}
	}
	
	public GaussianMixtureStandard(String filename_training, String filename_result, int f, int k, int d, String filename_testing) {
//		this.pathname = pn;
		this.filename_training = filename_training;
		this.filename_result = filename_result;
		this.FEATURE_NUM = f;
		this.CLUSTER_NUM = k;
		this.DATA_SIZE = d;
		this.filename_testing = filename_testing;
		clusters = new Cluster[CLUSTER_NUM];
		for (int c=0; c<CLUSTER_NUM; c++) {
			clusters[c] = new Cluster();
		}
	}
	
	//double[][] w = new double[DATA_SIZE][CLUSTER_NUM];
	//double[][] x = new double[DATA_SIZE][FEATURE_NUM];
	
	Cluster[] clusters;
	
	public void saveClusters(String filename) {
		try {
			FileWriter fw = new FileWriter(filename);
			
			for (int c=0; c<clusters.length; c++) {
				// alpha
				fw.write("" + clusters[c].alpha);
				
				// mu
				for (int f=0; f<FEATURE_NUM; f++) {
					if (f==0) {
						fw.write(";");
					} else {
						fw.write(",");
					}
					fw.write("" + clusters[c].mu[f]);
				}
				
				// sigma
				for (int f=0; f<FEATURE_NUM; f++) {
					for (int f2=0; f2<FEATURE_NUM; f2++) {
						if (f==0 && f2==0) {
							fw.write(";");
						} else {
							fw.write(",");
						}
						fw.write("" + clusters[c].sigma[f][f2]);
					}
				}
				
				fw.write("\n");
			}
			
			fw.close();
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public void loadClusters(String filename) {
		try {
			FileReader fr = new FileReader(filename);
			BufferedReader buf = new BufferedReader(fr);
			String s;
			
			int c=0;
			while ((s=buf.readLine())!=null) {
				String[] ss=s.split(";");
				
				// alpha
				clusters[c].alpha = Double.valueOf(ss[0]);
				
				// mu
				String[] sss1 = ss[1].split(",");
				for (int f=0; f<FEATURE_NUM; f++) {
					clusters[c].mu[f] = Double.valueOf(sss1[f]);
				}
				
				// sigma
				String[] sss2 = ss[2].split(",");
				for (int f=0; f<FEATURE_NUM; f++) {
					for (int f2=0; f2<FEATURE_NUM; f2++) {
						clusters[c].sigma[f][f2] = Double.valueOf(sss2[f+f2*FEATURE_NUM]);
					}
				}
				
				// distribution
				clusters[c].updateDistribution();
				
				c++;
			}
			
			fr.close();
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public void init() {
		for (int c=0; c<CLUSTER_NUM; c++) {
			clusters[c] = new Cluster();
		}
	}
	
	int count_sum0;
	
	public void EStep() {
		for (int c=0; c<CLUSTER_NUM; c++) {
			clusters[c].initWorkArea();
		}
		
		try {
			FileReader fr = new FileReader(filename_training);
			BufferedReader buf = new BufferedReader(fr);
			String s=buf.readLine();
			
			int count = 0;
			int count_sum0 = 0;
			while ((s=buf.readLine())!=null) {
				String[] ss=s.split("\t");
				
				double[] x = new double[FEATURE_NUM];
				for (int idx=0; idx<FEATURE_NUM; idx++) {
					x[idx] = Double.valueOf(ss[idx+1]);
				}
				
				EStep(x);
				
				count++;
				if (count%5000000==0)
					System.out.print("," + count);
			}
			
			if (count_sum0>0) {
				System.out.println(" [# of sum=0] is " + count_sum0);
			}
			
			fr.close();
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public void EStep(double[] x) {
		// compute w
		double[] density = new double[CLUSTER_NUM];
		double sum = 0;
		double[] w = new double[CLUSTER_NUM];
		
		for (int c=0; c<CLUSTER_NUM; c++) {
			double den = clusters[c].getDensity(x);
			density[c] = den * clusters[c].alpha;
			if (Double.isNaN(density[c])) {
				density[c] = 4.9E-324;
			}
			sum += density[c];
		}
		
		if (sum==0) {
			count_sum0++;
			for (int c=0; c<CLUSTER_NUM; c++) {
				w[c] = 1.0/CLUSTER_NUM;
			}
		} else {
			for (int c=0; c<CLUSTER_NUM; c++) {
				w[c] = density[c]/sum;
				if (Double.isNaN(w[c])) {
					System.out.println("NaN1!");
				}
			}
		}
		
		// update for M-step (Nk)
		for (int c=0; c<CLUSTER_NUM; c++) {
			clusters[c].Nk += w[c];
		}
		
		// update for M-step (sum mu)
		for (int c=0; c<CLUSTER_NUM; c++) {
			for (int f=0; f<FEATURE_NUM; f++) {
				clusters[c].sum_mu[f] += w[c]*x[f];
			}
		}
		
		// update for M-step (sum sigma)
		for (int c=0; c<CLUSTER_NUM; c++) {
			double[] x_mu = new double[FEATURE_NUM];
			for (int f=0; f<FEATURE_NUM; f++) {
				x_mu[f] = x[f] - clusters[c].mu[f];
			}
			for (int f=0; f<FEATURE_NUM; f++) {
				for (int f2=0; f2<FEATURE_NUM; f2++) {
					clusters[c].sum_sigma[f][f2] += w[c]*x_mu[f]*x_mu[f2];
				}
			}
		}
	}
	
	// update parameters
	public void MStep() {
		// compute alpha
		for (int c=0; c<CLUSTER_NUM; c++) {
			clusters[c].updateAlpha();
		}
		
		// compute mean
		for (int c=0; c<CLUSTER_NUM; c++) {
			clusters[c].updateMu();
		}
		
		// compute covariance
		for (int c=0; c<CLUSTER_NUM; c++) {
			clusters[c].updateSigma();
		}
		
		// combine mean and covariance
		for (int c=0; c<CLUSTER_NUM; c++) {
			clusters[c].updateDistribution();
		}
	}
	
	public void countSize() {
		int recordcount = 0;
		
		try {
			FileReader fr = new FileReader(filename_training);
			BufferedReader buf = new BufferedReader(fr);
//			String s;
			
			while ((buf.readLine())!=null) {
//				String[] ss = s.split(",");
				recordcount++;
			}
			
			fr.close();
		} catch (Exception e) {
			System.out.println(e);
		}
		
		DATA_SIZE = recordcount;
		System.out.println("User # is " + recordcount);
	}
	
	public void train(int r_begin, int r_end) {
		countSize();
		
		try {
			FileWriter fw = new FileWriter("loglikelihood.dmp");
			
			for (int r=r_begin; r<r_end; r++) {
				System.out.println("Round " + r + "...");
				EStep();
				MStep();
				saveClusters("Clusters_" + r + ".log");
				
				double ll = computeLogLikelihood();
				fw.write("" + ll + "\n");
			}
			
			fw.close();
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public void trainContinue(int round, String filename_testing) {
		loadClusters("Clusters_" + round + ".log");
		int r_begin=round+1;
		int r_end=1000;
		
		train(r_begin, r_end);
	}
	
	public void assignLabels() {
		try {
			FileReader fr = new FileReader(filename_training);
			BufferedReader buf = new BufferedReader(fr);
			String s;
			
			FileWriter fw = new FileWriter(filename_result);
			
			while ((s=buf.readLine())!=null) {
				String[] ss=s.split("\t");
				
				double[] x = new double[FEATURE_NUM];
				for (int idx=0; idx<FEATURE_NUM; idx++) {
					double d = Double.valueOf(ss[idx+1]);
					x[idx] = d;
				}
				
				int c_max = assignLabels(x);
				
				fw.write(""
						+ c_max + "\t" + s + "\n");
			}
			
			fw.close();
			fr.close();
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public int assignLabels(double[] x) {
		// compute w
		double[] density = new double[CLUSTER_NUM];
		double sum = 0;
		double[] w = new double[CLUSTER_NUM];
		
		for (int c=0; c<CLUSTER_NUM; c++) {
			double den = clusters[c].getDensity(x);
			density[c] = den * clusters[c].alpha;
			if (Double.isNaN(density[c])) {
				density[c] = 4.9E-324;
			}
			sum += density[c];
		}
		
		if (sum==0) {
			count_sum0++;
			for (int c=0; c<CLUSTER_NUM; c++) {
				w[c] = 1.0/CLUSTER_NUM;
			}
			
			return -1;
		} else {
			int c_max = 0;
			for (int c=0; c<CLUSTER_NUM; c++) {
				w[c] = density[c]/sum;
				if (Double.isNaN(w[c])) {
					System.out.println("NaN2!");
				}
				
				if (w[c] > w[c_max]) {
					c_max = c;
				}
			}
			
			return c_max;
		}
	}
	
	public double computeLogLikelihood(double[] x) {
		double l = 0;
		
		for (int c=0; c<CLUSTER_NUM; c++) {
			double component = clusters[c].alpha * clusters[c].normal.density(x);
			if (!Double.isInfinite(component) && !Double.isNaN(component)) {
				l += component;
			}
		}
		
		if (l<=0.0) {
			System.out.println("log (0)");
		}
		
		double ll = Math.log(l);
		
		return ll;
	}
	
	public double computeLogLikelihood() {
		double sum_ll = 0;
		int count = 0;
		
		try {
			FileReader fr = new FileReader(filename_testing);
			BufferedReader buf = new BufferedReader(fr);
			String s;
			
			while ((s=buf.readLine())!=null) {
				String[] ss=s.split("\t");
				
				double[] x = new double[FEATURE_NUM];
				for (int idx=0; idx<FEATURE_NUM; idx++) {
					x[idx] = Double.valueOf(ss[idx+1]);
				}
				
				double logLikelihood = computeLogLikelihood(x);
				sum_ll += logLikelihood;
				count++;
			}
			
			fr.close();
		} catch(Exception e) {
			System.out.println(e);
		}
		
		double ll = sum_ll/count;
		
		System.out.println("Log Likelihood: " + ll);
		
		return ll;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if (args.length!=7) {
			System.out.println("The parameter number is invalid!");
		}
		
		String flg = args[0];
		String filename_training = args[1];
		String filename_result = args[2];
		int FEATURE_NUM = Integer.valueOf(args[3]);
		int CLUSTER_NUM = Integer.valueOf(args[4]);
		String filename_model = args[5];
		String filename_testing = args[6];
		
		GaussianMixtureStandard gm = new GaussianMixtureStandard(filename_training, filename_result, FEATURE_NUM, CLUSTER_NUM, 0, filename_testing);
		
		if (flg.equalsIgnoreCase("testing")) {
			System.out.println("Compute log likelihood...");
			gm.loadClusters(filename_model);
			gm.computeLogLikelihood();
		} else {
			System.out.println("Run from the begining...");
			gm.train(0, 300);
			gm.assignLabels();
			
			System.out.println("Compute log likelihood...");
			gm.computeLogLikelihood();
		}
	}

}
