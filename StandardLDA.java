/*****************************
 * Parameters:
 *   - Input file name;
 *   - Output path name;
 *   - Topic number;
 *   - Thread number.
 */
package lda;

import java.io.*;
import java.util.*;
import org.apache.commons.math3.distribution.*;

public class StandardLDA {
	
	public class Profile {
		String pathnameInput;
		String pathnameOutput;
		String filenameInputTraining;
		String filenameInputTesting;
		
		int THREAD_NUM = 8;
		
		int TOPIC_NUM = 6;
		
		int ROUND_ITERATION = 300;
		
		int ROUND_START = -1;
		
		long ALPHA = 1;
		long BETA = 1;
		
		public boolean set(String[] args) {
			if (args.length < 5) {
				return false;
			}
			
			filenameInputTraining = args[0];
			System.out.println("filenameInputTraining: " + filenameInputTraining);
			
			pathnameOutput = args[1];
			if (pathnameOutput.charAt(pathnameOutput.length()-1) != File.separatorChar)
				pathnameOutput += File.separator;
			File dir = new File(pathnameOutput);
			if (!dir.exists()) {
				dir.mkdir();
			}
			System.out.println("pathnameOutput: " + pathnameOutput);
			
			TOPIC_NUM = Integer.valueOf(args[2]);
			System.out.println("TOPIC_NUM = " + TOPIC_NUM);
			
			THREAD_NUM = Integer.valueOf(args[3]);
			System.out.println("THREAD_NUM = " + THREAD_NUM);
			
			ROUND_ITERATION = Integer.valueOf(args[4]);
			System.out.println("ROUND_ITERATION = " + ROUND_ITERATION);
			
			ALPHA = Integer.valueOf(args[5]);
			System.out.println("ALPHA = " + ALPHA);
			
			BETA = Integer.valueOf(args[6]);
			System.out.println("BETA = " + BETA);
			
			if (args.length >= 8) {
				ROUND_START = Integer.valueOf(args[7]);
				System.out.println("ROUND_START = " + ROUND_START);
			}
			
			return true;
		}
	}
	
	public Profile profile = new Profile();
	
	int WORD_NUM;
	int DOC_NUM;
	
	public class Counts {
		long[] nz;
		long[][] nzw;
		long[] nd;
		long[][] ndz;
		
		public void init() {
			nz = new long[profile.TOPIC_NUM];
			nzw = new long[profile.TOPIC_NUM][WORD_NUM];
			nd = new long[DOC_NUM];
			ndz = new long[DOC_NUM][profile.TOPIC_NUM];
			System.out.println("allocation is finished!");
		}
		
		public void add(int z, int w, int d) {
			nz[z]++;
			nzw[z][w]++;
			nd[d]++;
			ndz[d][z]++;
		}
		
		public synchronized void update(int z_old, int z, int w, int d) {
			nz[z_old]--;
			nz[z]++;
			nzw[z_old][w]--;
			nzw[z][w]++;
			ndz[d][z_old]--;
			ndz[d][z]++;
		}
		
		public class NZW {
			int w;
			long count;
		}
		
		public class NZ {
			int z;
			long count;
		}
		
		public void output(int r) {
			String filename = profile.pathnameOutput + "dump" + r + ".dmp";
			try {
				FileWriter fw = new FileWriter(filename);
				ArrayList<NZ> listnz = new ArrayList<NZ>();
				for (int j=0; j<profile.TOPIC_NUM; j++) {
					NZ nz_new = new NZ();
					nz_new.z = j;
					nz_new.count = nz[j];
					listnz.add(nz_new);
				}
				Collections.sort(listnz, new Comparator<NZ>() {
					public int compare(NZ obj1, NZ obj2) {
						return (int)(obj2.count-obj1.count);
					}
				});
				
				for (Iterator<NZ> it=listnz.iterator(); it.hasNext(); ) {
					int j = it.next().z;
					fw.write("Topic " + j + ":" + nz[j]);
					ArrayList<NZW> listnzx = new ArrayList<NZW>();
					for (int i=0; i<WORD_NUM; i++) {
						NZW nzx_new = new NZW();
						nzx_new.w = i;
						nzx_new.count = nzw[j][i];
						listnzx.add(nzx_new);
					}
					Collections.sort(listnzx, new Comparator<NZW>() {
						public int compare(NZW obj1, NZW obj2) {
							return (int)(obj2.count-obj1.count);
						}
					});
					for (int i=0; i<listnzx.size(); i++) {
						//if (listnzx.get(i).count<10)
						//	break;
						if (i>100)
							break;
						fw.write(", " + mapWord.getName(listnzx.get(i).w) + ":" + listnzx.get(i).count);
					}
					fw.write("\n");
				}
				fw.close();
			} catch (Exception e) {
				System.out.println(e);
			}
		}
		
		public void dumptopics() {
			String filename = profile.pathnameOutput + "topics.dmp";
			try {
				FileWriter fw = new FileWriter(filename);
				
				// nzw = new long[profile.TOPIC_NUM][WORD_NUM];
				for (int j=0; j<profile.TOPIC_NUM; j++) {
					for (int i=0; i<WORD_NUM; i++) {
						long count = nzw[j][i];
						fw.write("" + j + "," + mapWord.getName(i) + "," + count + "\n");
					}
				}
				
				fw.close();
			} catch (Exception e) {
				System.out.println(e);
			}
		}
		
		public void check() {
			long sum1 = 0;
			for (int z=0; z<profile.TOPIC_NUM; z++) {
				long sumz = 0;
				for (int w=0; w<WORD_NUM; w++) {
					sumz += nzw[z][w];
				}
				if (sumz != nz[z]) {
					System.out.println("sumz not equal:" + z);
				}
				sum1 += sumz;
			}
			long sum2 = 0;
			for (int d=0; d<DOC_NUM; d++) {
				long sumd = 0;
				for(int z=0; z<profile.TOPIC_NUM; z++) {
					sumd += ndz[d][z];
				}
				if (sumd != nd[d]) {
					System.out.println("sumd not equal:" + d);
				}
				sum2 += sumd;
			}
			System.out.println("" + sum1 + ":" + sum2);
		}
	}
	
	Counts counts = new Counts();
	
	public class MapId {
		String filename;
		HashMap<Integer, String> mapIdtoName = new HashMap<Integer, String>();
		HashMap<String, Integer> mapNametoId = new HashMap<String, Integer>();
		
		public void set(String fn) {
			this.filename = fn;
		}
		
		public void insert(String name) {
			if (!mapNametoId.containsKey(name)) {
				int size = mapIdtoName.size();
				mapNametoId.put(name, size);
				mapIdtoName.put(size, name);
			}
		}
		
		public String getName(int id) {
			if (mapIdtoName.containsKey(id)) {
				return mapIdtoName.get(id);
			}
			
			return "";
		}
		
		public int getId(String name) {
			if (mapNametoId.containsKey(name)) {
				return mapNametoId.get(name);
			}
			
			return -1;
		}
		
		public int getSize() {
			return mapIdtoName.size();
		}
		
		public void readMap() {
			if (!mapIdtoName.isEmpty())
				return;
			
			mapIdtoName = new HashMap<Integer, String>();
			mapNametoId = new HashMap<String, Integer>();
			
			try {
				FileReader fr = new FileReader(filename);
				BufferedReader buf = new BufferedReader(fr);
				String s;
				
				while ((s=buf.readLine())!=null) {
					String[] ss=s.split("\t");
					int id = Integer.valueOf(ss[0]);
					String name = ss[1];
					mapIdtoName.put(id, name);
					mapNametoId.put(name, id);
				}
				
				fr.close();
			} catch (Exception e) {
				System.out.println(e);
			}
		}
		
		public void writeMap() {
			try {
				FileWriter fw = new FileWriter(filename);
				
				for (Iterator<Map.Entry<Integer, String>> it=mapIdtoName.entrySet().iterator(); it.hasNext(); ) {
					Map.Entry<Integer, String> entry = it.next();
					int id = entry.getKey();
					String name = entry.getValue();
					fw.write("" + id + "\t" + name + "\n");
				}
				
				fw.close();
			} catch (Exception e) {
				System.out.println(e);
			}
		}
	}
	
	// Mapping d and doc
	MapId mapDoc = new MapId();
	
	// Mapping w and word
	MapId mapWord = new MapId();
	
	public void split(int round) {
		System.out.println("Spliting data from iteration " + round);
		try {
			FileWriter[] fws = new FileWriter[profile.THREAD_NUM];
			for (int t=0; t<profile.THREAD_NUM; t++) {
				fws[t] = new FileWriter(profile.pathnameOutput + "work_" + round + "_" + t + ".dat");
			}
			
			FileReader fr = new FileReader(profile.filenameInputTraining);
			BufferedReader buf = new BufferedReader(fr);
			String s;
			
			int recordcount = 0;
			int t_new = 0;
			
			while ((s=buf.readLine())!=null) {
				String[] ss=s.split("\t");
				
				recordcount++;
				
				String docName = ss[0];
				mapDoc.insert(docName);
				
				String wordName = ss[1];
				mapWord.insert(wordName);
				
				int z;
				if (ss.length==2) {
					UniformIntegerDistribution uniform = new UniformIntegerDistribution(0, profile.TOPIC_NUM-1);
					z = uniform.sample();
				} else {
					z = Integer.valueOf(ss[2]);
				}
				fws[t_new].write("" + mapDoc.getId(docName) + "\t" + mapWord.getId(wordName) + "\t" + z + "\n");
				t_new = (t_new+1)%profile.THREAD_NUM;
			}
			
			fr.close();
			
			for (int t=0; t<profile.THREAD_NUM; t++) {
				fws[t].close();
			}
			
			DOC_NUM = mapDoc.getSize();
			WORD_NUM = mapWord.getSize();
			
			mapDoc.writeMap();
			mapWord.writeMap();
			
			System.out.println("Doc count: " + DOC_NUM);
			System.out.println("Word count: " + WORD_NUM);
			System.out.println("Record count:" + recordcount);
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public void construct(int round) {
		mapDoc.readMap();
		mapWord.readMap();
		
		DOC_NUM = mapDoc.getSize();
		WORD_NUM = mapWord.getSize();
		
		try {
			counts.init();
			
			int recordcount = 0;
			for (int t=0; t<profile.THREAD_NUM; t++) {
				FileReader fr = new FileReader(profile.pathnameOutput + "work_" + round + "_" + t + ".dat");
				BufferedReader buf = new BufferedReader(fr);
				String s;
				
				while ((s=buf.readLine())!=null) {
					String[] ss=s.split("\t");
					
					int d = Integer.valueOf(ss[0]);
					int w = Integer.valueOf(ss[1]);
					int z = Integer.valueOf(ss[2]);
					
					counts.add(z, w, d);
					
					recordcount++;
				}
				
				fr.close();
			}
			
			System.out.println("Doc count: " + DOC_NUM);
			System.out.println("Word count: " + WORD_NUM);
			System.out.println("Record count:" + recordcount);
		} catch (Exception e) {
			System.out.println(e);
		}
		
		counts.check();
	}
	
	public class Sampler implements Runnable {
		int round;
		String filename;
		String filename_new;
		
		public Sampler(int r, String file, String file_new) {
			round = r;
			filename = file;
			filename_new = file_new;
		}
		
		public void run() {
			System.out.println("[Round" + round + "] Thread for [" + filename + "] is started.");
			
			int changedNum = 0;
			try {
				FileReader fr = new FileReader(filename);
				BufferedReader buf = new BufferedReader(fr);
				String s;
				
				FileWriter fw = new FileWriter(filename_new);
				
				while ((s=buf.readLine())!=null) {
					String[] ss=s.split("\t");
					
					int d = Integer.valueOf(ss[0]);
					int w = Integer.valueOf(ss[1]);
					int z = Integer.valueOf(ss[2]);
					
					int z_new = sample(d, w, z);
					if (z_new != z) {
						changedNum++;
					}
					fw.write(ss[0] + "\t" + ss[1] + "\t" + z_new + "\n");
				}
				
				fw.close();
				
				fr.close();
			} catch (Exception e) {
				System.out.println(e);
			}
			
			try {
				File file = new File(filename);
			    if (file.exists()) {
			    	file.delete();
			    }
			} catch (Exception e) {
				System.out.println(e);
			}
			
			System.out.println("[Round" + round + "] Thread  for [" + filename + "] finished (" + changedNum + ").");
		}
		
		public int sample(int d, int w, int z) {
			int[] singletons = new int[profile.TOPIC_NUM];
			double[] probabilities = new double[profile.TOPIC_NUM];
			
			for (int j=0; j<profile.TOPIC_NUM; j++) {
				singletons[j] = j;
				if (j==z) {
					probabilities[j] = (double)(counts.nzw[j][w] - 1 + profile.BETA)/(counts.nz[j] - 1 + WORD_NUM*profile.BETA);
					probabilities[j] *= (double)(counts.ndz[d][j] - 1 + profile.ALPHA)/(counts.nd[d] - 1 + profile.TOPIC_NUM*profile.ALPHA);
				} else {
					probabilities[j] = (double)(counts.nzw[j][w] + profile.BETA)/(counts.nz[j] + WORD_NUM*profile.BETA);
					probabilities[j] *= (double)(counts.ndz[d][j] + profile.ALPHA)/(counts.nd[d] + profile.TOPIC_NUM*profile.ALPHA);
				}
			}
			
			EnumeratedIntegerDistribution discrete = new EnumeratedIntegerDistribution(singletons, probabilities);
			int z_new = discrete.sample();
			
			counts.update(z, z_new, w, d);
			
			return z_new;
		}
	}
	
	public void sampling(int round) {
		// run threads
		Thread[] threads = new Thread[profile.THREAD_NUM];
		
		for (int t=0; t<profile.THREAD_NUM; t++) {
			String filename = profile.pathnameOutput + "work_" + round + "_" + t + ".dat";
			String filename_new = profile.pathnameOutput + "work_" + (round+1) + "_" + t + ".dat";
			threads[t] = new Thread(new Sampler(round, filename, filename_new));
			threads[t].start();
		}
		
		try {
			for (int i=0; i<profile.THREAD_NUM; i++) {
				threads[i].join();
			}
		} catch (InterruptedException e) {
			System.out.println(e);
		}
		
		counts.check();
		
		if (round%1==0) {
			counts.output(round);
		}
	}
	
	public void combine(int round) {
		mapDoc.readMap();
		mapWord.readMap();
		
		DOC_NUM = mapDoc.getSize();
		WORD_NUM = mapWord.getSize();
		
		try {
			FileReader[] frs = new FileReader[profile.THREAD_NUM];
			BufferedReader[] bufs = new BufferedReader[profile.THREAD_NUM];
			String[] s = new String[profile.THREAD_NUM];
			
			for (int t=0; t<profile.THREAD_NUM; t++) {
				frs[t] = new FileReader(profile.pathnameOutput + "work_" + round + "_" + t + ".dat");
				bufs[t] = new BufferedReader(frs[t]);
			}
			
			String filenamew = profile.pathnameOutput + "result.dat";
			FileWriter fw = new FileWriter(filenamew);
			
			boolean bContinue = true;
			while (bContinue) {
				bContinue = false;
				for (int t=0; t<profile.THREAD_NUM; t++) {
					if ((s[t]=bufs[t].readLine())!=null) {
						String[] ss=s[t].split("\t");
						
						int d = Integer.valueOf(ss[0]);
						String docName = mapDoc.getName(d);
						int w = Integer.valueOf(ss[1]);
						String wordName = mapWord.getName(w);
						int z = Integer.valueOf(ss[2]);
						
						fw.write(docName + "\t" + wordName + "\t" + z + "\n");
						
						bContinue = true;
					}
				}
			}
			
			fw.close();
			
			for (int t=0; t<profile.THREAD_NUM; t++) {
				frs[t].close();
			}
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public void train() {
		int r;
		int r_begin=0;
		int r_end=profile.ROUND_ITERATION;
		
		split(r_begin);
		construct(r_begin);
		for (r=r_begin; r<r_end; r++) {
			System.out.println("Round " + r + "...");
			sampling(r);
		}
		combine(r_end);
		counts.dumptopics();
	}
	
	public void trainContinue() {
		int r;
		int r_begin=profile.ROUND_START;
		int r_end=profile.ROUND_START + profile.ROUND_ITERATION;
		
		construct(r_begin);
		for (r=r_begin; r<r_end; r++) {
			System.out.println("Round " + r + "...");
			sampling(r);
		}
		combine(r_end);
		counts.dumptopics();
	}
	
	public void init(String[] args) {
//		String[] argslist = {"E:\\IchibaFashion\\RawData\\DataFromJoe\\LDA_model_input_20150212.tsv",
//				"E:\\IchibaFashion\\OutputPhase1\\",
//				"60",
//				"8",
//				"35",
//				"E:\\IchibaFashion\\RawData\\DataFromJoe\\mapStyleToTopic.csv",
//				"E:\\IchibaFashion\\RawData\\DataFromJoe\\LadyFashionStyleModelCheckData_sadakane_round2_checked.csv"};
		
//		String[] argslist = {"E:\\IchibaFashion\\semiLDA\\lda_datamart.tsv",
//		"E:\\IchibaFashion\\semiLDA\\tmp\\",
//		"60",
//		"4",
//		"35",
//		"E:\\IchibaFashion\\semiLDA\\stl_mst_input.csv",
//		"E:\\IchibaFashion\\semiLDA\\super_input.csv"};
		
		String[] argslist = args;
		profile.set(argslist);
		mapDoc.set(profile.pathnameOutput + "mapDoc.dat");
		mapWord.set(profile.pathnameOutput + "mapWord.dat");
	}
	
	public void run(String[] args) {
		init(args);
		if (profile.ROUND_START == -1) {
			System.out.println("Starting training ...");
			train();
		} else {
			System.out.println("Continuing training from iteration " + profile.ROUND_START);
			trainContinue();
		}
	
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		StandardLDA lda = new StandardLDA();
		lda.run(args);
	}

}
