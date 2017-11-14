package readSimulator;
import java.io.File;

import genomeExtractor.GenomeSequenceExtractor;;


public class ReadSimulator {
	
	public static int length;
	public static int frlength;
	public static int SD;	
	public static String readcounts ;
	public static double mutationrate ;
	public static String fasta;
	public static String fidx;
	public static String gtf;
	public static String od;
	
	public static void main(String args[]) {
		
		String noInp = "Error while reading input instructions!!";
		
		for (int i = 0; i < args.length; i++) {
			switch (args[i]) {
			case "-length":
				int len = Integer.parseInt(args[i+1]);
				length = len;
				i++;
				break;
			case "-frlength":
				int flen = Integer.parseInt(args[i+1]);
				frlength = flen;
				i++;
			case "-SD":
				int sd = Integer.parseInt(args[i+1]);
				frlength = sd;
				i++;
			case "-readcounts":
				readcounts = args[i+1];
				i++;
			case "-mutationrate":
				double mutrate = Double.parseDouble(args[i+1]);
				mutationrate = mutrate;
				i++;
			case "-fasta":
				fasta = args[i + 1];
				i++;
			case "-fidx":
				fidx = args[i + 1];
				i++;
			case "-gtf":
				gtf = args[i + 1];
				i++;
			case "-od":
				od = args[i + 1];
				i++;
				break;
			default:
				System.out.println(noInp);
			}
		}
		
		File fasta = new File("C:\\Users\\Anja\\Desktop\\GoBi\\GoBi_2\\input\\Homo_sapiens.GRCh37.75.cdna.all.fa");
		File idx = new File("C:\\Users\\Anja\\Desktop\\GoBi\\GoBi_2\\input\\Homo_sapiens.GRCh37.75.dna.toplevel.fa.fai");
		
//		File fastaG = new File(fasta);
//		File idx = new File(fidx);
//		File GTF = new File(gtf);
		
	GenomeSequenceExtractor gse = new GenomeSequenceExtractor(fasta, idx);

		
		
		String simc = "C:\\Users\\Anja\\Desktop\\GoBi\\GoBi_2\\input\\readcounts.simulation";
		
		gse.readCounts(simc);
		
		
		
		
//		int start = 56;
//		int end = start + 117;
//		
//		for(int i = start; i<end ; i++) {
//			gse.getSequence("chr", i, 2000);			
//		}
		
	}

}
