package readSimulator;

import java.io.File;
import java.io.RandomAccessFile;

import genomeExtractor.GenomeSequenceExtractor;;

public class ReadSimulator {

	public static int length;
	public static int frlength;
	public static int SD;
	public static String readcounts;
	public static double mutationrate;
	public static String fasta;
	public static String fidx;
	public static String gtf;
	public static String od;

	public static void main(String args[]) {

		String noInp = "Error while reading input instructions!!";

		for (int i = 0; i < args.length; i++) {
			switch (args[i]) {
			case "-length":
				int len = Integer.parseInt(args[i + 1]);
				length = len;
				i++;
				System.out.println("length: " + length);
				break;

			case "-frlength":
				int flen = Integer.parseInt(args[i + 1]);
				frlength = flen;
				i++;
				System.out.println("frlength: " + frlength);
				break;

			case "-SD":
				int sd = Integer.parseInt(args[i + 1]);
				SD = sd;
				i++;
				System.out.println("SD: " + SD);
				break;

			case "-readcounts":
				readcounts = args[i + 1];
				i++;
				System.out.println("readcounts: " + readcounts);
				break;

			case "-mutationrate":
				double mutrate = Double.parseDouble(args[i + 1]);
				mutationrate = mutrate;
				i++;
				System.out.println("mutationrate: " + mutationrate);
				break;

			case "-fasta":
				fasta = args[i + 1];
				i++;
				System.out.println("fasta: " + fasta);
				break;

			case "-fidx":
				fidx = args[i + 1];
				i++;
				System.out.println("fidx: " + fidx);
				break;

			case "-gtf":
				gtf = args[i + 1];
				System.out.println("gtf: " + gtf);
				i++;
				break;

			case "-od":
				od = args[i + 1];
				i++;
				System.out.println("od: " + od);
				break;

			default:
				System.out.println(noInp);
			}
		}

		File thisFasta = new File(fasta);
		File fastaIndex = new File(fidx);

		GenomeSequenceExtractor gse = new GenomeSequenceExtractor(thisFasta, fastaIndex);
		System.out.println("Finished gse creation");
		
		gse.readIndex();
		System.out.println("Finished reading Index");
		
		gse.readCounts(readcounts);
		System.out.println("Finished reading counts");
		
		gse.readGTF(gtf);
		System.out.println("Finished reading GTF");
		
		gse.getAllSequences();
		System.out.println("Finished getting sequences");
		

//		for(String trans : gse.genes.get("ENSG00000131018").transcripts.keySet()) 
//		{
//			int len = gse.genes.get("ENSG00000131018").transcripts.get(trans).getLength();
//			System.out.println(trans + " " +len);
//		}
		
//		String seq = "";
//		try {
//			RandomAccessFile raffasta = new RandomAccessFile(fasta, "r");
//
//			raffasta.seek(4203);
//
//			for (int i = 0; i < 100 ; i++) {
//				int intChar = raffasta.read();
//				String character = Character.toString((char) intChar);
//				System.out.print(i + character);
//				seq += character;
//			}
//
//			raffasta.close();
//
//		} catch (Exception e) {
//			throw new RuntimeException("got error while reading input raf files", e);
//		}
//		// System.out.println("Sequence: " + seq);
//
//		seq = seq.replaceAll("\n", "");
//		
	}

}
