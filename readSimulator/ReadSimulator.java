package readSimulator;

import java.io.File;

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

		gse.readCounts(readcounts);
		gse.readGTF(gtf);

	}

}
