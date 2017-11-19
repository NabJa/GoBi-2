package readSimulator;

import java.io.File;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

import genomeExtractor.GenomeSequenceExtractor;
import genomicUtils.Tuple;;

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

		// File thisFasta = new File(fasta);
		// File fastaIndex = new File(fidx);
		//
		// GenomeSequenceExtractor gse = new GenomeSequenceExtractor(thisFasta,
		// fastaIndex);
		// System.out.println("Finished gse creation");
		//
		// gse.readIndex();
		// System.out.println("Finished reading Index");
		//
		// gse.readCounts(readcounts);
		// System.out.println("Finished reading counts");
		//
		// gse.readGTF(gtf);
		// System.out.println("Finished reading GTF");
		//
		//
		// long a = System.currentTimeMillis();
		//
		// gse.getAllSequences();
		// System.out.println("Finished getting sequences");
		//
		// System.out.println(System.currentTimeMillis() - a);

		// String seq = "";
		// try {
		// RandomAccessFile raffasta = new RandomAccessFile(fasta, "r");
		//
		// raffasta.seek(56);
		// int intChar = raffasta.read();
		// String character = Character.toString((char) intChar);
		// System.out.print(character);
		// seq += character;
		//
		// raffasta.close();
		//
		// } catch (Exception e) {
		// throw new RuntimeException("got error while reading input raf files", e);
		// }
		// // System.out.println("Sequence: " + seq);
		//
		// seq = seq.replaceAll("\n", "");

		String transcript = "ATGGCGGGGAGGAGGAGGAGAAGGCGGCGGCGGACCGAGCTGCGCTCTGTCAGTACCATTTGAGCCATTCGCTTCCTGACAAGGCCCGTGGCGAGGGGAGAGGAGCTGAAGGGGCCGTGGGGGATCAGTGTGACTGTGGGAAGATGGAGGAGTATGAGAAGTTCTGTGAAAAAAGTCTTGCAGAATACAAGAAGCATCACTATCCACAGAGAGCTTTCTCCCTGCTCAGTCTGAAAGTATCTCACTTATTCGCTTTCATGGAGTGGCTATCCTTTCTCCACTGAAAGCGAGGAGTTACTAAAAAGCAAGATGTTAGCTTTTGAAGAAATGCGGAAGAGACTAGAAGAACAGCACGCCCAGCAATTATCACTACTCATAGCTGAGCAGGAAAGGGAACAAGAAAGACTGCAAAAGGAAATAGAAGAGCAGGAGAAAATGTTAAAAGAGAAGAAGGCAATGACAGCGGAAGCCTCTGAGTTGGACATTAACAATGCAGTGGAATTAGAATGGAGAAAAATAAGTGACTCTAGTTTGCTGGAAACAATGCTGTCTCAAGCGGACTCACTCCATACTTCAAATTCAAATAGTTCTGGT";
		generateFragment(transcript);

		String mutante = "AAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCC";
		Tuple<StringBuilder, ArrayList<Integer>> mutierte = mutate(mutante);

	}

	public static String generateFragment(String trans) {

		Random rdm = new Random();
		double rand = rdm.nextGaussian() * SD + frlength;
		int fragmentLength1 = (int) Math.round(rand);
		int fragmentLength = fragmentLength1 < 0 ? -fragmentLength1 : fragmentLength1;

		// System.out.println();
		// System.out.println("Normal distrubution length: " + fragmentLength);

		if (fragmentLength >= length) // fragment length >= read Length
		{
			int possibleStarts = trans.length() - fragmentLength;
			int randomStart = ThreadLocalRandom.current().nextInt(possibleStarts + 1);
			String fragSeq = trans.substring(randomStart, randomStart + fragmentLength);
			return fragSeq;

			// System.out.println("Track 1");
			// System.out.println(fragSeq);
			// System.out.println("Fragment Length" + fragSeq.length());

		} else {
			int possibleStarts = trans.length() - length;
			int randomStart = ThreadLocalRandom.current().nextInt(possibleStarts + 1);
			String fragSeq = trans.substring(randomStart, randomStart + length);
			return fragSeq;

			// System.out.println("Track 2");
			// System.out.println(fragSeq);
			// System.out.println("Fragment Length" + fragSeq.length());
		}

	}

	public static String generateRead(String fragment) {

		int possibleStarts = fragment.length() - length;
		int randomStart = ThreadLocalRandom.current().nextInt(possibleStarts + 1);

		String read = fragment.substring(randomStart, randomStart + length);

		return read;

	}

	public static Tuple<StringBuilder, ArrayList<Integer>> mutate(String read) {

		Random rdm = new Random();
		StringBuilder readS = new StringBuilder(read);
		ArrayList<Integer> mutPos = new ArrayList<Integer>();

		for (int i = 0; i < read.length(); i++) {
			double rdmNumber = rdm.nextDouble() * 100;

			if (rdmNumber < mutationrate) {

				mutPos.add(i);
				int rdmChar = rdm.nextInt(2);

				switch (read.charAt(i)) {
				case 'A':
					char[] basesA = { 'T', 'G', 'C' };
					readS.setCharAt(i, basesA[rdmChar]);
					break;
				case 'T':
					char[] basesT = { 'A', 'G', 'C' };
					readS.setCharAt(i, basesT[rdmChar]);
					break;
				case 'G':
					char[] basesG = { 'T', 'A', 'C' };
					readS.setCharAt(i, basesG[rdmChar]);
					break;
				case 'C':
					char[] basesC = { 'T', 'G', 'A' };
					readS.setCharAt(i, basesC[rdmChar]);
					break;
				default:
					break;
				}
			}
		}
		
		Tuple<StringBuilder, ArrayList<Integer>>mutations  = new Tuple<StringBuilder, ArrayList<Integer>>(readS, mutPos);

		return mutations;

	}

}
