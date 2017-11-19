package readSimulator;

import java.io.File;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

import genomeExtractor.GenomeSequenceExtractor;
import genomicUtils.Triplet;
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

		File thisFasta = new File(fasta);
		File fastaIndex = new File(fidx);

		long zeit = System.currentTimeMillis();

		GenomeSequenceExtractor gse = new GenomeSequenceExtractor(thisFasta, fastaIndex);
		System.out.println("Finished gse creation");

		gse.readIndex();
		System.out.println("Finished reading Index");

		gse.readCounts(readcounts);
		System.out.println("Finished reading counts");

		gse.readGTF(gtf);
		System.out.println("Finished reading GTF");

		gse.getAllSequences();
		System.out.println("Finished reading Genome");

		long endTime = System.currentTimeMillis() - zeit;
		System.out.println("Alles eingelesen in: " + endTime);

		long length = 0;
		long amountTrans = 0;

		for (Triplet<String, String, Integer> readcount : gse.readcounts) {
			String transSeq = gse.sequences.get(readcount);

			length += transSeq.length();
			amountTrans++;

			for (int i = 0; i < readcount.getThird(); i++) {
				String transSeqFrag = generateFragment(transSeq);
				String fragmentRead = generateRead(transSeqFrag);
				Tuple<StringBuilder, ArrayList<Integer>> mutRead = mutate(fragmentRead);
			}
		}

		System.out.println(length / amountTrans);

		long endTime1 = System.currentTimeMillis() - zeit;
		System.out.println("Fertig in " + endTime1);

		// System.out.println(System.currentTimeMillis() - a);

	}

	public static String generateFragment(String trans) {

		Random rdm = new Random();
		int fragmentLength = Integer.MAX_VALUE;

		while (fragmentLength > trans.length()) {
			double rand = rdm.nextGaussian() * SD + frlength;
			int fragmentLength1 = (int) Math.round(rand);
			fragmentLength = fragmentLength1 < 0 ? -fragmentLength1 : fragmentLength1;
		}

		if (fragmentLength >= length) // fragment length >= read Length
		{
			int possibleStarts = trans.length() - fragmentLength;
			int randomStart = ThreadLocalRandom.current().nextInt(possibleStarts + 1);
			String fragSeq = trans.substring(randomStart, randomStart + fragmentLength);
			return fragSeq;

		} else {
			int possibleStarts = trans.length() - length;
			int randomStart = ThreadLocalRandom.current().nextInt(possibleStarts + 1);
			String fragSeq = trans.substring(randomStart, randomStart + length);
			return fragSeq;

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

		Tuple<StringBuilder, ArrayList<Integer>> mutations = new Tuple<StringBuilder, ArrayList<Integer>>(readS,
				mutPos);

		return mutations;

	}

}
