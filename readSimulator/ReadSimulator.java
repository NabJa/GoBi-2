package readSimulator;

import java.io.File;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

import fastqWriter.FastqWriter;
import genomeExtractor.GenomeSequenceExtractor;
import genomicUtils.DNAUtils;
import genomicUtils.Gene;
import genomicUtils.Region;
import genomicUtils.RegionVector;
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

		File rwod = new File(od, "rw.fastq");
		File fwod = new File(od, "fw.fastq");
		File map = new File(od, "read.mappinginfo");

		System.out.println("Output destinations: " + rwod + " " + fwod + " " + map);

		// Initialize Output wirter and Utils
		ReadInfoMap mapInfo = new ReadInfoMap(map);
		FastqWriter fastqRW = new FastqWriter(rwod);
		FastqWriter fastqFW = new FastqWriter(fwod);
		DNAUtils dnautil = new DNAUtils();

		// Generate String of length = read length with Score values. Here only I used
		StringBuffer qualityBuffer = new StringBuffer(length);
		for (int i = 0; i < length; i++) {
			qualityBuffer.append("I");
		}
		String qualityString = qualityBuffer.toString();

		// Generate all fragments and reads and write them directly in output
		long id = 0;

		for (Triplet<String, String, Integer> readcount : gse.readcounts) // For every transcript in readcounts ...
		{
			Gene gene = gse.genes.get(readcount.getFirst());
			String strand = gene.strand;
			String chr = gene.geneChr;
			String geneID = gene.geneID;

			if (strand.equals("-")) {

				String transSeq = gse.sequences.get(readcount);
				String revTransSeq = dnautil.revcomp(transSeq);
				String transID = readcount.getSecond();

				RegionVector transGenomicRegion = gse.transGenomicRegions.get(readcount);

				for (int i = 0; i < readcount.getThird(); i++) // ... Make count many reads
				{
					generateMinusReads(transSeq, fastqFW, fastqRW, id, qualityString, transGenomicRegion, mapInfo, chr,
							geneID, transID);
					id++;
				}
			} else {

				String transSeq = gse.sequences.get(readcount);
				String revTransSeq = dnautil.revcomp(transSeq);
				String transID = readcount.getSecond();

				RegionVector transGenomicRegion = gse.transGenomicRegions.get(readcount);
				for (int i = 0; i < readcount.getThird(); i++) // ... Make count many reads
				{
					generatePlusReads(transSeq, fastqFW, fastqRW, id, qualityString, transGenomicRegion, mapInfo, chr,
							geneID, transID);
					id++;
				}
			}
		}

		// Close all writer
		fastqRW.closeFastq();
		fastqFW.closeFastq();
		mapInfo.closeMapinfo();

		// Program Time
		long endTime1 = System.currentTimeMillis() - zeit;
		System.out.println("Fertig in " + endTime1);

	}

	public static void generateMinusReads(String transSeq, FastqWriter fastqFW, FastqWriter fastqRW, long id,
			String qualityString, RegionVector transGenomicRegion, ReadInfoMap mapInfo, String chr, String geneID,
			String transID) {

		Tuple<String, Integer> transSeqFrag = generateFragment(transSeq);

		String fragmentRwRead = generateRead(transSeqFrag.getFirst());
		String fragmentRead = generateRwRead(transSeqFrag.getFirst());

		Tuple<StringBuilder, ArrayList<Integer>> mutRead = mutate(fragmentRead);
		Tuple<StringBuilder, ArrayList<Integer>> mutRwRead = mutate(fragmentRwRead);

		fastqFW.writeFastq(id, mutRead.getFirst().toString(), id, qualityString);
		fastqRW.writeFastq(id, mutRwRead.getFirst().toString(), id, qualityString);

		RegionVector fw_regvec = getGenomicRegion(mutRwRead.getFirst().toString(), transGenomicRegion,
				(transSeqFrag.getSecond() + (transSeqFrag.getFirst().length() - length)));

		RegionVector rw_regvec = getGenomicRegion(mutRead.getFirst().toString(), transGenomicRegion,
				transSeqFrag.getSecond());

		int revStartTrans = transSeq.length() - transSeqFrag.getSecond() - length;
		int revEndTrans = revStartTrans + length;
		Region t_rw_regvec = new Region(revStartTrans, revEndTrans);

		int fwStartTrans = revEndTrans - transSeqFrag.getFirst().length() - 1;
		int fwEndTrans = fwStartTrans + length;
		Region t_fw_regvec = new Region(fwStartTrans, fwEndTrans);

		// if(!transSeq.substring(revStartTrans, revEndTrans).equals(fragmentRwRead)) {
		// System.out.println(transSeq.length());
		// System.out.println("Fragment: " + transSeqFrag);
		// System.out.println("TranSeq:\t\t" + transSeq.substring(revStartTrans,
		// revEndTrans));
		// System.out.println("RvRead:\t\t" + fragmentRwRead);
		// System.out.println();
		// }
		//
		// if(!transSeq.substring(fwStartTrans, fwEndTrans).equals(fragmentRead)) {
		// System.out.println(transSeq.length());
		// System.out.println("Fragment: " + transSeqFrag);
		// System.out.println("TranSeq:\t\t" + transSeq.substring(fwStartTrans,
		// fwEndTrans));
		// System.out.println("FwRead:\t\t" + fragmentRead);
		// System.out.println();
		// }

//		if (mutRead.getSecond().contains(0) || mutRwRead.getSecond().contains(0)) {
//			System.out.println(mutRead);
//			System.out.println(mutRwRead);
//		}

		mapInfo.writeMapinfo(id, chr, geneID, transID, rw_regvec, fw_regvec, t_fw_regvec, t_rw_regvec,
				mutRead.getSecond(), mutRwRead.getSecond());
	}

	public static void generatePlusReads(String transSeq, FastqWriter fastqFW, FastqWriter fastqRW, long id,
			String qualityString, RegionVector transGenomicRegion, ReadInfoMap mapInfo, String chr, String geneID,
			String transID) {

		Tuple<String, Integer> transSeqFrag = generateFragment(transSeq);

		// Get Fw and Rev read sequences:
		String fragmentRead = generateRead(transSeqFrag.getFirst());
		String fragmentRwRead = generateRwRead(transSeqFrag.getFirst());
		Tuple<StringBuilder, ArrayList<Integer>> mutRead = mutate(fragmentRead);
		Tuple<StringBuilder, ArrayList<Integer>> mutRwRead = mutate(fragmentRwRead);
		// Write reads in Fastq
		fastqFW.writeFastq(id, mutRead.getFirst().toString(), id, qualityString);
		fastqRW.writeFastq(id, mutRwRead.getFirst().toString(), id, qualityString);

		RegionVector rw_regvec = getGenomicRegion(mutRwRead.getFirst().toString(), transGenomicRegion,
				(transSeqFrag.getSecond() + (transSeqFrag.getFirst().length() - length)));

		RegionVector fw_regvec = getGenomicRegion(mutRead.getFirst().toString(), transGenomicRegion,
				transSeqFrag.getSecond());

		Region t_fw_regvec = new Region(transSeqFrag.getSecond(), (transSeqFrag.getSecond() + length));
		// System.out.println("t_fw_regvec: " + t_fw_regvec);

		Region t_rw_regvec = new Region((transSeqFrag.getSecond() + transSeqFrag.getFirst().length() - length),
				(transSeqFrag.getSecond() + transSeqFrag.getFirst().length()));
		// System.out.println("t_rw_regvec: " + t_rw_regvec);

//		if (mutRead.getSecond().contains(0) || mutRwRead.getSecond().contains(0)) {
//			System.out.println(mutRead);
//			System.out.println(mutRwRead);
//		}
		
		mapInfo.writeMapinfo(id, chr, geneID, transID, fw_regvec, rw_regvec, t_fw_regvec, t_rw_regvec,
				mutRead.getSecond(), mutRwRead.getSecond());

	}

	public static Tuple<String, Integer> generateFragment(String trans) {

		if (trans.length() < length) {
			throw new RuntimeException("Transcript is smaller then read length");
		}

		Random rdm = new Random();
		int fragmentLength = Integer.MAX_VALUE;

		while (fragmentLength > trans.length() || fragmentLength < length) {
			double rand = rdm.nextGaussian() * SD + frlength;
			int fragmentLength1 = (int) Math.round(rand);
			fragmentLength = fragmentLength1 < 0 ? -fragmentLength1 : fragmentLength1;
		}

		int randomStart = 1;
		if (trans.length() != fragmentLength) {
			int possibleStarts = trans.length() - fragmentLength;
			if (possibleStarts > 1) {
				try {

					randomStart = ThreadLocalRandom.current().nextInt(possibleStarts - 1); // MAYBE: possibleStarts +1
					if (randomStart == 0) {
						randomStart++;
					}

				} catch (Exception e) {
					throw new RuntimeException("Error with trans " + trans + " and rdmStart " + randomStart
							+ " possible starts: " + possibleStarts + " fragment length " + fragmentLength);
				}
			} else {
				randomStart = 1;
			}
		}
		// System.out.println(randomStart + " " + trans.length() + " " + trans);
		String fragSeq = trans.substring(randomStart, randomStart + fragmentLength);
		Tuple<String, Integer> fragment = new Tuple<String, Integer>(fragSeq, randomStart);
		return fragment;

	}

	public static String generateRwRead(String fragment) {

		DNAUtils reverser = new DNAUtils();
		String read;

		if (fragment.length() > length) {
			int start = fragment.length() - length;
			String fwRread = fragment.substring(start, fragment.length());
			read = reverser.revcomp(fwRread);
		} else {
			int start = fragment.length() - length;
			String fwRread = fragment.substring(start, fragment.length());
			read = reverser.revcomp(fwRread);
		}
		return read;

	}

	public static String generateRead(String fragment) {
		String read = fragment.substring(0, length);
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

	public static RegionVector getGenomicRegion(String fragment, RegionVector parent, int rdmStart) {
		RegionVector genomicRegions = new RegionVector();

		int FL = fragment.length();
		int distanceTravelled = 0;
		int i = 0;

		while (rdmStart > distanceTravelled) // skip regions if start isnt in it
		{
			distanceTravelled += parent.regions.get(i).getLength() + 1;
			i++;
		}

		int pos1 = 0;

		if (i - 1 == 0) // if start is in first region
		{
			pos1 = parent.regions.get(i - 1).getX1() + rdmStart - 1;
		} else if (i == 0) {
			pos1 = parent.regions.get(i).getX1();
		} else {
			pos1 = parent.regions.get(i - 1).getX1()
					+ (rdmStart - (distanceTravelled - parent.regions.get(i - 1).getLength()));
		}

		if (rdmStart + FL <= distanceTravelled + 1) // if end is in first region
		{
			int pos12 = pos1 + FL;
			Region onlyRegion = new Region(pos1, pos12);
			genomicRegions.addRegion(onlyRegion);

		} else {
			Region firstRegion = new Region(pos1, parent.regions.get(i - 1).getX2() + 1);
			genomicRegions.addRegion(firstRegion);

			try {
				while (rdmStart + FL >= distanceTravelled) // save regions inside of fragment
				{
					distanceTravelled += parent.regions.get(i).getLength();
					if (rdmStart + FL > distanceTravelled) // if FALSE: arrived at last region
					{
						genomicRegions.addRegion(parent.regions.get(i));
					}
					i++;
				}
			} catch (Exception e) {
				throw new RuntimeException("RandomStart:  " + rdmStart + " FragmentL: " + FL + " DT: "
						+ distanceTravelled + "\n" + parent + "\n" + fragment);
			}

			int pos12 = parent.regions.get(i - 1).getX2() - (distanceTravelled - (rdmStart + FL));
			Region lastRegion = new Region(parent.regions.get(i - 1).getX1(), pos12 - 1);

			genomicRegions.addRegion(lastRegion);

		}
		return genomicRegions;
	}
}
