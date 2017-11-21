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

		String rwod = od + "/rw.fastq";
		String fwod = od + "/fw.fastq";
		String map = od + "/read.mappinginfo";

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
			
			String transSeq = gse.sequences.get(readcount);
			String transID = readcount.getSecond();
			
			RegionVector transGenomicRegion = gse.transGenomicRegions.get(readcount);

			Gene gene = gse.genes.get(readcount.getFirst());
			String strand = gene.strand;
			String chr = gene.geneChr;
			String geneID = gene.geneID;
			
			for (int i = 0; i < readcount.getThird(); i++) // ... Make count many reads
			{
				String transSeqFrag = generateFragment(transSeq);
				String fragmentRead = generateRead(transSeqFrag);
				Tuple<StringBuilder, ArrayList<Integer>> mutRead = mutate(fragmentRead);
				
				RegionVector fw_regvec;
				RegionVector rw_regvec;
				Region t_fw_regvec;
				Region t_rw_regvec;
				
				
				if (strand.equals("-")) {
					fastqRW.writeFastq(id, mutRead.getFirst().toString(), id, qualityString);
					fastqFW.writeFastq(id, dnautil.revcomp(mutRead.getFirst().toString()), id, qualityString);
					
					RegionVector fw_regvec = getGenomicRegion(mutRead.getFirst().toString(), transGenomicRegion, 3);
					RegionVector rw_regvec = getGenomicRegion();
					Region t_fw_regvec = getGenomicRegion();
					Region t_rw_regvec = getGenomicRegion();
					
				} else {
					fastqFW.writeFastq(id, mutRead.getFirst().toString(), id, qualityString);
					fastqRW.writeFastq(id, dnautil.revcomp(mutRead.getFirst().toString()), id, qualityString);
				
					RegionVector fw_regvec = getGenomicRegion();
					RegionVector rw_regvec = getGenomicRegion();
					Region t_fw_regvec = getGenomicRegion();
					Region t_rw_regvec = getGenomicRegion();
					
				}
				mapInfo.writeMapinfo(id, chr, geneID, transID, fw_regvec, rw_regvec, t_fw_regvec, t_rw_regvec, mutRead.getSecond(), mutRead.getSecond());
				id++;
			}
		}

		//Close all writer
		fastqRW.closeFastq();
		fastqFW.closeFastq();
		mapInfo.closeMapinfo();

		// Program Time
		long endTime1 = System.currentTimeMillis() - zeit;
		System.out.println("Fertig in " + endTime1);

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

	public RegionVector getGenomicRegions(RegionVector genomic, String query, int randStart) {

		RegionVector genomicRegions = new RegionVector();

		int queryLen = query.length();
		int genomicStart = genomic.getX1() + randStart;

		if (genomic.regions.get(0).getLength() >= queryLen) // Fragment is in first Region (=exon)
		{
			Region genomRegion = new Region(genomicStart, genomicStart + queryLen);
			genomicRegions.addRegion(genomRegion);
		} else {
			int i = 0;
			int subTrans = genomic.regions.get(0).getLength() - randStart;
			while (subTrans < queryLen) {
				genomicRegions.addRegion(genomic.regions.get(i));
				subTrans += genomic.regions.get(i).getLength();
				i++;
			}
			subTrans += genomic.regions.get(i).getLength();
			int endPos = genomic.regions.get(i).getX2() - subTrans - queryLen;
			Region genomRegion2 = new Region(genomic.regions.get(i).getX1(), endPos);
			genomicRegions.addRegion(genomRegion2);
		}
		return genomicRegions;
	}
	
	public static RegionVector getGenomicRegion(String fragment, RegionVector parent, int rdmStart) {
		RegionVector genomicRegions = new RegionVector();

		int FL = fragment.length();
		int distanceTravelled = 0;
		int i = 0;

		while (rdmStart > distanceTravelled) // skip regions if start isnt in it
		{
			distanceTravelled += parent.regions.get(i).getLength();
			i++;
		}

		int pos1;
		if (i - 1 == 0) // if start is in first region
		{
			pos1 = parent.regions.get(i - 1).getX1() + rdmStart;
		} else {
			pos1 = parent.regions.get(i - 1).getX1()
					+ (rdmStart - (distanceTravelled - parent.regions.get(i - 1).getLength()) - 1);
		}

		if (rdmStart + FL < distanceTravelled) // if end ist in first region
		{
			int pos12 = pos1 + FL;
			Region onlyRegion = new Region(pos1, pos12);
			genomicRegions.addRegion(onlyRegion);

		} else {
			Region firstRegion = new Region(pos1, parent.regions.get(i - 1).getX2());
			genomicRegions.addRegion(firstRegion);

			while (rdmStart + FL > distanceTravelled) // save regions inside of fragment
			{
				distanceTravelled += parent.regions.get(i).getLength();
				if (rdmStart + FL > distanceTravelled) // if FALSE: arrived at last region
				{
					genomicRegions.addRegion(parent.regions.get(i));
				}
				i++;
			}
			int pos12 = parent.regions.get(i - 1).getX2() - (distanceTravelled - (rdmStart + FL));
			Region lastRegion = new Region(parent.regions.get(i - 1).getX1(), pos12);
			genomicRegions.addRegion(lastRegion);
		}
		return genomicRegions;
	}
	
}
