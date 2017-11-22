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
				Tuple<String, Integer> transSeqFrag = generateFragment(transSeq);
				
				Tuple<String, Integer> fragmentRead = generateRead(transSeqFrag.getFirst());
				Tuple<String, Integer> fragmentRwRead = generateRwRead(transSeqFrag.getFirst());
				
				Tuple<StringBuilder, ArrayList<Integer>> mutRead = mutate(fragmentRead.getFirst());
				Tuple<StringBuilder, ArrayList<Integer>> mutRwRead = mutate(fragmentRwRead.getFirst());

				int readRelativeRandomStart = transSeqFrag.getSecond() + fragmentRead.getSecond();

				RegionVector fw_regvec;
				RegionVector rw_regvec;
				Region t_fw_regvec = new Region(0,0);
				Region t_rw_regvec = new Region(0,0);
				
				if (strand.equals("-")) {
					
					fastqRW.writeFastq(id, mutRead.getFirst().toString(), id, qualityString);
					fastqFW.writeFastq(id, mutRwRead.getFirst().toString(), id, qualityString);
					
					fw_regvec = getGenomicRegion(mutRwRead.getFirst().toString(), transGenomicRegion, (transSeqFrag.getSecond() + fragmentRwRead.getSecond()));
					rw_regvec = getGenomicRegion(mutRead.getFirst().toString(), transGenomicRegion, transSeqFrag.getSecond());
					
					 t_fw_regvec.setRegions( (transSeqFrag.getFirst().length() -length) , (transSeqFrag.getSecond() + transSeqFrag.getFirst().length()) );
					 t_rw_regvec.setRegions(transSeqFrag.getSecond(), (transSeqFrag.getSecond() + length));

				} else {

					fastqFW.writeFastq(id, mutRead.getFirst().toString(), id, qualityString);
					fastqRW.writeFastq(id, dnautil.revcomp(mutRead.getFirst().toString()), id, qualityString);

					rw_regvec = getGenomicRegion(mutRwRead.getFirst().toString(), transGenomicRegion, (transSeqFrag.getSecond() + fragmentRwRead.getSecond()));
					fw_regvec = getGenomicRegion(mutRead.getFirst().toString(), transGenomicRegion, transSeqFrag.getSecond());
					
					t_fw_regvec.setRegions( transSeqFrag.getSecond(), (transSeqFrag.getSecond() + length) );
					t_rw_regvec.setRegions(  (transSeqFrag.getFirst().length() -length) , (transSeqFrag.getSecond() + transSeqFrag.getFirst().length()) );

				}
				mapInfo.writeMapinfo(id, chr, geneID, transID, fw_regvec, rw_regvec, t_fw_regvec, t_rw_regvec, mutRead.getSecond(), mutRead.getSecond());
				id++;
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

	public static Tuple<String, Integer> generateFragment(String trans) {

		Random rdm = new Random();
		int fragmentLength = Integer.MAX_VALUE;

		while (fragmentLength > trans.length() ||  fragmentLength < length ) {
			double rand = rdm.nextGaussian() * SD + frlength;
			int fragmentLength1 = (int) Math.round(rand);
			fragmentLength = fragmentLength1 < 0 ? -fragmentLength1 : fragmentLength1; 
		}

//		if (fragmentLength >= length) // fragment length >= read Length
//		{
//		} 
			int possibleStarts = trans.length() - fragmentLength;
			System.out.println("TransLen: " + trans.length() + " " + fragmentLength + " "+possibleStarts);
			int randomStart = ThreadLocalRandom.current().nextInt(possibleStarts);  //MAYBE: possibleStarts +1
			String fragSeq = trans.substring(randomStart, randomStart + fragmentLength);

			Tuple<String, Integer> fragment = new Tuple<String, Integer>(fragSeq, randomStart);
			
			System.out.println("Fragment > Read");
			System.out.println("RandomStart: " + randomStart + "TransLength: " + trans.length() + "FragLength: " + fragSeq.length());

			return fragment;
		
//		else {
//			int possibleStarts = trans.length() - length;
//			int randomStart = ThreadLocalRandom.current().nextInt(0, possibleStarts);  //MAYBE: possibleStarts +1
//			String fragSeq = trans.substring(randomStart, randomStart + length);
//
//			Tuple<String, Integer> fragment = new Tuple<String, Integer>(fragSeq, randomStart);
//
//			System.out.println("Fragment < Read");
//			System.out.println("RandomStart: " + randomStart + "TransLength: " + trans.length() + "FragLength: " + fragSeq.length());
//			
//			return fragment;
//
//		}

	}

	public static Tuple<String, Integer> generateRwRead(String fragment) {

		DNAUtils reverser = new DNAUtils();

		int randomStart = fragment.length() - length;
		String fwRread = fragment.substring(randomStart, randomStart + length);
		String read = reverser.revcomp(fwRread);
		Tuple<String, Integer> Read = new Tuple<String, Integer>(read, randomStart);
		return Read;

	}

	public static Tuple<String, Integer> generateRead(String fragment) {
		// int possibleStarts = fragment.length() - length;
		// int randomStart = ThreadLocalRandom.current().nextInt(possibleStarts + 1);
		// String read = fragment.substring(randomStart, randomStart + length);
		int randomStart = 1;
		String read = fragment.substring(0, length - 1);

		Tuple<String, Integer> Read = new Tuple<String, Integer>(read, randomStart);

		return Read;

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

		System.out.println();
		System.out.println(parent.regions);
		System.out.println(fragment);
		System.out.println(parent.regions.size());
		System.out.println("Parent RegionL: " + parent.getRegionLength());
		System.out.println("Rondom Start: " + rdmStart);

		
		int FL = fragment.length();
		int distanceTravelled = 0;
		int i = 0;

		while (rdmStart >= distanceTravelled) // skip regions if start isnt in it
		{
			distanceTravelled += parent.regions.get(i).getLength();
			i++;
		}

		int pos1 = 0;

		if (i - 1 == 0) // if start is in first region
		{
			pos1 = parent.regions.get(i - 1).getX1() + rdmStart;
		} else if(i == 0) {
			pos1 = parent.regions.get(i).getX1();
		} else {
			pos1 = parent.regions.get(i - 1).getX1()
					+ (rdmStart - (distanceTravelled - parent.regions.get(i - 1).getLength()));
		}

		if (rdmStart + FL < distanceTravelled) // if end is in first region
		{
			int pos12 = pos1 + FL;
			Region onlyRegion = new Region(pos1, pos12);
			genomicRegions.addRegion(onlyRegion);

		} else {
			Region firstRegion = new Region(pos1, parent.regions.get(i - 1).getX2());
			genomicRegions.addRegion(firstRegion);

			System.out.println(i);
			
			while (rdmStart + FL > distanceTravelled + i) // save regions inside of fragment
			{
				distanceTravelled += parent.regions.get(i).getLength();
				if (rdmStart + FL > distanceTravelled) // if FALSE: arrived at last region
				{
					genomicRegions.addRegion(parent.regions.get(i));
				}
				i++;
			}
			if (i != parent.regions.size()) {

				int pos12;
				Region lastRegion = new Region();
				if (distanceTravelled != (rdmStart + FL)) {
					pos12 = parent.regions.get(i - 1).getX2() - (distanceTravelled - (rdmStart + FL));
					lastRegion = new Region(parent.regions.get(i - 1).getX1(), pos12 + 1);
				} else {
					pos12 = parent.regions.get(i).getX2() - distanceTravelled;
					lastRegion = new Region(parent.regions.get(i).getX1(), pos12 + 1);
				}
				genomicRegions.addRegion(lastRegion);
			}
		}
		return genomicRegions;
	}

}
