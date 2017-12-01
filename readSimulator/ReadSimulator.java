package readSimulator;

import java.io.File;
import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import fastqWriter.FastqWriter;
import genomeExtractor.GenomeSequenceExtractor;
import genomicUtils.DNAUtils;
import genomicUtils.Gene;
import genomicUtils.MutationDistribution;
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

//		String filePath = "C:\\Users\\Anja\\Desktop\\GoBi\\GoBi_2\\input\\Homo_sapiens.GRCh37.75.cdna.all.fa";
//		gse.readTranscriptom(filePath);
//		System.out.println("Finished reading Transcriptome Fasta");
		
		gse.readIndex();
		System.out.println("Finished reading Index");
		
		gse.readCounts(readcounts);
		System.out.println("Finished reading counts");

		gse.readGTF(gtf);
		System.out.println("Finished reading GTF");

		gse.getAllSequences();
		System.out.println("Finished reading Genome");
		

//		int check =  0;
//		for(Triplet<String,String, Integer> trpl : gse.sequences.keySet()) {
//			check++;
//			String trans = gse.sequences.get(trpl);
//			String transID =trpl.getSecond();
//			gse.checkTrans(transID, trans);	
//		}
//		System.out.println("Checked all " + check + " transcripts");	
//		System.out.println(gse.sequences.size());
		
		long endTime = System.currentTimeMillis() - zeit;
		System.out.println("Alles eingelesen in: " + endTime);

		File rwod = new File(od, "rw.fastq");
		File fwod = new File(od, "fw.fastq");
		File map = new File(od, "read.mappinginfo");
		
//		File mut = new File(od, "mutDistribution.txt");

		System.out.println("Output destinations: " + rwod + " " + fwod + " " + map);

		// Initialize Output wirter and Utils
		ReadInfoMap mapInfo = new ReadInfoMap(map);
		FastqWriter fastqRW = new FastqWriter(rwod);
		FastqWriter fastqFW = new FastqWriter(fwod);
		
//		DNAUtils dnautil = new DNAUtils();
//		MutationDistribution mutDis = new MutationDistribution(mut);

		// Generate String of length = read length with Score values. Here only I used
		StringBuffer qualityBuffer = new StringBuffer(length);
		for (int i = 0; i < length; i++) {
			qualityBuffer.append("I");
		}
		String qualityString = qualityBuffer.toString();

		// Generate all fragments and reads and write them directly in output
		long id = 0;
		
		int count = 0;
		for (Triplet<String, String, Integer> readcount : gse.readcounts) // For every transcript in readcounts ...
		{
			count ++;
			if(count % 10 == 0) {
				System.out.println("Prcessed "+ count + " transcripts!");
			}
			
			Gene gene = gse.genes.get(readcount.getFirst());
			String strand = gene.strand;
			String chr = gene.geneChr;
			String geneID = gene.geneID;

			if (strand.equals("-")) {

				String transSeq = gse.sequences.get(readcount);
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
//		mutDis.close();

		// Program Time
		long endTime1 = System.currentTimeMillis() - zeit;
		System.out.println("Fertig in " + endTime1);

	}

	public static void generateMinusReads(String transSeq, FastqWriter fastqFW, FastqWriter fastqRW, long id,
			String qualityString, RegionVector transGenomicRegion, ReadInfoMap mapInfo, String chr, String geneID,
			String transID) {

		Tuple<String, Integer> transSeqFrag = generateFragment(transSeq, "-");

		String fragmentRwRead = generateRead(transSeqFrag.getFirst());
		String fragmentRead = generateRwRead(transSeqFrag.getFirst());

		Tuple<StringBuilder, ArrayList<Integer>> mutRead = mutate(fragmentRead);
		Tuple<StringBuilder, ArrayList<Integer>> mutRwRead = mutate(fragmentRwRead);

		fastqFW.writeFastq(id, mutRead.getFirst().toString(), id, qualityString);
		fastqRW.writeFastq(id, mutRwRead.getFirst().toString(), id, qualityString);

		int revStartTrans = transSeq.length() - transSeqFrag.getSecond() - length;
		int revEndTrans = revStartTrans + length;

		Region t_rw_regvec = new Region(revStartTrans, revEndTrans);

		int fwStartTrans = revEndTrans - transSeqFrag.getFirst().length();
		int fwEndTrans = fwStartTrans + length;

		Region t_fw_regvec = new Region(fwStartTrans, fwEndTrans);
		
		RegionVector rw_regvec = getGenomicRegion(mutRwRead.getFirst().toString(), transGenomicRegion,
				transSeqFrag.getSecond() + (transSeqFrag.getFirst().length() - length));
		rw_regvec.id = "rw_regvec";

		RegionVector fw_regvec = getGenomicRegion(mutRead.getFirst().toString(), transGenomicRegion,
				transSeqFrag.getSecond());
		fw_regvec.id = "fw_regvec";

		mapInfo.writeMapinfo(id, chr, geneID, transID, rw_regvec, fw_regvec, t_fw_regvec, t_rw_regvec,
				mutRead.getSecond(), mutRwRead.getSecond());
	}

	public static void generatePlusReads(String transSeq, FastqWriter fastqFW, FastqWriter fastqRW, long id,
			String qualityString, RegionVector transGenomicRegion, ReadInfoMap mapInfo, String chr, String geneID,
			String transID) {

		Tuple<String, Integer> transSeqFrag = generateFragment(transSeq, "+");

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
		Region t_rw_regvec = new Region((transSeqFrag.getSecond() + transSeqFrag.getFirst().length() - length),
				(transSeqFrag.getSecond() + transSeqFrag.getFirst().length()));

		mapInfo.writeMapinfo(id, chr, geneID, transID, fw_regvec, rw_regvec, t_fw_regvec, t_rw_regvec,
				mutRead.getSecond(), mutRwRead.getSecond());
		 
	}

	public static Tuple<String, Integer> generateFragment(String trans, String strand) {

		if (trans.length() < length) {
			throw new RuntimeException("Transcript is smaller then read length");
		}

		Random rdm = new Random();
		int fragmentLength = Integer.MAX_VALUE;

		while (fragmentLength >= trans.length() || fragmentLength <= length) {
			double rand = rdm.nextGaussian() * SD + frlength;
			fragmentLength = (int) Math.floor(rand);
		}

		int randomStart;
		try {
			randomStart = rdm.nextInt(trans.length() - fragmentLength);			
		} catch (Exception e) {
			throw new RuntimeException("Error in fragment length: \n" + "fragmentLength: "  +fragmentLength +"trans: " + trans );
		}

		String fragSeq = trans.substring(randomStart, randomStart + fragmentLength);

		Tuple<String, Integer> fragment = new Tuple<String, Integer>(fragSeq, randomStart);
		return fragment;

	}

	public static String generateRwRead(String fragment) {

		DNAUtils reverser = new DNAUtils();
		String read;

		int start = fragment.length() - length;
		String fwRread = fragment.substring(start, fragment.length());
		read = reverser.revcomp(fwRread);

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
				 mutPos.add(i); int rdmChar = rdm.nextInt(2); 
				 
				 switch (read.charAt(i)) {
				 
				 case 'A': char[] basesA = { 'T', 'G', 'C' }; 
				 readS.setCharAt(i, basesA[rdmChar]);
				 break; 
				 
				 case 'T': char[] basesT = { 'A', 'G', 'C' };
				 readS.setCharAt(i, basesT[rdmChar]);
				 break;
				 
				 case 'G': char[] basesG = { 'T', 'A', 'C' }; 
				 readS.setCharAt(i, basesG[rdmChar]);
				 break; 
				 
				 case 'C': char[] basesC = { 'T', 'G', 'A' }; 
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
		int remainingDist = rdmStart;
		boolean inEnd = false;

		for (Region r : parent.regions) {
			if (remainingDist > r.getLength()) {
				remainingDist -= r.getLength();
				if (inEnd) {
					genomicRegions.addRegion(r);
				}
				continue;
			}
			if (inEnd) {
				Region region = new Region(r.getX1(), r.getX1() + remainingDist);
				if (region.getX1() != region.getX2()) {
					genomicRegions.addRegion(region);
				}
				break;
			}
			if (r.getX2() >= r.getX1() + remainingDist + FL) {
				Region region = new Region(r.getX1() + remainingDist, r.getX1() + remainingDist + FL);
				if (region.getX1() != region.getX2()) {
					genomicRegions.addRegion(region);
				}
				break; // found start/end in same exon
			} else {
				Region region = new Region(r.getX1() + remainingDist, r.getX2());
				if (region.getX1() != region.getX2()) {
					genomicRegions.addRegion(region);
				}
				remainingDist = FL - (r.getX2() - (r.getX1() + remainingDist));
				inEnd = true;
			}
		}
		return genomicRegions;
	}
}
