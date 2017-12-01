package genomeExtractor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;

import genomicUtils.DNAUtils;
import genomicUtils.GTFReader;
import genomicUtils.Gene;
import genomicUtils.Region;
import genomicUtils.RegionVector;
import genomicUtils.Triplet;

public class GenomeSequenceExtractor {

	// InputStream fastaInputStream;
	File idx;
	File fasta;

	BufferedReader indexReader = null;
	ArrayList<IndexLine> genomicIndex = new ArrayList<IndexLine>();

	BufferedReader countReader = null;
	public ArrayList<Triplet<String, String, Integer>> readcounts = new ArrayList<Triplet<String, String, Integer>>();

	BufferedReader gtfReader = null;
	public HashMap<String, Gene> genes = new HashMap<String, Gene>();

	RandomAccessFile raffasta = null;

	public HashMap<Triplet<String, String, Integer>, String> sequences = new HashMap<Triplet<String, String, Integer>, String>();
	public HashMap<Triplet<String, String, Integer>, RegionVector> transGenomicRegions = new HashMap<Triplet<String, String, Integer>, RegionVector>();

	public HashMap<String, String> transcripts = new HashMap<String, String>();

	public GenomeSequenceExtractor() {
	}

	/**
	 * @param fasta
	 * @param idx
	 */
	public GenomeSequenceExtractor(File fasta, File idx) {

		this.idx = idx;
		this.fasta = fasta;

		try {
			this.raffasta = new RandomAccessFile(fasta, "r");
		} catch (Exception e) {
			throw new RuntimeException("RAF build failed", e);
		}
	}

	/**
	 * 
	 */
	public void readIndex() {
		String line = "";

		try {
			indexReader = new BufferedReader(new FileReader(idx));
			while ((line = indexReader.readLine()) != null) {

				String[] Line = line.split("\t");

				String chr = Line[0];
				long length = Long.valueOf(Line[1]);
				long start = Long.valueOf(Line[2]);
				int lineLength = Integer.parseInt(Line[3]);
				int absoluteLineLength = Integer.parseInt(Line[4]);

				IndexLine indexline = new IndexLine(chr, length, start, lineLength, absoluteLineLength);

				genomicIndex.add(indexline);

			}
		} catch (Exception e) {
			throw new RuntimeException("got error while reading input index files", e);
		}
	}

	/**
	 * 
	 * @param simr
	 */
	public void readCounts(String simr) {
		int lineNumber = 0;

		try {
			countReader = new BufferedReader(new FileReader(simr));
			String line = "";
			while ((line = countReader.readLine()) != null) {
				if (lineNumber != 0) {
					String[] Line = line.split("\t");

					String geneID = Line[0];
					String transID = Line[1];
					Integer count = Integer.parseInt(Line[2]);

					Triplet<String, String, Integer> readcount = new Triplet<String, String, Integer>(geneID, transID,
							count);

					readcounts.add(readcount);
				}
				lineNumber++;
			}
		} catch (Exception e) {
			throw new RuntimeException("got error while reading readcounts in line: " + lineNumber, e);
		}
	}

	/**
	 * @param gtf
	 */
	public void readGTF(String gtf) {
		GTFReader gtfreader = new GTFReader();
		gtfreader.readCDS(gtf);
		genes = gtfreader.getGenes();
	}

	public void getAllSequences() {

		DNAUtils dnau = new DNAUtils();
			
		for (Triplet<String, String, Integer> readcount : readcounts) // readcount = geneID, transID, count
		{
			String splicedTrans = "";
			RegionVector genomicRegions = new RegionVector(readcount.getSecond());

			String chr = genes.get(readcount.getFirst()).geneChr;
			String geneStrand = genes.get(readcount.getFirst()).strand;

			// Gets sequence and genomic position of every region
			for (Region r : genes.get(readcount.getFirst()).transcripts.get(readcount.getSecond()).regions) {

				String seq = "";
				seq = getSequence(chr, r.getX1(), r.getX2(), geneStrand);
				splicedTrans += seq;
				genomicRegions.addRegion(r);
			}
			transGenomicRegions.put(readcount, genomicRegions);
			sequences.put(readcount, splicedTrans);

			
//Only for transcriptom check;
//			if(geneStrand.equals("-")) {
//				System.out.println(readcount.getSecond());
//				splicedTrans = dnau.revcomp(splicedTrans);
//			}
		
		}
		try {
			raffasta.close();
		} catch (Exception e) {
			throw new RuntimeException("got error while closing RAF", e);
		}
	}

	/**
	 * 
	 * @param chr
	 * @param start
	 * @param end
	 */
	public String getSequence(String chr, int start, int end, String strand) {

		IndexLine index = null;
		for (IndexLine idx : genomicIndex) {
			if (idx.getChr().equals(chr)) {
				index = idx;
				break;
			}
		}

		long indexStart = index.getStart();
		int indexLineLen = index.lineLength();
		int indexLineTotal = index.lineTotalLength();
		int breakLen = indexLineTotal - indexLineLen;
		int seqLen = end - start;
		int skipLines = (int) start / indexLineLen;
		int readStart = start + skipLines;

		String seq = "";
		try {

			if (readStart % 61 == 0) {
				raffasta.seek(indexStart + readStart - 2);
			} else {
				raffasta.seek(indexStart + readStart - 1);
			}
	
			byte b;
			while (seq.length() < seqLen) {

				b = raffasta.readByte();

				if ((char) b != '\n') {
					seq += (char) b;
				}

			}

		} catch (Exception e) {
			throw new RuntimeException("got error while reading input raf files", e);
		}
		return seq;
	}

	public void readTranscriptom(String transcriptom) {

		File transFasta = new File(transcriptom);
		BufferedReader bfreader = null;

		try {
			bfreader = new BufferedReader(new FileReader(transFasta));

			String rline = "";
			String transID = "";
			String sequence = "";
			while ((rline = bfreader.readLine()) != null) {
				if (rline.startsWith(">")) {
					transcripts.put(transID, sequence);
					sequence = "";
					int idLen = rline.indexOf(" ");
					transID = rline.substring(1, idLen);
					if(transID.length()>20) {
						System.out.println("Very long transID: " + transID);
					}
				} else {
					sequence += rline;
				}
			}
			transcripts.remove("");
			transcripts.put(transID, sequence);
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			try {
				bfreader.close();
			} catch (Exception e) {
				throw new RuntimeException("got error while closing fastaReader.", e);
			}
		}
	}

	public void checkTrans(String transID, String trans) {
		String correctSeq = transcripts.get(transID);

		if (correctSeq == null) {
			System.out.println("Transcipt not found: "+ transID+ "\n" + trans);
		} else {
			if (!correctSeq.equals(trans)) {
				System.out.println("Correct Sequence: " + transID  +"\n" + correctSeq + "\n" + "Youre Seq:" + "\n" + trans);
			}
		}
	}
}
