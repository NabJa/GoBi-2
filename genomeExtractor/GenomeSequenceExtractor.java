package genomeExtractor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;
import genomicUtils.*;

public class GenomeSequenceExtractor {

	// InputStream fastaInputStream;
	File idx;
	File fasta;

	BufferedReader indexReader = null;
	ArrayList<IndexLine> genomicIndex = new ArrayList<IndexLine>();

	BufferedReader countReader = null;
	ArrayList<Triplet<String, String, Integer>> readcounts = new ArrayList<Triplet<String, String, Integer>>();

	BufferedReader gtfReader = null;
	public HashMap<String, Gene> genes = new HashMap<String, Gene>();

	RandomAccessFile raffasta = null;

	// ArrayList<String> sequences = new ArrayList<String>();

	HashMap<Triplet<String, String, Integer>, String> sequences = new HashMap<Triplet<String, String, Integer>, String>();

	/**
	 * 
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
	 * 
	 * @param gtf
	 */
	public void readGTF(String gtf) {
		GTFReader gtfreader = new GTFReader();
		gtfreader.readCDS(gtf);
		genes = gtfreader.getGenes();
	}

	/**
	 * @param geneID
	 * @param transcriptID
	 * @return
	 */
	public Triplet<String, Integer, Integer> getTranscriptPos(String geneID, String transcriptID) {
		Gene gene = genes.get(geneID);
		RegionVector trans = gene.transcripts.get(transcriptID);

		Integer start = trans.getX1();
		Integer end = trans.getX2();
		String chr = gene.geneChr;

		Triplet<String, Integer, Integer> position = new Triplet<String, Integer, Integer>(chr, start, end);
		return position;
	}

	public void getAllSequences() {

		for (Triplet<String, String, Integer> readcount : readcounts) // readcount = geneID, transID, count
		{
			String splicedTrans = "";
			String chr = genes.get(readcount.getFirst()).geneChr;

			for (Region r : genes.get(readcount.getFirst()).transcripts.get(readcount.getSecond()).regions) {
				String seq = getSequence(chr, r.getX1(), r.getX2());
				splicedTrans += seq;
			}
			splicedTrans.trim();

			// if(genes.get(readcount.getFirst()).strand.equals("-")) {
			// DNAUtils reverser = new DNAUtils();
			// splicedTrans = reverser.revcomp(splicedTrans);
			// }

			// System.out.println(readcount.getSecond() + " " + splicedTrans.length() +" "+
			// splicedTrans);
			// if (readcount.getSecond().equals("ENST00000562083")) {
			// System.out.println();
			// System.out.println("Chromosom: " + genes.get(readcount.getFirst()).geneChr);
			// System.out.println("Strand: " + genes.get(readcount.getFirst()).strand);
			// System.out.println(readcount.getSecond() + " " + splicedTrans.length() + " "
			// + splicedTrans);
			// }

			sequences.put(readcount, splicedTrans);
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
	public String getSequence(String chr, int start, int end) {

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
		int seqLen = (end - start + 1);
		int skipLines = (int) Math.floor(start / indexLineLen);
		int readStart = start + skipLines * breakLen;

		String seq = "";
		try {

			raffasta.seek(indexStart + readStart - 1);
			// raffasta.seek(indexStart + start);

			// byte[] buffer = new byte[seqLen];
			// raffasta.read(buffer);
			// String exon = new String(buffer);
			// seq += exon;

			int len = 0;
			while (len < seqLen) {
				int read = raffasta.read();
				char readChar = (char) read;
				if (readChar != '\n') {
					seq += readChar;
					len++;
				}
			}

		} catch (Exception e) {
			throw new RuntimeException("got error while reading input raf files", e);
		}

		seq = seq.replaceAll("\n", "");

		return seq;
	}
}
