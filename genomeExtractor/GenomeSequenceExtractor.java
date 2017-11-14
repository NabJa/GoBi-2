package genomeExtractor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

import genomicUtils.Triplet;
import readSimulator.IndexLine;

public class GenomeSequenceExtractor {

	// InputStream fastaInputStream;
	File idx;
	File fasta;
	
	BufferedReader reader = null;
	HashMap<String, IndexLine> genomicIndex = new HashMap<String, IndexLine>();

	BufferedReader countReader = null;
	ArrayList<Triplet<String, String, Integer>> readcounts = new ArrayList<Triplet<String, String, Integer>>();

	/**
	 * 
	 * @param fasta
	 * @param idx
	 */
	public GenomeSequenceExtractor(File fasta, File idx) {

		this.idx = idx;
		this.fasta = fasta;

		/*
		 * Reads Fasta Index
		 */
		String rline = "";
		try {
			reader = new BufferedReader(new FileReader(idx));
			while ((rline = reader.readLine()) != null) {

				String[] indexLine = rline.split("\t");

				String chr = indexLine[0];
				long length = Long.valueOf(indexLine[1]);
				long start = Long.valueOf(indexLine[2]);
				int lineLength = Integer.parseInt(indexLine[3]);
				int absoluteLineLength = Integer.parseInt(indexLine[4]);

				IndexLine indexline = new IndexLine(length, start, lineLength, absoluteLineLength);

				genomicIndex.put(chr, indexline);

			}
		} catch (Exception e) {
			throw new RuntimeException("got error while reading input index files", e);
		}

	}

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

	public void getSequence(String chr, int start, int end) {

		/**
		 * RAF for big FASTA file
		 */
		try {
			RandomAccessFile raffasta = new RandomAccessFile(fasta, "r");

			raffasta.seek(start);

			int intChar = raffasta.read();
			String character = Character.toString((char) intChar);

			System.out.println(character);

		} catch (Exception e) {
			throw new RuntimeException("got error while reading input raf files", e);
		}

	}

}
