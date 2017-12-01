package fastqWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

public class FastqWriter {

	String id;
	String seq;
	int score;
	String quality;
	
	String outputDestination;	
	FileWriter file;
	BufferedWriter writer;

	public FastqWriter() {
	}

	public FastqWriter(String id, String seq, int score, String quality, String outputDestination) {
		this.id = id;
		this.seq = seq;
		this.score = score;
		this.quality = quality;
		
		this.outputDestination = outputDestination;
		try {
			this.file = new FileWriter(outputDestination); 
			this.writer = new BufferedWriter(file);
		} catch (Exception e) {
			throw new RuntimeException("got error while creating BufferedWriter.", e);
		}
	}

	public FastqWriter(String id, String seq, int score, String outputDestination) {
		this.id = id;
		this.seq = seq;
		this.score = score;
		
		this.outputDestination = outputDestination;
		try {
			this.file = new FileWriter(outputDestination); 
			this.writer = new BufferedWriter(file);
		} catch (Exception e) {
			throw new RuntimeException("got error while creating BufferedWriter.", e);
		}
	}

	public FastqWriter(String id, String seq, String outputDestination) {
		this.id = id;
		this.seq = seq;
		
		this.outputDestination = outputDestination;
		try {
			this.file = new FileWriter(outputDestination); 
			this.writer = new BufferedWriter(file);
		} catch (Exception e) {
			throw new RuntimeException("got error while /.", e);
		} 
	}
	
	public FastqWriter(String outputDestination) {

		this.outputDestination = outputDestination;
		try {
			this.file = new FileWriter(outputDestination); 
			this.writer = new BufferedWriter(file);
		} catch (Exception e) {
			throw new RuntimeException("got error while /.", e);
		}
	}	
	
	public FastqWriter(File outputDestination) {
		try {
			this.file = new FileWriter(outputDestination); 
			this.writer = new BufferedWriter(file);
		} catch (Exception e) {
			throw new RuntimeException("got error while /.", e);
		}
	}	

	public void writeFastq(long id, String seq, long score, String quality) {
		try {
			
			writer.write("@" + id);
			writer.newLine();
			writer.write(seq);
			writer.newLine();
			writer.write("+" + score);
			writer.newLine();
			writer.write(quality);
			writer.newLine();
			
		} catch (Exception e) {
			throw new RuntimeException("Got error while writing fastq." + e);
		}
	}

	public void closeFastq () {
		try {
			writer.close();
		} catch (Exception e) {
			throw new RuntimeException("Got error while closing fastq writer", e);
		}
	}
	
}
