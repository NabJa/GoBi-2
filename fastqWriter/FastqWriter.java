package fastqWriter;

import java.io.BufferedWriter;

public class FastqWriter {

	String id;
	String seq;
	int score;
	String quality;
	
	Writer fw = new FileWriter( filename );
	BufferedWriter wrtr = new BufferedWriter(wrtr);

	public FastqWriter() {
	}

	public FastqWriter(String id, String seq, int score, String quality) {
		this.id = id;
		this.seq = seq;
		this.score = score;
		this.quality = quality;
	}

	public FastqWriter(String id, String seq, int score) {
		this.id = id;
		this.seq = seq;
		this.score = score;
	}

	public FastqWriter(String id, String seq) {
		this.id = id;
		this.seq = seq;
	}

	public void writeFastq(String id, String seq) {
		
	}

}
