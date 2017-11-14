package readSimulator;

public class IndexLine {

	long start;
	long length;
	int lineLength;
	int lineTotalLength;

	public IndexLine(long start, long length, int lineLength, int lineTotalLenth) {
		this.start = start;
		this.length = length;
		this.lineLength = lineLength;
		this.lineTotalLength = lineTotalLenth;
	}

	public long getStart() {
		return start;
	}

	public long getLength() {
		return length;
	}

	public int lineLength() {
		return lineLength;
	}

	public int lineTotalLength() {
		return lineTotalLength;
	}

		@Override 
		public String toString() {
			return "" + start + "\t" + length + "\t" + lineLength + "\t" + lineTotalLength;
		}
}
